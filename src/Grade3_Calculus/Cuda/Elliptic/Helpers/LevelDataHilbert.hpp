#pragma once

#include "FunctorsLibrary.hpp"
#include "LevelData.H"
#include "isHilbert.hpp"
namespace Elliptic
{
namespace Helpers
{
template <typename T>
class LevelDataHilbert : public LevelData<T>
{
  static_assert(isHilbert<T>::value,
                "T must be a Hilbert Space. See isHilbert.hpp for the members "
                "that must be provided by T.");

public:
  typedef typename T::MyMemory MyMemory;

public:
  // standard constructors
  LevelDataHilbert() : LevelData<T>() {}
  LevelDataHilbert(const DisjointBoxLayout grids, int ncomp,
                   const IntVect ghosts)
      : LevelData<T>(grids, ncomp, ghosts) {}
  LevelDataHilbert(const DisjointBoxLayout grids, int ncomp)
      : LevelData<T>(grids, ncomp) {}
  // destructor for good measure
   virtual ~LevelDataHilbert() {}

  // these constructors are used to transfer a LD from GPU to CPU or viceversa.
  // Note that this is an exact deep copy.
  template <class S,
            typename std::enable_if<
                std::is_same<typename S::MyMemory, cusp::host_memory>::value &&
                    !std::is_same<S, T>::value,
                S>::type * = nullptr>

  LevelDataHilbert(const S &a_src)
      : LevelDataHilbert<CudaFArrayBox<cusp::device_memory>>(
            a_src.getBoxes(), a_src.nComp(), a_src.ghostVect())
  {
    auto dit = this->getBoxes().dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      thrust::copy(a_src[dit].begin(), a_src[dit].end(), (*this)[dit].begin());
    }
  }
  template <class S,
            typename std::enable_if<std::is_same<typename S::MyMemory,
                                                 cusp::device_memory>::value &&
                                        !std::is_same<S, T>::value,
                                    S>::type * = nullptr>
  LevelDataHilbert(const S &a_src)
      : LevelDataHilbert<CudaFArrayBox<cusp::host_memory>>(
            a_src.getBoxes(), a_src.nComp(), a_src.ghostVect())
  {
    auto dit = this->getBoxes().dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      thrust::copy(a_src[dit].begin(), a_src[dit].end(), (*this)[dit].begin());
    }
  }
  //  template <class S, bool B>
  //  LevelDataHilbert(const S &a_src) { static_assert("LevelDataHilbert copy
  //  constructor called with wrong type"); }
  // template <class S, typename std::enable_if<std::is_same<S, T>::value,
  // S>::type* = nullptr>
  //    LevelDataHilbert(const S &a_src) {
  //    static_assert(!std::is_same<S,T>::value,"LevelDataHilbert copy
  //    constructor called with wrong type"); }

  Real dotProduct(const LevelDataHilbert<T> &a1) const
  {
    const DisjointBoxLayout &grids = a1.getBoxes();
    DataIterator dit = grids.dataIterator();
    Real val = 0.0;
    CH_assert(this->getBoxes() == a1.getBoxes());
    for (dit.reset(); dit.ok(); ++dit)
    {

      val += a1[dit].dotProduct((*this)[dit]);
    }

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalVal = 0.0;
    int result = MPI_Allreduce(&val, &globalVal, 1, MPI_CH_REAL, MPI_SUM,
                               Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in "
                    "ZeroAvgConstInterpPS::prolongIncrement");
    }

#else
    Real globalVal = val;
#endif

    return globalVal;
  }

  Real norm() const
  {
    const DisjointBoxLayout &grids = this->getBoxes();
    return sqrt(this->dotProduct(*this) / grids.numCells());
  }

  void incr(const LevelDataHilbert<T> &a1, const Real scale)
  {
    CH_assert(this->getBoxes() == a1.getBoxes());
    const DisjointBoxLayout &grids = a1.getBoxes();
    DataIterator dit = grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {

      (*this)[dit].plus(a1[dit], scale);
    }
  }
  void scale(const Real scale)
  {
    const DisjointBoxLayout &grids = this->getBoxes();
    DataIterator dit = grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      (*this)[dit] *= scale;
    }
  }
  void setToZero()
  {
    const DisjointBoxLayout &grids = this->getBoxes();
    DataIterator dit = grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      (*this)[dit].setVal(0.0);
    }
  }
  void axpy(const LevelDataHilbert<T> &a_x, const LevelDataHilbert<T> &a_y,
            const Real a, const Real b)
  {
    CH_assert(a_x.getBoxes() == a_y.getBoxes());
    CH_assert(a_y.getBoxes() == (*this).getBoxes());
    CH_assert(a_x.nComp() == this->nComp());
    CH_assert(a_y.nComp() == this->nComp());

    DataIterator dit = this->dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      //  FABAlgebra::axby(a_lhs[dit], a_x[dit], a_y[dit], a_a, a_b);
      Box box = a_x[dit].box();
      box &= a_y[dit].box();
      (*this)[dit].axby(a_x[dit], a_y[dit], box, a, b);
    }
  }
  void create(const LevelData<T> &a_x)
  {
    this->define(a_x.getBoxes(), a_x.nComp(), a_x.ghostVect());
  }

  void assign(const LevelData<T> &a_x)
  {
    CH_assert(a_x.getBoxes() == this->getBoxes());
    CH_assert(a_x.nComp() == this->nComp());

    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      (*this)[dit].copy(a_x[dit]);
    }
  }
  virtual void fillRandom(Real min, Real max)
  {
    Functors::PseudoNumberGenerator<Real> rnd(min, max);
    const DisjointBoxLayout &grids = this->getBoxes();
    DataIterator dit = this->dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      thrust::counting_iterator<unsigned int> index_sequence_begin(
          grids.index(dit()));
      thrust::transform(index_sequence_begin,
                        index_sequence_begin + (*this)[dit].box().numPts(),
                        (*this)[dit].m_view.begin(), rnd);
    }
  }
  Real removeAverage(Real dV)
  {
    LevelDataHilbert<T> VolEl(this->getBoxes(),1);
    DataIterator dit=this->getBoxes().dataIterator();
    for (dit.reset();dit.ok(); ++dit)
    {VolEl[dit].setVal(dV);}
    Real r = this->dotProduct(VolEl);
    
    for (dit.reset();dit.ok(); ++dit)
    {VolEl[dit].setVal(r);}
    
    this->incr(VolEl,-1.0);

    return r;
  }
};
} // namespace Helpers
} // namespace Elliptic
