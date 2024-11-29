/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2024 Thomas Jefferson University and Arizona State University
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/MUON-CFD/somar.
 ******************************************************************************/

//#define GlobalOps
#ifndef __ELLIKIT_H__
#define __ELLIKIT_H__
// #ifdef NDEBUG
// #undef NDEBUG
// #endif
#include <assert.h>
#include <cusp/array1d.h>
#include <cusp/blas/blas.h>
#include <cusp/csr_matrix.h>
#include <cusp/elementwise.h>
#include <cusp/functional.h>
#include <cusp/linear_operator.h>
#include <cusp/multiply.h>
#include <cusp/print.h>
#include <cusp/transpose.h>
//#include <cusp/gallery/poisson.h>
// include cusp relaxators header file
#include <cusp/relaxation/gauss_seidel.h>
#include <cusp/relaxation/jacobi.h>
#include <cusp/relaxation/sor.h>
#include <cusp/relaxation/polynomial.h>
//#include <cusp/krylov/cg.h>
//#include <cusp/krylov/bicgstab.h>
//#include <cusp/krylov/bicg.h>
#ifdef USE_CUSPARSE
#include "cusparse_v2.h"
#endif
#include <cuda_runtime.h>
//#include <chrono>
#include <cusp/precond/diagonal.h>
#include <typeinfo>
//#include <cusp/precond/aggregation/smoothed_aggregation.h>
//#include <cusp/precond/ainv.h>
#include <map>
#include <string>
#include "BCTools.H"
#include "CudaFArrayBox.hpp"
#include "IO.H"
#include "Operator.hpp"
#include "VectorCalculusKit.hpp"
#include "icecream.hpp"
#if CH_USE_PYTHON
#   include "PyGlue.H"
#endif

//// This class uses ElliKit.py to generate the sparse vector representation of
/// the elliptic operators This class is intended to be used as part of a
/// LayoutData, so we need to keep the constructor trivial and do the real work
/// on define

namespace Elliptic
{
namespace PythonInterface
{
/*!  The class provide access to three structures: Operator, Relaxator and
Preconditioner -Operator encapsulates the Linear operators, defined via a sparse
matrix. It provides the operator (x,y) which returns Lx in y. -Relaxator
encapsulates the relaxator (JA, GS or SOR) associated to the operator L
 - Preconditioner encapsulates the diagonal preconditioner associated to the
operator L

 It depends on 4 template parameters:
 - IntType, the type used to label rows and columns (usually int)
 - ValueType, the type of the entries in the matrix (usually double)
 - InputsMemory, where the inputs are stored, which can be host_memory or
device_memory
 - OwnMemory,  where the representation of the matrix operator is stored,
host_memory or device_memory

 Note that if OwnMemory is different than InputsMemory, the members will copy
data in and out of GPU/CPU as required.

 The sparse representation is obtained invoking the appropriate Python module.


 The python functions are called by the ElliKit::generate(bool bc) function. Its
job is to pass the box, domain and boundary conditions to the python scripts,
and collect the vectors that specify the CSR representation of the operators. It
also defines the appropriate relaxator and preconditioner.

 Note that the Python script generates several operators, all of which are
harvested by this class. They are -Laplacian: This is the laplacian. On a given
boundary, it can have three BC: Neumann, Dirichlet or None. If None, the field
is supposed to have ghost points, and the Laplacian extends to those points. If
Dirichlet or Neumann, the ghost point is virtually filled by symmetry or
antisymmetry as appropriate, and so the field does not need to have ghost points
in that direction. This matrix is not generally square unless all BCs are
Neumann or Dirichlet. -LaplacianNG: This is the Laplacian restricted to the
interior points. This is a square matrix -LaplacianOG: If PhiNG is the
restriction of phi to the interior points, and Phi is the full field on the
interior+ghost (where needed by the None Bc), then Laplacian Phi=LaplacianOG
Phi+LaplacianNG PhiNG. Note that LaplacianOG is not square, and for the most
part is a matrix of zeros
- The divergence operators: these are really one sided differences in the three
directions. It also defines the corresponding DMinusOne operators. The Python
Script in turn is an interface to the ElliKit.py class in Python, where the
heavy lifting is done.
*/
template <typename IntType, typename ValueType, typename InputsMemory,
          typename OwnMemory, typename BCType>
class ElliKit
{
  // first some typedefs to make like a little easier
  // public:
  //   typedef cusp::host_memory OwnMemory;
public:
  typedef typename cusp::csr_matrix<IntType, ValueType, OwnMemory> Matrix;
  typedef typename Elliptic::CudaInterface::Operator<IntType, ValueType,
                                                     InputsMemory, OwnMemory>
      Operator;

private:
  typedef typename cusp::relaxation::gauss_seidel<ValueType, OwnMemory> GS;
  typedef typename cusp::relaxation::jacobi<ValueType, OwnMemory> JA;
  typedef typename cusp::relaxation::sor<ValueType, OwnMemory> SOR;
  typedef typename cusp::relaxation::polynomial<ValueType, OwnMemory> PR;
  typedef typename CudaBaseFab<ValueType, InputsMemory>::Array InputArr;
  typedef typename CudaBaseFab<ValueType, InputsMemory>::View InputView;
  typedef typename cusp::array1d<ValueType, OwnMemory> MyArr;
  typedef BaseFab<ValueType> BaseFAB;
  typedef CudaFArrayBox<InputsMemory> DiagType;
  // typedef typename CudaBaseFab<ValueType, cusp::host_memory> GeoArray; //
  // they reside in CPU
  typedef typename std::vector<IntType> PyIntArray;
  typedef typename std::vector<ValueType> PyRealArray;

#ifndef USE_CUSPARSE
  typedef int cudaStream_t;
#endif
public:
  // the following is a hack so that arrays of the same type as MyArr can be
  // inferred using declval
  MyArr myarr;
  //! The generic interface for relaxators
  struct Relaxator
  {
    Relaxator(const Matrix *Mptr) : m_MatrixPtr(Mptr) {}
    virtual ~Relaxator() {}
    //! Applies the relaxator to the problem Ax=b. Purely virtual, must be overriden.
    virtual void operator()(const MyArr &b, MyArr &x) = 0;

  protected:
    const Matrix *m_MatrixPtr; /*!< Stores a pointer to the Operator matrix*/
  };

  /*!
Implementation of the Jacobi Relaxator
*/
  struct JacobiRelax : public Relaxator
  {
    //! Constructor. Requires a pointer to the matrix and the relaxation parameter Omega
    JacobiRelax(const Matrix *Mptr, const Real Omega)
        : Relaxator(Mptr), m_Relax(*Mptr, Omega) {}

    ~JacobiRelax() {}

    void operator()(const MyArr &b, MyArr &x) { m_Relax(*(this->m_MatrixPtr), b, x); }

  private:
    JA m_Relax;
  };
  /*!
Implementation of the Successive Over Relaxation relaxer
*/
  struct SORRelax : public Relaxator
  { //! Constructor. Requires a pointer to the matrix and the relaxation parameter Omega
    SORRelax(const Matrix *Mptr, const Real Omega)
        : Relaxator(Mptr), m_Relax(*Mptr, Omega) {}
    ~SORRelax() {}

    void operator()(const MyArr &b, MyArr &x) { m_Relax(*(this->m_MatrixPtr), b, x); }

  private:
    SOR m_Relax;
  };
  /*!
Implementation of the Gauss Seidel relaxer
*/
  struct GSRelax : public Relaxator
  { //! Constructor. Requires a pointer to the matrix. Note that Omega is not needed.
    GSRelax(const Matrix *Mptr)
        : Relaxator(Mptr), m_Relax(*Mptr) {}
    ~GSRelax() {}

    void operator()(const MyArr &b, MyArr &x) { m_Relax(*(this->m_MatrixPtr), b, x); }

  private:
    GS m_Relax;
  };
  /*!
Implementation of a polynomial relaxer. Experimental
*/
  struct PolynomialRelax : public Relaxator
  { //! Constructor. Requires a pointer to the matrix.
    /*!
    rho is the ritz spectral radius. If negative (default is -1)
    the constructor will calculate it using 8 iterations.
    */
    PolynomialRelax(const Matrix *Mptr, float Rho = -1.0) : Relaxator(Mptr)
    {
      float rho = Rho < 0.0 ? cusp::eigen::ritz_spectral_radius(*(this->m_MatrixPtr), 8, true) : Rho;
      cusp::relaxation::detail::chebyshev_polynomial_coefficients(rho, this->m_coefficients);
      m_RelaxPtr = new PR(*(this->m_MatrixPtr), this->m_coefficients);
    }
    ~PolynomialRelax()
    {
      delete m_RelaxPtr;
      this->m_coefficients.resize(0);
    }
    void operator()(const MyArr &b, MyArr &x) { m_RelaxPtr->operator()(*(this->m_MatrixPtr), b, x); }

  private:
    cusp::array1d<float, cusp::host_memory> m_coefficients;

    PR *m_RelaxPtr;
  };
  //! Generic interface for preconditioner
  struct Preconditioner
  {
    Preconditioner(const Matrix *Mptr) : m_MatrixPtr(Mptr) {}
    virtual ~Preconditioner() {}
    virtual void operator()(const MyArr &b, MyArr &x) = 0;

  protected:
    const Matrix *m_MatrixPtr;
  };
  //! Jacobi Preconditioner
  /*!
  Preconditions with the inverse of the diagonal elements.
  */
  struct JacobiPreconditioner : public Preconditioner
  {
    JacobiPreconditioner(const Matrix *Mptr)
        : Preconditioner(Mptr), m_Precond(*Mptr) {}
    ~JacobiPreconditioner() {}
    void operator()(const MyArr &b, MyArr &x) { m_Precond(b, x); }

  private:
    cusp::precond::diagonal<ValueType, OwnMemory> m_Precond;
  };

public:
  /*!
Standard constructor, required by LayoutData. Does nothing. All the heavy lifting is done in define
*/
  ElliKit() : m_isDefined(false) {}

  ~ElliKit()
  {
    for (auto &M : m_Operators)
    {
      if (M.second != nullptr)
        delete M.second;
    }
    m_Operators.clear();
    for (auto &M : m_Relaxators)
    {
      if (M.second != nullptr)
        delete M.second;
    }
    m_Relaxators.clear();
    for (auto &M : m_Preconditioners)
    {
      if (M.second != nullptr)
        delete M.second;
    }
    m_Preconditioners.clear();
    for (auto &M : m_Diagonals)
    {
      if (M.second != nullptr)
        delete M.second;
    }

    m_Diagonals.clear();
    m_bcDesc.clear();
    m_opNames.clear();
  }

  Box getBox() const { return this->m_box; }

  Box getBoxIntOnly() const { return this->m_boxIntOnly; }

  void operatorsNeeded(const std::vector<std::string> &opNames)
  {
    this->m_opNames = opNames;
  }
  std::vector<std::string> operatorsProvided() const
  {
    return this->m_opNames;
  }

  void define(const RealVect &dx, const BCType *a_BCD, const Box &a_box,
              const DisjointBoxLayout &a_grids, const Real Omega,
              const ProblemDomain &domain,
              const IntVect a_refRatio = IntVect::Unit,
              const int level = 0,
              const cudaStream_t StreamId = 0,
              const std::string PyModuleName = "ElliKitWrapper");

  inline bool isDefined() const { return m_isDefined; }

  void Relax(std::string Key, const CudaBaseFab<ValueType, InputsMemory> &a_b,
             CudaBaseFab<ValueType, InputsMemory> &a_x) const
  {
    if (a_b.box() != this->m_boxIntOnly) // in case they are not on the same
                                         // box, we restrict it.
    {
      CudaBaseFab<ValueType, InputsMemory> b_local(m_boxIntOnly, a_b.nComp());
      b_local.copy(a_b);
      this->Relax(Key, b_local, a_x);
    }
    assert(a_b.nComp() == a_x.nComp());
    assert(a_b.nComp() == 1);
    bool useViews = !std::is_same<OwnMemory, InputsMemory>::value ||
                    a_b.isAliased() || a_x.isAliased();
    useViews ? this->Relax(Key, a_b.m_view, a_x.m_view)
             : this->Relax(Key, a_b.m_varr, a_x.m_varr);
  }

  void Relax(std::string Key, const InputView &a_b, InputView &a_x) const
  {
    assert(m_Relaxators.at(Key) != nullptr);
    Relaxator &M = *(m_Relaxators.at(Key));

    const MyArr b(a_b.begin(), a_b.end());
    MyArr x(a_x.begin(), a_x.end());
    M(b, x);
    thrust::copy(x.begin(), x.end(), a_x.begin());
  }

  void Relax(std::string Key, const InputArr &a_b, InputArr &a_x) const
  {
    assert(m_Relaxators.at(Key) != nullptr);
    MyArr b(a_b.begin(), a_b.end());
    MyArr x(a_x.begin(), a_x.end());
    (*(m_Relaxators.at(Key)))(b, x);
    thrust::copy(x.begin(), x.end(), a_x.begin());
  }

  template <typename C_, typename M_ = OwnMemory>
  auto ApplyOperatorAndAccumulate(std::string Key, const C_ &a_x, const C_ &a_y,
                                  Real alpha, Real beta) ->
      typename std::enable_if<std::is_same<M_, cusp::device_memory>::value &&
                                  std::is_same<Real, double>::value,
                              void>::type
  {
    (*(m_Operators.at(Key)))(a_x, a_y, alpha, beta);
  }
  void ApplyOperator(std::string Key, const MyArr &a_x, MyArr &a_y)
  {
    // first check on my own ops
    if (m_Operators.count(Key) == 1)
    {
      assert(m_Operators.at(Key) != nullptr);
      (*(m_Operators.at(Key)))(a_x, a_y);
    }
    else // it must be in Vector Calculus
    {
      this->m_VCK.ApplyOperator(Key, a_x, a_y);
    }
  }
  void ApplyOperator(std ::string Key, const DiagType &a_x,
                     DiagType &a_y) const
  {
    Interval I(0, a_x.nComp() - 1);

    this->ApplyOperator(Key, a_x, a_y, I);
  }

  void ApplyOperator(std::string Key, const DiagType &a_x, DiagType &a_y,
                     Interval I) const
  {
    if (m_Operators.count(Key) == 1)
    {
      assert(m_Operators.at(Key) != nullptr);
      (*(m_Operators.at(Key)))(a_x, a_y, I);
    }
    else
    {
      this->m_VCK.ApplyOperator(Key, a_x, a_y, I);
    }
  }

  void Precond(std::string Key, const DiagType &a_x, DiagType &a_y) const
  {
    assert(m_Diagonals.at(Key) != nullptr);
    DiagType x(m_boxIntOnly, a_x.nComp());
    x.copy(a_x);
    x.divide(*(m_Diagonals.at(Key)));
    a_y.copy(x);
  }
  Matrix &retrieveOperatorMatrix(std::string Key) const
  {
    return *(m_Operators.at(Key)->getMatrixPtr());
  }
  Operator &retrieveOperator(std::string Key) const
  {
    return *(m_Operators.at(Key));
  }
  void deleteOperator(std::string Key)
  {
    delete m_Operators.at(Key);
    m_Operators.erase(Key);
  }

  // data members
public:
  std::map<std::string, DiagType *> m_Diagonals;

private:
  void generate(
      bool bc = false); // the actual wrapper around the Python script

  std::map<std::string, Operator *> m_Operators;             /*!< contains pointers to the operator matrices*/
  RealVect m_dx;                                             ///< Grid spacing
  IntVect m_refRatio;                                        ///< Refinement ratio used to define the Restriction operator
  Box m_box, m_boxIntOnly;                                   ///< The over which the operator input (m_box) and output (m_boxIntOnly) are defined.
  ProblemDomain m_domain;                                    ///< The global domain which contains m_box.
  std::map<std::string, Relaxator *> m_Relaxators;           ///< Map of pointers to the relaxators
  std::map<std::string, Preconditioner *> m_Preconditioners; ///<Map of pointers to the preconditioners
  const BCType *m_BCD;                                       ///< pointer to the class that contains ::getBCType
  DisjointBoxLayout m_grids;                                 ///<  the DisjointBoxLayout to which m_box belongs.
  cudaStream_t m_StreamId;                                   ///< only used if cuda multiply is used instead of cusp::multiply
  bool m_isDefined;
  Real m_Omega;             ///< relaxation parameter used by SOR and JA relaxators
  std::string m_MyPyModule; ///< the name of the python module that we use
  int m_level = 0;          ///< The level of this solver.

public:
  std::vector<int> m_bcDesc; ///< stores the boundary conditions to be passed to the python script.

private:
  /*! stores the names of the operators we need. Must be a subset of the list of
   names for the Operator CSRData in ElliKit.py note that for the moment the
   python script generates all of them. Can be overwritten */
  std::vector<std::string> m_opNames{
      "Laplacian", "LaplacianDMinusOne", "LaplacianNG",
      "LaplacianOG", "LaplacianNGDMinusOne", "LaplacianOGDMinusOne",
      "Restriction", "Prolongation",
      "GridToDomain",
      "DomainToGrid"};
  /*! stores the names of the basic calculus we need. Must be a subset of the list of
   names for the Operator CSRData in ElliKit.py note that for the moment the
   python script generates all of them. */
  std::vector<std::string> m_VectorCalculusOperations
  {
    "Divergence x", "Divergence y", "Gradient x",
        "Gradient y"
#if CH_SPACEDIM > 2
        ,
        "Divergence z", "Gradient z"
#endif
  };
  /*! holds the class that handles vector calculus operators. Other operators
  are built on top of these.
  */
  Elliptic::PythonInterface::VectorCalculusKit<IntType, ValueType, InputsMemory,
                                               OwnMemory, BCType>
      m_VCK;
  // GeoArray m_... elements of the metric...
};
/// defines the class
template <typename IntType, typename ValueType, typename InputsMemory,
          typename OwnMemory, typename BCType>
void ElliKit<IntType, ValueType, InputsMemory, OwnMemory, BCType>::define(
    const RealVect &a_dx, const BCType *a_BCD, const Box &a_grid,
    const DisjointBoxLayout &a_grids, const Real a_Omega,
    const ProblemDomain &a_domain, const IntVect a_refRatio, const int level,
    const cudaStream_t StreamId, const std::string PyModuleName)
{
  m_dx = a_dx;
  m_refRatio = a_refRatio;
  m_box = a_grid;
  m_boxIntOnly = a_grid;
  m_MyPyModule = PyModuleName;
  m_Omega = a_Omega;
  m_domain = a_domain;
  m_BCD = a_BCD;
  m_grids = a_grids;
  m_StreamId = StreamId;
  m_level = level;
  // a_bcPtr->getBCDescriptor();
  this->m_VCK.define(a_dx, a_BCD, a_grid, a_grids, a_domain, a_refRatio, false,
                     this->m_VectorCalculusOperations, level, StreamId);
  this->generate(true);
  m_isDefined = true;
}

template <typename IntType, typename ValueType, typename InputsMemory,
          typename OwnMemory, typename BC>
void ElliKit<IntType, ValueType, InputsMemory, OwnMemory, BC>::generate(
    bool bc)
{
  std::vector<std::string> *OPN;
  OPN = &m_opNames;
  if (!bc)
  {
    m_box = m_boxIntOnly;
    m_box.grow(1);
    m_bcDesc.resize(SpaceDim);
    for (auto &i : m_bcDesc)
    {
      i = 0;
    }
  }
  else
  {
    m_bcDesc.resize(0);
    m_box = m_boxIntOnly;
    m_BCD->setBoxAndBCDescription(m_grids, m_boxIntOnly, m_domain, m_level, m_box, m_bcDesc);
  }

  if (!this->m_VCK.isDefined() || this->m_VCK.isPythonMatrixCleared())
    Py::PythonFunction(m_MyPyModule, "GenerateMatrix", m_dx, m_box,
                          m_refRatio, m_domain.domainBox(), m_bcDesc);

  for (const auto &OperatorName : *OPN)
  {
    int Nrows;
    int Ncolumns;
    int nnz;
    if (OperatorName != "Laplacian")
    {
      nnz = Py::PythonReturnFunction<int>(m_MyPyModule, "GetMatrixnnz",
                                             OperatorName);
      Nrows = Py::PythonReturnFunction<int>(
          m_MyPyModule, "GetMatrixRows", OperatorName);
      Ncolumns = Py::PythonReturnFunction<int>(
          m_MyPyModule, "GetMatrixColumns", OperatorName);
      {
        PyIntArray IndPtr(Nrows + 1);
        PyIntArray Col(nnz);
        PyRealArray Data(nnz);
        auto N = Py::PythonReturnFunction<int>(
            m_MyPyModule, "GetMatrixData", OperatorName, IndPtr, Col, Data);
        // Allocate space for a new matrix, and store the pointer in the map

        m_Operators.insert(std::pair<std::string, Operator *>(
            OperatorName, new Operator(Nrows, Ncolumns, nnz, m_StreamId)));
        Operator &A = *m_Operators.at(OperatorName); // get a local reference

        A.setOperator(IndPtr, Col, Data);

        if (OperatorName == "Restriction")
        {
          Box m_crseBox = m_boxIntOnly;
          m_crseBox.coarsen(m_refRatio);
          A.setGrids(m_boxIntOnly, m_crseBox);
        }
        else if (OperatorName == "Prolongation")
        {
          Box m_crseBox = m_boxIntOnly;
          m_crseBox.coarsen(m_refRatio);

          A.setGrids(m_crseBox, m_boxIntOnly);
        }
        else if (OperatorName == "DomainToGrid")
        {
          A.setGrids(m_domain.domainBox(), m_box);
        }
        else if (OperatorName == "GridToDomain")
        {
          A.setGrids(m_boxIntOnly, m_domain.domainBox());
        }
        else if (OperatorName == "LaplacianNG")
        {
          A.setGrids(m_boxIntOnly, m_boxIntOnly);
        }
        //  else if (OperatorName == "Divergence x") {
        //   Box edgeBox(surroundingNodes(m_boxIntOnly, 0));
        //   A.setGrids(edgeBox, m_boxIntOnly);
        // } else if (OperatorName == "Divergence y") {
        //   Box edgeBox(surroundingNodes(m_boxIntOnly, 1));
        //   A.setGrids(edgeBox, m_boxIntOnly);
        // } else if (OperatorName == "Divergence z") {
        //   Box edgeBox(surroundingNodes(m_boxIntOnly, 2));
        //   A.setGrids(edgeBox, m_boxIntOnly);
        // } else if (OperatorName == "Gradient x") {
        //   Box edgeBox(surroundingNodes(m_boxIntOnly, 0));
        //   A.setGrids(m_box, edgeBox);
        // } else if (OperatorName == "Gradient y") {
        //   Box edgeBox(surroundingNodes(m_boxIntOnly, 1));
        //   A.setGrids(m_box, edgeBox);
        // } else if (OperatorName == "Gradient z") {
        //   Box edgeBox(surroundingNodes(m_boxIntOnly, 2));
        //   A.setGrids(m_box, edgeBox);
        // }
        else
        {
          A.setGrids(m_box, m_boxIntOnly);
        }
      }
    }
    else
    {
      Nrows = m_boxIntOnly.numPts();
      Ncolumns = m_box.numPts();
      nnz = 0;
      this->m_Operators.insert(std::pair<std::string, Operator *>(
          OperatorName, new Operator(Nrows, Ncolumns, nnz, m_StreamId)));
      auto &Op = *(this->m_Operators.at(OperatorName));
      auto &Laplacian = *(Op.getMatrixPtr());
      Op.setGrids(m_box, m_boxIntOnly);
      switch (SpaceDim)
      {
      case 3:
      {
        auto &D = this->m_VCK.retrieveOperatorMatrix("Divergence z");
        auto &G = this->m_VCK.retrieveOperatorMatrix("Gradient z");
        cusp::multiply(D, G, Laplacian);
      }
      case 2:
      {
        auto &D = this->m_VCK.retrieveOperatorMatrix("Divergence y");
        auto &G = this->m_VCK.retrieveOperatorMatrix("Gradient y");

        if (SpaceDim == 3)
        {
          Matrix dummy(Nrows, Ncolumns, nnz);
          cusp::multiply(D, G, dummy);
          cusp::add(Laplacian, dummy, Laplacian);
        }
        else
        {
          cusp::multiply(D, G, Laplacian);
        }
      }
      default:
      {
        auto &D = this->m_VCK.retrieveOperatorMatrix("Divergence x");
        auto &G = this->m_VCK.retrieveOperatorMatrix("Gradient x");
        Matrix dummy(Nrows, Ncolumns, nnz);
        cusp::multiply(D, G, dummy);
        cusp::add(Laplacian, dummy, Laplacian);
      }
      }
    }

    // generate the relaxator that goes with this operator, but only for square
    // matrices
    if (Nrows != Ncolumns)
    {
      m_Relaxators.insert(
          std::pair<std::string, Relaxator *>(OperatorName, nullptr));
    }
    else
    {
      // m_Relaxators.insert(std::pair<std::string, Relaxator *>(OperatorName,
      // new GSRelax(m_Operators[OperatorName]->getMatrixPtr(), m_Omega)));
      m_Relaxators.insert(std::pair<std::string, Relaxator *>(
          OperatorName,
          new SORRelax(m_Operators[OperatorName]->getMatrixPtr(), m_Omega)));
      // m_Relaxators.insert(std::pair<std::string, Relaxator *>(OperatorName,
      // new JacobiRelax(m_Operators[OperatorName]->getMatrixPtr(), .666)));
    }

    // collect the main diagonal
    // reduce the size of the box to the interior points
    if (OperatorName == "Laplacian" ||
        OperatorName == "LaplacianNG") // the diagonals are common to NG and full
    {
      m_Diagonals.insert(std::pair<std::string, DiagType *>(
          OperatorName, new DiagType(m_boxIntOnly, 1)));
      PyRealArray Data(Nrows);

      Py::PythonFunction(m_MyPyModule, "GetMatrixDiagonal", OperatorName,
                            Data);
      thrust::copy(Data.begin(), Data.end(),
                   m_Diagonals[OperatorName]->begin());
    }
  }
  // KEEP THIS LAST!
  Py::PythonFunction(m_MyPyModule, "ClearMatrix");
} // end of generate()
} // namespace PythonInterface
} // namespace Elliptic

#endif
