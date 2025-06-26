#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REALVECT_H_
#define _REALVECT_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <array>

#include "SPACE.H"
#include "Misc.H"
#include "REAL.H"
#include "SPACE.H"
#include "IntVect.H"
#include "Vector.H"
#include "SPMD.H"
#include "IndexTM.H"
#ifdef CH_USE_PYTHON
#include "Python.h"
#endif
#include "NamespaceHeader.H"


//template<typename T, int n>
//class IndexTM;

/// A Real vector in SpaceDim-dimensional space
/**
  The class RealVect is an implementation of a Real vector in a
  SpaceDim-dimensional space.
  RealVect values are accessed using the operator[] function, as for a normal
  C++ array.  In addition, the basic arithmetic operators have been overloaded
  to implement scaling and translation operations.
*/

class RealVect
{
public:

  /**
     \name Constructors and Accessors
  */
  /*@{*/

    // These make RealVects trivially copyable.
    constexpr RealVect(const RealVect& rhs)  noexcept = default;
    RealVect(RealVect&& rhs)                 noexcept = default;
    RealVect& operator=(const RealVect& rhs) noexcept = default;
    RealVect& operator=(RealVect&& rhs)      noexcept = default;
    ~RealVect()                              noexcept = default;


    ///
    /**
       Construct a RealVect whose components are zero.
    */
    RealVect() noexcept;

    explicit RealVect(const Vector<Real>&);

    ///
    /**
       Construct a RealVect given the specific values for its
       coordinates.  D_DECL6 is a macro that sets the constructor to
       take CH_SPACEDIM arguments.
    */
    constexpr RealVect(
        D_DECL6(Real i, Real j, Real k, Real l, Real m, Real n)) noexcept
    : vect{ D_DECL6(i, j, k, l, m, n) }
    {
  }

  ///
  /**
     Construct a RealVect from a C-style array.
  */
  explicit RealVect(const Real* r)
  : vect{ D_DECL6(r[0], r[1], r[2], r[3], r[4], r[5]) }
  {
  }

  ///
  /**
     Construct a RealVect from a std::array.
  */
  explicit RealVect(const std::array<Real, CH_SPACEDIM>& r) noexcept
  : vect{ D_DECL6(r[0], r[1], r[2], r[3], r[4], r[5]) }
  {
  }

#ifdef CH_USE_PYTHON
   explicit RealVect(PyObject *a_pin);
   static constexpr bool isConvertibleToNumpy = true;
   static std::string numpyName() {return std::to_string(SpaceDim)+"DRealVect";}

#endif

 ///
  /**
     Construct a RealVect from an IndexTM<Real,SpaceDim>.
  */
  RealVect(const IndexTM<Real,CH_SPACEDIM>& a_tm) noexcept;

  ///
  /**
     Construct a RealVect from an IntVect by coercing each component
     from <tt>int</tt> to Real.
  */
  RealVect(const IntVect& iv) noexcept
  {
      for (int d = 0; d < SpaceDim; ++d) {
          vect[d] = (Real)iv[d];
      }
  }

  ///
  /**
     Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
     RealVect.
  */
  inline
  Real& operator[] (int i);

  ///
  /**
     Returns the <i>i</i>'th coordinate of the RealVect.
  */
  inline
  constexpr Real operator[] (int i) const;

  inline constexpr Real
  getConst(const int i) const
  {
      CH_assert(i >= 0 && i < SpaceDim);
      return vect[i];
  }

  /*@}*/

  /**
     \name Comparison Operators
  */
  /*@{*/

  ///
  /**
     Returns true if this RealVect is equivalent to argument RealVect.  All
     comparisons between analogous components must be satisfied.
  */
  bool operator== (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is different from argument RealVect.
     All comparisons between analogous components must be satisfied.
  */
  bool operator!= (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is less than argument RealVect.  All
     comparisons between analogous components must be satisfied.  Note
     that, since the comparison is component-wise, it is possible for
     an RealVect to be neither greater than, less than, nor equal to
     another.
  */
  bool operator< (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is less than or equal to argument
     RealVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RealVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */
  bool operator<= (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is greater than argument RealVect.
     All comparisons between analogous components must be satisfied.
     Note that, since the comparison is component-wise, it is possible
     for an RealVect to be neither greater than, less than, nor equal
     to another.
  */
  bool operator> (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is greater than or equal to argument
     RealVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RealVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */

  bool operator>= (const RealVect& p) const noexcept;

  /*@}*/

  /**
     \name Arithmetic Operators
  */
  /*@{*/

  ///
  /**
     Modifies this RealVect by addition of a scalar to each component.
  */
  RealVect& operator+= (Real s) noexcept;

  ///
  /**
     Returns a RealVect that is this RealVect with a scalar s added
     to each component.
  */
  inline constexpr RealVect operator+ (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by component-wise addition by argument.
  */
  RealVect& operator+= (const RealVect& p) noexcept;

  ///
  /**
     Modifies this RealVect by subtraction of a scalar from each
     component.
  */
  RealVect& operator-= (Real s) noexcept;

  ///
  /**
     Modifies this RealVect by component-wise subtraction by argument.
  */
  RealVect& operator-= (const RealVect& p) noexcept;

  ///
  /**
     Returns a RealVect that is this RealVect with a scalar s subtracted
     from each component.
  */
  inline constexpr RealVect operator- (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by multiplying each component by a scalar.
  */
  RealVect& operator*= (Real s) noexcept;

  ///
  /**

  */
  Real dotProduct(const RealVect& a_rhs) const noexcept;

  ///
  /**
     Modifies this RealVect by component-wise multiplication by argument.
  */
  RealVect& operator*= (const RealVect& p) noexcept;

//XXX  ///
//XXX  /**
//XXX     Returns component-wise product of this RealVect and argument.
//XXX  */
//XXX  RealVect operator* (const RealVect& p) const;

  ///
  /**
     Returns a RealVect that is this RealVect with each component
     multiplied by a scalar.
  */
  inline constexpr RealVect operator* (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by dividing each component by a scalar.
  */
  RealVect& operator/= (Real s);

  ///
  /**
     Modifies this RealVect by component-wise division by argument.
  */
  RealVect& operator/= (const RealVect& p);

//XXX  ///
//XXX  /**
//XXX     Returns component-wise quotient of this RealVect by argument.
//XXX  */
//XXX  RealVect operator/ (const RealVect& p) const;

  ///
  /**
     Returns a RealVect that is this RealVect with each component
     divided by a scalar.
  */
  inline constexpr RealVect operator/ (Real s) const;

  ///
  /**
     Modifies this RealVect by multiplying each component by a scalar.
  */
  RealVect& scale (Real s) noexcept;

  /*@}*/

  /**
     \name Other arithmetic operators
  */
  /*@{*/

  ///
  /**
     Modifies this RealVect by taking component-wise min with RealVect
     argument.
  */
  RealVect& min (const RealVect& p);

  ///
  /**
     Returns the RealVect that is the component-wise minimum of two
     argument RealVects.
  */
  friend inline RealVect min (const RealVect& p1,
                              const RealVect& p2);

  ///
  /**
     Modifies this RealVect by taking component-wise max with RealVect
     argument.
  */
  RealVect& max (const RealVect& p);

  ///
  /**
     Returns the RealVect that is the component-wise maximum of two
     argument RealVects.
  */
  friend inline RealVect max (const RealVect& p1,
                              const RealVect& p2);

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  inline constexpr RealVect operator+ () const noexcept;

  ///
  /**
     Unary minus -- negates all components of this RealVect.
  */
  inline constexpr RealVect operator- () const noexcept;

  ///
  /**
     Sum of all components of this RealVect.
  */
  Real sum () const noexcept;

  ///
  /**
     sqrt(sum squares)
  */
  Real vectorLength() const noexcept;

  ///
  /**
     sum squares. No sqrt.
  */
  Real vectorLengthSq() const noexcept;

  ///
  /**
     Product of all components of this RealVect.
  */
  Real product () const noexcept;

  ///
  /**
     Component with the minimum value of this RealVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int minDir(const bool& a_doAbs) const;

  ///
  /**
     Component with the maximum value of this RealVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int maxDir(const bool& a_doAbs) const;

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Only for sending stuff to Fortran
   */
  const Real* dataPtr() const;

  ///
  /**
     Only for sending stuff to Fortran
   */
  Real* dataPtr() ;

  /*@}*/

  /**
     \name Constants
  */
  /*@{*/

  ///
  /**
     Returns a basis vector in the given coordinate direction.<br>
     In 2-D:<br>
     BASISREALV(0) == (1.,0.);
     BASISREALV(1) == (0.,1.).<br>
     In 3-D:<br>
     BASISREALV(0) == (1.,0.,0.);
     BASISREALV(1) == (0.,1.,0.);
     BASISREALV(2) == (0.,0.,1.).<br>
     Note that the coordinate directions are based at zero.
  */
  friend RealVect BASISREALV(int dir);

  ///
  /**
     This is a RealVect all of whose components are equal to zero.
  */
  static const RealVect Zero;

  ///
  /**
     This is a RealVect all of whose components are equal to one.
  */
  static const RealVect Unit;

  /*@}*/

  /**
     \name Arithmetic friend functions
  */
  /*@{*/

  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with
     a scalar <i>s</i> added to each component.
  */
  friend inline constexpr RealVect
  operator+(Real s, const RealVect& p) noexcept;

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend inline constexpr RealVect
  operator-(Real s, const RealVect& p) noexcept;

  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with each component
     multiplied by a scalar <i>s</i>.
  */
  friend inline constexpr RealVect
  operator*(Real s, const RealVect& p) noexcept;
  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with each component
     divided by a scalar <i>s</i>.
  */
  friend inline constexpr RealVect
  operator/(Real s, const RealVect& p);

  ///
  /**
     Returns component-wise sum of RealVects <i>s</i> and <i>p</i>.
  */
  friend inline constexpr RealVect
  operator+(const RealVect& s, const RealVect& p) noexcept;

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend inline constexpr RealVect
  operator-(const RealVect& s, const RealVect& p) noexcept;

  ///
  /**
     Returns component-wise product of <i>s</i> and <i>p</i>.
  */
  friend inline constexpr RealVect
  operator*(const RealVect& s, const RealVect& p) noexcept;
  ///
  /**
     Returns component-wise quotient <i>p / s</i>.
  */
  friend inline constexpr RealVect
  operator/(const RealVect& s, const RealVect& p);


  ///
  /**
     Returns a RealVect obtained by multiplying each of the components
     of the given RealVect by a scalar.
  */
  friend inline RealVect scale (const RealVect& p,
                                Real            s) noexcept;


  /*@}*/

  ///
  /**
     Print to the given output stream in ASCII.
  */
  friend std::ostream& operator<< (std::ostream& ostr,
                                   const RealVect& p);

  friend class HDF5Handle;

  static size_t io_offset;
#ifdef CH_USE_PYTHON
   PyObject *pack() const;
#endif
protected:

  /**
     The individual components of this RealVect.
   */
  Real vect[SpaceDim];
};

#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"

#include "NamespaceVar.H"
// template specialization of linearIn and LinearOut for RealVect
// RealVect spcializations of linearization
template < >
int linearSize(const CH_XDIR::RealVect& vindex);

//VolIndex specialization of linearIn
template < >
void linearIn(CH_XDIR::RealVect& a_outputT, const void* const inBuf);

//VolIndex specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const CH_XDIR::RealVect& a_inputT);

#include "BaseNamespaceFooter.H"

#include "NamespaceHeader.H"

inline Real& RealVect::operator[] (int i)
{
  CH_assert(i>=0 && i < SpaceDim);
  return vect[i];
}

inline constexpr Real RealVect::operator[] (int i) const
{
  CH_assert(i>=0 && i < SpaceDim);
  return vect[i];
}

inline
RealVect&
RealVect::operator-= (Real s) noexcept
{
  D_EXPR6(vect[0] -= s, vect[1] -= s, vect[2] -= s,
          vect[3] -= s, vect[4] -= s, vect[5] -= s);
  return *this;
}

inline
RealVect&
RealVect::operator-= (const RealVect& p) noexcept
{
  D_EXPR6(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2],
          vect[3] -= p[3], vect[4] -= p[4], vect[5] -= p[5]);
  return *this;
}

inline
constexpr RealVect
RealVect::operator+ () const noexcept
{
  return RealVect(*this);
}

inline
constexpr RealVect
RealVect::operator- () const noexcept
{
  return RealVect(D_DECL6(-vect[0], -vect[1], -vect[2],
                          -vect[3], -vect[4], -vect[5] ));
}

inline
RealVect&
RealVect::scale (Real s) noexcept
{
  D_EXPR6(vect[0] *= s, vect[1] *= s, vect[2] *= s,
          vect[3] *= s, vect[4] *= s, vect[5] *= s);
  return *this;
}

inline
Real
RealVect::sum () const noexcept
{
  return D_TERM6(vect[0], + vect[1], + vect[2], +
                 vect[3], + vect[4], + vect[5]);
}

inline
Real
RealVect::vectorLength () const noexcept
{
  Real len = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      len = len + vect[idir]*vect[idir];
    }
  len = sqrt(len);

  return len;
}

inline
Real
RealVect::vectorLengthSq () const noexcept
{
    return D_TERM(vect[0] * vect[0], + vect[1] * vect[1], + vect[2] * vect[2]);
}

inline
Real
RealVect::product () const noexcept
{
  return D_TERM6(vect[0], * vect[1], * vect[2], *
                 vect[3], * vect[4], * vect[5]);
}

inline
RealVect
scale (const RealVect& p,
       Real            s) noexcept
{
  return RealVect(D_DECL6(s * p[0], s * p[1], s * p[2],
                          s * p[3], s * p[4], s * p[5]));
}

inline
bool
RealVect::operator< (const RealVect& p) const noexcept
{
  return D_TERM6(vect[0] < p[0], && vect[1] < p[1], && vect[2] < p[2],
              && vect[3] < p[3], && vect[4] < p[4], && vect[5] < p[5]);
}

inline
bool
RealVect::operator<= (const RealVect& p) const noexcept
{
  return D_TERM6(vect[0] <= p[0], && vect[1] <= p[1], && vect[2] <= p[2],
              && vect[3] <= p[3], && vect[4] <= p[4], && vect[5] <= p[5]);
}


inline
bool
RealVect::operator> (const RealVect& p) const noexcept
{
  return D_TERM6(vect[0] > p[0], && vect[1] > p[1], && vect[2] > p[2],
              && vect[3] > p[3], && vect[4] > p[4], && vect[5] > p[5]);
}

inline
bool
RealVect::operator>= (const RealVect& p) const noexcept
{
  return D_TERM6(vect[0] >= p[0], && vect[1] >= p[1], && vect[2] >= p[2],
              && vect[3] >= p[3], && vect[4] >= p[4], && vect[5] >= p[5]);
}

inline
RealVect&
RealVect::min (const RealVect& p)
{
  D_EXPR6(vect[0] = Min(vect[0], p.vect[0]),
          vect[1] = Min(vect[1], p.vect[1]),
          vect[2] = Min(vect[2], p.vect[2]),
          vect[3] = Min(vect[3], p.vect[3]),
          vect[4] = Min(vect[4], p.vect[4]),
          vect[5] = Min(vect[5], p.vect[5]));
  return *this;
}

inline
RealVect&
RealVect::max (const RealVect& p)
{
  D_EXPR6(vect[0] = Max(vect[0], p.vect[0]),
          vect[1] = Max(vect[1], p.vect[1]),
          vect[2] = Max(vect[2], p.vect[2]),
          vect[3] = Max(vect[3], p.vect[3]),
          vect[4] = Max(vect[4], p.vect[4]),
          vect[5] = Max(vect[5], p.vect[5]));
  return *this;
}

inline
RealVect
min (const RealVect& p1,
     const RealVect& p2)
{
  RealVect p(p1);
  return p.min(p2);
}

inline
RealVect
max (const RealVect& p1,
     const RealVect& p2)
{
  RealVect p(p1);
  return p.max(p2);
}

extern RealVect BASISREALV(int idir);


inline constexpr RealVect
operator/ (Real            s,
           const RealVect& p)
{
  return RealVect(D_DECL6(s/p[0], s/p[1], s/p[2],
                          s/p[3], s/p[4], s/p[5]));
}
inline constexpr RealVect
operator+ (Real            s,
           const RealVect& p) noexcept
{
  return RealVect(D_DECL6(p[0] + s, p[1] + s, p[2] + s,
                          p[3] + s, p[4] + s, p[5] + s));
}

inline constexpr RealVect
operator- (Real            s,
           const RealVect& p) noexcept
{
  return RealVect(D_DECL6(s - p[0], s - p[1], s - p[2],
                          s - p[3], s - p[4], s - p[5]));
}

inline constexpr RealVect
operator* (Real            s,
           const RealVect& p) noexcept
{
  return RealVect(D_DECL6(s * p[0], s * p[1], s * p[2],
                          s * p[3], s * p[4], s * p[5]));
}

inline constexpr RealVect
operator/ (const RealVect& s,
           const RealVect& p)
{
  return RealVect(D_DECL6(s[0] / p[0], s[1] /p[1], s[2] / p[2],
                          s[3] / p[3], s[4] /p[4], s[5] / p[5]));
}

inline constexpr RealVect
operator+ (const RealVect& s,
           const RealVect& p) noexcept
{
  return RealVect(D_DECL6(p[0] + s[0], p[1] +s[1], p[2] + s[2],
                          p[3] + s[3], p[4] +s[4], p[5] + s[5]));
}

inline constexpr RealVect
operator- (const RealVect& s,
           const RealVect& p) noexcept
{
  return RealVect(D_DECL6(s[0] - p[0], s[1] - p[1], s[2] - p[2],
                          s[3] - p[3], s[4] - p[4], s[5] - p[5]));
}

inline constexpr RealVect
operator* (const RealVect& s,
           const RealVect& p) noexcept
{
  return RealVect(D_DECL6(p[0] * s[0], p[1] *s[1], p[2] * s[2],
                          p[3] * s[3], p[4] *s[4], p[5] * s[5]));
}

inline constexpr RealVect
RealVect::operator* (Real s) const noexcept
{
  RealVect v(D_DECL6(vect[0]*s, vect[1]*s, vect[2]*s,
                     vect[3]*s, vect[4]*s, vect[5]*s));
  return v;
}

inline constexpr RealVect
RealVect::operator- (Real s) const noexcept
{
  RealVect v(D_DECL6(vect[0]-s, vect[1]-s, vect[2]-s,
                     vect[3]-s, vect[4]-s, vect[5]-s));
  return v;
}


inline constexpr RealVect
RealVect::operator+ (Real s) const noexcept
{
  RealVect v(D_DECL6(vect[0]+s, vect[1]+s, vect[2]+s,
                     vect[3]+s, vect[4]+s, vect[5]+s));
  return v;
}

inline constexpr RealVect
RealVect::operator/ (Real s) const
{
  RealVect result( D_DECL6( vect[0] / s, vect[1] / s, vect[2] / s,
                            vect[3] / s, vect[4] / s, vect[5] / s) );
  return result ;
}

#include "NamespaceFooter.H"
#endif
