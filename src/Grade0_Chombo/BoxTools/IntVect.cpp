#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SPACE.H"
#include <iostream>
using std::ostream;
using std::istream;
using std::ws;

#include "MayDay.H"
#include "Misc.H"
#include "IntVect.H"
#include "parstream.H"
#include "IndexTM.H"

#include "NamespaceHeader.H"

//const size_t  IntVect::IntVectSize = SpaceDim*sizeof(int);

//
// Returns IntVect which is the componentwise integer projection
// of IntVect p1 by IntVect p2.
//

ostream&
operator<< (ostream&       os,
            const IntVect& p)
{
  os << D_TERM6( '(' << p[0] , <<
                 ',' << p[1] , <<
                 ',' << p[2] , <<
                 ',' << p[3] , <<
                 ',' << p[4] , <<
                 ',' << p[5])  << ')';
    if (os.fail())
        MayDay::Error("operator<<(ostream&,IntVect&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define CH_IGNORE_MAX 100000

istream&
operator>> (istream& is,
            IntVect& p)
{
    is >> ws;
    char c;
    is >> c;
    is.putback(c);
    if (c == '(')
    {
        D_EXPR6(is.ignore(CH_IGNORE_MAX, '(') >> p[0],
                is.ignore(CH_IGNORE_MAX, ',') >> p[1],
                is.ignore(CH_IGNORE_MAX, ',') >> p[2],
                is.ignore(CH_IGNORE_MAX, '(') >> p[3],
                is.ignore(CH_IGNORE_MAX, ',') >> p[4],
                is.ignore(CH_IGNORE_MAX, ',') >> p[5]);

        is.ignore(CH_IGNORE_MAX, ')');
    }
    else if (c == '<')
    {
        D_EXPR6(is.ignore(CH_IGNORE_MAX, '<') >> p[0],
                is.ignore(CH_IGNORE_MAX, ',') >> p[1],
                is.ignore(CH_IGNORE_MAX, ',') >> p[2],
                is.ignore(CH_IGNORE_MAX, '<') >> p[3],
                is.ignore(CH_IGNORE_MAX, ',') >> p[4],
                is.ignore(CH_IGNORE_MAX, ',') >> p[5]);
        is.ignore(CH_IGNORE_MAX, '>');
    }
    else
        MayDay::Error("operator>>(istream&,IntVect&): expected \'(\' or \'<\'");

    if (is.fail())
        MayDay::Error("operator>>(istream&,IntVect&) failed");

    return is;
}

void
IntVect::printOn (ostream& os) const
{
    os << "IntVect: " << *this << '\n';
}

void
IntVect::p() const
{
    pout() << *this << '\n';
}

IntVect::IntVect(const IndexTM<int, CH_SPACEDIM>& a_tm) noexcept
{
  D_EXPR6(vect[0] = a_tm[0], vect[1] = a_tm[1], vect[2] = a_tm[2],
          vect[3] = a_tm[3], vect[4] = a_tm[4], vect[5] = a_tm[5]);
}

void
IntVect::dumpOn (ostream& os) const
{
    os << "IntVect " << *this << '\n';
}

//#ifdef __INTEL_COMPILER
const IntVect IntVect::Zero(D_DECL6(0,0,0,0,0,0));
const IntVect IntVect::Unit(D_DECL6(1,1,1,1,1,1));
//#else
//constexpr IntVect IntVect::Zero(D_DECL6(0,0,0,0,0,0));
//constexpr IntVect IntVect::Unit(D_DECL6(1,1,1,1,1,1));
//#endif
 
 

#ifdef CH_USE_PYTHON
IntVect::IntVect(PyObject *a_pin)
{
    
  
  for (auto i = 0; i < SpaceDim; ++i) {
    vect[i] = static_cast<int>(PyLong_AsLong(PyTuple_GetItem(a_pin, i)));
  }
}

PyObject *IntVect::pack() const {
  PyObject *pIntVect = PyTuple_New(SpaceDim+1);
  for (auto i = 0; i < SpaceDim; ++i) {
    PyTuple_SetItem(pIntVect, i, PyLong_FromLong(vect[i]));
  }
  std::string Label = "IntVect";
  PyObject *pLabel = PyUnicode_FromString(Label.c_str());
  PyTuple_SetItem(pIntVect, SpaceDim, pLabel);
  return pIntVect;
}
#endif
#include "NamespaceFooter.H"
