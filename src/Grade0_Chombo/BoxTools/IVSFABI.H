#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IVSFABI_H_
#define _IVSFABI_H_

#include "MayDay.H"
#include "IntVectSet.H"
#include "parstream.H"
#include "SPMD.H"
#include "DebugOut.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

template <class T>
bool IVSFAB<T>::s_verbose = false;

template <class T>
void
IVSFAB<T>::setVerbose(bool a_verbose)
{
  s_verbose = a_verbose;
}
/*************/
template <class T> inline
const IntVectSet&
IVSFAB<T>::getIVS() const
{
  return m_ivs;
}
/******************/
template <class T> inline
IVSFAB<T>::IVSFAB()
{
  setDefaultValues();
}
/******************/
template <class T> inline
IVSFAB<T>::~IVSFAB()
{
  clear();
}
/******************/
template <class T> inline
IVSFAB<T>::IVSFAB(const IntVectSet& a_ivsin,
                        const int&        a_nvarin)
{
  setDefaultValues();
  define(a_ivsin, a_nvarin);
}
/******************/
template <class T> inline
void
IVSFAB<T>::define(const IntVectSet& a_ivsin,
                     const int& a_nvarin)
{
  clear();
  m_isDefined = true;
  CH_assert(a_nvarin > 0);
  m_ivs = a_ivsin;
  m_nComp = a_nvarin;
  m_nIvs = 0;
  if (!a_ivsin.isEmpty())
    {
      Box minbox = a_ivsin.minBox();
      m_ivmap.resize(minbox, 1);
      IVSIterator ivsit(a_ivsin);
      for (ivsit.reset(); ivsit.ok(); ++ivsit)
        {
          const IntVect& iv = ivsit();
          m_ivmap(iv,0) = m_nIvs;
          m_nIvs++;
        }

      if ( (m_nIvs > 0) && m_nComp > 0)
        {
          m_dataPtr = new T[m_nComp*m_nIvs];
        }
      else
        {
          m_dataPtr = NULL;
        }
    }
  else
    {
      m_dataPtr = NULL;
    }
  if (m_nIvs > 0)
    {
      //set up face low and high vectors
      for (int dir = 0; dir < CH_SPACEDIM; ++dir)
        {
          m_loVect[dir] = 1;
          m_hiVect[dir] = 1;
        }
      m_hiVect[0] = m_nIvs;

    }
  else
    {
      m_loVect = IntVect::Unit;
      m_hiVect = IntVect::Zero;
      m_dataPtr = NULL;
    }

}
/******************/
template <class T> inline
void
IVSFAB<T>::setVal(const T& a_value)
{
  CH_assert(isDefined());
  for (int ivec = 0; ivec < m_nIvs*m_nComp; ivec++)
    m_dataPtr[ivec] = a_value;
}
/******************/
template <class T>
inline void
IVSFAB<T>::copy(const Box&                  a_fromBox,
                const Interval&             a_dstInterval,
                [[maybe_unused]] const Box& a_toBox,
                const IVSFAB<T>&            a_src,
                const Interval&             a_srcInterval)
{
    CH_assert(a_fromBox == a_toBox);
    CH_assert(isDefined());
    CH_assert(a_src.isDefined());
    CH_assert(a_srcInterval.size() == a_dstInterval.size());
    CH_assert(a_dstInterval.begin() >= 0);
    CH_assert(a_srcInterval.begin() >= 0);
    CH_assert(a_dstInterval.end() < m_nComp);
    CH_assert(a_srcInterval.end() < a_src.m_nComp);

    if ((!m_ivs.isEmpty()) && (!a_src.m_ivs.isEmpty())) {
        Box        intBox       = a_fromBox;
        IntVectSet ivsIntersect = m_ivs;
        ivsIntersect &= a_src.m_ivs;
        ivsIntersect &= intBox;
        IVSFAB<T>& thisFAB  = *this;
        int        compSize = a_srcInterval.size();
        for (IVSIterator ivit(ivsIntersect); ivit.ok(); ++ivit) {
            const IntVect& iv = ivit();
            for (int icomp = 0; icomp < compSize; icomp++) {
                int isrccomp          = a_srcInterval.begin() + icomp;
                int idstcomp          = a_dstInterval.begin() + icomp;
                thisFAB(iv, idstcomp) = a_src(iv, isrccomp);
            }
        }
    }
}
/******************/

template <> inline
int IVSFAB<Real>::size(const Box& a_region,
                       const Interval& a_comps) const
{
  // size of actual data plus size of count ..
  // specialization for T=Real.
  IntVectSet subset = m_ivs & a_region;
  int numberThings = subset.numPts();
  return numberThings*sizeof(Real)*a_comps.size() + sizeof(int);
}

template <class T>
inline int
IVSFAB<T>::size(const Box&                       a_region,
                [[maybe_unused]] const Interval& a_comps) const
{
    // Note:  it is assumed that there is one component.
    //        and, that all "T" have equal size.
    // message includes count ...
    CH_assert(a_comps.begin() == 0 && a_comps.end() == 0);
    IntVectSet subset       = m_ivs & a_region;
    int        numberThings = subset.numPts();
    T          dummy;
    int        sizeofThing = dummy.linearSize();
    int        retval      = numberThings * sizeofThing + sizeof(int);
    return retval;
}


/********************/
template <> inline void IVSFAB<Real>::linearOut(void* a_buf,
                             const Box& a_region,
                             const Interval& a_comps) const
{
  IntVectSet subset = m_ivs & a_region;
  int count = subset.numPts();
  int* iptr = (int*)a_buf;
  *iptr++ = count;
  Real* rptr = (Real*)iptr;
  for ( IVSIterator ivsit(subset); ivsit.ok(); ++ivsit )
  {
    const IntVect& iv = ivsit();
    for ( int c=a_comps.begin(); c<=a_comps.end(); c++ )
    {
      int iloc = getIndex(iv, c);
      *rptr++ = m_dataPtr[iloc];
    }
  }
  CH_assert( ((char*)rptr) - ((char*)a_buf) == size(a_region,a_comps) );
}


template <class T>
inline void
IVSFAB<T>::linearOut(void*                            a_buf,
                     const Box&                       a_region,
                     [[maybe_unused]] const Interval& a_comps) const
{
    // Note:  it is assumed that there is one component.
    //        and all "T" have equal size.
    // message includes count
    IntVectSet subset = m_ivs & a_region;
    int        count  = subset.numPts();
    int*       iptr   = (int*)a_buf;
    *iptr++           = count;
    T     dummy;
    int   Tsz  = dummy.linearSize();
    char* cptr = (char*)iptr;
    for (IVSIterator ivsit(subset); ivsit.ok(); ++ivsit) {
        const IntVect& iv   = ivsit();
        int            iloc = getIndex(iv, 0);  // assume 1 component ...
        m_dataPtr[iloc].linearOut((void*)cptr);
        cptr += Tsz;
    }
    CH_assert(a_comps.begin() == 0 && a_comps.end() == 0);
    CH_assert(cptr - ((char*)a_buf) == size(a_region, a_comps));
}


/********************/
template <> inline void IVSFAB<Real>::linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  // message includes count
  IntVectSet subset = m_ivs & a_region;
  int count_expected = subset.numPts();
  int *iptr = (int*)a_buf;
  int count_received = *iptr++;
  if ( count_expected != count_received )
    {
      MayDay::Error("IVSFAB::linearin: counts do not match");
    }
  Real* rptr = (Real*)iptr;
  for ( IVSIterator ivsit(subset); ivsit.ok(); ++ivsit )
  {
    const IntVect& iv = ivsit();
    for ( int c=a_comps.begin(); c<=a_comps.end(); c++ )
    {
      int iloc = getIndex(iv, c);
      m_dataPtr[iloc] = *rptr++;
    }
  }
  CH_assert( ((char*)rptr) - ((char*)a_buf) == size(a_region,a_comps) );
}

template <class T>
inline void
IVSFAB<T>::linearIn(void*                            a_buf,
                    const Box&                       a_region,
                    [[maybe_unused]] const Interval& a_comps)
{
    // Note:  it is assumed that there is one component.
    //        and all "T" have same size
    // message includes count
    IntVectSet           subset         = m_ivs & a_region;
    [[maybe_unused]] int count_expected = subset.numPts();
    int*                 iptr           = (int*)a_buf;
    [[maybe_unused]] int count_received = *iptr++;
    CH_assert(count_expected == count_received);
    T     dummy;
    int   Tsz  = dummy.linearSize();
    char* cptr = (char*)iptr;
    for (IVSIterator ivsit(subset); ivsit.ok(); ++ivsit) {
        const IntVect& iv   = ivsit();
        int            iloc = getIndex(iv, 0);
        m_dataPtr[iloc].linearIn((void*)cptr);
        cptr += Tsz;
    }
    CH_assert(a_comps.begin() == 0 && a_comps.end() == 0);
    CH_assert(cptr - ((char*)a_buf) == size(a_region, a_comps));
}

template <class T> inline
int
IVSFAB<T>::getIndex(const IntVect& a_iv, const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(m_ivs.contains(a_iv));
  CH_assert((a_comp >= 0) && (a_comp < m_nComp));

  int ioffset = m_ivmap(a_iv, 0);
  CH_assert(ioffset >= 0);
  CH_assert(ioffset < m_nIvs);
  //now add offset from componentnitude
  ioffset += m_nIvs*a_comp;
  return ioffset;
}
/********************/
template <class T> inline
void
IVSFAB<T>::clear()
{
  m_nComp = 0;
  m_nIvs = 0;
  m_ivs.makeEmpty();
  m_ivmap.clear();
  if (m_dataPtr != NULL)
    {
      delete[] m_dataPtr;
      m_dataPtr = NULL;
    }
  m_isDefined = false;
}
/*************************/
template <class T> inline
bool
IVSFAB<T>::isDefined() const
{
  return (m_isDefined);
}
/*************************/
template <class T> inline
int
IVSFAB<T>::numIvs() const
{
  return m_nIvs;
}
/*************************/
template <class T> inline
int
IVSFAB<T>::nComp() const
{
  return m_nComp;
}
/*************************/
template <class T> inline
T&
IVSFAB<T>::operator() (const IntVect& a_ndin,
                       const int& a_comp)
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  int iloc = getIndex(a_ndin, a_comp);
  return(m_dataPtr[iloc]);
}
/**************************/
template <class T> inline
const T&
IVSFAB<T>::operator() (const IntVect& a_ndin,
                       const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  int iloc = getIndex(a_ndin, a_comp);
  return(m_dataPtr[iloc]);
}
/******************/
template <class T> inline
const T*
IVSFAB<T>::dataPtr(const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  return m_dataPtr + a_comp*m_nIvs;
}
/******************/
template <class T> inline
T*
IVSFAB<T>::dataPtr(const int& a_comp)
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  return m_dataPtr + a_comp*m_nIvs;
}
/******************/
template <class T> inline
const int*
IVSFAB<T>::loVect() const
{
  return m_loVect.getVect();
}
/******************/
template <class T> inline
const int*
IVSFAB<T>::hiVect() const
{
  return m_hiVect.getVect();
}
/******************/
template <class T> inline
void
IVSFAB<T>::setDefaultValues()
{
  m_isDefined = false;
  m_dataPtr = NULL;
  m_nIvs = 0;
  m_nComp = 0;
  m_loVect = IntVect::Unit;
  m_hiVect = IntVect::Zero;
}
/******************/
template <class T> inline
IVSFAB<T>&
IVSFAB<T>::operator= (const IVSFAB<T>& /*a_input*/)
{
  MayDay::Error("IVSFAB operator = not defined");
  return *this;
}
/******************/
template <class T> inline
IVSFAB<T>::IVSFAB (const IVSFAB<T>& /*a_input*/)
{
  MayDay::Error("IVSFAB copy constructor not defined");
}

#include "NamespaceFooter.H"
#endif
