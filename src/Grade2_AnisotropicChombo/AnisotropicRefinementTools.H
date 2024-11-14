/**
 * \file   AnisotropicRefinementTools.H
 * \brief  Tools to deal with anisotropic refinement and coarsening.
 * \date   Mon 03 Jun 2013 12:59:01 AM EDT
 * \author Edward.Santilli@jefferson.edu
 */
#ifndef __AnisotropicRefinementTools_H__INCLUDED__
#define __AnisotropicRefinementTools_H__INCLUDED__

class Box;
class BoxLayout;
class DisjointBoxLayout;
class IntVect;
class IntVectSet;
class ProblemDomain;
class Copier;
class CFRegion;
class CFIVS;
class RealVect;


// -----------------------------------------------------------------------------
/// \name Box tools
// -----------------------------------------------------------------------------
/// \{

/// Returns true if a_box can be coarsened by a_refRatio and
/// return back to the original Box when refined by a_refRatio.
bool coarsenable (const Box&     a_box,
                  const IntVect& a_refRatio);

/// A very simple utility to calculate how much refinement is needed to send
/// data from a_src to a_dest. Returns a negative value if coarsening is
/// needed and zero if the two are incompatible.
IntVect calculateRefinementRatio (const Box& a_src,
                                  const Box& a_dest);
/// \}


// -----------------------------------------------------------------------------
/// \name ProblemDomain tools
// -----------------------------------------------------------------------------
/// \{

/// Use this instead of Chombo's coarsen.
void coarsen (ProblemDomain&       a_crseDomain,
              const ProblemDomain& a_fineDomain,
              const IntVect&       a_refRatio);

/// Use this instead of Chombo's refine.
void refine (ProblemDomain&       a_fineDomain,
             const ProblemDomain& a_crseDomain,
             const IntVect&       a_refRatio);

/// Computes how much we need to coarsen to get *close* to isotropic cells.
void
coarsenToIsotropy(ProblemDomain&       a_isoDomain,
                  IntVect&             a_refToIsotropy,
                  const ProblemDomain& a_domain,
                  const RealVect&      a_L);
/// \}


// -----------------------------------------------------------------------------
/// \name BoxLayout tools
// -----------------------------------------------------------------------------
/// \{

/// Returns true if every Box in the BoxLayout can be coarsened by a_refRatio and
/// return back to the original Box when refined by a_refRatio.
bool coarsenable (const BoxLayout& a_layout,
                  const IntVect&   a_refRatio);

///
void coarsen (BoxLayout&       a_output,
              const BoxLayout& a_input,
              const IntVect&   a_refRatio);

///
void refine (BoxLayout&       a_output,
             const BoxLayout& a_input,
             const IntVect&   a_refinement);

///
void coarsen (DisjointBoxLayout&       a_output,
              const DisjointBoxLayout& a_input,
              const IntVect&           a_refRatio);

///
void refine (DisjointBoxLayout&       a_output,
             const DisjointBoxLayout& a_input,
             const IntVect&           a_refinement);
/// \}


// -----------------------------------------------------------------------------
/// \name Copier tools
// -----------------------------------------------------------------------------
/// \{

///
void coarsen (Copier&        a_copier,
              const IntVect& a_refRatio);
/// \}


// -----------------------------------------------------------------------------
/// \name IntVectSet tools
// -----------------------------------------------------------------------------
/// \{

/// \bug
///  The IVS design makes it difficult to coarsen quickly in parallel. I've
///  avoided the issue by making thus a serial operation. We should revisit
///  the issue someday!
void coarsen (IntVectSet&    a_ivs,
              const IntVect& a_refRatio);

///
void refine (IntVectSet&    a_ivs,
             const IntVect& a_refRatio);
/// \}


// -----------------------------------------------------------------------------
/// \name CFRegion tools
// -----------------------------------------------------------------------------
/// \{

///
void coarsen (CFRegion&      a_cfRegion,
              const IntVect& a_refRatio);
/// \}


// -----------------------------------------------------------------------------
/// \name CFIVS tools
// -----------------------------------------------------------------------------
/// \{

///
void coarsen (CFIVS&         a_cfivs,
              const IntVect& a_refRatio);
/// \}


#endif //!__AnisotropicRefinementTools_H__INCLUDED__


