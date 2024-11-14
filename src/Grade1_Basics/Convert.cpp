#include "Convert.H"
#include "ConvertF_F.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// The Swiss Army knife. Performs basic, 2nd-order averages.
// Used to convert a CC FAB to FC FAB and so on.
// -----------------------------------------------------------------------------
void
Convert::Simple (FArrayBox&       a_dest,
                 const Interval&  a_destIvl,
                 const Box&       a_destBox,
                 const FArrayBox& a_src,
                 const Interval&  a_srcIvl)
{
    // Gather some needed info.
    const int ncomp = a_destIvl.size();
    const Box& srcBox = a_src.box();
    const IntVect& srcBoxType = srcBox.type();
    const IntVect& destBoxType = a_destBox.type();

#ifndef NDEBUG
    // Sanity checks
    CH_assert(a_destIvl.size() == a_srcIvl.size());

    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (srcBoxType[dir] == 0 && destBoxType[dir] == 1) {
            CH_assert(srcBox.smallEnd(dir) <  a_destBox.smallEnd(dir));
            CH_assert(srcBox.  bigEnd(dir) >= a_destBox.  bigEnd(dir));
        } else if (srcBoxType[dir] == 1 && destBoxType[dir] == 0) {
            CH_assert(srcBox.smallEnd(dir) <= a_destBox.smallEnd(dir));
            CH_assert(srcBox.  bigEnd(dir) >  a_destBox.  bigEnd(dir));
        } else {
            CH_assert(srcBox.smallEnd(dir) <= a_destBox.smallEnd(dir));
            CH_assert(srcBox.  bigEnd(dir) >= a_destBox.  bigEnd(dir));
        }
    }
#endif // !NDEBUG

    // Loop over FAB components and recenter each.
    for (int comp = 0; comp < ncomp; ++comp) {
        const int destComp = a_destIvl.begin() + comp;
        const int  srcComp =  a_srcIvl.begin() + comp;

        CH_assert(a_dest.interval().contains(destComp));
        CH_assert( a_src.interval().contains( srcComp));

        FORT_CONVERT_SIMPLE (
            CHF_FRA1(a_dest, destComp),
            CHF_BOX(a_destBox),
            CHF_CONST_INTVECT(destBoxType),
            CHF_CONST_FRA1(a_src, srcComp),
            CHF_CONST_INTVECT(srcBoxType));
    }
}


// -----------------------------------------------------------------------------
// The Swiss Army knife. Performs basic, 2nd-order averages.
// Used to convert a CC FAB to FC FAB and so on.
// -----------------------------------------------------------------------------
void
Convert::Simple (FArrayBox&       a_dest,
                 const int        a_destComp,
                 const Box&       a_destBox,
                 const FArrayBox& a_src,
                 const int        a_srcComp)
{
    Interval destIvl(a_destComp, a_destComp);
    Interval  srcIvl( a_srcComp,  a_srcComp);
    Simple(a_dest, destIvl, a_destBox, a_src, srcIvl);
}


// -----------------------------------------------------------------------------
// The Swiss Army knife. Performs basic, 2nd-order averages.
// Used to convert a CC FAB to FC FAB and so on.
// -----------------------------------------------------------------------------
void
Convert::Simple (FArrayBox&       a_dest,
                 const Box&       a_destBox,
                 const FArrayBox& a_src)
{
    Simple(a_dest, a_dest.interval(), a_destBox, a_src, a_src.interval());
}


// -----------------------------------------------------------------------------
// The Swiss Army knife. Performs basic, 2nd-order averages.
// Used to convert a CC FAB to FC FAB and so on.
// -----------------------------------------------------------------------------
void
Convert::Simple (FArrayBox&       a_dest,
                 const FArrayBox& a_src)
{
    Simple(a_dest, a_dest.interval(), a_dest.box(), a_src, a_src.interval());
}


// -----------------------------------------------------------------------------
// This version requires a_src.nComp() == SpaceDim * a_dest.nComp().
// Comp c, FC in the d direction of a_dest is interpolated from comp
// d + SpaceDim*c of a_dest.
// -----------------------------------------------------------------------------
void
Convert::CellsToFaces(FluxBox&         a_destFlub,
                      const FArrayBox& a_srcFAB,
                      const Box&       a_ccDestBox)
{
    CH_assert(a_srcFAB.nComp() == SpaceDim * a_destFlub.nComp());

    const int numDestComps = a_destFlub.nComp();

    for (int destComp = 0; destComp < numDestComps; ++destComp) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            const int srcComp = SpaceDim * destComp + fcDir;
            Convert::Simple(
                a_destFlub[fcDir], destComp, a_ccDestBox, a_srcFAB, srcComp);
        }  // end loop over fcDir
    } // end loop over destComp
}


// -----------------------------------------------------------------------------
// This version requires a_src.nComp() == SpaceDim * a_dest.nComp().
// Comp c, FC in the d direction of a_dest is interpolated from comp
// d + SpaceDim*c of a_dest.
// -----------------------------------------------------------------------------
void
Convert::CellsToFaces (LevelData<FluxBox>&         a_dest,
                       const LevelData<FArrayBox>& a_src,
                       const bool                  a_validOnly)
{
    CH_assert(a_src.nComp() == SpaceDim * a_dest.nComp());
    CH_assert(a_dest.getBoxes().compatible(a_src.getBoxes()));

    const DisjointBoxLayout& grids        = a_dest.getBoxes();
    DataIterator             dit          = a_dest.dataIterator();
    const int                numDestComps = a_dest.nComp();

    for (dit.reset(); dit.ok(); ++dit) {
        for (int destComp = 0; destComp < numDestComps; ++destComp) {
            for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                FArrayBox&       destFAB = a_dest[dit][fcDir];
                const FArrayBox& srcFAB  = a_src[dit];
                const int        srcComp = SpaceDim * destComp + fcDir;

                Box destBox = destFAB.box();
                if (a_validOnly) destBox &= surroundingNodes(grids[dit], fcDir);

                Convert::Simple(destFAB, destComp, destBox, srcFAB, srcComp);
            } // end loop over fcDir
        } // end loop over destComp
    } // end loop over grids (dit)
}


// -----------------------------------------------------------------------------
// This version takes each comp of CC a_src and sends it to a comp of FC
// a_dest. For example, if a_src has 2 comps, then a_dest must also have
// 2 comps. Comp c of a_src will be sent to comp c of a_dest and will be
// interpolated to each face of the original cell.
// -----------------------------------------------------------------------------
void
Convert::CellsToAllFaces (FluxBox&         a_destFlub,
                          const FArrayBox& a_srcFAB)
{
    CH_assert(a_destFlub.nComp() == a_srcFAB.nComp());

    const Interval& ivl = a_srcFAB.interval();

    D_TERM(
    Convert::Simple(a_destFlub[0], ivl, a_destFlub[0].box(), a_srcFAB, ivl);,
    Convert::Simple(a_destFlub[1], ivl, a_destFlub[1].box(), a_srcFAB, ivl);,
    Convert::Simple(a_destFlub[2], ivl, a_destFlub[2].box(), a_srcFAB, ivl);)
}


// -----------------------------------------------------------------------------
// This version takes each comp of CC a_src and sends it to a comp of FC
// a_dest. For example, if a_src has 2 comps, then a_dest must also have
// 2 comps. Comp c of a_src will be sent to comp c of a_dest and will be
// interpolated to each face of the original cell.
// -----------------------------------------------------------------------------
void
Convert::CellsToAllFaces (LevelData<FluxBox>&         a_dest,
                          const LevelData<FArrayBox>& a_src)
{
    CH_assert(a_dest.getBoxes().compatible(a_src.getBoxes()));

    DataIterator dit = a_dest.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        CellsToAllFaces(a_dest[dit], a_src[dit]);
    }
}


// -----------------------------------------------------------------------------
// This version takes each comp and centering of FC a_src and sends it
// to a comp of CC a_dest. Suppose a_src has 2 comps, then a_dest
// must have 2*SpaceDim comps. If a_fcIsFastest == true, then comp c,
// FC in the d direction of a_src, will be sent to comp d + SpaceDim*c of
// a_dest. Otherwise, it wil be sent to comp c + a_src.nCpmp()*d.
// This works on the overlap of the src and dest boxes.
// -----------------------------------------------------------------------------
void
Convert::FacesToCells (LevelData<FArrayBox>&     a_dest,
                       const LevelData<FluxBox>& a_src,
                       const bool                a_fcIsFastest)
{
    CH_assert(a_dest.nComp() == SpaceDim * a_src.nComp());
    CH_assert(a_dest.getBoxes().compatible(a_src.getBoxes()));

    // I think this will fail if there is no ghost layer in the src.
    // Without a ghost, we will not fill the lower valid cells.
    // CH_assert(a_src.ghostVect() >= IntVect::Unit);

    const int numSrcComps = a_src.nComp();
    DataIterator dit = a_dest.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& destFAB = a_dest[dit];
        const FluxBox& srcFB = a_src[dit];
        int destComp = 0;

        if (a_fcIsFastest) {
            for (int srcComp = 0; srcComp < numSrcComps; ++srcComp) {
                for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                    Box b(srcFB.box());
                    b &= destFAB.box();

                    Convert::Simple(
                        destFAB, destComp++, b, srcFB[fcDir], srcComp);
                } // fcdir
            } // srcComp
        } else {
            for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                for (int srcComp = 0; srcComp < numSrcComps; ++srcComp) {
                    Box b(srcFB.box());
                    b &= destFAB.box();

                    Convert::Simple(
                        destFAB, destComp++, b, srcFB[fcDir], srcComp);
                } // srcComp
            } // fcDir
        } // if / if not fcIsFastest
    } // dit
}

