#include "StaggeredCopier.H"
#include "NeighborIterator.H"
#include "DataIterator.H"
// #include "MiscUtils.H"
#include "Debug.H"
// #include <algorithm>


// -----------------------------------------------------------------------------
StaggeredCopier::StaggeredCopier()
: Copier()
, m_ghostVect(D_DECL(-1, -1,- 1))
, m_fcDir(-1)
{
}


// -----------------------------------------------------------------------------
///
void
StaggeredCopier::define(const DisjointBoxLayout& /*a_srcLayout*/,
                        const BoxLayout&         /*a_destLayout*/,
                        bool                     /*a_exchange*/,
                        IntVect                  /*a_shift*/)
{
    UNDEFINED_FUNCTION();
}


// -----------------------------------------------------------------------------
/// contains support for periodic BCs
void
StaggeredCopier::define(const DisjointBoxLayout& /*a_srcLayout*/,
                        const BoxLayout&         /*a_destLayout*/,
                        const ProblemDomain&     /*a_domain*/,
                        bool                     /*a_exchange*/,
                        IntVect                  /*a_shift*/)
{
    UNDEFINED_FUNCTION();
}


// -----------------------------------------------------------------------------
///
void
StaggeredCopier::define(const DisjointBoxLayout& /*a_srcLayout*/,
                        const BoxLayout&         /*a_destLayout*/,
                        const IntVect&           /*a_destGhost*/,
                        bool                     /*a_exchange*/,
                        IntVect                  /*a_shift*/)
{
    UNDEFINED_FUNCTION();
}


// -----------------------------------------------------------------------------
/// contains support for periodic BCs
void
StaggeredCopier::define(const BoxLayout&     /*a_srcLayout*/,
                        const BoxLayout&     /*a_destLayout*/,
                        const ProblemDomain& /*a_domain*/,
                        const IntVect&       /*a_destGhost*/,
                        bool                 /*a_exchange*/,
                        IntVect              /*a_shift*/)
{
    UNDEFINED_FUNCTION();
}


// -----------------------------------------------------------------------------
///  alternative version of define for exchange copiers that uses new
///  optimized neighborIterator
void
StaggeredCopier::exchangeDefine(const DisjointBoxLayout& /*a_grids*/,
                                const IntVect&           /*a_ghost*/)
{
    UNDEFINED_FUNCTION();
}


// -----------------------------------------------------------------------------
void
StaggeredCopier::define(const BoxLayout&     a_srcLayout,
                        const BoxLayout&     a_destLayout,
                        const ProblemDomain& a_domain,
                        const IntVect&       a_destGhost,
                        const int            a_fcDir,
                        const IntVect        a_shift)
{
    CH_assert(a_srcLayout.isClosed());
    CH_assert(a_destLayout.isClosed());
    CH_assert(0 <= a_fcDir && a_fcDir < SpaceDim);

    this->clear();
    buffersAllocated = false;

    m_isDefined = true;
    m_ghostVect = a_destGhost;
    m_fcDir     = a_fcDir;

    const Box fcDomainBox = surroundingNodes(a_domain.domainBox(), a_fcDir);
    const unsigned int myProcID = (unsigned int)procID();

    // The following 4 for loops are the result of a performance optimization.
    // When increasing the size of the problem, we found that the code was
    // looping over every destination box for every source box which was N1*N2
    // loop iterations (essentially an N-squared approach).
    // The following code attempts to simply reduce N1 and N2 by first
    // separating the boxes (or LayoutIndexes to boxes) that reside on the
    // current processor. Then the loop to determine which boxes of the first
    // list intersect with which boxes of the second list can be done in N1' *
    // N2' iterations, where N1' is the reduced N1 and N2' is the reduced N2. We
    // have to break up the assigning of MotionItems into two separate loops and
    // be careful about the local copies.  These 4 loops are significantly
    // faster than the original for loop -- _especially_ for large problems.
    // (ndk)

#ifdef CH_MPI  // don't need to do this in serial
    // make a vector of boxes (or LayoutIndexes to boxes) from dest layout
    // that are known to reside on this processor.
    std::vector<DataIndex> vectorDestDI;
    std::vector<DataIndex> vectorDestOnProcDI;
    for (LayoutIterator to(a_destLayout.layoutIterator()); to.ok(); ++to) {
        vectorDestDI.push_back(DataIndex(to()));
        if (myProcID == a_destLayout.procID(to())) {
            vectorDestOnProcDI.push_back(DataIndex(to()));
        }
    }

    // make a vector of boxes (or LayoutIndexes to boxes) from src layout
    // that are known to reside on this processor.
    std::vector<DataIndex> vectorSrcDI;
    std::vector<DataIndex> vectorSrcOnProcDI;
    for (LayoutIterator from(a_srcLayout.layoutIterator()); from.ok(); ++from) {
        vectorSrcDI.push_back(DataIndex(from()));
        if (myProcID == a_srcLayout.procID(from())) {
            vectorSrcOnProcDI.push_back(DataIndex(from()));
        }
    }
#else
    // in serial, it's not very interesting as it's all of them.
    std::vector<DataIndex> vectorDestOnProcDI;
    std::vector<DataIndex> vectorSrcDI;
    for (LayoutIterator to(a_destLayout.layoutIterator()); to.ok(); ++to) {
        vectorDestOnProcDI.push_back(DataIndex(to()));
    }
    for (LayoutIterator from(a_srcLayout.layoutIterator()); from.ok(); ++from) {
        vectorSrcDI.push_back(DataIndex(from()));
    }
#endif

    const bool isSorted = (a_srcLayout.isSorted() && a_destLayout.isSorted());

    // Create MotionItems for (local or remote) srcValid --> (local) destFABBox.
    for (const auto& destDI : vectorDestOnProcDI) {
        const Box destFABBox = Box(a_destLayout[destDI] - a_shift)
                                   .grow(a_destGhost)
                                   .surroundingNodes(a_fcDir);

        for (const auto& srcDI : vectorSrcDI) {
            const unsigned int srcProcID = a_srcLayout.procID(srcDI);
            const Box srcValid = surroundingNodes(a_srcLayout[srcDI], a_fcDir);

            if (isSorted && (srcValid.bigEnd(0) < destFABBox.smallEnd(0)))
                continue;

            if (destFABBox.intersectsNotEmpty(srcValid)) {
                const Box srcBox  = srcValid & destFABBox;
                const Box destBox = srcBox + a_shift;
                CH_assert(srcBox.type() == BASISV(a_fcDir));

                MotionItem* item = new (s_motionItemPool.getPtr())
                    MotionItem(srcDI, destDI, srcBox, destBox);
                if (item == nullptr) {
                    MayDay::Error("Out of Memory in copier::define");
                }
                if (srcProcID == myProcID) {
                    // if (a_exchange && srcDI == destDI) {
                    //     s_motionItemPool.returnPtr(item);
                    // } else {
                        m_localMotionPlan.push_back(item);
                    // }
                } else {
                    item->procID = srcProcID;
                    m_toMotionPlan.push_back(item);
                }
            }

            if (isSorted && (srcValid.smallEnd(0) > destFABBox.bigEnd(0)))
                break;
        }
    }

#ifdef CH_MPI
    // Create MotionItems for (local) srcValid --> (remote) destFABBox.
    // We already took care of the local copy, so skip this in serial.
    for (const auto& destDI : vectorDestDI) {
        const unsigned int destProcID = a_destLayout.procID(destDI);
        if (destProcID == myProcID) continue; // We already did the local move.

        const Box destFABBox = Box(a_destLayout[destDI] - a_shift)
                                   .grow(a_destGhost)
                                   .surroundingNodes(a_fcDir);

        for (const auto& srcDI : vectorSrcOnProcDI) {
            // srcProcID = myProcID
            const Box srcValid = surroundingNodes(a_srcLayout[srcDI], a_fcDir);

            if (isSorted && (srcValid.bigEnd(0) < destFABBox.smallEnd(0)))
                continue;

            if (destFABBox.intersectsNotEmpty(srcValid)) {
                const Box srcBox  = srcValid & destFABBox;
                const Box destBox = srcBox + a_shift;
                CH_assert(srcBox.type() == BASISV(a_fcDir));

                MotionItem* item = new (s_motionItemPool.getPtr())
                    MotionItem(srcDI, destDI, srcBox, destBox);
                if (item == nullptr) {
                    MayDay::Error("Out of Memory in copier::define");
                }

                item->procID = destProcID;
                m_fromMotionPlan.push_back(item);
            }

            if (isSorted && (srcValid.smallEnd(0) > destFABBox.bigEnd(0)))
                break;
        }
    }
#endif

    // Everything below this point deals with FABs that extend beyond a periodic
    // domain. Skip it all if not periodic.
    {
        const bool isPeriodic =
            (fcDomainBox.isEmpty() ? false : a_domain.isPeriodic());

        if (!isPeriodic) {
            this->sort();
            return;
        }
    }

    // (dfm -- 9/13/05) as currently written, the Copier won't correctly
    // handle periodic cases where the number of ghost cells is greater
    // than the width of the domain.  We _should_ do multiple wraparounds,
    // but we don't. So, put in this assertion. We can revisit this if it
    // becomes an issue
    D_TERM(
    CH_assert(!a_domain.isPeriodic(0) || a_destGhost[0] <= fcDomainBox.size(0));,
    CH_assert(!a_domain.isPeriodic(1) || a_destGhost[1] <= fcDomainBox.size(1));,
    CH_assert(!a_domain.isPeriodic(2) || a_destGhost[2] <= fcDomainBox.size(2));)

    if (a_shift != IntVect::Zero) {
        MayDay::Error(
            "StaggeredCopier::define - domain periodic and a non-zero shift in "
            "the copy is not implemented");
    }

    // Set up vector of DIs to keep track of which dest boxes are not completely
    // contained within the primary domain. These boxes are then candidates for
    // filling by periodic images of the src data.
    std::vector<DataIndex> periodicDestVect;
    Box grownDestDomainCheckBox = grow(fcDomainBox, 1);
    int periodicDestCheckRadius = 1;

    const Box loDomainBdry = bdryBox(a_domain.domainBox(), a_fcDir, Side::Lo);
    const Box hiDomainBdry = bdryBox(a_domain.domainBox(), a_fcDir, Side::Hi);

    for (LayoutIterator lit(a_destLayout.layoutIterator()); lit.ok(); ++lit) {
        const LayoutIndex& destLI = lit();
        const Box destFABBox = surroundingNodes(a_destLayout[destLI], a_fcDir)
                              .grow(a_destGhost);

        if (!fcDomainBox.contains(destFABBox)) {
            // If dest hangs over the domain's edge, add it to the list.
            periodicDestVect.push_back(DataIndex(destLI));

            // Also, grow check radius to include all ghost faces.
            if (!grownDestDomainCheckBox.contains(destFABBox)) {
                while (!grownDestDomainCheckBox.contains(destFABBox)) {
                    grownDestDomainCheckBox.grow(1);
                    ++periodicDestCheckRadius;
                }
            }
        } else if (a_domain.isPeriodic(a_fcDir) &&
                   loDomainBdry.intersectsNotEmpty(destFABBox)) {
            // Add destDIs that lie AT the lower FC domain boundary.
            periodicDestVect.push_back(DataIndex(destLI));

        } else if (a_domain.isPeriodic(a_fcDir) &&
                   hiDomainBdry.intersectsNotEmpty(destFABBox)) {
            // Add destDIs that lie AT the upper FC domain boundary.
            periodicDestVect.push_back(DataIndex(destLI));
        }
    } // lit, destLI

    // The only src boxes we will need to check will be those within
    // periodicDestCheckRadius of the domain boundary. So, create a box to
    // screen out those which we will need to check.
    Box shrunkDomainBox = grow(fcDomainBox, -periodicDestCheckRadius);
    shrunkDomainBox.grow(a_fcDir, -1); // Added by ES.

    // Valid regions of the src DBL may also be outside the primary domain,
    // so we need to keep track of these too.
    std::vector<DataIndex> periodicSrcVect;
    Box grownSrcDomainCheckBox = grow(fcDomainBox, 1); // Grown by ES.
    int periodicSrcCheckRadius = 1;

    for (LayoutIterator lit(a_srcLayout.layoutIterator()); lit.ok(); ++lit) {
        const LayoutIndex& srcLI = lit();
        const Box srcValid = surroundingNodes(a_srcLayout[srcLI], a_fcDir);

        if (!shrunkDomainBox.contains(srcValid)) {
            if (!fcDomainBox.contains(srcValid)) {
                // If src hangs over the domain's edge, add it to the list.
                periodicSrcVect.push_back(DataIndex(srcLI));

                // Also, grow check radius.
                if (!grownSrcDomainCheckBox.contains(srcValid)) {
                    while (!grownSrcDomainCheckBox.contains(srcValid)) {
                        grownSrcDomainCheckBox.grow(1);
                        ++periodicSrcCheckRadius;
                    }
                }
            } else if (a_domain.isPeriodic(a_fcDir) &&
                       loDomainBdry.intersectsNotEmpty(srcValid)) {
                // Add srcDIs that lie AT the lower FC domain boundary.
                periodicSrcVect.push_back(DataIndex(srcLI));

            } else if (a_domain.isPeriodic(a_fcDir) &&
                       hiDomainBdry.intersectsNotEmpty(srcValid)) {
                // Add srcDIs that lie AT the upper FC domain boundary.
                periodicSrcVect.push_back(DataIndex(srcLI));
            }
        }
    } // lit, srcLI

    // At this point, periodicDestVect and periodicSrcVect contain a complete
    // list of candidates for data transfer.
    // Now, we actually do the periodic checking.

    for (LayoutIterator lit(a_srcLayout.layoutIterator()); lit.ok(); ++lit) {
        const LayoutIndex& srcLI = lit();
        const unsigned int srcProcID = a_srcLayout.procID(srcLI);
        const Box srcValid = surroundingNodes(a_srcLayout[srcLI], a_fcDir);

        if (shrunkDomainBox.contains(srcValid)) continue;

        // Loop over those dest boxes which were not contained by the domain.
        for (const auto& destDI : periodicDestVect) {
            const unsigned int destProcID = a_destLayout.procID(destDI);
            if (destProcID != myProcID && srcProcID != myProcID) continue;

            Box destFABBox = Box(a_destLayout[destDI])
                                 .grow(a_destGhost)
                                 .surroundingNodes(a_fcDir);

            // destFABBox is temporarily shifted to overlap srcValid.
            // destBox is the *unshifted* overlap region.
            for (auto sit = a_domain.shiftIterator(); sit.ok(); ++sit) {
                const IntVect shiftVect(sit() * a_domain.size());
                destFABBox.shift(shiftVect);

                if (destFABBox.intersectsNotEmpty(srcValid)) {
                    const Box srcBox(destFABBox & srcValid);
                    const Box destBox = Box(srcBox).shift(-shiftVect);

                    CH_assert(srcBox.type() == BASISV(a_fcDir));
                    CH_assert(destBox.type() == BASISV(a_fcDir));

                    MotionItem* item = new (s_motionItemPool.getPtr())
                        MotionItem(DataIndex(srcLI),
                                    DataIndex(destDI),
                                    srcBox,
                                    destBox);
                    if (item == nullptr) {
                        MayDay::Error("Out of Memory in copier::define");
                    }

                    if (destProcID == srcProcID) {
                        m_localMotionPlan.push_back(item);
                    } else if (srcProcID == myProcID) {
                        item->procID = destProcID;
                        m_fromMotionPlan.push_back(item);
                    } else {
                        item->procID = srcProcID;
                        m_toMotionPlan.push_back(item);
                    }
                } // end if shifted box intersects
                destFABBox.shift(-shiftVect);
            }  // end loop over shift vectors, sit
        } // destDI
    } // lit, srcLI

    // Loop through the candidate src boxes.
    if (periodicSrcVect.size() != 0) {
        // The only dest boxes we will need to check are those within
        // periodicSrcCheckRadius of the domain boundary. So, create a box to
        // screen out those which we will need to check.
        shrunkDomainBox = grow(fcDomainBox, -periodicSrcCheckRadius);
        shrunkDomainBox.grow(a_fcDir, -1); // Added by ES.

        for (LayoutIterator lit(a_destLayout.layoutIterator()); lit.ok(); ++lit) {
            const LayoutIndex& destLI = lit();
            const unsigned int destProcID = a_destLayout.procID(destLI);
            Box destFABBox = surroundingNodes(a_destLayout[destLI], a_fcDir)
                                 .grow(a_destGhost);

            if (shrunkDomainBox.contains(destFABBox)) continue;

            // Loop over those src boxes which are not contained by the domain.
            for (const auto& srcDI : periodicSrcVect) {
                const unsigned int srcProcID = a_srcLayout.procID(srcDI);
                if (destProcID != myProcID && srcProcID != myProcID) continue;

                const Box srcValid = surroundingNodes(a_srcLayout[srcDI], a_fcDir);

                // Again, we temporarily shift destFABBox to overlap srcValid.
                for (auto sit = a_domain.shiftIterator(); sit.ok(); ++sit) {
                    const IntVect shiftVect(sit() * a_domain.size());
                    destFABBox.shift(shiftVect);

                    if (destFABBox.intersectsNotEmpty(srcValid)) {
                        const Box srcBox(destFABBox & srcValid);
                        const Box destBox = Box(srcBox).shift(-shiftVect);

                        CH_assert(srcBox.type() == BASISV(a_fcDir));
                        CH_assert(destBox.type() == BASISV(a_fcDir));

                        MotionItem* item = new (s_motionItemPool.getPtr())
                            MotionItem(DataIndex(srcDI),
                                       DataIndex(destLI),
                                       srcBox,
                                       destBox);
                        if (item == nullptr) {
                            MayDay::Error("Out of Memory in copier::define");
                        }

                        if (destProcID == srcProcID) {
                            m_localMotionPlan.push_back(item);
                        } else if (srcProcID == myProcID) {
                            item->procID = destProcID;
                            m_fromMotionPlan.push_back(item);
                        } else {
                            item->procID = srcProcID;
                            m_toMotionPlan.push_back(item);
                        }
                    }  // end if shifted box intersects

                    destFABBox.shift(-shiftVect);
                }  // end loop over shift vectors, sit
            } // srcDI
        } // lit, destLI
    } // end if any periodicSrcVects

    this->sort();
}


// -----------------------------------------------------------------------------
void
StaggeredCopier::ghostDefine(const DisjointBoxLayout& /*a_srcGrids*/,
                             const DisjointBoxLayout& /*a_destGrids*/,
                             const ProblemDomain&     /*a_domain*/,
                             const IntVect&           /*a_srcGhost*/)
{
    // UNDEFINED_FUNCTION();
    MAYDAYERROR("You attempted to call StaggeredCopier::ghostDefine. "
                "Just so that I know which define() function to use, please "
                "call either StaggeredCopier::standardGhostDefine, or "
                "StaggeredCopier::staggeredGhostDefine.");
}


// -----------------------------------------------------------------------------
void
StaggeredCopier::ghostDefine(const DisjointBoxLayout& a_srcGrids,
                             const DisjointBoxLayout& a_destGrids,
                             const ProblemDomain&     a_domain,
                             const IntVect&           a_srcGhost,
                             const int                a_fcDir)
{
    // first, define a regular copier operation
    this->define(a_destGrids, a_srcGrids, a_domain, a_srcGhost, a_fcDir);

    // now, reverse the direction of the operation.
    this->reverse();

    this->sort();
}


// -----------------------------------------------------------------------------
void
StaggeredCopier::defineValidExchange(const DisjointBoxLayout& a_grids,
                                     const IntVect&           a_ghostVect,
                                     const int                a_fcDir,
                                     const bool               a_doValidCorners)
{
    CH_assert(0 <= a_fcDir && a_fcDir < SpaceDim);

    this->clear();

    m_isDefined = true;
    m_ghostVect = a_ghostVect;
    m_fcDir     = a_fcDir;

    if (a_ghostVect.product() == 0) return;

    DataIterator         dit    = a_grids.dataIterator();
    NeighborIterator     nit(a_grids);
    const int            myRank = procID();

    for (dit.begin(); dit.ok(); ++dit) {
        // 1. Valid neighbors --> our ghosts
        Box interior = a_grids[dit];
        for (int boffset = 0; boffset < SpaceDim; ++boffset) {
            const int bdryDir = (m_fcDir + boffset) % SpaceDim;
            const int ighost  = m_ghostVect[bdryDir];

            if (ighost == 0) continue;

            for (SideIterator bdrySit; bdrySit.ok(); ++bdrySit) {
                const Side::LoHiSide& side  = bdrySit();
                const int             isign = sign(side);

                const Box myGhost = [&] {
                    Box bx = adjCellBox(interior, bdryDir, side, ighost);

                    if (bdryDir == m_fcDir) {
                        bx.shiftHalf(m_fcDir, isign);
                    } else {
                        bx.surroundingNodes(m_fcDir);
                    }

                    return bx;
                }();

                for (nit.begin(dit()); nit.ok(); ++nit) {
                    const int nbrRank  = a_grids.procID(nit());
                    const Box nbrValid = Box(nit.box()).surroundingNodes(m_fcDir);

                    if (nbrValid.intersectsNotEmpty(myGhost)) {
                        const Box overlap(nbrValid & myGhost);

                        MotionItem* item = new (s_motionItemPool.getPtr())
                            MotionItem(DataIndex(nit()),
                                       dit(),
                                       nit.unshift(overlap),
                                       overlap);

                        if (nbrRank == myRank) {  // local move
                            m_localMotionPlan.push_back(item);
                        } else {
                            item->procID = nbrRank;
                            m_toMotionPlan.push_back(item);
                        }
                    }
                } // nit
            } // bdrySit

            if (a_doValidCorners) {
                interior.grow(bdryDir, ighost);
            }
        } // bdryDir

        // 2. Our valid --> neighbor ghosts
        for (nit.begin(dit()); nit.ok(); ++nit) {
            const int nbrRank = a_grids.procID(nit());

            // We already created the local plan.
            if (nbrRank == myRank) continue;

            interior = nit.box();
            for (int boffset = 0; boffset < SpaceDim; ++boffset) {
                const int bdryDir = (m_fcDir + boffset) % SpaceDim;
                const Box myValid = Box(a_grids[dit]).surroundingNodes(m_fcDir);
                const int ighost  = m_ghostVect[bdryDir];

                if (ighost == 0) continue;

                for (SideIterator bdrySit; bdrySit.ok(); ++bdrySit) {
                    const Side::LoHiSide& side  = bdrySit();
                    const int             isign = sign(side);

                    const Box nbrGhost = [&] {
                        Box bx = adjCellBox(interior, bdryDir, side, ighost);

                        if (bdryDir == m_fcDir) {
                            bx.shiftHalf(m_fcDir, isign);
                        } else {
                            bx.surroundingNodes(m_fcDir);
                        }

                        return bx;
                    }();

                    if (myValid.intersectsNotEmpty(nbrGhost)) {
                        const Box overlap(myValid & nbrGhost);

                        MotionItem* item = new (s_motionItemPool.getPtr())
                            MotionItem(dit(),
                                       DataIndex(nit()),
                                       overlap,
                                       nit.unshift(overlap));

                        item->procID = nbrRank;
                        m_fromMotionPlan.push_back(item);
                    }
                } // nit

                if (a_doValidCorners) {
                    interior.grow(bdryDir, ighost);
                }
            } // bdrySit
        } // bdryDir
    } // dit

    this->sort();
}


// -----------------------------------------------------------------------------
void
StaggeredCopier::defineInvalidCornerExchange1(const DisjointBoxLayout& a_grids,
                                              const IntVect& a_ghostVect,
                                              const int      a_fcDir)
{
    CH_assert(0 <= a_fcDir && a_fcDir < SpaceDim);

    this->clear();

    m_isDefined = true;
    m_ghostVect = a_ghostVect;
    m_fcDir     = a_fcDir;

    if (a_ghostVect.product() == 0) return;

    DataIterator     dit = a_grids.dataIterator();
    NeighborIterator nit(a_grids);
    const int        myRank = procID();

    // The ghosts that should have already been filled/validated using
    // StaggeredCopier::defineValidExchange will be reffered to as
    // "valid" data thoughout this function.

    for (dit.begin(); dit.ok(); ++dit) {
        // 1. Neighbors' "valid" faces --> our ghosts
        for (int dir0off = 0; dir0off < SpaceDim; ++dir0off) {
            const int dir0 = (m_fcDir + dir0off) % SpaceDim;
            const int dir1 = (dir0 + 1) % SpaceDim;  // 1st bdry dir
            const int dir2 = (dir0 + 2) % SpaceDim;  // 2nd bdry dir in 3D

            const int ghost1 = m_ghostVect[dir1];
            const int ghost2 = m_ghostVect[dir2];
            if (ghost1 == 0 || ghost2 == 0) continue;

            for (SideIterator sit1; sit1.ok(); ++sit1) {
                const Side::LoHiSide& side1 = sit1();
                const int             sign1 = sign(side1);

                for (SideIterator sit2; sit2.ok(); ++sit2) {
                    const Side::LoHiSide& side2 = sit2();
                    const int             sign2 = sign(side2);

                    const Box myGhost = [&]{
                        Box bx = adjCellBox(a_grids[dit], dir1, side1, ghost1);
                        bx     = adjCellBox(bx, dir2, side2, ghost2);

                        if (dir0 == m_fcDir) {
                            bx.surroundingNodes(m_fcDir);
                        } else if (dir1 == m_fcDir) {
                            bx.shiftHalf(dir1, sign1);
                        } else {
                            bx.shiftHalf(dir2, sign2);
                        }

                        return bx;
                    }();

                    for (nit.begin(dit()); nit.ok(); ++nit) {
                        const int nbrRank = a_grids.procID(nit());

                        Box nbrInterior = nit.box();
                        for (int nbrDirOff = 0; nbrDirOff < SpaceDim; ++nbrDirOff) {
                            const int nbrDir = (m_fcDir + nbrDirOff) % SpaceDim;
                            const int nbrGhost  = m_ghostVect[nbrDir];

                            if (nbrGhost == 0) continue;

                            for (SideIterator nbrSit; nbrSit.ok(); ++nbrSit) {
                                const Side::LoHiSide& nbrSide = nbrSit();
                                const int             nbrSign = sign(nbrSide);

                                Box nbrValid = adjCellBox(nbrInterior, nbrDir, nbrSide, nbrGhost);
                                if (nbrDir == m_fcDir) {
                                    nbrValid.shiftHalf(m_fcDir, nbrSign);
                                } else {
                                    nbrValid.surroundingNodes(m_fcDir);
                                }

                                // TODO: Remove the truly valid faces from nbrValid.

                                if (nbrValid.intersectsNotEmpty(myGhost)) {
                                    const Box overlap(nbrValid & myGhost);

                                    MotionItem* item =
                                        new (s_motionItemPool.getPtr())
                                            MotionItem(DataIndex(nit()),
                                                       dit(),
                                                       nit.unshift(overlap),
                                                       overlap);

                                    if (nbrRank == myRank) {  // local move
                                        m_localMotionPlan.push_back(item);
                                    } else {
                                        item->procID = nbrRank;
                                        m_toMotionPlan.push_back(item);
                                    }
                                }
                            } // nbrSit
                        } // nbrDir
                    }  // nit
                } // sit2
            } // sit1
        } // directions

        // 2. Our "valid" faces --> neighbor ghosts
        for (nit.begin(dit()); nit.ok(); ++nit) {
            // We already made the local motion plan.
            const int nbrRank = a_grids.procID(nit());
            if (myRank == nbrRank) continue;

            for (int dir0off = 0; dir0off < SpaceDim; ++dir0off) {
                const int dir0 = (m_fcDir + dir0off) % SpaceDim;
                const int dir1 = (dir0 + 1) % SpaceDim;  // 1st bdry dir
                const int dir2 = (dir0 + 2) % SpaceDim;  // 2nd bdry dir in 3D

                const int ghost1 = m_ghostVect[dir1];
                const int ghost2 = m_ghostVect[dir2];
                if (ghost1 == 0 || ghost2 == 0) continue;

                for (SideIterator sit1; sit1.ok(); ++sit1) {
                    const Side::LoHiSide& side1 = sit1();
                    const int             sign1 = sign(side1);

                    for (SideIterator sit2; sit2.ok(); ++sit2) {
                        const Side::LoHiSide& side2 = sit2();
                        const int             sign2 = sign(side2);

                        const Box nbrGhost = [&]{
                            Box bx = adjCellBox(nit.box(), dir1, side1, ghost1);
                            bx     = adjCellBox(bx, dir2, side2, ghost2);

                            if (dir0 == m_fcDir) {
                                bx.surroundingNodes(m_fcDir);
                            } else if (dir1 == m_fcDir) {
                                bx.shiftHalf(dir1, sign1);
                            } else {
                                bx.shiftHalf(dir2, sign2);
                            }

                            return bx;
                        }();

                        Box myInterior = a_grids[dit];
                        for (int myDirOff = 0; myDirOff < SpaceDim; ++myDirOff) {
                            const int myDir = (m_fcDir + myDirOff) % SpaceDim;
                            const int myGhost  = m_ghostVect[myDir];

                            if (myGhost == 0) continue;

                            for (SideIterator mySit; mySit.ok(); ++mySit) {
                                const Side::LoHiSide& mySide  = mySit();
                                const int             mySign = sign(mySide);

                                Box myValid = adjCellBox(myInterior, myDir, mySide, myGhost);
                                if (myDir == m_fcDir) {
                                    myValid.shiftHalf(m_fcDir, mySign);
                                } else {
                                    myValid.surroundingNodes(m_fcDir);
                                }

                                // TODO: Remove the truly valid faces from myValid.

                                if (myValid.intersectsNotEmpty(nbrGhost)) {
                                    const Box overlap(myValid & nbrGhost);

                                    MotionItem* item =
                                        new (s_motionItemPool.getPtr())
                                            MotionItem(dit(),
                                                       DataIndex(nit()),
                                                       overlap,
                                                       nit.unshift(overlap));

                                    item->procID = nbrRank;
                                    m_fromMotionPlan.push_back(item);
                                }
                            } // mySit
                        } // myDirOff
                    } // sit2
                } // sit1
            } // directions
        } // nit
    } // dit

    this->sort();
}


// -----------------------------------------------------------------------------
void
StaggeredCopier::defineInvalidCornerExchange2(const DisjointBoxLayout& a_grids,
                                              const IntVect& a_ghostVect,
                                              const int      a_fcDir)
{
    CH_assert(0 <= a_fcDir && a_fcDir < SpaceDim);

    this->clear();

    m_isDefined = true;
    m_ghostVect = a_ghostVect;
    m_fcDir     = a_fcDir;

    if (SpaceDim == 2 || m_ghostVect.product() == 0) return;

    DataIterator     dit = a_grids.dataIterator();
    NeighborIterator nit(a_grids);
    const int        myRank = procID();

    // The ghosts that should have already been filled/validated using
    // StaggeredCopier::defineValidExchange will be reffered to as
    // "valid" data thoughout this function.

    for (dit.begin(); dit.ok(); ++dit) {
        // 3. Neighbor's non-valid edges* --> our non-valid vertices
        // *Stay within the domain.
        for (SideIterator mySit0; mySit0.ok(); ++mySit0) {
        for (SideIterator mySit1; mySit1.ok(); ++mySit1) {
        for (SideIterator mySit2; mySit2.ok(); ++mySit2) {
            const int myDir0 = m_fcDir;
            const int myDir1 = (m_fcDir + 1) % SpaceDim;
            const int myDir2 = (m_fcDir + 2) % SpaceDim;

            const Box myVertex = [&]{
                Box bx = a_grids[dit];
                bx = adjCellBox(bx, myDir0, mySit0(), m_ghostVect[myDir0]);
                bx = adjCellBox(bx, myDir1, mySit1(), m_ghostVect[myDir1]);
                bx = adjCellBox(bx, myDir2, mySit2(), m_ghostVect[myDir2]);
                bx.shiftHalf(m_fcDir, sign(mySit0()));
                return bx;
            }();

            // // Skip valid faces.
            // {
            //     bool isValid = false;
            //     for (nit.begin(dit()); nit.ok(); ++nit) {
            //         const Box nbrInterior = surroundingNodes(nit.box(), m_fcDir);
            //         if (nbrInterior.intersectsNotEmpty(myVertex)) {
            //             isValid = true;
            //             break;
            //         }
            //     } // nit
            //     if (isValid) continue;
            // }

            for (nit.begin(dit()); nit.ok(); ++nit) {
                const int nbrRank  = a_grids.procID(nit());

                for (int nbrDir0off = 0; nbrDir0off < SpaceDim; ++nbrDir0off) {
                    const int nbrDir0 = (m_fcDir + nbrDir0off) % SpaceDim;
                    const int nbrDir1 = (nbrDir0 + 1) % SpaceDim;  // 1st bdry dir
                    const int nbrDir2 = (nbrDir0 + 2) % SpaceDim;  // 2nd bdry dir in 3D
                    const int nbrGhost1 = m_ghostVect[nbrDir1];
                    const int nbrGhost2 = m_ghostVect[nbrDir2];

                    if (nbrGhost1 == 0 || nbrGhost2 == 0) continue;

                    for (SideIterator nbrSit1; nbrSit1.ok(); ++nbrSit1) {
                    for (SideIterator nbrSit2; nbrSit2.ok(); ++nbrSit2) {
                        const Side::LoHiSide& nbrSide1 = nbrSit1();
                        const int             nbrSign1 = sign(nbrSide1);
                        const Side::LoHiSide& nbrSide2 = nbrSit2();
                        const int             nbrSign2 = sign(nbrSide2);

                        const Box nbrEdge = [&] {
                            Box bx = nit.box();
                            bx = adjCellBox(bx, nbrDir1, nbrSide1, nbrGhost1);
                            bx = adjCellBox(bx, nbrDir2, nbrSide2, nbrGhost2);

                            if (nbrDir0 == m_fcDir) {
                                bx.surroundingNodes(m_fcDir);
                            } else if (nbrDir1 == m_fcDir) {
                                bx.shiftHalf(nbrDir1, nbrSign1);
                            } else {
                                bx.shiftHalf(nbrDir2, nbrSign2);
                            }

                            return bx;
                        }();

                        if (nbrEdge.intersectsNotEmpty(myVertex)) {
                            const Box overlap(nbrEdge & myVertex);

                            MotionItem* item = new (s_motionItemPool.getPtr())
                                MotionItem(DataIndex(nit()),
                                           dit(),
                                           nit.unshift(overlap),
                                           overlap);

                            if (nbrRank == myRank) {  // local move
                                m_localMotionPlan.push_back(item);
                            } else {
                                item->procID = nbrRank;
                                m_toMotionPlan.push_back(item);
                            }
                        }
                    } // nbrSit2
                    } // nbrSit1
                } // nbrDirs
            } // nit
        } // mySit2
        } // mySit1
        } // mySit0

        // 4. Our non-valid edges* --> neighbor's non-valid vertices
        // *Stay within the domain.
        for (nit.begin(dit()); nit.ok(); ++nit) {
            // We already created the local plan.
            const int nbrRank  = a_grids.procID(nit());
            if (nbrRank == myRank) continue;

            for (SideIterator nbrSit0; nbrSit0.ok(); ++nbrSit0) {
            for (SideIterator nbrSit1; nbrSit1.ok(); ++nbrSit1) {
            for (SideIterator nbrSit2; nbrSit2.ok(); ++nbrSit2) {
                const int nbrDir0 = m_fcDir;
                const int nbrDir1 = (m_fcDir + 1) % SpaceDim;
                const int nbrDir2 = (m_fcDir + 2) % SpaceDim;

                const Box nbrVertex = [&]{
                    Box bx = nit.box();
                    bx = adjCellBox(bx, nbrDir0, nbrSit0(), m_ghostVect[nbrDir0]);
                    bx = adjCellBox(bx, nbrDir1, nbrSit1(), m_ghostVect[nbrDir1]);
                    bx = adjCellBox(bx, nbrDir2, nbrSit2(), m_ghostVect[nbrDir2]);
                    bx.shiftHalf(m_fcDir, sign(nbrSit0()));
                    return bx;
                }();

                // // Skip valid faces.
                // {
                //     bool isValid = false;
                //     NeighborIterator nit2(a_grids);
                //     for (nit2.begin(DataIndex(nit())); nit2.ok(); ++nit2) {
                //         const Box nInterior = surroundingNodes(nit2.box(), m_fcDir);
                //         if (nInterior.intersectsNotEmpty(nbrVertex)) {
                //             isValid = true;
                //             break;
                //         }
                //     } // nit2
                //     if (isValid) continue;
                // }

                for (int myDir0off = 0; myDir0off < SpaceDim; ++myDir0off) {
                    const int myDir0 = (m_fcDir + myDir0off) % SpaceDim;
                    const int myDir1 = (myDir0 + 1) % SpaceDim;  // 1st bdry dir
                    const int myDir2 = (myDir0 + 2) % SpaceDim;  // 2nd bdry dir in 3D
                    const int myGhost1 = m_ghostVect[myDir1];
                    const int myGhost2 = m_ghostVect[myDir2];

                    if (myGhost1 == 0 || myGhost2 == 0) continue;

                    for (SideIterator mySit1; mySit1.ok(); ++mySit1) {
                    for (SideIterator mySit2; mySit2.ok(); ++mySit2) {
                        const Side::LoHiSide& mySide1 = mySit1();
                        const int             mySign1 = sign(mySide1);
                        const Side::LoHiSide& mySide2 = mySit2();
                        const int             mySign2 = sign(mySide2);

                        const Box myEdge = [&] {
                            Box bx = a_grids[dit];
                            bx = adjCellBox(bx, myDir1, mySide1, myGhost1);
                            bx = adjCellBox(bx, myDir2, mySide2, myGhost2);

                            if (myDir0 == m_fcDir) {
                                bx.surroundingNodes(m_fcDir);
                            } else if (myDir1 == m_fcDir) {
                                bx.shiftHalf(myDir1, mySign1);
                            } else {
                                bx.shiftHalf(myDir2, mySign2);
                            }

                            return bx;
                        }();

                        if (myEdge.intersectsNotEmpty(nbrVertex)) {
                            const Box overlap(myEdge & nbrVertex);

                            MotionItem* item = new (s_motionItemPool.getPtr())
                                MotionItem(dit(),
                                           DataIndex(nit()),
                                           overlap,
                                           nit.unshift(overlap));

                            item->procID = nbrRank;
                            m_fromMotionPlan.push_back(item);
                        }
                    } // mySit2
                    } // mySit1
                } // myDirs
            } // nbrSit2
            } // nbrSit1
            } // nbrSit0
        } // nit
    } // dit

    this->sort();
}

