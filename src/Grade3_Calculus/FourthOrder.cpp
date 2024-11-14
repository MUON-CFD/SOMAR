#include "FourthOrder.H"
#include "FourthOrderF_F.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
void
FourthOrder::nodalToAvg(FArrayBox& a_dataFAB, const Box& a_conversionBox)
{
    const auto& boxType = a_dataFAB.box().type();
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (boxType[dir] == IndexType::NODE) continue;
        FourthOrder::nodalToAvg(a_dataFAB, a_conversionBox, dir);
    }
}


// -----------------------------------------------------------------------------
void
FourthOrder::avgToNodal(FArrayBox& a_dataFAB, const Box& a_conversionBox)
{
    const auto& boxType = a_dataFAB.box().type();
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (boxType[dir] == IndexType::NODE) continue;
        FourthOrder::avgToNodal(a_dataFAB, a_conversionBox, dir);
    }
}


// -----------------------------------------------------------------------------
void
FourthOrder::nodalToAvg(FArrayBox& a_dataFAB,
                        const Box& a_conversionBox,
                        const int  a_conversionDir,
                        const Real a_scale)
{
    CH_assert(a_dataFAB.box().type() == a_conversionBox.type());
    CH_assert(a_dataFAB.contains(a_conversionBox));
    CH_assert(0 <= a_conversionDir && a_conversionDir < SpaceDim);

    Box loBox, hiBox;
    Box midBox = a_conversionBox;
    if (a_conversionBox.type(a_conversionDir) == IndexType::CELL) {
        if (midBox.smallEnd(a_conversionDir) ==
            a_dataFAB.box().smallEnd(a_conversionDir)) {
            loBox =
                adjCellLo(midBox, a_conversionDir).shift(a_conversionDir, 1);
            midBox.growLo(a_conversionDir, -1);
        }
        if (midBox.bigEnd(a_conversionDir) ==
            a_dataFAB.box().bigEnd(a_conversionDir)) {
            hiBox =
                adjCellHi(midBox, a_conversionDir).shift(a_conversionDir, -1);
            midBox.growHi(a_conversionDir, -1);
        }
    } else {
        if (midBox.smallEnd(a_conversionDir) ==
            a_dataFAB.box().smallEnd(a_conversionDir)) {
            loBox = bdryLo(midBox, a_conversionDir);
            midBox.growLo(a_conversionDir, -1);
        }
        if (midBox.bigEnd(a_conversionDir) ==
            a_dataFAB.box().bigEnd(a_conversionDir)) {
            hiBox = bdryHi(midBox, a_conversionDir);
            midBox.growHi(a_conversionDir, -1);
        }
    }

    FArrayBox srcFAB(a_dataFAB.box(), a_dataFAB.nComp());
    srcFAB.copy(a_dataFAB);

    FORT_FOURTHORDER_NODALTOAVG(CHF_FRA(a_dataFAB),
                                CHF_FRA(srcFAB),
                                CHF_BOX(loBox),
                                CHF_BOX(midBox),
                                CHF_BOX(hiBox),
                                CHF_CONST_INT(a_conversionDir),
                                CHF_CONST_REAL(a_scale));
}


// -----------------------------------------------------------------------------
void
FourthOrder::avgToNodal(LevelData<FluxBox>& a_data)
{
    const auto& grids = a_data.getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& dataFAB = a_data[dit][fcDir];
            const Box  fcValid = grids[dit].surroundingNodes(fcDir);

            D_TERM(,
            FourthOrder::avgToNodal(dataFAB, fcValid, (fcDir + 1) % SpaceDim);,
            FourthOrder::avgToNodal(dataFAB, fcValid, (fcDir + 2) % SpaceDim);)
        }
    }
}


// -----------------------------------------------------------------------------
void
FourthOrder::nodalToAvg(LevelData<FluxBox>& a_data)
{
    const auto& grids = a_data.getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& dataFAB = a_data[dit][fcDir];
            const Box  fcValid = grids[dit].surroundingNodes(fcDir);

            D_TERM(,
            FourthOrder::nodalToAvg(dataFAB, fcValid, (fcDir + 1) % SpaceDim);,
            FourthOrder::nodalToAvg(dataFAB, fcValid, (fcDir + 2) % SpaceDim);)
        }
    }
}


// -----------------------------------------------------------------------------
void
FourthOrder::conservativeInterp(FArrayBox&       a_destFAB,
                                const int        a_destFABComp,
                                const Box&       a_destBox,
                                const FArrayBox& a_srcFAB,
                                const int        a_srcFABComp,
                                const int        a_interpDir)
{
    CH_assert(a_destFAB.box().sameType(a_destBox));
    CH_assert(a_destFAB.box().contains(a_destBox));

    CH_assert(0 <= a_destFABComp && a_destFABComp < a_destFAB.nComp());
    CH_assert(0 <= a_srcFABComp && a_srcFABComp < a_srcFAB.nComp());

    Box lolo, lo, mid, hi, hihi;
    FourthOrder::getDestBoxes(
        lolo, lo, mid, hi, hihi, a_srcFAB.box(), a_destBox, a_interpDir);

    if (a_destBox.type(a_interpDir) == IndexType::NODE) {
        FORT_FOURTHORDER_CONSERVATIVEINTERP1D_CELLTONODE(
            CHF_FRA1(a_destFAB, a_destFABComp),
            CHF_CONST_FRA1(a_srcFAB, a_srcFABComp),
            CHF_BOX(lolo),
            CHF_BOX(lo),
            CHF_BOX(mid),
            CHF_BOX(hi),
            CHF_BOX(hihi),
            CHF_CONST_INT(a_interpDir));

    } else {
        CH_assert(lolo.isEmpty());
        CH_assert(hihi.isEmpty());

        FORT_FOURTHORDER_CONSERVATIVEINTERP1D_NODETOCELL(
            CHF_FRA1(a_destFAB, a_destFABComp),
            CHF_CONST_FRA1(a_srcFAB, a_srcFABComp),
            CHF_BOX(lo),
            CHF_BOX(mid),
            CHF_BOX(hi),
            CHF_CONST_INT(a_interpDir));
    }
}


// -----------------------------------------------------------------------------
void
FourthOrder::interp(FArrayBox&       a_destFAB,
                    const int        a_destFABComp,
                    const Box&       a_destBox,
                    const FArrayBox& a_srcFAB,
                    const int        a_srcFABComp,
                    const int        a_interpDir)
{
    CH_assert(a_destFAB.box().sameType(a_destBox));
    CH_assert(a_destFAB.box().contains(a_destBox));

    CH_assert(0 <= a_destFABComp && a_destFABComp < a_destFAB.nComp());
    CH_assert(0 <= a_srcFABComp && a_srcFABComp < a_srcFAB.nComp());

    Box lolo, lo, mid, hi, hihi;
    FourthOrder::getDestBoxes(
        lolo, lo, mid, hi, hihi, a_srcFAB.box(), a_destBox, a_interpDir);
    CH_assert(lolo.isEmpty());
    CH_assert(hihi.isEmpty());

    if (a_destBox.type(a_interpDir) == IndexType::NODE) {
        FORT_FOURTHORDER_INTERP1D_CELLTONODE(
            CHF_FRA1(a_destFAB, a_destFABComp),
            CHF_CONST_FRA1(a_srcFAB, a_srcFABComp),
            CHF_BOX(lo),
            CHF_BOX(mid),
            CHF_BOX(hi),
            CHF_CONST_INT(a_interpDir));

    } else {
        FORT_FOURTHORDER_INTERP1D_NODETOCELL(
            CHF_FRA1(a_destFAB, a_destFABComp),
            CHF_CONST_FRA1(a_srcFAB, a_srcFABComp),
            CHF_BOX(lo),
            CHF_BOX(mid),
            CHF_BOX(hi),
            CHF_CONST_INT(a_interpDir));
    }
}


// // -----------------------------------------------------------------------------
// void
// FiniteDiff::multiplyFourthOrderAverages(FArrayBox&       a_fg,
//                                         const int        a_fgComp,
//                                         const FArrayBox& a_f,
//                                         const int        a_fComp,
//                                         const FArrayBox& a_g,
//                                         const int        a_gComp,
//                                         const Box&       a_validBox)
// {
//     // Sanity checks.
//     CH_assert(0 <= a_fgComp && a_fgComp < a_fg.nComp());
//     CH_assert(0 <= a_fComp && a_fComp < a_f.nComp());
//     CH_assert(0 <= a_gComp && a_gComp < a_g.nComp());

//     CH_assert(a_fg.box().type() == a_f.box().type());
//     CH_assert(a_fg.box().type() == a_g.box().type());
//     CH_assert(a_fg.box().type() == a_validBox.type());

//     CH_assert(a_validBox.contains(a_fg.box()));
//     CH_assert(a_validBox.contains(a_f.box()));
//     CH_assert(a_validBox.contains(a_g.box()));

//     a_fg.copy(a_f, a_validBox, a_fComp, a_validBox, a_fgComp, 1);
//     a_fg.mult(a_g, a_validBox, a_gComp, a_fgComp, 1);

//     for (int avgDir = 0; avgDir < SpaceDim; ++avgDir) {
//         if (a_validBox.type(avgDir) == IndexType::NODE) continue;

//         const Box loBox     = adjCellLo(a_validBox, avgDir, -1) & a_f.box();
//         const Box centerBox = grow(a_validBox, -BASISV(avgDir)) & a_f.box();
//         const Box hiBox     = adjCellHi(a_validBox, avgDir, -1) & a_f.box();

//         CH_assert(loBox.isEmpty() || a_fg.box().contains(loBox));
//         CH_assert(a_fg.box().contains(centerBox));
//         CH_assert(hiBox.isEmpty() || a_fg.box().contains(hiBox));

//         CH_assert(loBox.isEmpty() || a_f.box().contains(loBox));
//         CH_assert(a_f.box().contains(centerBox));
//         CH_assert(hiBox.isEmpty() || a_f.box().contains(hiBox));

//         CH_assert(loBox.isEmpty() || a_g.box().contains(loBox));
//         CH_assert(a_g.box().contains(centerBox));
//         CH_assert(hiBox.isEmpty() || a_g.box().contains(hiBox));

//         FORT_FOURTHORDER_MULTIPLYAVERAGES_DIR(
//             CHF_FRA1(a_fg, a_fgComp),
//             CHF_CONST_FRA1(a_f, a_fComp),
//             CHF_CONST_FRA1(a_g, a_gComp),
//             CHF_CONST_INT(avgDir),
//             CHF_BOX(loBox),
//             CHF_BOX(centerBox),
//             CHF_BOX(hiBox)
//         );
//     }
// }


// -----------------------------------------------------------------------------
bool
FourthOrder::transverseSameType(const Box& a_box1,
                                const Box& a_box2,
                                const int  a_dir)
{
    const int t1 = (a_dir + 1) % SpaceDim;
    bool retval = a_box1.type(t1) == a_box2.type(t1);

    if constexpr (SpaceDim > 2) {
        const int t2 = (a_dir + 2) % SpaceDim;
        retval &= a_box1.type(t2) == a_box2.type(t2);
    }

    return retval;
}


// -----------------------------------------------------------------------------
void
FourthOrder::getDestBoxes(Box&       a_lolo,
                          Box&       a_lo,
                          Box&       a_mid,
                          Box&       a_hi,
                          Box&       a_hihi,
                          const Box& a_srcBox,
                          const Box& a_destBox,
                          const int  a_interpDir)
{
    CH_assert(FourthOrder::transverseSameType(a_destBox, a_srcBox, a_interpDir));
    CH_assert(FourthOrder::oppositeType(a_destBox, a_srcBox, a_interpDir));

    if (a_destBox.type(a_interpDir) == IndexType::NODE) {
        CH_verify(a_srcBox.size(a_interpDir) >= 4);

        a_lolo = a_destBox & bdryLo(a_srcBox, a_interpDir);
        a_lo   = a_destBox & Box(a_lolo).shift(a_interpDir, 1);
        a_mid  = a_destBox & surroundingNodes(a_srcBox, a_interpDir).grow(a_interpDir, -2);
        a_hihi = a_destBox & bdryHi(a_srcBox, a_interpDir);
        a_hi   = a_destBox & Box(a_hihi).shift(a_interpDir, -1);

    } else {
        CH_verify(a_srcBox.size(a_interpDir) >= 4);

        a_lolo = Box();
        a_lo   = a_destBox & adjCellLo(a_srcBox, a_interpDir, -1);
        a_mid  = a_destBox & enclosedCells(a_srcBox, a_interpDir).grow(a_interpDir, -1);
        a_hi   = a_destBox & adjCellHi(a_srcBox, a_interpDir, -1);
        a_hihi = Box();
    }
}
