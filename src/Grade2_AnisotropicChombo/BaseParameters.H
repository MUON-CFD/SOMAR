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
#ifndef ___BaseParameters_H__INCLUDED___
#define ___BaseParameters_H__INCLUDED___

#include <memory>
#include "RealVect.H"
#include "ProblemDomain.H"


class BaseParameters
{
public:
    // Constructor
    BaseParameters();

    RealVect      L;
    IntVect       nx;
    IntVect       nxOffset;
    bool          isPeriodic[CH_SPACEDIM];
    ProblemDomain domain;
    IntVect       splitDirs;

    IntVect maxBaseGridSize;
    int     blockFactor;

    // You shouldn't need to call this. AnisotropicAMR will do it for you.
    static void
    freeMemory();

    // It's nice to be able to see these parameters in pout.*.
    void
    dump() const;


private:
    static void
    createDefaults();

    //
    static void
    chooseMaxBaseGridSize(IntVect&             a_maxBaseGridSize,
                          const ProblemDomain& a_domain,
                          const RealVect&      a_L);

    //
    static void
    powOf2Strategy(IntVect&             a_maxGridSize,
                   int&                 a_idealNumProc,
                   const ProblemDomain& a_domain,
                   const int            a_nproc);

    static std::unique_ptr<BaseParameters> s_defPtr;
    static bool                            s_constructorLock;
};


#endif //!___BaseParameters_H__INCLUDED___
