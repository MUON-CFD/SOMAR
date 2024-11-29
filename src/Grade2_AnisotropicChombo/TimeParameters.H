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
#ifndef ___TimeParameters_H__INCLUDED___
#define ___TimeParameters_H__INCLUDED___

#include <array>
#include <memory>
#include <string>
#include "REAL.H"
#include "RealVect.H"


class TimeParameters
{
public:
    // The constructor sets default / reads parameters.
    TimeParameters();

    Real        stopTime;
    int         maxSteps;
    Real        fixedDt;
    Real        maxDt;
    Real        initDtMult;
    Real        dtMult;
    Real        maxDtGrow;

    bool        useElementaryController;
    bool        usePIController;
    bool        usePIDController;
    Real        absTol;
    Real        relTol;

    bool        isRestart;
    std::string restartFile;

    // You shouldn't need to call this. AnisotropicAMR will do it for you.
    static void
    freeMemory();

    // It's nice to be able to see these parameters in pout.*.
    void
    dump() const;

private:
    static void
    createDefaults();

    static std::unique_ptr<TimeParameters> s_defPtr;
    static bool                            s_constructorLock;
};


#endif //!___TimeParameters_H__INCLUDED___
