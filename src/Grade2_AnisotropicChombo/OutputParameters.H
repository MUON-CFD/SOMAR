/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2022
 *    Thomas Jefferson University and
 *    University of North Carolina at Chapel Hill
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
 *  https://github.com/UNC-CFD/LES-SOMAR.
 ******************************************************************************/
#ifndef ___OutputParameters_H__INCLUDED___
#define ___OutputParameters_H__INCLUDED___

#include <memory>
#include <string>
#include "REAL.H"


class OutputParameters
{
public:
    // Constructor
    OutputParameters();

    int         verbosity;
    bool        doFlowchart;

    int         plotInterval;
    Real        plotPeriod;
    std::string plotPrefix;

    int         checkpointInterval;
    std::string checkpointPrefix;

    // You shouldn't need to call this. AnisotropicAMR will do it for you.
    static void
    freeMemory();

    // It's nice to be able to see these parameters in pout.*.
    void
    dump() const;

private:
    static void
    createDefaults();

    static std::unique_ptr<OutputParameters> s_defPtr;
    static bool                              s_constructorLock;
};


#endif //!___OutputParameters_H__INCLUDED___
