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
 *  https://github.com/somarhub.
 ******************************************************************************/
#ifndef ___AMRParameters_H__INCLUDED___
#define ___AMRParameters_H__INCLUDED___

#include <memory>
#include <vector>
#include "BaseParameters.H"


class AMRParameters
{
public:
    // The constructor sets default / reads parameters.
    AMRParameters(const BaseParameters& a_baseParams);

    int                  maxLevel;
    int                  numLevels;
    std::vector<IntVect> refRatios;
    IntVect              maxGridSize;
    int                  bufferSize;
    Real                 fillRatio;

    std::vector<int> regridIntervals;
    bool             useSubcycling;

    bool              tagIB;
    Real              velTagTol;
    Real              bTagTol;
    Real              TTagTol;
    Real              STagTol;
    Real              bpertTagTol;
    Real              TpertTagTol;
    Real              SpertTagTol;
    std::vector<Real> scalarsTagTol;
    int               growTags;

    // You shouldn't need to call this. AnisotropicAMR will do it for you.
    static void
    freeMemory();

    // It's nice to be able to see these parameters in pout.*.
    void
    dump() const;

private:
    static void
    createDefaults(const BaseParameters& a_baseParams);

    //
    static void
    chooseMaxGridSize(IntVect&                    a_maxGridSize,
                      const std::vector<IntVect>& a_refRatio,
                      const BaseParameters&       a_baseParams);

    //
    static void
    chooseMaxGridSize(IntVect&              a_maxGridSize,
                      const IntVect&        a_refRatio,
                      const BaseParameters& a_baseParams);

    static std::unique_ptr<AMRParameters> s_defPtr;
    static bool                           s_constructorLock;
};


#endif  //!___AMRParameters_H__INCLUDED___
