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
 *  https://github.com/MUON-CFD/SOMAR.
 ******************************************************************************/
#ifndef ___AMRNSLevelFactory_H__INCLUDED___
#define ___AMRNSLevelFactory_H__INCLUDED___
#if defined(CH_HAVE_THRUST) || defined(__NVCC__)
#include "CudaPoissonOp.hpp"
#endif
#include "AnisotropicAMRLevelFactory.H"

class GeoSourceInterface;


class AMRNSLevelFactory: public AnisotropicAMRLevelFactory
{
public:
    // Constructor.
    // This is where you set things once in the beginning of the simulation.
    AMRNSLevelFactory();

    // Destructor
    // This is where you will perform whatever cleanup is needed.
    virtual ~AMRNSLevelFactory();

    // The level factory.
    // This is where you create new user-defined physics classes that will be
    // used as a level solver.
    virtual AnisotropicAMRLevel*
    new_amrlevel() const;

protected:
    GeoSourceInterface*  m_geoSrcPtr;
};


#endif  //!___AMRNSLevelFactory_H__INCLUDED___
