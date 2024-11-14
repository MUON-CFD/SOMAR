/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2018
 *    Jefferson (Philadelphia University + Thomas Jefferson University) and
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
 *  https://github.com/somarhub.
 ******************************************************************************/
#include "AMRNSLevelFactory.H"
#include "AMRNSLevel.H"


class UserMain
{
public:
    static AMRNSLevel*
    createPhysics();

    static GeoSourceInterface*
    oneTimeSetup();
};


// -----------------------------------------------------------------------------
// Constructor.
// This is where you set things once in the beginning of the simulation.
// -----------------------------------------------------------------------------
AMRNSLevelFactory::AMRNSLevelFactory()
: m_geoSrcPtr(nullptr)
{
    m_geoSrcPtr = UserMain::oneTimeSetup();
}


// -----------------------------------------------------------------------------
// Destructor
// This is where you will perform whatever cleanup is needed.
// -----------------------------------------------------------------------------
AMRNSLevelFactory::~AMRNSLevelFactory()
{
    delete m_geoSrcPtr;
    m_geoSrcPtr = nullptr;
}


// -----------------------------------------------------------------------------
// The level factory.
// This is where you create new user-defined physics classes that will be
// used as a level solver.
// -----------------------------------------------------------------------------
AnisotropicAMRLevel*
AMRNSLevelFactory::new_amrlevel() const
{
    // Let the user create their own new level.
    AMRNSLevel* newPhysPtr = UserMain::createPhysics();

    // Set some basics for the user.
    newPhysPtr->m_geoSrcPtr = m_geoSrcPtr;
    newPhysPtr->m_isFactoryDefined = true;

    // Return a slice of the new object.
    return (static_cast<AnisotropicAMRLevel*>(newPhysPtr));
}
