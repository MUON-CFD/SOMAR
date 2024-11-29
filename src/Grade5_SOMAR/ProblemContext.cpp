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
#include "ProblemContext.H"
#include "Debug.H"


std::unique_ptr<ProblemContext> ProblemContext::s_ctx;


// -----------------------------------------------------------------------------
// Private constructor
// -----------------------------------------------------------------------------
ProblemContext::ProblemContext ()
: base()
, time()
, output()
, amr(base)
, rhs(base, amr.numLevels)
, proj()
{
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
ProblemContext::~ProblemContext ()
{
}


// -----------------------------------------------------------------------------
// This returns the single ProblemContext object. (const version)
// -----------------------------------------------------------------------------
const ProblemContext* ProblemContext::getInstance ()
{
    if (!s_ctx) {
        s_ctx.reset(new ProblemContext());
    }
    return &*s_ctx;
}


// -----------------------------------------------------------------------------
// This returns the single ProblemContext object. (non-const version)
// -----------------------------------------------------------------------------
ProblemContext*
ProblemContext::getNonConstInstance ()
{
    if (!s_ctx) {
        s_ctx.reset(new ProblemContext());
    }
    return &*s_ctx;
}


// -----------------------------------------------------------------------------
// This deletes the static pointer in hopes valgrind will take notice.
// -----------------------------------------------------------------------------
void ProblemContext::freeMemory ()
{
    s_ctx.reset();
    BUG("Are all *Parameters getting freed?");
}

