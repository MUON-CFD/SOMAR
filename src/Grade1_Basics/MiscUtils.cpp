/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
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
#include "MiscUtils.H"


// -----------------------------------------------------------------------------
bool
compare_insensitive(const std::string& a, const std::string& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (std::tolower(a[i]) != std::tolower(b[i])) return false;
    return true;
}


// -----------------------------------------------------------------------------
// Rounds up to the nearest multiple.
// So, roundUp(4, 3) = 6
// This code was pilfered from
// https://stackoverflow.com/questions/3407012/c-rounding-up-to-the-nearest-multiple-of-a-number
// -----------------------------------------------------------------------------
int
roundUp(int numToRound, int factor)
{
    if (factor == 0)
        return numToRound;

    int remainder = abs(numToRound) % factor;
    if (remainder == 0)
        return numToRound;

    if (numToRound < 0)
        return -(abs(numToRound) - remainder);
    else
        return numToRound + factor - remainder;
}


// -----------------------------------------------------------------------------
// This produces a_n random numbers on the interval 0.0 to 1.0.
// All MPI ranks will receive the same set of numbers.
// This function is blocking.
// -----------------------------------------------------------------------------
void generateRandomNumbers (Vector<Real>& a_randos)
{
    static const int myRank = procID();
    const int numRandos = a_randos.size();

    // Generate the random numbers on a single rank.
    if (myRank == 0) {
        for (int idx = 0; idx < numRandos; ++idx) {
            a_randos[idx] = Real(rand()) / Real(RAND_MAX);
        }
    }

#ifdef CH_MPI
    // Broadcast the numbers.
    MPI_Bcast(&a_randos[0], numRandos, MPI_CH_REAL, 0, Chombo_MPI::comm);
#endif
}
