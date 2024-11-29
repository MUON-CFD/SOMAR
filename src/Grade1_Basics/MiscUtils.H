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

// This file contains random tools that have no home.
#ifndef __MiscUtils_H__INCLUDED__
#define __MiscUtils_H__INCLUDED__

#include "RefCountedPtr.H"
#include "Vector.H"
#include "BoxLayoutData.H"
#include "FluxBox.H"
#include "NodeFArrayBox.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "IntVect.H"
#include "SOMAR_Constants.H"
#include <cfenv>
#include <cmath>


/**
 * @brief Case-insensitive string comparison.
 *
 * @param a
 * @param b
 * @return true if a and b are equal.
 * @return false if a and b are not equal.
 */
bool compare_insensitive(const std::string& a, const std::string& b);


/**
 * @brief      Takes a Vector of RefCountedPtrs and converts it to a Vector of C
 * pointers.
 *
 * @warning    I claim no responsibility for irresponsible use.
 */
template <class T>
Vector<T*>
stripRCPs(Vector<RefCountedPtr<T> >& a_src)
{
    const int size = a_src.size();
    Vector<T*> dest(size, NULL);
    for (int idx = 0; idx < size; ++idx) {
        dest[idx] = &(*(a_src[idx]));
    }
    return dest;
}

/**
 * @brief      Converts Vector<T*> to Vector<const T*>.
 */
template <class T>
Vector<const T*>
makeConstVector(Vector<T*>& a_in)
{
    Vector<const T*> out(a_in.size());
    for (size_t idx = 0; idx < a_in.size(); ++idx) {
        out[idx] = a_in[idx];
    }
    return out;
}


/**
 * @brief Scales a vector so that it has a sum of 1.0.
 */
template <typename T>
void
normalizeSum(std::vector<T>& a_vec)
{
    const T sum = std::accumulate(a_vec.begin(), a_vec.end(), 0.0);
    for (auto& vi : a_vec) {
        vi /= sum;
    }
}


/**
 * @brief Removes a set of elements from a std::vector.
 *
 * @tparam T     The element type in std::vector<T>.
 * @param a_vec  The vector whose elements will be removed.
 * @param a_list The indices of the elements we want to remove.
 *
 * @details For example, suppose a_vec = {'a', 'b', 'c', 'd', 'e'}
 *          and a_list = {1, 4}. This function will remove a_vec[1] = 'b'
 *          and a_vec[4] = '3', leaving a_vec = {'a', 'c', 'd'}.
 */
template <typename T>
void
eraseElements(std::vector<T>& a_vec, const std::vector<size_t>& a_list)
{
    if (a_list.size() == 0) return;

    std::vector<T> newVec;
    newVec.reserve(a_vec.size());

    for (size_t next = 0, e = 0; e < a_vec.size(); ++e) {
        if (e == a_list[next]) {
            if (next < a_list.size() - 1) {
                ++next;
            } else {
                next = 0;
            }
        } else {
            newVec.push_back(a_vec[e]);
        }
    }

    std::swap(a_vec, newVec);
}

/**
 * @brief Creates a plan to sort vec. To use the plan, call applyPermutation.
 * @details Shamelessly borrowed from
 * https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 *
 */
template <typename T, typename Compare>
std::vector<std::size_t>
sortPermutation(const std::vector<T>& vec, Compare comp)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(
        p.begin(),
        p.end(),
        [&](const std::size_t i, const std::size_t j) {
            return comp(vec[i], vec[j]);
        });
    return p;
}

/// This version uses a trivial Compare function.
template <typename T>
std::vector<std::size_t>
sortPermutation(const std::vector<T>& vec)
{
    auto comp = [](const T a, const T b) { return a < b; };
    return sortPermutation(vec, comp);
}


/**
 * @brief Sorts a vector using the provided permutation.
 * @details Shamelessly borrowed from
 * https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 *
 */
template <typename T>
void
applyPermutation(std::vector<T>& vec, const std::vector<std::size_t>& p)
{
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (done[i]) continue;

        done[i]            = true;
        std::size_t prev_j = i;
        std::size_t j      = p[i];

        while (i != j) {
            std::swap(vec[prev_j], vec[j]);

            done[j] = true;
            prev_j  = j;
            j       = p[j];
        }
    }
}


/**
 * @brief      DataFactory for creating FABs with arbitrary centering.
 */
template<class T>
class NCDataFactory: public DataFactory<T>
{
public:
    NCDataFactory(const IntVect& a_type)
    {
        for (int dir = 0; dir < SpaceDim; ++dir)
            m_type[dir] = (a_type[dir] == 0) ? 0 : 1;
    }

    virtual ~NCDataFactory()
    {}

    // factory function.  creates a new 'T' object
    // creates a new 'T' object and returns a pointer to it.
    virtual T*
    create(const Box& a_box, int a_ncomps, const DataIndex& /*a_datInd*/) const
    {
        Box fcBox = a_box;
        fcBox.convert(m_type);
        return new T(fcBox, a_ncomps);
    }

protected:
    IntVect m_type;
};


template<>
class NCDataFactory<FArrayBox>: public DataFactory<FArrayBox>
{
public:
    NCDataFactory(const IntVect& a_type)
    {
        for (int dir = 0; dir < SpaceDim; ++dir)
            m_type[dir] = (a_type[dir] == 0) ? 0 : 1;
    }

    virtual ~NCDataFactory()
    {}

    // factory function.  creates a new 'T' object
    // creates a new 'T' object and returns a pointer to it.
    virtual FArrayBox*
    create(const Box& a_box, int a_ncomps, const DataIndex& /*a_datInd*/) const
    {
        Box fcBox = a_box;
        fcBox.convert(m_type);
        return new FArrayBox(fcBox, a_ncomps);
    }

protected:
    IntVect m_type;
};



/**
 * @brief      Rounds up to the nearest multiple.
 *
 * @param[in]  numToRound  The number to round
 * @param[in]  factor      Result will be a multiple of this number.
 *
 * @details
 *  For example, roundUp(4, 3) = 6
 *  This code was pilfered from
 *  https://stackoverflow.com/questions/3407012/c-rounding-up-to-the-nearest-multiple-of-a-number
 */
int
roundUp(int numToRound, int factor);


/**
 * @brief Rounds each component up to the nearest integer.
 */
constexpr IntVect
roundUpToIntVect(const RealVect& a_rv)
{
    return IntVect{ D_DECL(
        static_cast<int>(std::ceil(a_rv.getConst(0) + smallReal)),
        static_cast<int>(std::ceil(a_rv.getConst(1) + smallReal)),
        static_cast<int>(std::ceil(a_rv.getConst(2) + smallReal))) };
}


/**
 * @brief      This produces a_n random numbers on the interval 0.0 to 1.0.
 *
 * @param[out] a_randos  The Vector of random numbers
 *
 * @details
 *  All MPI ranks will receive the same set of numbers.
 *  This function is blocking.
 */
void
generateRandomNumbers(Vector<Real>& a_randos);


#endif //!__MiscUtils_H__INCLUDED__
