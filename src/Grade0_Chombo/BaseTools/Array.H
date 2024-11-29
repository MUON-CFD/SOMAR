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
#ifndef ___Array_H__INCLUDED___
#define ___Array_H__INCLUDED___

#include <cstddef>
#include <stdexcept>
#include "MayDay.H"


/// Creates a 2D array using one block of contiguous memory.
/// I pilfered this code from
/// https://stackoverflow.com/questions/21943621/how-to-create-a-contiguous-2d-array-in-c
template <typename T>
T**
create2DArray(const size_t a_N,
              const size_t a_M);


/// De-allocates a 2D array made with create2DArray.
template <typename T>
void
delete2DArray(T** a_array);



/// Might as well have a 1D version.
template <typename T>
T*
create1DArray(const size_t a_N);


/// De-allocates a 2D array made with create1DArray.
template <typename T>
void
delete1DArray(T* a_array);


// =============================================================================
#include <array>
#include <memory>
#include <numeric>
#include "CH_assert.H"

/// @brief Allocates and maintains an N-dimensional, contiguous chunk of memory.
/// @tparam T      The underlying data type.
/// @tparam ND     The number of dimensions of this array.
/// @details
///  Performance:
///  I used 5 storage methods and compared thier access times (in a double loop)
///  in QuickBench. The storage methods were:
///   1. double** data = create2DArray<double>(Nx, Ny)
///   2. double* data = new double[Nx*Ny]
///   3. std::shared_ptr<double[]> data(new double[Nx*Ny])
///   4. std::vector<double> data(Nx*Ny)
///   5. MultiDimData<double, 2> data({Nx,Ny})
///  Of these, MultiDimData access times were significantly faster than options
///  1 and 2 and equal to options 3 and 4 when using GCC 12.2. MultiDimData was
///  the fastest when using GCC 8.5.
///
///  Memory management:
///  This class will automatically clean up after itself unless you define it
///  as a reference of a slice.
template <typename T, size_t ND>
class MultiDimData
{
public:
    using                   Value   = T;
    static constexpr size_t NumDims = ND;

    /// @name Construction / destruction
    /// @{

    /// @brief Create an empty object. User must call define() before use.
    MultiDimData()
    : m_dims{0}
    , m_dataPtr(nullptr)
    , m_dataPtrMM()
    {
    }

    /// @brief Creates an array of given size.
    /// @param a_dims
    /// @param a_ptr    If nullptr, this object will allocate memory and take
    /// care of memory managment. If not nullptr, this function will behave like
    /// a placement new and will NOT manage memory.
    MultiDimData(const std::array<size_t, NumDims>& a_dims,
                 Value*                             a_ptr = nullptr)
    : m_dims{0}
    , m_dataPtr(nullptr)
    , m_dataPtrMM()
    {
        this->define(a_dims, a_ptr);
    }

    /// @brief Creates an array of given size.
    /// @param a_dims
    /// @param a_ptr    If nullptr, this object will allocate memory and take
    /// care of memory managment. If not nullptr, this function will behave like
    /// a placement new and will NOT manage memory.
    void
    define(const std::array<size_t, NumDims>& a_dims,
           Value*                             a_ptr = nullptr)
    {
        const size_t numPts = std::accumulate(a_dims.data(),
                                              a_dims.data() + NumDims,
                                              1.0,
                                              std::multiplies<size_t>());
        if (numPts == 0) {
            MayDay::Error("MultiDimData received a dimension of zero.");
        }

        if (a_ptr) {
            m_dims = a_dims;
            m_dataPtr = a_ptr;
            m_dataPtrMM.reset();
        } else {
            m_dims = a_dims;
            m_dataPtr = new Value[numPts];
            m_dataPtrMM.reset(m_dataPtr);
        }
    }

    /// @brief References a slice of another MultiDimData.
    ///
    /// @param a_src
    /// @param a_sliceRange  The beggining and ending indexes of the slice.
    ///  For example, if a_sliceRange = {2, 5}, then
    ///  (*this)(i, j, k, 0) = a_src(i, j, k, 2)
    ///  (*this)(i, j, k, 1) = a_src(i, j, k, 3)
    ///  (*this)(i, j, k, 2) = a_src(i, j, k, 4)
    ///  This follows iterator rules. The end is actually 1 past the end!
    ///
    /// @details Memory of the original data will persist at least as long as
    ///  this slice persists.
    ///
    ///  Slicing can only occur on the slowest (last) index.
    MultiDimData(MultiDimData<Value, NumDims>& a_src,
                 const std::array<size_t, 2>&  a_sliceRange)
    : m_dims{0}
    , m_dataPtr(nullptr)
    , m_dataPtrMM()
    {
        this->define(a_src, a_sliceRange);
    }

    /// @brief References a slice of another MultiDimData.
    ///
    /// @param a_src
    /// @param a_sliceRange  The beggining and ending indexes of the slice.
    ///  For example, if a_sliceRange = {2, 5}, then
    ///  (*this)(i, j, k, 0) = a_src(i, j, k, 2)
    ///  (*this)(i, j, k, 1) = a_src(i, j, k, 3)
    ///  (*this)(i, j, k, 2) = a_src(i, j, k, 4)
    ///  This follows iterator rules. The end is actually 1 past the end!
    ///
    /// @details Memory of the original data will persist at least as long as
    ///  this slice persists.
    ///
    ///  Slicing can only occur on the slowest (last) index.
    void
    define(MultiDimData<Value, NumDims>& a_src,
           const std::array<size_t, 2>&  a_sliceRange)
    {
        if (a_sliceRange[0] >= a_src.size(NumDims - 1) ||
            a_sliceRange[1] >  a_src.size(NumDims - 1) ||
            a_sliceRange[0] > a_sliceRange[1]) {
            MayDay::Error(
                "MultiDimData::define received an invalid sliceRange.");
        }

        m_dims = a_src.m_dims;
        m_dims[NumDims - 1] = a_sliceRange[1] - a_sliceRange[0] + 1;

        const size_t jump = std::accumulate(m_dims.data(),
                                            m_dims.data() + NumDims - 1,
                                            1.0,
                                            std::multiplies<size_t>());
        m_dataPtr = a_src.m_dataPtr + jump * a_sliceRange[0];

        m_dataPtrMM = a_src.m_dataPtrMM;
    }

    /// @brief Free memory. Leaves object unusable.
    void
    clear()
    {
        m_dataPtr = nullptr;
        m_dataPtrMM.reset();
    }
    /// @}


    /// @name Basic accessors
    /// @{

    template <typename... Idxs>
    Value&
    operator()(const Idxs... a_idxs)
    {
        static_assert(sizeof...(a_idxs) == NumDims);
        return m_dataPtr[this->linearIndex(a_idxs...)];
    }

    template <typename... Idxs>
    const Value&
    operator()(const Idxs... a_idxs) const
    {
        static_assert(sizeof...(a_idxs) == NumDims);
        return m_dataPtr[this->linearIndex(a_idxs...)];
    }
    /// @}

    /// @name Inquiries
    /// @{

    const std::array<size_t, NumDims>&
    size() const noexcept
    {
        return m_dims;
    }

    size_t
    size(const size_t a_axis) const
    {
        CH_assert(a_axis < NumDims);
        return m_dims[a_axis];
    }

    /// @brief Returns the product of the sizes in the specified dimensions.
    /// @param a_begin The starting index.
    /// @param a_end   The end index. Uses iterator rules. For example,
    ///                sizeProduct(2, 4) produces m_dims[2] * m_dims[3].
    size_t
    sizeProduct(const size_t a_begin = 0, const size_t a_end = NumDims) const {
        return std::accumulate(m_dims.data() + a_begin,
                               m_dims.data() + a_end,
                               1.0,
                               std::multiplies<size_t>());
    }
    /// @}

    template <typename... Ts>
    size_t
    linearIndex(const Ts... a_ts) const
    {
        return _linearIndex(0, a_ts...);
    }

    Value*
    dataPtr(const size_t a_linearOffset = 0) noexcept
    {
        return m_dataPtr + a_linearOffset;
    }

    /// \brief Linearization
    /// \{
    int
    linearSize() const
    {
        return sizeof(m_dims) + this->sizeProduct() * sizeof(Value);
    }

    void
    linearOut(void* const a_outBuf) const
    {
        char* charBuf = static_cast<char*>(a_outBuf);
        ::linearOut(charBuf, m_dims);
        charBuf += ::linearSize(m_dims);
        std::memcpy(charBuf, m_dataPtr, this->sizeProduct() * sizeof(Value));
    }

    void
    linearIn(const void* const a_inBuf)
    {
        const char* charBuf = static_cast<const char*>(a_inBuf);

        decltype(m_dims) incomingDims;
        ::linearIn(incomingDims, charBuf);
        if(incomingDims != m_dims) MayDay::Error("MultiDimData has incorrect dimensions.");
        charBuf += ::linearSize(m_dims);

        std::memcpy(m_dataPtr, charBuf, this->sizeProduct() * sizeof(Value));
    }
    /// \}

private:
    size_t
    _linearIndex(const size_t d, const size_t i) const
    {
        CH_assert(i < m_dims[d]);
        return i;
    }

    template <typename... Ts>
    size_t
    _linearIndex(const size_t d, const size_t i, const Ts... a_ts) const
    {
        CH_assert(i < m_dims[d]);
        return i + m_dims[d] * _linearIndex(d + 1, a_ts...);
    }

    std::array<size_t, NumDims> m_dims;
    Value*                      m_dataPtr;
    std::shared_ptr<Value[]>    m_dataPtrMM;
};



#define M4b858ba808f4a6a391f8ea74ce227cdf
#include "ArrayI.H"
#undef M4b858ba808f4a6a391f8ea74ce227cdf

#endif //!___Array_H__INCLUDED___
