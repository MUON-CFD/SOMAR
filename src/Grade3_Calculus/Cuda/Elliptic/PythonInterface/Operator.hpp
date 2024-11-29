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

//#define GlobalOps
#ifndef __OPERATOR_H__
#define __OPERATOR_H__
// #ifdef NDEBUG
// #undef NDEBUG
// #endif
#include <assert.h>
#include <cusp/array1d.h>
#include <cusp/blas/blas.h>
#include <cusp/csr_matrix.h>
#include <cusp/elementwise.h>
#include <cusp/functional.h>
#include <cusp/linear_operator.h>
#include <cusp/multiply.h>
#include <cusp/print.h>
#include <cusp/transpose.h>

#ifdef USE_CUSPARSE
#include "cusparse_v2.h"
#endif
#include <cuda_runtime.h>
//#include <chrono>

#include <typeinfo>
//#include <cusp/precond/aggregation/smoothed_aggregation.h>
//#include <cusp/precond/ainv.h>
#include <map>
#include <string>

#include "CudaFArrayBox.hpp"
#include "IO.H"

#include "icecream.hpp"
namespace Elliptic {
namespace CudaInterface {
template <typename IntType, typename ValueType, typename InputsMemory,
          typename OwnMemory>
struct Operator {
  typedef typename cusp::csr_matrix<IntType, ValueType, OwnMemory> Matrix;
  typedef typename CudaBaseFab<ValueType, InputsMemory>::Array InputArr;
  typedef typename CudaBaseFab<ValueType, InputsMemory>::View InputView;
  typedef typename cusp::array1d<ValueType, OwnMemory> MyArr;
  typedef BaseFab<ValueType> BaseFAB;
  typedef CudaFArrayBox<InputsMemory> DiagType;
  // typedef typename CudaBaseFab<ValueType, cusp::host_memory> GeoArray; //
  // they reside in CPU
  typedef typename std::vector<IntType> PyIntArray;
  typedef typename std::vector<ValueType> PyRealArray;

#ifndef USE_CUSPARSE
  typedef int cudaStream_t;
#endif
  // data members
 private:
  mutable Matrix m_Matrix;
  int m_Nrows;
  int m_Ncolumns;
  int m_nnz;
  Box m_gridInput;
  Box m_gridOutput;
  bool m_isMatrixDefined;
  bool m_areGridsDefined;
#ifdef USE_CUSPARSE
  cusparseHandle_t m_handle = 0;
  cusparseMatDescr_t m_descr = 0;
  cudaError_t m_cudaStat;
  cusparseStatus_t m_status;
#endif
public:
  ~Operator() {
#ifdef USE_CUSPARSE
    m_status = cusparseDestroy(m_handle);
#endif
  }
  Operator() {}

  Operator(const int Nrows, const int Ncolumns, const int nnz,
           cudaStream_t StreamId = 0)
      : m_Matrix(Nrows, Ncolumns, nnz),
        m_Nrows(Nrows),
        m_Ncolumns(Ncolumns),
        m_nnz(nnz) {
#ifdef USE_CUSPARSE
    m_status = cusparseCreate(&m_handle);
    m_status = cusparseSetStream(m_handle, StreamId);
    /* create and setup matrix descriptor */
    m_status = cusparseCreateMatDescr(&m_descr);
    cusparseSetMatType(m_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(m_descr, CUSPARSE_INDEX_BASE_ZERO);
#endif
  }

  // void define(const int Nrows, const int Ncolumns, const int nnz)
  // {
  //     m_Nrows = Nrows;
  //     m_Ncolumns = Ncolumns;
  //     m_nnz = nnz;
  //     m_Matrix(Nrows, Ncolumns, nnz);
  //     /* initialize cusparse library */
  //     m_status = cusparseCreate(&m_handle);
  //     /* create and setup matrix descriptor */
  //     m_status = cusparseCreateMatDescr(&m_descr);
  //     cusparseSetMatType(m_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  //     cusparseSetMatIndexBase(m_descr, CUSPARSE_INDEX_BASE_ZERO);
  // }
  void setOperator(const PyIntArray &IndPtr, const PyIntArray &Col,
                   const PyRealArray &Data) {
    assert(IndPtr.size() == m_Nrows + 1);
    assert(Col.size() == m_nnz);
    assert(Data.size() == m_nnz);

    thrust::copy(IndPtr.begin(), IndPtr.end(), m_Matrix.row_offsets.begin());

    thrust::copy(Col.begin(), Col.end(), m_Matrix.column_indices.begin());
    thrust::copy(Data.begin(), Data.end(), m_Matrix.values.begin());
    m_isMatrixDefined = true;
  }

  bool isDefined() const { return m_isMatrixDefined && m_areGridsDefined; }

  void setGrids(const Box &gridInput, const Box &gridOutput) {
    m_gridInput = gridInput;
    m_gridOutput = gridOutput;
    m_areGridsDefined = true;
  }

  void fillTheseGrids(Box &a_gridInput, Box &a_gridOutput) const {
    a_gridInput = m_gridInput;
    a_gridOutput = m_gridOutput;
  }
  auto &getMyInputGrid() const { return m_gridInput; }
  auto &getMyOutputGrid() const { return m_gridOutput; }
  auto getMatrixPtr() const { return &m_Matrix; }

  void operator()(const DiagType &a_x, DiagType &a_y, Interval I) {
    if (!a_x.box().contains(this->m_gridInput)) {
      IO::tout(0) << " x.box() " << a_x.box() << "; grid Input " << m_gridInput
                  << endl;

      MayDay::Error("Operator needs larger input array");
    }
    if (!a_y.box().contains(this->m_gridOutput)) {
      IO::tout(0) << " y.box() " << a_y.box() << "; grid Output "
                  << m_gridOutput << endl;
      MayDay::Error("Operator needs larger output array");
    }
    if (a_y.box() != this->m_gridOutput) {
      DiagType Y(m_gridOutput, a_y.nComp());
      Y.copy(a_y);
      this->operator()(a_x, Y, I);

      a_y.copy(Y);
      return;
    }
    if (a_x.box() != this->m_gridInput) {
      DiagType X(m_gridInput, a_x.nComp());
      X.copy(a_x);
      this->operator()(X, a_y, I);
      return;
    }

    bool useViews = !std::is_same<OwnMemory, InputsMemory>::value ||
                    a_x.isAliased() || a_y.isAliased();
    useViews ? this->operator()(a_x.m_view, a_y.m_view, I)
             : this->operator()(a_x.m_varr, a_y.m_varr, I);
  }


 public:
  template <
      typename C_, typename M_ = OwnMemory,
      typename std::enable_if<std::is_same<M_, cusp::device_memory>::value &&
                                  std::is_same<Real, double>::value,
                              int>::type * = nullptr>
  void operator()(const C_ &a_x, C_ &a_y, Real alpha = 1.0, Real beta = 0.0) {
    assert(this->isDefined());
    assert(a_x.size() == m_gridInput.numPts());
    assert(a_y.size() == m_gridOutput.numPts());
#ifdef USE_CUSPARSE
    cusparseStatus_t status;
    // apply

    status = cusparseDcsrmv(
        m_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_Matrix.num_rows,
        m_Matrix.num_cols, m_Matrix.num_entries, &alpha, m_descr,
        (&m_Matrix.values[0]).get(), (&m_Matrix.row_offsets[0]).get(),
        (&m_Matrix.column_indices[0]).get(), (&a_x[0]).get(), &beta,
        (&a_y[0]).get());
    checkErrorStatus(status);
#else
    assert(m_Matrix.num_cols = a_x.size());
    assert(m_Matrix.num_rows = a_y.size());
    cusp::multiply(m_Matrix, a_x, a_y);
#endif
  }

  template <
      typename C_, typename M_ = OwnMemory,
      typename std::enable_if<std::is_same<M_, cusp::device_memory>::value &&
                                  std::is_same<Real, float>::value,
                              int>::type * = nullptr>
  void operator()(const C_ &a_x, C_ &a_y, Real alpha = 1.0, Real beta = 0.0) {
    assert(this->isDefined());
    assert(a_x.size() == m_gridInput.numPts());
    assert(a_y.size() == m_gridOutput.numPts());
//
#ifdef USE_CUSPARSE
    cusparseStatus_t status;
    // apply

    status = cusparseScsrmv(
        m_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_Matrix.num_rows,
        m_Matrix.num_cols, m_Matrix.num_entries, &alpha, m_descr,
        (&m_Matrix.values[0]).get(), (&m_Matrix.row_offsets[0]).get(),
        (&m_Matrix.column_indices[0]).get(), (&a_x[0]).get(), &beta,
        (&a_y[0]).get());
    checkErrorStatus(status);
#else
    assert(m_Matrix.num_cols = a_x.size());
    assert(m_Matrix.num_rows = a_y.size());
    cusp::multiply(m_Matrix, a_x, a_y);
#endif
  }

  template <typename C_, typename M_ = OwnMemory,
            typename std::enable_if<std::is_same<M_, cusp::host_memory>::value,
                                    int>::type * = nullptr>
  auto operator()(const C_ &a_x, C_ &a_y, Real alpha = 1.0, Real beta = 1.0) {
    assert(this->isDefined());
    assert(a_x.size() == m_gridInput.numPts());
    assert(a_y.size() == m_gridOutput.numPts());
    assert(m_Matrix.num_cols = a_x.size());
    assert(m_Matrix.num_rows = a_y.size());
    cusp::multiply(m_Matrix, a_x, a_y);
  }

  void operator_with_copy(const InputArr &a_x, InputArr &a_y) {
    assert(this->isDefined());
    assert(a_x.size() == m_gridInput.numPts());
    assert(a_y.size() == m_gridOutput.numPts());
    MyArr x(a_x.begin(), a_x.end());
    MyArr y(a_y.size());
    assert(m_Matrix.num_cols = x.size());
    assert(m_Matrix.num_rows = y.size());
    this->operator()(x, y);
    thrust::copy(y.begin(), y.end(), a_y.begin());
  }
  private:
  void operator()(const InputArr &a_x, InputArr &a_y, Interval I) {
    if (I.size() == 1 && I.begin() == 0)  // saves a copy
    {
      if (std::is_same<MyArr, InputArr>::value) {
        this->operator()(a_x, a_y);
      } else {
        this->operator_with_copy(a_x, a_y);
      }
    } else {
      auto XSize = a_x.size() / I.size();
      auto YSize = a_y.size() / I.size();
      for (auto i = I.begin(); i < I.end() + 1; ++i) {
        MyArr x(a_x.subarray(i * XSize, XSize));
        MyArr y(YSize);
        this->operator()(x, y);
        thrust::copy(y.begin(), y.end(), a_y.begin() + i * YSize);
      }
    }
  }
  void operator()(const InputView &a_x, InputView &a_y, Interval I) {
    auto XSize = a_x.size() / I.size();
    auto YSize = a_y.size() / I.size();

    for (auto i = I.begin(); i < I.end() + 1; ++i) {
      MyArr x(a_x.begin() + i * XSize, a_x.begin() + (i + 1) * XSize);
      MyArr y(YSize);
      this->operator()(x, y);
      thrust::copy(y.begin(), y.end(), a_y.begin() + i * YSize);
    }
  }
#ifdef USE_CUSPARSE
  void checkErrorStatus(cusparseStatus_t status) {
    switch (status) {
      case CUSPARSE_STATUS_NOT_INITIALIZED:
        MayDay::Error("CUSPARSE_STATUS_NOT_INITIALIZED ");
        break;
      case CUSPARSE_STATUS_ALLOC_FAILED:
        MayDay::Error("CUSPARSE_STATUS_ALLOC_FAILED ");
        break;
      case CUSPARSE_STATUS_INVALID_VALUE:
        MayDay::Error("CUSPARSE_STATUS_INVALID_VALUE");
        break;
      case CUSPARSE_STATUS_ARCH_MISMATCH:
        MayDay::Error("CUSPARSE_STATUS_ARCH_MISMATCH");
        break;
      case CUSPARSE_STATUS_MAPPING_ERROR:
        MayDay::Error("CUSPARSE_STATUS_MAPPING_ERROR");
        break;
      case CUSPARSE_STATUS_EXECUTION_FAILED:
        MayDay::Error("CUSPARSE_STATUS_EXECUTION_FAILED");
        break;
      case CUSPARSE_STATUS_INTERNAL_ERROR:
        MayDay::Error("CUSPARSE_STATUS_INTERNAL_ERROR");
        break;
      case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
        MayDay::Error("CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED");
        break;
      case CUSPARSE_STATUS_ZERO_PIVOT:
        MayDay::Error("CUSPARSE_STATUS_ZERO_PIVOT");
        break;

      default:
        break;
    }
  }
#endif
};
}  // namespace CudaInterface
}  // namespace Elliptic
#endif  // __OPERATOR_H__
