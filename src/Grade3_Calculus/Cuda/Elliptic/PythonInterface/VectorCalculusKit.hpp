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
#pragma once
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
//#include <cusp/gallery/poisson.h>
// include cusp relaxators header file
#include <cusp/relaxation/gauss_seidel.h>
#include <cusp/relaxation/jacobi.h>
#include <cusp/relaxation/sor.h>
//#include <cusp/krylov/cg.h>
//#include <cusp/krylov/bicgstab.h>
//#include <cusp/krylov/bicg.h>
#ifdef USE_CUSPARSE
#include "cusparse_v2.h"
#endif
#include <cuda_runtime.h>
//#include <chrono>
#include <cusp/precond/diagonal.h>
#include <typeinfo>
//#include <cusp/precond/aggregation/smoothed_aggregation.h>
//#include <cusp/precond/ainv.h>
#include <map>
#include <string>
#include "BCTools.H"
#include "CudaFArrayBox.hpp"
#include "IO.H"
#include "Operator.hpp"
#if CH_USE_PYTHON
#   include "PyGlue.H"
#endif
#include "icecream.hpp"
namespace Elliptic
{
namespace PythonInterface
{
template <typename IntType, typename ValueType, typename InputsMemory,
          typename OwnMemory, typename BCType>
class VectorCalculusKit
{
private:
  typedef typename cusp::csr_matrix<IntType, ValueType, OwnMemory> Matrix;
  typedef typename Elliptic::CudaInterface::Operator<IntType, ValueType,
                                                     InputsMemory, OwnMemory>
      Operator;
  typedef BaseFab<ValueType> BaseFAB;
  typedef typename std::vector<IntType> PyIntArray;
  typedef typename std::vector<ValueType> PyRealArray;
  typedef typename CudaBaseFab<ValueType, InputsMemory>::Array InputArr;
  typedef typename CudaBaseFab<ValueType, InputsMemory>::View InputView;
  typedef typename cusp::array1d<ValueType, OwnMemory> MyArr;
  typedef CudaFArrayBox<InputsMemory> DiagType;
#ifndef USE_CUSPARSE
  typedef int cudaStream_t;
#endif
  // Where we store the data
  std::map<std::string, Operator *> m_Operators;
  std::vector<std::string> m_OpNames;
  RealVect m_dx;
  IntVect m_refRatio;
  Box m_box, m_boxIntOnly;
  ProblemDomain m_domain;

  const BCType *m_BCD;
  DisjointBoxLayout m_grids;
  cudaStream_t m_StreamId;
  bool m_isDefined;

  std::string m_MyPyModule;
  // stores the boundary conditions
  bool m_ClearPythonAtTheEnd = false;
  bool m_isPythonMatrixCleared = true;
  int m_level;

public:
  std::vector<int> m_bcDesc;

private:
  // stores the names of the operators we need. Must be a subset of the list of
  // names for the Operator CSRData in ElliKit.py note that for the moment the
  // python script generates all of them.
  std::vector<std::string> m_opNames;

public:
  VectorCalculusKit() : m_isDefined(false), m_isPythonMatrixCleared(true) {} // as required by LayoutData

  ~VectorCalculusKit()
  {
    for (auto &M : m_Operators)
    {
      if (M.second != nullptr)
        delete M.second;
    }
    m_Operators.clear();
    m_bcDesc.clear();
    m_opNames.clear();
  }
  void define(const RealVect &a_dx, const BCType *a_BCD, const Box &a_box,
              const DisjointBoxLayout &a_grids, const ProblemDomain &domain, const IntVect refRatio,
              const bool ClearPythonAtTheEnd,
              const std::vector<std::string> OpNames,
              const int level,
              const cudaStream_t StreamId = 0,
              const std::string PyModuleName = "ElliKitWrapper")
  {
    m_dx = a_dx;
    m_box = a_box;
    m_boxIntOnly = a_box;
    m_MyPyModule = PyModuleName;
    m_domain = domain;
    m_BCD = a_BCD;
    m_grids = a_grids;
    m_StreamId = StreamId;
    m_OpNames = OpNames;
    m_ClearPythonAtTheEnd = ClearPythonAtTheEnd;
    m_isDefined = true;
    m_isPythonMatrixCleared = ClearPythonAtTheEnd;
    m_refRatio = refRatio;
    m_level = level;
    this->generate();
  }

  inline bool isDefined() { return m_isDefined; }
  inline bool isPythonMatrixCleared() { return m_isPythonMatrixCleared; }

  void ApplyOperator(std::string Key, const MyArr &a_x, MyArr &a_y)
  {
    assert(m_Operators.at(Key) != nullptr);
    (*(m_Operators.at(Key)))(a_x, a_y);
  }
  void ApplyOperator(std ::string Key, const DiagType &a_x,
                     DiagType &a_y) const
  {
    Interval I(0, a_x.nComp() - 1);

    this->ApplyOperator(Key, a_x, a_y, I);
  }

  void ApplyOperator(std::string Key, const DiagType &a_x, DiagType &a_y,
                     Interval I) const
  {
    assert(m_Operators.at(Key) != nullptr);
    (*(m_Operators.at(Key)))(a_x, a_y, I);
  }
  Matrix &retrieveOperatorMatrix(std::string Key) const
  {
    return *(m_Operators.at(Key)->getMatrixPtr());
  }
  Operator &retrieveOperator(std::string Key) const
  {
    return *(m_Operators.at(Key));
  }
  void deleteOperator(std::string Key)
  {
    delete m_Operators.at(Key);
    m_Operators.erase(Key);
  }

private:
  void generate();
};
template <typename IntType, typename ValueType, typename InputsMemory,
          typename OwnMemory, typename BC>
void VectorCalculusKit<IntType, ValueType, InputsMemory, OwnMemory, BC>::generate()
{
  m_bcDesc.resize(0);
  m_box = m_boxIntOnly;
  m_BCD->setBoxAndBCDescription(m_grids, m_boxIntOnly, m_domain, m_level, m_box, m_bcDesc);

  Py::PythonFunction(this->m_MyPyModule, "GenerateMatrix", this->m_dx,
                        this->m_box, this->m_refRatio, this->m_domain.domainBox(),
                        this->m_bcDesc);

  for (const auto &OperatorName : this->m_OpNames)
  {
    auto nnz = Py::PythonReturnFunction<int>(m_MyPyModule, "GetMatrixnnz",
                                                OperatorName);
    auto Nrows = Py::PythonReturnFunction<int>(m_MyPyModule, "GetMatrixRows",
                                                  OperatorName);
    auto Ncolumns = Py::PythonReturnFunction<int>(
        m_MyPyModule, "GetMatrixColumns", OperatorName);
    {
      PyIntArray IndPtr(Nrows + 1);
      PyIntArray Col(nnz);
      PyRealArray Data(nnz);
      auto N = Py::PythonReturnFunction<int>(
          m_MyPyModule, "GetMatrixData", OperatorName, IndPtr, Col, Data);
      // Allocate space for a new matrix, and store the pointer in the map

      m_Operators.insert(std::pair<std::string, Operator *>(
          OperatorName, new Operator(Nrows, Ncolumns, nnz, m_StreamId)));
      Operator &A = *m_Operators[OperatorName]; // get a local reference

      A.setOperator(IndPtr, Col, Data);

      if (OperatorName == "Divergence x")
      {
        Box edgeBox(surroundingNodes(m_boxIntOnly, 0));
        A.setGrids(edgeBox, m_boxIntOnly);
      }
      else if (OperatorName == "Divergence y")
      {
        Box edgeBox(surroundingNodes(m_boxIntOnly, 1));
        A.setGrids(edgeBox, m_boxIntOnly);
      }
      else if (OperatorName == "Divergence z")
      {
        Box edgeBox(surroundingNodes(m_boxIntOnly, 2));
        A.setGrids(edgeBox, m_boxIntOnly);
      }
      else if (OperatorName == "Gradient x")
      {
        Box edgeBox(surroundingNodes(m_boxIntOnly, 0));
        A.setGrids(m_box, edgeBox);
      }
      else if (OperatorName == "Gradient y")
      {
        Box edgeBox(surroundingNodes(m_boxIntOnly, 1));
        A.setGrids(m_box, edgeBox);
      }
      else if (OperatorName == "Gradient z")
      {
        Box edgeBox(surroundingNodes(m_boxIntOnly, 2));
        A.setGrids(m_box, edgeBox);
      }
    }
  }
  // KEEP THIS LAST!
  if (m_ClearPythonAtTheEnd)
  {
    Py::PythonFunction(m_MyPyModule, "ClearMatrix");
    m_isPythonMatrixCleared = true;
  }
}
} // namespace PythonInterface
} // namespace Elliptic
