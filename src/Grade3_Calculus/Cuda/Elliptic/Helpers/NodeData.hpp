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
// a class to represent LevelData objects as single sparse matrix over the whole
// domain
// This cannot have virtual members, or else it will cause problems with PoissonHelpers.
#include "cusp/array1d.h"

#include "FunctorsLibrary.hpp"
namespace Elliptic {
namespace Helpers {
template <typename ValueType, typename Memory,
          typename Container = cusp::array1d<ValueType, Memory>>
class NodeData {

public:
  typedef Memory MyMemory;
  NodeData() : m_isDefined(false), m_container(0) {} //leaves object unusable
  NodeData(const NodeData &x) : m_container(x()), m_isDefined(true) {}
  NodeData(int i) : m_container(i), m_isDefined(true) {}

  ~NodeData() {}

  cusp::array1d<ValueType, Memory> &operator()() { return m_container; }
  const cusp::array1d<ValueType, Memory> &operator()() const {
    return m_container;
  }
  ValueType dotProduct(const NodeData &y) const {

    assert(y.isCompatible(*this));

    return thrust::inner_product(y().begin(), y().end(),
                                 this->m_container.begin(), 0.0);
  }
  ValueType norm() const {
    return sqrt(this->dotProduct(*this) / this->m_container.size());
  }
  bool isCompatible(const NodeData &y) const {
    return this->m_container.size() == y().size();
  }
  NodeData &operator*=(const ValueType r) {
    this->scale(r);
    return *this;
  }
  NodeData &scale(const ValueType r) {

    scale_only<ValueType> M(r);
    thrust::for_each(this->m_container.begin(), this->m_container.end(), M);
    return *this;
  }
  void setVal(const ValueType r) {
    thrust::fill(this->m_container.begin(), this->m_container.end(), r);
  }
  void setToZero() { this->setVal(0.0); }
  NodeData &incr(const NodeData &y, const ValueType r) {
    assert(y.isCompatible(*this));
    scale_and_add<ValueType> SaD(r);
    thrust::transform(y().begin(), y().end(), this->m_container.begin(),
                      this->m_container.begin(), SaD);
    return *this;
  }
  NodeData &plus(const NodeData &y, const ValueType r) {
    assert(y.isCompatible(*this));
    this->incr(y, r);
    return *this;
  }
  NodeData &axpy(const NodeData &x, const NodeData &y, const ValueType a,
                 const ValueType b) {
    assert(y.isCompatible(*this));
    assert(x.isCompatible(*this));
    saxpby_functor<ValueType> SaD(a, b);
    thrust::transform(x().begin(), x().end(), y().begin(),
                      this->m_container.begin(), SaD);
    return *this;
  }
  NodeData &axby(const NodeData &x, const NodeData &y, const ValueType a,
                 const ValueType b) {
    this->axpy(x, y, a, b);
    return *this;
  }
  void create(const NodeData &x) {
    m_isDefined = true;
    this->m_container.resize(x().size());
  }
  void resize(int i)
  {
    m_isDefined = true;
    this->m_container.resize(i);
  }

  void assign(const NodeData &x) {
    this->m_container = x();
    m_isDefined = true;
  }
  void clear() {
    this->m_container.resize(0);
    m_isDefined = false;
  }
  void fillRandom(const ValueType min, const ValueType max) {
    assert(m_isDefined);
    Functors::PseudoNumberGenerator<ValueType> rnd(min, max);
    thrust::counting_iterator<unsigned int> index_sequence_begin(procID());
    thrust::transform(index_sequence_begin,
                      index_sequence_begin + m_container.size(),
                      m_container.begin(), rnd);
  }
  Real removeAverage(Real dV)
  {
    NodeData VolEl(this->m_container.size());
    VolEl.setVal(dV);
    Real r = this->dotProduct(VolEl);
    VolEl.setVal(r);


    this->incr(VolEl,-1.0);

    return r;
  }
private:
  bool m_isDefined = false;
  cusp::array1d<ValueType, Memory> m_container;
  NodeData &operator=(const NodeData &that){};
};

} // namespace Helpers
} // namespace Elliptic
