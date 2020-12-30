/*
 * slice_diagnostic.hpp
 *
 *  Created on: 26 Dec 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_DIAGNOSTIC_SLICE_DIAGNOSTIC_HPP_
#define HUERTO_DIAGNOSTIC_SLICE_DIAGNOSTIC_HPP_

#include "../../huerto/simulation/simulation_context.hpp"
#include "../../huerto/simulation/task.hpp"

#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

template<typename GridType, typename GridPtrType, typename DiagnosticType>
class GridSliceDiagnostic :
        public schnek::HDFGridDiagnostic<GridType, GridPtrType, DiagnosticType>,
        public SimulationEntity,
        public SimulationTask
{
  protected:
    typedef typename schnek::HDFGridDiagnostic<GridType, GridPtrType, DiagnosticType>::IndexType IndexType;
  private:
    int dim;
    int pos;

    bool outside;
    IndexType localMin;
    IndexType localMax;
    schnek::Range<int, DIMENSION> srcRange;
    int count;

    pField sourceField;
    GridPtrType field;

    void expandIndex(const IndexType &src, Index &dest);
  protected:
    IndexType getGlobalMin() override;
    IndexType getGlobalMax() override;
    bool isDerived() override { return true; }
    void initParameters(schnek::BlockParameters &parameters) override;
    void init() override;
    void write() override;
  public:
    std::string getPhase() override;
    void execute() override;
};

template<typename GridType, typename GridPtrType, typename DiagnosticType>
inline void GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::expandIndex(const IndexType &src, Index &dest) {
  for (int i=0; i<dim; ++i) {
    dest[i] = src[i];
  }
  dest[dim] = pos;

  for (size_t i=dim+1; i<DIMENSION; ++i) {
    dest[i] = src[i-1];
  }
}


template<typename GridType, typename GridPtrType, typename DiagnosticType>
typename GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::IndexType
  GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::getGlobalMin() {
  IndexType lo =  getContext().getSubdivision().getGlobalDomain().getLo();
  lo[dim] = 0;
  return lo;
}

template<typename GridType, typename GridPtrType, typename DiagnosticType>
typename GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::IndexType
  GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::getGlobalMax() {
  IndexType hi = getContext().getSubdivision().getGlobalDomain().getHi();
  hi[dim] = this->getInterval() - 1;
  return hi;
}

template<typename GridType, typename GridPtrType, typename DiagnosticType>
void GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::initParameters(schnek::BlockParameters& parameters) {
  schnek::HDFGridDiagnostic<GridType, GridPtrType, DiagnosticType>::initParameters(parameters);
  parameters.addParameter("dim", &dim);
  parameters.addParameter("pos", &pos);
}

template<typename GridType, typename GridPtrType, typename DiagnosticType>
void GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::init() {
  SimulationEntity::init(this);
  count = 0;

  this->retrieveData(this->getFieldName(), sourceField);
  Index gLocalMin = sourceField->getInnerLo();
  Index gLocalMax = sourceField->getInnerHi();

  outside = (pos < gLocalMin[dim]) || (pos > gLocalMax[dim]);

  if (!outside) {
    localMin = gLocalMin;
    localMin[dim] = 0;
    localMax = gLocalMax;
    localMax[dim] = this->getInterval() - 1;

    field = std::make_shared<GridType>(localMin, localMax);
    (*field) = 0.0;

    this->container.grid = field.get();
    this->container.local_min = localMin;
    this->container.local_max = localMax;
    this->container.global_min = this->getGlobalMin();
    this->container.global_max = this->getGlobalMax();

    gLocalMin[dim] = pos;
    gLocalMax[dim] = pos;

    srcRange = schnek::Range<int, DIMENSION>(gLocalMin, gLocalMax);
  }

  // deliberately not calling HDFGridDiagnostic::init() because we have done the initialisation
  // already if it is required
  schnek::Block::init();
}

template<typename GridType, typename GridPtrType, typename DiagnosticType>
void GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::write() {
  schnek::HDFGridDiagnostic<GridType, GridPtrType, DiagnosticType>::write();
  count = 0;
  (*field) = 0.0;
}

template<typename GridType, typename GridPtrType, typename DiagnosticType>
std::string GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::getPhase() {
  return "pre-diagnostic";
}
template<typename GridType, typename GridPtrType, typename DiagnosticType>
void GridSliceDiagnostic<GridType, GridPtrType, DiagnosticType>::execute() {
  GridType &destField = *field;
  Field &srcField = *sourceField;
  if (!outside && count < this->getInterval()) {
    for (Index src: srcRange) {
      Index dest = src;
      dest[dim] = count;
      destField[dest] = srcField[src];
    }
    ++count;
  }
}

#endif /* HUERTO_DIAGNOSTIC_SLICE_DIAGNOSTIC_HPP_ */
