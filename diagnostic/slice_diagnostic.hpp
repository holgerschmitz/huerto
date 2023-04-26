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

template<typename GridType, typename DiagnosticType>
class GridSliceDiagnostic :
        public schnek::HDFGridDiagnostic<GridType, DiagnosticType>,
        public SimulationEntity,
        public SimulationTask
{
  protected:
    typedef typename schnek::HDFGridDiagnostic<GridType, DiagnosticType>::IndexType IndexType;
  private:
    int dim;
    int pos;

    bool outside;
    IndexType localMin;
    IndexType localMax;
    schnek::Range<int, DIMENSION> srcRange;
    int count;

    Field sourceField;
    GridType field;

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

template<typename GridType, typename DiagnosticType>
inline void GridSliceDiagnostic<GridType, DiagnosticType>::expandIndex(const IndexType &src, Index &dest) {
  for (int i=0; i<dim; ++i) {
    dest[i] = src[i];
  }
  dest[dim] = pos;

  for (size_t i=dim+1; i<DIMENSION; ++i) {
    dest[i] = src[i-1];
  }
}


template<typename GridType, typename DiagnosticType>
typename GridSliceDiagnostic<GridType, DiagnosticType>::IndexType
  GridSliceDiagnostic<GridType, DiagnosticType>::getGlobalMin() {
  IndexType lo =  getContext().getSubdivision().getGlobalDomain().getLo();
  lo[dim] = 0;
  return lo;
}

template<typename GridType, typename DiagnosticType>
typename GridSliceDiagnostic<GridType, DiagnosticType>::IndexType
  GridSliceDiagnostic<GridType, DiagnosticType>::getGlobalMax() {
  IndexType hi = getContext().getSubdivision().getGlobalDomain().getHi();
  hi[dim] = this->getInterval() - 1;
  return hi;
}

template<typename GridType, typename DiagnosticType>
void GridSliceDiagnostic<GridType, DiagnosticType>::initParameters(schnek::BlockParameters& parameters) {
  schnek::HDFGridDiagnostic<GridType, DiagnosticType>::initParameters(parameters);
  parameters.addParameter("dim", &dim);
  parameters.addParameter("pos", &pos);
}

template<typename GridType, typename DiagnosticType>
void GridSliceDiagnostic<GridType, DiagnosticType>::init() {
  SimulationEntity::init(this);
  count = 0;

  this->retrieveData(this->getFieldName(), sourceField);
  Index gLocalMin = sourceField.getInnerLo();
  Index gLocalMax = sourceField.getInnerHi();

  outside = (pos < gLocalMin[dim]) || (pos > gLocalMax[dim]);

  if (!outside) {
    localMin = gLocalMin;
    localMin[dim] = 0;
    localMax = gLocalMax;
    localMax[dim] = this->getInterval() - 1;
    this->output.setActive(true);
  } else {
    localMin = Index(-1);
    localMax = Index(-1);
    this->output.setActive(false);
  }

  field.resize(localMin, localMax);
  field = 0.0;

  this->container.grid = field.get();
  this->container.local_min = localMin;
  this->container.local_max = localMax;
  this->container.global_min = this->getGlobalMin();
  this->container.global_max = this->getGlobalMax();

  if (!outside) {
    gLocalMin[dim] = pos;
    gLocalMax[dim] = pos;

    srcRange = schnek::Range<int, DIMENSION>(gLocalMin, gLocalMax);
  }

  // deliberately not calling HDFGridDiagnostic::init() because we have done the initialisation
  // already if it is required
  schnek::Block::init();
}

template<typename GridType, typename DiagnosticType>
void GridSliceDiagnostic<GridType, DiagnosticType>::write() {
  schnek::HDFGridDiagnostic<GridType, DiagnosticType>::write();
  count = 0;
  field = 0.0;
}

template<typename GridType, typename DiagnosticType>
std::string GridSliceDiagnostic<GridType, DiagnosticType>::getPhase() {
  return "pre-diagnostic";
}
template<typename GridType, typename DiagnosticType>
void GridSliceDiagnostic<GridType, DiagnosticType>::execute() {
  if (outside || count >= this->getInterval()) {
    return;
  }

  for (Index src: srcRange) {
    Index dest = src;
    dest[dim] = count;
    field[dest] = sourceField[src];
  }
  ++count;
}

#endif /* HUERTO_DIAGNOSTIC_SLICE_DIAGNOSTIC_HPP_ */
