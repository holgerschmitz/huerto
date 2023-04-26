/*
 * field_diagnostic.hpp
 *
 *  Created on: 26 Dec 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_DIAGNOSTIC_FIELD_DIAGNOSTIC_HPP_
#define HUERTO_DIAGNOSTIC_FIELD_DIAGNOSTIC_HPP_

#include "../../huerto/simulation/simulation_context.hpp"

#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

template<typename FieldType, typename DiagnosticType>
class FieldDiagnostic :
        public schnek::HDFGridDiagnostic<FieldType, DiagnosticType>,
        public SimulationEntity
{
  protected:
    typedef typename schnek::HDFGridDiagnostic<FieldType, DiagnosticType>::IndexType IndexType;
    IndexType getGlobalMin();
    IndexType getGlobalMax();
    void init();
};


template<typename FieldType, typename DiagnosticType>
typename FieldDiagnostic<FieldType, DiagnosticType>::IndexType
  FieldDiagnostic<FieldType, DiagnosticType>::getGlobalMin() {
  return getContext().getSubdivision().getGlobalDomain().getLo();
}

template<typename FieldType, typename DiagnosticType>
typename FieldDiagnostic<FieldType, DiagnosticType>::IndexType
  FieldDiagnostic<FieldType, DiagnosticType>::getGlobalMax() {
  return getContext().getSubdivision().getGlobalDomain().getHi();
}

template<typename FieldType, typename DiagnosticType>
void FieldDiagnostic<FieldType, DiagnosticType>::init() {
  SimulationEntity::init(this);
  schnek::HDFGridDiagnostic<FieldType, DiagnosticType>::init();
}

#endif /* HUERTO_DIAGNOSTIC_FIELD_DIAGNOSTIC_HPP_ */
