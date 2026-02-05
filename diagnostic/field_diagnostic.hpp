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
        public schnek::HDFGridRegistrationDiagnostic<FieldType, HuertoDecomposition, DiagnosticType>,
        public SimulationEntity
{
  protected:
    typedef typename schnek::HDFGridRegistrationDiagnostic<FieldType, HuertoDecomposition, DiagnosticType>::IndexType IndexType;
    IndexType getGlobalMin();
    IndexType getGlobalMax();
    void init() override;
    HuertoDecomposition &getDecomposition() override;
};


template<typename FieldType, typename DiagnosticType>
typename FieldDiagnostic<FieldType, DiagnosticType>::IndexType
  FieldDiagnostic<FieldType, DiagnosticType>::getGlobalMin() {
  return getContext().getDecomposition().getGlobalRange().getLo();
}

template<typename FieldType, typename DiagnosticType>
typename FieldDiagnostic<FieldType, DiagnosticType>::IndexType
  FieldDiagnostic<FieldType, DiagnosticType>::getGlobalMax() {
  return getContext().getDecomposition().getGlobalRange().getHi();
}

template<typename FieldType, typename DiagnosticType>
void FieldDiagnostic<FieldType, DiagnosticType>::init() {
  SimulationEntity::init(this);
  schnek::HDFGridRegistrationDiagnostic<FieldType, HuertoDecomposition, DiagnosticType>::init();
}

template<typename FieldType, typename DiagnosticType>
HuertoDecomposition &FieldDiagnostic<FieldType, DiagnosticType>::getDecomposition() {
  return getContext().getDecomposition();
}

#endif /* HUERTO_DIAGNOSTIC_FIELD_DIAGNOSTIC_HPP_ */
