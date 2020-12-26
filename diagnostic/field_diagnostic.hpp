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

template<typename FieldType, typename FieldPtrType, typename DiagnosticType>
class FieldDiagnostic :
        public schnek::HDFGridDiagnostic<FieldType, FieldPtrType, DiagnosticType>,
        public SimulationEntity
{
  protected:
    typedef typename schnek::HDFGridDiagnostic<FieldType, FieldPtrType, DiagnosticType>::IndexType IndexType;
    IndexType getGlobalMin();
    IndexType getGlobalMax();
    void init();
};


template<typename FieldType, typename FieldPtrType, typename DiagnosticType>
typename FieldDiagnostic<FieldType, FieldPtrType, DiagnosticType>::IndexType
  FieldDiagnostic<FieldType, FieldPtrType, DiagnosticType>::getGlobalMin() {
  return getContext().getSubdivision().getGlobalDomain().getLo();
}

template<typename FieldType, typename FieldPtrType, typename DiagnosticType>
typename FieldDiagnostic<FieldType, FieldPtrType, DiagnosticType>::IndexType
  FieldDiagnostic<FieldType, FieldPtrType, DiagnosticType>::getGlobalMax() {
  return getContext().getSubdivision().getGlobalDomain().getHi();
}

template<typename FieldType, typename FieldPtrType, typename DiagnosticType>
void FieldDiagnostic<FieldType, FieldPtrType, DiagnosticType>::init() {
  SimulationEntity::init(this);
  schnek::HDFGridDiagnostic<FieldType, FieldPtrType, DiagnosticType>::init();
}

#endif /* HUERTO_DIAGNOSTIC_FIELD_DIAGNOSTIC_HPP_ */
