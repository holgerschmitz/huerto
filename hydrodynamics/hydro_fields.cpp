/*
 * em_fields.cpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */

#include "hydro_fields.hpp"
#include "../constants.hpp"
#include <schnek/grid/domainsubdivision.hpp>
#include <schnek/tools/fieldtools.hpp>

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include <string>
#include <iostream>

void HydroFields::initParameters(schnek::BlockParameters &parameters)
{
  Rho.parameter = parameters.addParameter("Rho", &Rho.value, 1.0);
  for (size_t i=0; i<DIMENSION; ++i)
  {
    M[i].parameter = parameters.addParameter(indexToCoord(i, "M"), &M[i].value , 0.0);
  }

  E.parameter = parameters.addParameter("E", &E.value, 1.0);
}

void HydroFields::registerData()
{
  Rho.field = std::make_shared<Field>();
  addData("Rho", Rho.field);

  for (size_t i=0; i<DIMENSION; ++i)
  {
    M[i].field = std::make_shared<Field>();
    addData(indexToCoord(i, "M"), M[i].field);
  }

  E.field = std::make_shared<Field>();
  addData("E", E.field);
}

void HydroFields::fillValues()
{
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  Vector &x = getContext().getX();
  schnek::Array<schnek::pParameter, DIMENSION> x_parameters = getContext().getXParameter();
  updater.addIndependentArray(x_parameters);

  schnek::fill_field(*Rho.field, x, Rho.value, updater, Rho.parameter);

  for (size_t i=0; i<DIMENSION; ++i)
  {
    schnek::fill_field(*M[i].field, x, M[i].value, updater, M[i].parameter);
  }

  schnek::fill_field(*E.field, x, E.value, updater, E.parameter);
}

void HydroFields::init()
{
  SimulationEntity::init(this);
  const schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();
  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  schnek::Range<double, DIMENSION> domainSize = subdivision.getInnerExtent(getContext().getSize());

  schnek::Array<bool, DIMENSION> stagger;
  stagger = true;

  Rho.field->resize(lowIn, highIn, domainSize, stagger, 2);
  for (size_t i=0; i<DIMENSION; ++i)
  {
    M[i].field->resize(lowIn, highIn, domainSize, stagger, 2);
  }
  E.field->resize(lowIn, highIn, domainSize, stagger, 2);
  fillValues();
}
