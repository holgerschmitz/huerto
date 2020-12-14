/*
 * em_fields.cpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */

#include "em_fields.hpp"
#include "../types.hpp"
#include "../constants.hpp"

#include <schnek/util/logger.hpp>
#include <schnek/grid/domainsubdivision.hpp>
#include <schnek/tools/fieldtools.hpp>

#include <string>
#include <iostream>


void EMFields::initParameters(schnek::BlockParameters &parameters)
{
  for (size_t i=0; i<3; ++i)
  {
    E[i].parameter = parameters.addParameter(indexToCoord(i, "E"), &E[i].value , 0.0);
    B[i].parameter = parameters.addParameter(indexToCoord(i, "B"), &B[i].value , 0.0);
  }
}

void EMFields::registerData() {
  for (size_t i=0; i<3; ++i) {
    E[i].field = std::make_shared<Field>();
    B[i].field = std::make_shared<Field>();
    addData(indexToCoord(i, "E"), E[i].field);
    addData(indexToCoord(i, "B"), B[i].field);
  }
}

void EMFields::fillValues() {
  std::cout << "Filling fields" << std::endl;
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  Vector &x = getContext().getX();
  schnek::Array<schnek::pParameter, DIMENSION> x_parameters = getContext().getXParameter();
  updater.addIndependentArray(x_parameters);

  for (size_t i=0; i<3; ++i) {
    schnek::fill_field(*E[i].field, x, E[i].value, updater, E[i].parameter);
    schnek::fill_field(*B[i].field, x, B[i].value, updater, B[i].parameter);
  }
}

void EMFields::preInit() {
  schnek::ChildBlock<EMFields>::preInit();
}

void EMFields::init() {
  schnek::ChildBlock<EMFields>::init();
  SimulationEntity::init(this);
  const schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();
  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  schnek::Range<double, DIMENSION> domainSize = subdivision.getInnerExtent(getContext().getSize());

  E[0].field->resize(lowIn, highIn, domainSize, exStaggerYee, 2);
  E[1].field->resize(lowIn, highIn, domainSize, eyStaggerYee, 2);
  E[2].field->resize(lowIn, highIn, domainSize, ezStaggerYee, 2);

  B[0].field->resize(lowIn, highIn, domainSize, bxStaggerYee, 2);
  B[1].field->resize(lowIn, highIn, domainSize, byStaggerYee, 2);
  B[2].field->resize(lowIn, highIn, domainSize, bzStaggerYee, 2);
  fillValues();
}
// end of main
