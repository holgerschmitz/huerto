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
  auto &decomposition = getContext().getDecomposition();
  E[0].field = decomposition.registerField(schnek::GridFactory<Field>{exStaggerYee, 2});
  E[1].field = decomposition.registerField(schnek::GridFactory<Field>{eyStaggerYee, 2});
  E[2].field = decomposition.registerField(schnek::GridFactory<Field>{ezStaggerYee, 2});

  B[0].field = decomposition.registerField(schnek::GridFactory<Field>{bxStaggerYee, 2});
  B[1].field = decomposition.registerField(schnek::GridFactory<Field>{byStaggerYee, 2});
  B[2].field = decomposition.registerField(schnek::GridFactory<Field>{bzStaggerYee, 2});

  for (size_t i=0; i<3; ++i) {
    addData(indexToCoord(i, "E"), E[i].field);
    addData(indexToCoord(i, "B"), B[i].field);
  }
}

void EMFields::fillValues() {
  std::cout << "Filling fields" << std::endl;
  auto &decomposition = getContext().getDecomposition();
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  Vector &x = getContext().getX();
  schnek::Array<schnek::pParameter, DIMENSION> x_parameters = getContext().getXParameter();
  updater.addIndependentArray(x_parameters);

  for (size_t i=0; i<3; ++i) {
    auto gridContext = decomposition.getGridContext({E[i].field, B[i].field});
    gridContext.forEach([&](Range& /* range */, Field &Efield, Field &Bfield) {
      schnek::fill_field(Efield, x, E[i].value, updater, E[i].parameter);
      schnek::fill_field(Bfield, x, B[i].value, updater, B[i].parameter);
    });
  }
}

void EMFields::preInit() {
  schnek::ChildBlock<EMFields>::preInit();
}

void EMFields::init() {
  schnek::ChildBlock<EMFields>::init();
  SimulationEntity::init(this);

  fillValues();
}
// end of main
