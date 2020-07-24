/*
 * constants.cpp
 *
 *  Created on: 20 Jul 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#include "constants.hpp"

void initConstantParameters(schnek::BlockParameters &parameters) {
  parameters.addConstant("pi", PI);
  parameters.addConstant("clight", clight);
  parameters.addConstant("me", mass_e);
  parameters.addConstant("mp", mass_p);
  parameters.addConstant("e", unit_charge);
  parameters.addConstant("mu0", mu_0);
  parameters.addConstant("eps0", eps_0);
  parameters.addConstant("planck_h", planck_h);
  parameters.addConstant("planck_hbar", planck_hbar);
}


