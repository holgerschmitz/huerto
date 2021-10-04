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

  parameters.addConstant("centi", 1e-2);
  parameters.addConstant("milli", 1e-3);
  parameters.addConstant("micro", 1e-6);
  parameters.addConstant("nano", 1e-9);
  parameters.addConstant("pico", 1e-12);
  parameters.addConstant("femto", 1e-15);
  parameters.addConstant("atto", 1e-18);

  parameters.addConstant("kilo", 1e3);
  parameters.addConstant("mega", 1e6);
  parameters.addConstant("giga", 1e9);
  parameters.addConstant("tera", 1e12);
  parameters.addConstant("peta", 1e15);
  parameters.addConstant("exa", 1e18);

  parameters.addConstant("true", int(true));
  parameters.addConstant("yes", int(true));
  parameters.addConstant("false", int(false));
  parameters.addConstant("no", int(false));
}


int huerto_debug_flag = 0;