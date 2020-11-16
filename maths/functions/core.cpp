/*
 * core.cpp
 *
 *  Created on: 16 Nov 2020
 *      Author: Holger Schmitz
 */

#include "core.hpp"
#include "../../constants.hpp"

double step(double x, double x0)
{
  return (x>=x0)?1.0:0.0;
}

double stepi(double x, double x0)
{
  return (x>=x0)?0.0:1.0;
}

double box(double x, double xmin, double xmax)
{
  return ((x>=xmin)&&(x<xmax))?1.0:0.0;
}

double diagf(std::string name, double value)
{
  std::cout << "LOG " << name << ": " << value << std::endl;
  return value;
}

void registerCoreFunctions(schnek::FunctionRegistry & registry) {
  registry.registerFunction("step", step);
  registry.registerFunction("stepi", stepi);
  registry.registerFunction("box", box);
  registry.registerFunction("diagf", diagf);
}

void registerConstants(schnek::BlockParameters &parameters) {
  parameters.addConstant("pi", PI);
  parameters.addConstant("clight", clight);
  parameters.addConstant("me", mass_e);
  parameters.addConstant("mp", mass_p);
  parameters.addConstant("e", unit_charge);
  parameters.addConstant("mu0", mu_0);
  parameters.addConstant("eps0", eps_0);
  parameters.addConstant("false", 0);
  parameters.addConstant("true", 1);
}
