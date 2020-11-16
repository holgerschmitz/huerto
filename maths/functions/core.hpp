/*
 * core.hpp
 *
 *  Created on: 16 Nov 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_FUNCTIONS_CORE_HPP_
#define HUERTO_FUNCTIONS_CORE_HPP_

#include <schnek/parser/parser.hpp>
#include <schnek/variables/blockparameters.hpp>

#include <string>

double step(double x, double x0);

double stepi(double x, double x0);

double box(double x, double xmin, double xmax);

double diagf(std::string name, double value);

void registerCoreFunctions(schnek::FunctionRegistry & registry);

void registerConstants(schnek::BlockParameters &parameters);

#endif // HUERTO_FUNCTIONS_CORE_HPP_