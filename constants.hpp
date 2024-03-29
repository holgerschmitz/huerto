/*
 * constants.hpp
 *
 *  Created on: 15 Jan 2020
 *  Author: Holger Schmitz
 */

#ifndef HUERTO_CONSTANTS_HPP_
#define HUERTO_CONSTANTS_HPP_

#include <schnek/variables/blockparameters.hpp>

#include <string>

static const double PI          = 3.141592653589793238462643383279502884L;
static const double TWO_PI      = 6.283185307179586476925286766559005768L;
static const double clight      = 299792458;
static const double clight2     = clight*clight;
static const double mass_e      = 9.10938291e-31;
static const double mass_p      = 1.672621777e-27;
static const double unit_charge = 1.602176565e-19;
static const double mu_0        = 4e-7*PI;
static const double eps_0       = 1/(mu_0*clight2);
static const double eps_0_inv   = (mu_0*clight2);
static const double planck_h    = 6.62607015e-34;
static const double planck_hbar = planck_h/TWO_PI;


static const std::string indexToCoordMapping[] = {"x", "y", "z", "u", "v", "w"};

inline const std::string indexToCoord(int index, std::string prefix = "", std::string postfix = "") {
  return prefix+indexToCoordMapping[index] + postfix;
}

void initConstantParameters(schnek::BlockParameters &parameters);

#endif /* HUERTO_CONSTANTS_HPP_ */
