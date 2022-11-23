/*
 * random.cpp
 *
 *  Created on: 23 Jan 2020
 *      Author: Holger Schmitz
 */

#include "random.hpp"

boost::random::mt19937 rng;
boost::random::uniform_real_distribution<> random_unit_interval(0, 1);
boost::random::normal_distribution<> random_normal(0, 1);
boost::random::uniform_on_sphere<> random_sphere(3);
