/*
 * random.hpp
 *
 *  Created on: 23 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_MATHS_RANDOM_HPP_
#define HUERTO_MATHS_RANDOM_HPP_

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_on_sphere.hpp>

extern boost::random::mt19937 rng;
extern boost::random::uniform_real_distribution<> random_unit_interval;
extern boost::random::uniform_on_sphere<> random_sphere;

#endif /* HUERTO_MATHS_RANDOM_HPP_ */
