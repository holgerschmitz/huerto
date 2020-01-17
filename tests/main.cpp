/*
 * main.cpp
 *
 *  Created on: 17 Jan 2020
 *      Author: Holger Schmitz
 */

#define BOOST_TEST_MODULE "Unit Tests for SchnAR"
#include <boost/test/included/unit_test.hpp>

#include <cmath>

bool is_equal(double a, double b)
{
  return ((a==0.0) && (b==0.0)) ||
      ( fabs(a-b)/(fabs(a)+fabs(b)) < 100*std::numeric_limits<double>::epsilon() );
}
