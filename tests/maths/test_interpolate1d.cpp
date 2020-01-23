/*
 * interpolate1d.cpp
 *
 *  Created on: 23 Jan 2020
 *      Author: Holger Schmitz
 */

#include "../test_types.hpp"

#include "../../maths/interpolate/interpolate1d.hpp"
#include "../../constants.hpp"

#include <boost/test/unit_test.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE( maths )

BOOST_AUTO_TEST_SUITE( interpolate1d )

BOOST_AUTO_TEST_CASE( findInsertIndex_regular )
{
  Grid1d v(Index1d(0),Index1d(100));
  for (int i=0; i<=100; i++)
  {
    v(i) = i/100.0;
  }

  for (int i=0; i<100; i++)
  {
    BOOST_CHECK_EQUAL(findInsertIndex(v, (i+0.5)/100.0), i);
  }

  // corner cases
  BOOST_CHECK_EQUAL(findInsertIndex(v, 0.0), 0);
  BOOST_CHECK_EQUAL(findInsertIndex(v, 1.0), 100);

  // this works because 1/2^n does not have rounding errors
  BOOST_CHECK_EQUAL(findInsertIndex(v, 0.25), 25);
  BOOST_CHECK_EQUAL(findInsertIndex(v, 0.5), 50);
  BOOST_CHECK_EQUAL(findInsertIndex(v, 0.75), 75);

  // out of range
  BOOST_CHECK_EQUAL(findInsertIndex(v, -1.0), -1);
  BOOST_CHECK_EQUAL(findInsertIndex(v, 1.5), 100);


  BOOST_CHECK_EQUAL(findInsertIndex(v, 0.5 - std::numeric_limits<double>::epsilon()), 49);
}

BOOST_AUTO_TEST_CASE( linearInterpolate_regular )
{
  Grid1d x(Index1d(0),Index1d(100));
  Grid1d y(Index1d(0),Index1d(100));

  double tolerance = 5e-3;

  for (int i=0; i<=100; i++)
  {
    x(i) = i/100.0;
    y(i) = sin(PI*i/200.);
  }

  for (int i=0; i<400; i++)
  {
    double v = i/400.0;
    BOOST_CHECK_CLOSE(linearInterpolate(x, y, v), sin(PI*v/2.0), tolerance);
  }

  // corner cases
  BOOST_CHECK_EQUAL(linearInterpolate(x, y, 0.0), 0.0);
  BOOST_CHECK_EQUAL(linearInterpolate(x, y, 1.0), 1.0);

  // this works because 1/2^n does not have rounding errors
  BOOST_CHECK_CLOSE(linearInterpolate(x, y, 0.25), 0.5*sqrt(2-sqrt(2)), tolerance);
  BOOST_CHECK_CLOSE(linearInterpolate(x, y, 0.5), sqrt(0.5), tolerance);
  BOOST_CHECK_CLOSE(linearInterpolate(x, y, 0.75), 0.5*sqrt(2+sqrt(2)), tolerance);

  // out of range
  BOOST_CHECK_EQUAL(linearInterpolate(x, y, -1.0), 0.0);
  BOOST_CHECK_EQUAL(linearInterpolate(x, y, 1.5), 1.0);

}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
