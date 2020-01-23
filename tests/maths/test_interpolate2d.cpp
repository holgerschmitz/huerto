/*
 * interpolate1d.cpp
 *
 *  Created on: 23 Jan 2020
 *      Author: Holger Schmitz
 */

#include "../test_types.hpp"

#include "../../maths/interpolate/interpolate2d.hpp"
#include "../../constants.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE( maths )

BOOST_AUTO_TEST_SUITE( interpolate2d )


BOOST_AUTO_TEST_CASE( linearInterpolate2d_regular )
{
  boost::progress_display show_progress(402);
  Grid1d x(Index1d(0), Index1d(100));
  Grid1d y(Index1d(0), Index1d(100));
  Grid2d z(Index2d(0, 0), Index2d(100, 100));

  double tolerance = 5e-3;

  for (int i=0; i<=100; i++)
    for (int j=0; j<=100; j++)
    {
      double xv = i/100.0;
      double yv = j/100.0;
      x(i) = xv;
      y(j) = yv;
      z(i, j) = sin(0.5*PI*xv)*cos(0.5*PI*yv);
    }

  for (int i=0; i<=400; i++)
  {
    for (int j=0; j<=400; j++)
    {
      double xv = i/400.0;
      double yv = j/400.0;
      BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, xv, yv) + 1.0, sin(0.5*PI*xv)*cos(0.5*PI*yv) + 1.0, tolerance);
    }
    ++show_progress;
  }

  for (int i=0; i<400; i++)
  {
    double xv = i/400.0;
    // edge cases
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, xv, 0.0) + 1.0, sin(0.5*PI*xv) + 1.0, tolerance);
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, xv, 1.0) + 1.0, 1.0, tolerance);
    BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, 0.0, xv), 0.0);
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, 1.0, xv) + 1.0, cos(0.5*PI*xv) + 1.0, tolerance);

    // out of range
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, xv, -1.0) + 1.0, sin(0.5*PI*xv) + 1.0, tolerance);
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, xv, -std::numeric_limits<double>::min()) + 1.0, sin(0.5*PI*xv) + 1.0, tolerance);
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, xv, 2.0) + 1.0, 1.0, tolerance);
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, xv, 1.0 + std::numeric_limits<double>::epsilon()) + 1.0, 1.0, tolerance);

    BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, -1.0, xv), 0.0);
    BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, -std::numeric_limits<double>::min(), xv), 0.0);
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, 2.0, xv) + 1.0, cos(0.5*PI*xv) + 1.0, tolerance);
    BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, 1.0 + std::numeric_limits<double>::epsilon(), xv) + 1.0, cos(0.5*PI*xv) + 1.0, tolerance);
  }

  ++show_progress;
  // corner cases
  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, 0.0, 0.0), 0.0);
  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, 1.0, 0.0), 1.0);
  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, 0.0, 1.0), 0.0);
  BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, 1.0, 1.0) + 1.0, 1.0, tolerance);

  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, -1.0, -1.0), 0.0);
  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, -std::numeric_limits<double>::min(), -std::numeric_limits<double>::min()), 0.0);
  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, 2.0, -1.0), 1.0);
  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, 1.0 + std::numeric_limits<double>::epsilon(), -std::numeric_limits<double>::min()), 1.0);
  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, -1.0, 2.0), 0.0);
  BOOST_CHECK_EQUAL(linearInterpolate2d(x, y, z, -std::numeric_limits<double>::min(), 1.0 + std::numeric_limits<double>::epsilon()), 0.0);
  BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, 2.0, 2.0) + 1.0, 1.0, tolerance);
  BOOST_CHECK_CLOSE(linearInterpolate2d(x, y, z, 1.0 + std::numeric_limits<double>::epsilon(), 1.0 + std::numeric_limits<double>::epsilon()) + 1.0, 1.0, tolerance);

}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
