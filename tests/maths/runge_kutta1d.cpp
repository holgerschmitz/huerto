/*
 * runge_kutta1d.cpp
 *
 *  Created on: 24 Jul 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */


#include "../test_types.hpp"

#include "../../maths/integrate/runge_kutta.hpp"
#include "../../constants.hpp"

#include <boost/test/unit_test.hpp>

#include <cmath>

void noopBoundary() {}

BOOST_AUTO_TEST_SUITE( maths )

BOOST_AUTO_TEST_SUITE( runge_kutta_1d )

BOOST_AUTO_TEST_CASE( oscillator ){
  struct Oscillator {
      Field1d &om, &x, &v;
      Oscillator(Field1d &om_, Field1d &x_, Field1d &v_): om(om_), x(x_), v(v_) {}

      void operator()(Index1d p, Vector2d &dudt) const {
        double w = om[p];
        dudt[0] = v[p];
        dudt[1] = -w*w*x[p];
      }
  };

  Domain1d range(Vector1d(0.0), Vector1d(1.0));
  Stagger1d stagger(false);
  Field1d om(Index1d(0),Index1d(100), range, stagger, 1);
  Field1d x(Index1d(0),Index1d(100), range, stagger, 1);
  Field1d v(Index1d(0),Index1d(100), range, stagger, 1);

  for (int i=0; i<=100; i++)
  {
    om(i) = 1.0/(i+1.0);
    x(i) = 1.0;
    v(i) = 0.0;
  }

  FieldRungeKutta4<1,2> rk4;
  rk4.setField(0, x);
  rk4.setField(1, v);

  Oscillator osc(om, x, v);

  for (int i=0; i<=100; i++)
  {
    rk4.integrateStep(PI, osc, noopBoundary);
    double err = fabs(1.0 + x(i));
    BOOST_CHECK(err <= 12.5*pow(1.0/(i+1.0), 3));
    std::cout << "Result: " << (i+1) << " " << err << " " << 12.5*pow(1.0/(i+1.0), 3) << std::endl;
  }
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
