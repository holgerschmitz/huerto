/*
 * vector.cpp
 *
 *  Created on: 17 Jan 2020
 *      Author: Holger Schmitz
 */

#include <boost/test/unit_test.hpp>

#include "../../maths/vector/vector.hpp"

#include "../test_types.hpp"


BOOST_AUTO_TEST_SUITE( maths )

BOOST_AUTO_TEST_SUITE( vector )

BOOST_AUTO_TEST_CASE( norm_2d )
{
  Vector2d v1(0, 0),
      v2(1.0, 0),
      v3(-1.0, 0),
      v4(0, 1.0),
      v5(0, -1.0),
      v6(3.0, 4.0),
      v7(3.0, 7.0),
      v8(-3.0, 4.0),
      v9(3.0, -7.0),
      v10(-3.0, 7.0),
      v11(-1e-42, 1.0),
      v12(1.0, 1e-42);

  BOOST_CHECK(is_equal(norm(v1), 0.0));
  BOOST_CHECK(is_equal(norm(v2), 1.0));
  BOOST_CHECK(is_equal(norm(v3), 1.0));
  BOOST_CHECK(is_equal(norm(v4), 1.0));
  BOOST_CHECK(is_equal(norm(v5), 1.0));
  BOOST_CHECK(is_equal(norm(v6), 5.0));
  BOOST_CHECK(is_equal(norm(v7), 7.61577310586391));
  BOOST_CHECK(is_equal(norm(v8), 5.0));
  BOOST_CHECK(is_equal(norm(v9), 7.61577310586391));
  BOOST_CHECK(is_equal(norm(v10), 7.61577310586391));
  BOOST_CHECK(is_equal(norm(v11), 1.0));
  BOOST_CHECK(is_equal(norm(v12), 1.0));
}

BOOST_AUTO_TEST_CASE( dot_2d )
{
  Vector2d v1(0, 0),
      v2(1.0, 0),
      v3(-1.0, 0),
      v4(0, 1.0),
      v5(0, -1.0),
      v6(3.0, 4.0),
      v7(3.0, 7.0),
      v8(-3.0, 4.0),
      v9(-4.0, 3.0),
      v10(7.0, -3.0),
      v11(-1e-42, 1.0),
      v12(1.0, 1e-42);

  BOOST_CHECK(is_equal(dot(v1, v2), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v3), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v4), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v5), 0.0));

  BOOST_CHECK(is_equal(dot(v2, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v2),  1.0));
  BOOST_CHECK(is_equal(dot(v2, v3), -1.0));
  BOOST_CHECK(is_equal(dot(v2, v4),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v5),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v6),  3.0));
  BOOST_CHECK(is_equal(dot(v2, v11), -1e-42));

  BOOST_CHECK(is_equal(dot(v3, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v2), -1.0));
  BOOST_CHECK(is_equal(dot(v3, v3),  1.0));
  BOOST_CHECK(is_equal(dot(v3, v4),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v5),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v6), -3.0));
  BOOST_CHECK(is_equal(dot(v3, v11), 1e-42));

  BOOST_CHECK(is_equal(dot(v4, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v2),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v3),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v4),  1.0));
  BOOST_CHECK(is_equal(dot(v4, v5), -1.0));
  BOOST_CHECK(is_equal(dot(v4, v6),  4.0));
  BOOST_CHECK(is_equal(dot(v4, v12), 1e-42));

  BOOST_CHECK(is_equal(dot(v6, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v6, v2),  3.0));
  BOOST_CHECK(is_equal(dot(v6, v3), -3.0));
  BOOST_CHECK(is_equal(dot(v6, v4),  4.0));
  BOOST_CHECK(is_equal(dot(v6, v5), -4.0));
  BOOST_CHECK(is_equal(dot(v6, v6),  25.0));
  BOOST_CHECK(is_equal(dot(v6, v9),  0.0));
}

BOOST_AUTO_TEST_CASE( absMaxIndex_2d )
{
  Vector2d v1(0, 0),
      v2(1.0, 0),
      v3(-1.0, 0),
      v4(0, 1.0),
      v5(0, -1.0),
      v6(3.0, 4.0),
      v7(3.0, 7.0),
      v8(-3.0, 4.0),
      v9(-4.0, 3.0),
      v10(7.0, -3.0),
      v11(-1e-42, 1.0),
      v12(1.0, 1e-42);

  BOOST_CHECK_EQUAL(absMaxIndex(v1), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v2), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v3), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v4), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v5), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v6), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v7), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v8), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v9), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v10), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v11), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v12), 0);
}

BOOST_AUTO_TEST_CASE( perpendicular_2d )
{
  Vector2d v1(0, 0),
      v2(1.0, 0),
      v3(-1.0, 0),
      v4(0, 1.0),
      v5(0, -1.0),
      v6(3.0, 4.0),
      v7(3.0, 7.0),
      v8(-3.0, 4.0),
      v9(-4.0, 3.0),
      v10(7.0, -3.0);

  BOOST_CHECK(is_equal(perpendicular(v1), v1));
  BOOST_CHECK(is_equal(perpendicular(v2), v5));
  BOOST_CHECK(is_equal(perpendicular(v3), v4));
  BOOST_CHECK(is_equal(perpendicular(v4), v2));
  BOOST_CHECK(is_equal(perpendicular(v5), v3));
  BOOST_CHECK(is_equal(perpendicular(v6), Vector2d(4.0, -3.0)));
  BOOST_CHECK(is_equal(perpendicular(v7), v10));
  BOOST_CHECK(is_equal(perpendicular(v8), Vector2d(4.0, 3.0)));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()


