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

BOOST_AUTO_TEST_CASE( norm_1d )
{
  Vector1d v1(0), v2(1.0), v3(-1.0),
      v4(1e-42), v5(-1e-42), v6(1e42), v7(-1e42);

  BOOST_CHECK(is_equal(norm(v1), 0.0));
  BOOST_CHECK(is_equal(norm(v2), 1.0));
  BOOST_CHECK(is_equal(norm(v3), 1.0));
  BOOST_CHECK(is_equal(norm(v4), 1e-42));
  BOOST_CHECK(is_equal(norm(v5), 1e-42));
  BOOST_CHECK(is_equal(norm(v6), 1e42));
  BOOST_CHECK(is_equal(norm(v7), 1e42));
}

BOOST_AUTO_TEST_CASE( dot_1d )
{
  Vector1d v1(0), v2(1.0), v3(-1.0),
      v4(1e-42), v5(-1e-42), v6(1e42), v7(-1e42);

  BOOST_CHECK(is_equal(dot(v1, v2), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v3), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v4), 0.0));

  BOOST_CHECK(is_equal(dot(v2, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v2),  1.0));
  BOOST_CHECK(is_equal(dot(v2, v3), -1.0));
  BOOST_CHECK(is_equal(dot(v2, v4), 1e-42));
  BOOST_CHECK(is_equal(dot(v2, v5), -1e-42));
  BOOST_CHECK(is_equal(dot(v2, v6), 1e42));
  BOOST_CHECK(is_equal(dot(v2, v7), -1e42));

  BOOST_CHECK(is_equal(dot(v3, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v2), -1.0));
  BOOST_CHECK(is_equal(dot(v3, v3),  1.0));
  BOOST_CHECK(is_equal(dot(v3, v4), -1e-42));
  BOOST_CHECK(is_equal(dot(v3, v5), 1e-42));
  BOOST_CHECK(is_equal(dot(v3, v6), -1e42));
  BOOST_CHECK(is_equal(dot(v3, v7), 1e42));

  BOOST_CHECK(is_equal(dot(v4, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v2),  1e-42));
  BOOST_CHECK(is_equal(dot(v4, v3), -1e-42));
  BOOST_CHECK(is_equal(dot(v4, v4),  1e-84));
  BOOST_CHECK(is_equal(dot(v4, v5), -1e-84));
  BOOST_CHECK(is_equal(dot(v4, v6),  1.0));
  BOOST_CHECK(is_equal(dot(v4, v7), -1.0));
}

BOOST_AUTO_TEST_CASE( absMaxIndex_1d )
{
  Vector1d v1(0), v2(1.0), v3(-1.0),
      v4(1e-42), v5(-1e-42), v6(1e42), v7(-1e42);

  BOOST_CHECK_EQUAL(absMaxIndex(v1), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v2), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v3), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v4), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v5), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v6), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v7), 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()


