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

BOOST_AUTO_TEST_CASE( norm_3d )
{
  Vector3d v1(0, 0, 0),
      v2(1.0, 0, 0),
      v3(-1.0, 0, 0),
      v4(0, 1.0, 0),
      v5(0, -1.0, 0),
      v6(0, 0, 1.0),
      v7(0, 0, -1.0),
      v8(3.0, 2.0, 6.0),
      v9(3.0, 7.0, 9.0),
      v10(-3.0, 2.0, 6.0),
      v11(3.0, -7.0, 0.0),
      v12(0.0, 3.0, -7.0),
      v13(-1e-42, 3.0, 4.0),
      v14(4.0, 1e-42, 3.0);

  BOOST_CHECK(is_equal(norm(v1), 0.0));
  BOOST_CHECK(is_equal(norm(v2), 1.0));
  BOOST_CHECK(is_equal(norm(v3), 1.0));
  BOOST_CHECK(is_equal(norm(v4), 1.0));
  BOOST_CHECK(is_equal(norm(v5), 1.0));
  BOOST_CHECK(is_equal(norm(v6), 1.0));
  BOOST_CHECK(is_equal(norm(v7), 1.0));
  BOOST_CHECK(is_equal(norm(v8), 7.0));
  BOOST_CHECK(is_equal(norm(v9), 11.7898261225516));
  BOOST_CHECK(is_equal(norm(v10), 7.0));
  BOOST_CHECK(is_equal(norm(v11), 7.61577310586391));
  BOOST_CHECK(is_equal(norm(v12), 7.61577310586391));
  BOOST_CHECK(is_equal(norm(v13), 5.0));
  BOOST_CHECK(is_equal(norm(v14), 5.0));
}

BOOST_AUTO_TEST_CASE( dot_3d )
{
  Vector3d v1(0, 0, 0),
      v2(1.0, 0, 0),
      v3(-1.0, 0, 0),
      v4(0, 1.0, 0),
      v5(0, -1.0, 0),
      v6(0, 0, 1.0),
      v7(0, 0, -1.0),
      v8(3.0, 2.0, 6.0),
      v9(3.0, 7.0, 9.0),
      v10(-3.0, 2.0, 6.0),
      v11(7.0, -3.0, 0.0),
      v12(0.0, 3.0, -7.0),
      v13(-1e-42, 3.0, 4.0),
      v14(4.0, 1e-42, 3.0),
      v15(3.0, 4.0, 1e-42);

  BOOST_CHECK(is_equal(dot(v1, v2), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v3), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v4), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v5), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v6), 0.0));
  BOOST_CHECK(is_equal(dot(v1, v7), 0.0));

  BOOST_CHECK(is_equal(dot(v2, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v2),  1.0));
  BOOST_CHECK(is_equal(dot(v2, v3), -1.0));
  BOOST_CHECK(is_equal(dot(v2, v4),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v5),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v6),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v7),  0.0));
  BOOST_CHECK(is_equal(dot(v2, v8),  3.0));
  BOOST_CHECK(is_equal(dot(v2, v13), -1e-42));

  BOOST_CHECK(is_equal(dot(v3, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v2), -1.0));
  BOOST_CHECK(is_equal(dot(v3, v3),  1.0));
  BOOST_CHECK(is_equal(dot(v3, v4),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v5),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v6),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v7),  0.0));
  BOOST_CHECK(is_equal(dot(v3, v8), -3.0));
  BOOST_CHECK(is_equal(dot(v3, v13), 1e-42));

  BOOST_CHECK(is_equal(dot(v4, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v2),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v3),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v4),  1.0));
  BOOST_CHECK(is_equal(dot(v4, v5), -1.0));
  BOOST_CHECK(is_equal(dot(v4, v6),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v7),  0.0));
  BOOST_CHECK(is_equal(dot(v4, v8),  2.0));
  BOOST_CHECK(is_equal(dot(v4, v14), 1e-42));

  BOOST_CHECK(is_equal(dot(v5, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v5, v2),  0.0));
  BOOST_CHECK(is_equal(dot(v5, v3),  0.0));
  BOOST_CHECK(is_equal(dot(v5, v4), -1.0));
  BOOST_CHECK(is_equal(dot(v5, v5),  1.0));
  BOOST_CHECK(is_equal(dot(v5, v6),  0.0));
  BOOST_CHECK(is_equal(dot(v5, v7),  0.0));
  BOOST_CHECK(is_equal(dot(v5, v8), -2.0));
  BOOST_CHECK(is_equal(dot(v5, v14), -1e-42));

  BOOST_CHECK(is_equal(dot(v6, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v6, v2),  0.0));
  BOOST_CHECK(is_equal(dot(v6, v3),  0.0));
  BOOST_CHECK(is_equal(dot(v6, v4),  0.0));
  BOOST_CHECK(is_equal(dot(v6, v5),  0.0));
  BOOST_CHECK(is_equal(dot(v6, v6),  1.0));
  BOOST_CHECK(is_equal(dot(v6, v7), -1.0));
  BOOST_CHECK(is_equal(dot(v6, v8),  6.0));
  BOOST_CHECK(is_equal(dot(v6, v15), 1e-42));

  BOOST_CHECK(is_equal(dot(v7, v1),  0.0));
  BOOST_CHECK(is_equal(dot(v7, v2),  0.0));
  BOOST_CHECK(is_equal(dot(v7, v3),  0.0));
  BOOST_CHECK(is_equal(dot(v7, v4),  0.0));
  BOOST_CHECK(is_equal(dot(v7, v5),  0.0));
  BOOST_CHECK(is_equal(dot(v7, v6), -1.0));
  BOOST_CHECK(is_equal(dot(v7, v7),  1.0));
  BOOST_CHECK(is_equal(dot(v7, v8), -6.0));
  BOOST_CHECK(is_equal(dot(v7, v15), -1e-42));

  BOOST_CHECK(is_equal(dot(v8, v8),  49.0));
  BOOST_CHECK(is_equal(dot(v8, v9),  77.0));
  BOOST_CHECK(is_equal(dot(v8, v10), 31.0));
  BOOST_CHECK(is_equal(dot(v9, v11),  0.0));
  BOOST_CHECK(is_equal(dot(v13, v14), 12.0));
  BOOST_CHECK(is_equal(dot(v13, v15), 12.0));
  BOOST_CHECK(is_equal(dot(v14, v15), 12.0));
}

BOOST_AUTO_TEST_CASE( cross_3d )
{
  Vector3d v1(0, 0, 0),
      v2(1.0, 0, 0),
      v3(-1.0, 0, 0),
      v4(0, 1.0, 0),
      v5(0, -1.0, 0),
      v6(0, 0, 1.0),
      v7(0, 0, -1.0),
      v8(3.0, 2.0, 6.0),
      v9(3.0, 7.0, 9.0),
      v10(-3.0, 2.0, 6.0),
      v11(7.0, -3.0, 0.0),
      v12(0.0, 3.0, -7.0),
      v13(-1e-42, 3.0, 4.0),
      v14(4.0, 1e-42, 3.0),
      v15(3.0, 4.0, 1e-42);

  BOOST_CHECK(is_equal(cross(v1, v2), v1));
  BOOST_CHECK(is_equal(cross(v1, v3), v1));
  BOOST_CHECK(is_equal(cross(v1, v4), v1));
  BOOST_CHECK(is_equal(cross(v1, v5), v1));
  BOOST_CHECK(is_equal(cross(v1, v6), v1));
  BOOST_CHECK(is_equal(cross(v1, v7), v1));

  BOOST_CHECK(is_equal(cross(v2, v1),  v1));
  BOOST_CHECK(is_equal(cross(v2, v2),  v1));
  BOOST_CHECK(is_equal(cross(v2, v3),  v1));
  BOOST_CHECK(is_equal(cross(v2, v4),  v6));
  BOOST_CHECK(is_equal(cross(v2, v5),  v7));
  BOOST_CHECK(is_equal(cross(v2, v6),  v5));
  BOOST_CHECK(is_equal(cross(v2, v7),  v4));
  BOOST_CHECK(is_equal(cross(v2, v8),  3.0));
  BOOST_CHECK(is_equal(cross(v2, v13), -1e-42));

  BOOST_CHECK(is_equal(cross(v3, v1),  v1));
  BOOST_CHECK(is_equal(cross(v3, v2),  v1));
  BOOST_CHECK(is_equal(cross(v3, v3),  v1));
  BOOST_CHECK(is_equal(cross(v3, v4),  0.0));
  BOOST_CHECK(is_equal(cross(v3, v5),  0.0));
  BOOST_CHECK(is_equal(cross(v3, v6),  0.0));
  BOOST_CHECK(is_equal(cross(v3, v7),  0.0));
  BOOST_CHECK(is_equal(cross(v3, v8), -3.0));
  BOOST_CHECK(is_equal(cross(v3, v13), 1e-42));

  BOOST_CHECK(is_equal(cross(v4, v1),  v1));
  BOOST_CHECK(is_equal(cross(v4, v2),  0.0));
  BOOST_CHECK(is_equal(cross(v4, v3),  0.0));
  BOOST_CHECK(is_equal(cross(v4, v4),  v1));
  BOOST_CHECK(is_equal(cross(v4, v5),  v1));
  BOOST_CHECK(is_equal(cross(v4, v6),  0.0));
  BOOST_CHECK(is_equal(cross(v4, v7),  0.0));
  BOOST_CHECK(is_equal(cross(v4, v8),  2.0));
  BOOST_CHECK(is_equal(cross(v4, v14), 1e-42));

  BOOST_CHECK(is_equal(cross(v5, v1),  v1));
  BOOST_CHECK(is_equal(cross(v5, v2),  0.0));
  BOOST_CHECK(is_equal(cross(v5, v3),  0.0));
  BOOST_CHECK(is_equal(cross(v5, v4),  v1));
  BOOST_CHECK(is_equal(cross(v5, v5),  v1));
  BOOST_CHECK(is_equal(cross(v5, v6),  0.0));
  BOOST_CHECK(is_equal(cross(v5, v7),  0.0));
  BOOST_CHECK(is_equal(cross(v5, v8), -2.0));
  BOOST_CHECK(is_equal(cross(v5, v14), -1e-42));

  BOOST_CHECK(is_equal(cross(v6, v1),  v1));
  BOOST_CHECK(is_equal(cross(v6, v2),  0.0));
  BOOST_CHECK(is_equal(cross(v6, v3),  0.0));
  BOOST_CHECK(is_equal(cross(v6, v4),  0.0));
  BOOST_CHECK(is_equal(cross(v6, v5),  0.0));
  BOOST_CHECK(is_equal(cross(v6, v6),  v1));
  BOOST_CHECK(is_equal(cross(v6, v7),  v1));
  BOOST_CHECK(is_equal(cross(v6, v8),  6.0));
  BOOST_CHECK(is_equal(cross(v6, v15), 1e-42));

  BOOST_CHECK(is_equal(cross(v7, v1),  v1));
  BOOST_CHECK(is_equal(cross(v7, v2),  0.0));
  BOOST_CHECK(is_equal(cross(v7, v3),  0.0));
  BOOST_CHECK(is_equal(cross(v7, v4),  0.0));
  BOOST_CHECK(is_equal(cross(v7, v5),  0.0));
  BOOST_CHECK(is_equal(cross(v7, v6),  v1));
  BOOST_CHECK(is_equal(cross(v7, v7),  v1));
  BOOST_CHECK(is_equal(cross(v7, v8), -6.0));
  BOOST_CHECK(is_equal(cross(v7, v15), -1e-42));

  BOOST_CHECK(is_equal(cross(v8, v8),  v1));
  BOOST_CHECK(is_equal(cross(v8, v9),  77.0));
  BOOST_CHECK(is_equal(cross(v8, v10), 31.0));
  BOOST_CHECK(is_equal(cross(v9, v11),  0.0));
  BOOST_CHECK(is_equal(cross(v13, v14), 12.0));
  BOOST_CHECK(is_equal(cross(v13, v15), 12.0));
  BOOST_CHECK(is_equal(cross(v14, v15), 12.0));
}

BOOST_AUTO_TEST_CASE( absMaxIndex_3d )
{
  Vector3d v1(0, 0, 0),
      v2(1.0, 0, 0),
      v3(-1.0, 0, 0),
      v4(0, 1.0, 0),
      v5(0, -1.0, 0),
      v6(0, 0, 1.0),
      v7(0, 0, -1.0),
      v8(3.0, 2.0, 6.0),
      v9(3.0, 7.0, 9.0),
      v10(-3.0, 2.0, 6.0),
      v11(3.0, -7.0, 0.0),
      v12(0.0, 3.0, -7.0),
      v13(-1e-42, 3.0, 4.0),
      v14(4.0, 1e-42, 3.0),
      v15(3.0, 4.0, 1e-42);

  BOOST_CHECK_EQUAL(absMaxIndex(v1), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v2), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v3), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v4), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v5), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v6), 2);
  BOOST_CHECK_EQUAL(absMaxIndex(v7), 2);
  BOOST_CHECK_EQUAL(absMaxIndex(v8), 2);
  BOOST_CHECK_EQUAL(absMaxIndex(v9), 2);
  BOOST_CHECK_EQUAL(absMaxIndex(v10), 2);
  BOOST_CHECK_EQUAL(absMaxIndex(v11), 1);
  BOOST_CHECK_EQUAL(absMaxIndex(v12), 2);
  BOOST_CHECK_EQUAL(absMaxIndex(v13), 2);
  BOOST_CHECK_EQUAL(absMaxIndex(v14), 0);
  BOOST_CHECK_EQUAL(absMaxIndex(v15), 1);
}

BOOST_AUTO_TEST_CASE( perpendicular_3d )
{
  Vector3d v1(0, 0, 0),
      v2(1.0, 0, 0),
      v3(-1.0, 0, 0),
      v4(0, 1.0, 0),
      v5(0, -1.0, 0),
      v6(0, 0, 1.0),
      v7(0, 0, -1.0),
      v8(3.0, 2.0, 6.0),
      v9(3.0, 9.0, 5.0),
      v10(-3.0, 2.0, 6.0),
      v11(3.0, -7.0, 0.0),
      v12(0.0, 3.0, -7.0);

  BOOST_CHECK(is_equal(perpendicular(v1), v1));
  BOOST_CHECK(is_equal(perpendicular(v2), v5));
  BOOST_CHECK(is_equal(perpendicular(v3), v4));
  BOOST_CHECK(is_equal(perpendicular(v4), v7));
  BOOST_CHECK(is_equal(perpendicular(v5), v6));
  BOOST_CHECK(is_equal(perpendicular(v6), v3));
  BOOST_CHECK(is_equal(perpendicular(v7), v2));
  BOOST_CHECK(is_equal(perpendicular(v8), Vector3d(-6.0, 0.0, 3.0)));
  BOOST_CHECK(is_equal(perpendicular(v9), Vector3d( 0.0, 5.0,-9.0)));
  BOOST_CHECK(is_equal(perpendicular(v10), Vector3d(-6.0, 0.0,-3.0)));
  BOOST_CHECK(is_equal(perpendicular(v11), Vector3d( 0.0, 0.0, 7.0)));
  BOOST_CHECK(is_equal(perpendicular(v12), Vector3d( 7.0, 0.0, 0.0)));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()


