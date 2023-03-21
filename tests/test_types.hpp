/*
 * test_types.hpp
 *
 *  Created on: 17 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_TESTS_TEST_TYPES_HPP_
#define HUERTO_TESTS_TEST_TYPES_HPP_

#include "../types.hpp"

bool is_equal(double a, double b);

template<size_t Rank>
bool is_equal(const schnek::Array<double, Rank> &a, const schnek::Array<double, Rank> &b)
{
  for (size_t d=0; d<Rank; ++d)
  {
    if (!is_equal(a[d], b[d])) return false;
  }
  return true;
}

#endif /* HUERTO_TESTS_TEST_TYPES_HPP_ */
