/*
 * fieldsum.hpp
 *
 *  Created on: 09 Jan 2026
 *      Author: Holger Schmitz
 */
#include <schnek/macros.hpp>

#include "../types.hpp"

template<typename AccType, typename FieldType, typename IndexType = Index>
struct FieldSum {
  AccType accumulator;
  FieldType field;
  SCHNEK_INLINE void operator()(IndexType pos) {
    accumulator[pos] += field[pos];
  }
};