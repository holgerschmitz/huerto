/*
 * generic_defs.hpp
 *
 *  Created on: 30 Jan 2020
 *      Author: Holger Schmitz
 */

#include <schnek/grid/arrayexpression.hpp>

#define HUERTO_FUNC_EXPR_EXPR(type, func)                                       \
  template<class exp1, class exp2, size_t rank>                                    \
  inline type func (                                                            \
    const schnek::ArrayExpression<exp1, rank> &A,                             \
    const schnek::ArrayExpression<exp2, rank> &B)                             \
  {                                                                             \
    return func(schnek::Array<type, rank>(A), schnek::Array<type, rank>(B));    \
  }

#define HUERTO_FUNC_ARR_EXPR(type, func)                                        \
  template<class exp, size_t rank, template<size_t> class CheckPolicy>                              \
  inline type func (                                                            \
    const schnek::ArrayExpression<exp, rank> &A,                              \
    const schnek::Array<type, rank, CheckPolicy> &B)                          \
  {                                                                             \
    return func(schnek::Array<type, rank>(A), B);                               \
  }

#define HUERTO_FUNC_EXPR_ARR(type, func)                                        \
  template<class exp, size_t rank, template<size_t> class CheckPolicy>                              \
  inline type func (                                                            \
    const schnek::Array<type, rank, CheckPolicy> &A,                           \
    const schnek::ArrayExpression<exp, rank> &B)                              \
  {                                                                             \
    return func(A, schnek::Array<type, rank>(B));                               \
  }
