/*
 * vector.hpp
 *
 *  Created on: 16 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_MATHS_VECTOR_VECTOR_HPP_
#define HUERTO_MATHS_VECTOR_VECTOR_HPP_

#include <schnek/grid.hpp>
#include <cmath>

#include "../../util/generic_defs.hpp"

template<size_t Rank, template<size_t> class CheckingPolicy>
inline double norm(const schnek::Array<double, Rank, CheckingPolicy> &v)
{
    return sqrt(v.sqr());
}

template<typename T, size_t Rank, template<size_t> class CheckingPolicy>
inline double min(const schnek::Array<T, Rank, CheckingPolicy> &v)
{
    T m = v[0];
    for (size_t d=0; d<Rank; ++d)
    {
      m = std::min(m, v[d]);
    }
    return m;
}

template<
  size_t Rank,
  template<size_t> class CheckingPolicyA,
  template<size_t> class CheckingPolicyB
>
inline double dot(const schnek::Array<double, Rank, CheckingPolicyA> &a,
                  const schnek::Array<double, Rank, CheckingPolicyB> &b)
{
    double sum = 0;
    for (size_t d=0; d<Rank; ++d)
    {
      sum += a[d]*b[d];
    }
    return sum;
}

HUERTO_FUNC_ARR_EXPR(double, dot)
HUERTO_FUNC_EXPR_ARR(double, dot)
HUERTO_FUNC_EXPR_EXPR(double, dot)

template<
  template<size_t> class CheckingPolicyA,
  template<size_t> class CheckingPolicyB
>
inline schnek::Array<double, 3, CheckingPolicyA>
    cross(const schnek::Array<double, 3, CheckingPolicyA> &a,
          const schnek::Array<double, 3, CheckingPolicyB> &b)
{
    return schnek::Array<double, 3, CheckingPolicyA>(
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}

template<class exp1, class exp2>
inline schnek::Array<double, 3> cross (
  const schnek::ArrayExpression<exp1, 3> &A,
  const schnek::ArrayExpression<exp2, 3> &B)
{
  return cross(schnek::Array<double, 3>(A), schnek::Array<double, 3>(B));
}

template<class exp, template<size_t> class CheckPolicy>
inline schnek::Array<double, 3, CheckPolicy> cross (
  const schnek::ArrayExpression<exp, 3> &A,
  const schnek::Array<double, 3, CheckPolicy> &B)
{
  return cross(schnek::Array<double, 3, CheckPolicy>(A), B);
}

template<class exp, template<size_t> class CheckPolicy>
inline schnek::Array<double, 3, CheckPolicy> cross (
  const schnek::Array<double, 3, CheckPolicy> &A,
  const schnek::ArrayExpression<exp, 3> &B)
{
  return cross(A, schnek::Array<double, 3, CheckPolicy>(B));
}

template<size_t Rank, template<size_t> class CheckingPolicy>
inline int absMaxIndex(const schnek::Array<double, Rank, CheckingPolicy> &x)
{
  int m = 0;
  double M = fabs(x[0]);
  for (size_t i=1; i<Rank; ++i)
    if (fabs(x[i])>M)
    {
      m = i;
      M = fabs(x[i]);
    }

  return m;
}


template<template<size_t> class CheckingPolicy>
inline schnek::Array<double, 2, CheckingPolicy>
    perpendicular(const schnek::Array<double, 2, CheckingPolicy> &v)
{
    return schnek::Array<double, 2, CheckingPolicy>(v[1], -v[0]);
}


template<template<size_t> class CheckingPolicy>
inline schnek::Array<double, 3, CheckingPolicy>
    perpendicular(const schnek::Array<double, 3, CheckingPolicy> &v)
{
    int xM = absMaxIndex(v);
    int xP = (xM+1) % 3;
    schnek::Array<double, 3, CheckingPolicy> perp(0.0);
    perp[xM] =  v[xP];
    perp[xP] = -v[xM];
    return perp;
}

#include "../../util/generic_undefs.hpp"


#endif /* HUERTO_MATHS_VECTOR_VECTOR_HPP_ */
