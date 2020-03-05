/*
 * vector.hpp
 *
 *  Created on: 16 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef SCHNAR_MATHS_VECTOR_VECTOR_HPP_
#define SCHNAR_MATHS_VECTOR_VECTOR_HPP_

#include <schnek/grid.hpp>
#include <cmath>

#include "../../util/generic_defs.hpp"

template<int Rank, template<int> class CheckingPolicy>
inline double norm(const schnek::Array<double, Rank, CheckingPolicy> &v)
{
    return sqrt(v.sqr());
}

template<typename T, int Rank, template<int> class CheckingPolicy>
inline double min(const schnek::Array<T, Rank, CheckingPolicy> &v)
{
    T m = v[0];
    for (int d=0; d<Rank; ++d)
    {
      m = std::min(m, v[d]);
    }
    return m;
}

template<
  int Rank,
  template<int> class CheckingPolicyA,
  template<int> class CheckingPolicyB
>
inline double dot(const schnek::Array<double, Rank, CheckingPolicyA> &a,
                  const schnek::Array<double, Rank, CheckingPolicyB> &b)
{
    double sum = 0;
    for (int d=0; d<Rank; ++d)
    {
      sum += a[d]*b[d];
    }
    return sum;
}

SCHNAR_FUNC_ARR_EXPR(double, dot)
SCHNAR_FUNC_EXPR_ARR(double, dot)
SCHNAR_FUNC_EXPR_EXPR(double, dot)

template<
  template<int> class CheckingPolicyA,
  template<int> class CheckingPolicyB
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

template<int Rank, template<int> class CheckingPolicy>
inline int absMaxIndex(const schnek::Array<double, Rank, CheckingPolicy> &x)
{
  int m = 0;
  double M = fabs(x[0]);
  for (int i=1; i<Rank; ++i)
    if (fabs(x[i])>M)
    {
      m = i;
      M = fabs(x[i]);
    }

  return m;
}


template<template<int> class CheckingPolicy>
inline schnek::Array<double, 2, CheckingPolicy>
    perpendicular(const schnek::Array<double, 2, CheckingPolicy> &v)
{
    return schnek::Array<double, 2, CheckingPolicy>(v[1], -v[0]);
}


template<template<int> class CheckingPolicy>
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


#endif /* SCHNAR_MATHS_VECTOR_VECTOR_HPP_ */
