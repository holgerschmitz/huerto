/*
 * interpolate1d.hpp
 *
 *  Created on: 22 Jan 2020
 *  Author: Holger Schmitz
 */

#ifndef SCHNAR_MATHS_INTERPOLATE_INTERPOLATE1D_HPP_
#define SCHNAR_MATHS_INTERPOLATE_INTERPOLATE1D_HPP_

#include "../../types.hpp"

inline int findInsertIndex(const Grid1d &X, double x)
{
  int lo = X.getLo(0);
  int hi = X.getHi(0);

  while (lo <= hi) {
    int mid = (hi + lo) / 2;

    if (x < X(mid))
    {
      hi = mid - 1;
    }
    else
    {
      lo = mid + 1;
    }
  }
  return hi;
}

inline bool checkStrictlyAscending(const Grid1d &X)
{
  int lo = X.getLo(0);
  int hi = X.getHi(0);

  double last = std::numeric_limits<double>::lowest();
  for (int i=lo; i<=hi; ++i)
  {
    if (X(i) <= last) return false;
    last = X(i);
  }
  return true;
}

inline double linearInterpolate(const Grid1d &X, const Grid1d &Y, double x)
{
  if (x<=X(X.getLo(0))) return Y(X.getLo(0));
  if (x>=X(X.getHi(0))) return Y(X.getHi(0));
  int p = findInsertIndex(X, x);
  return (Y(p+1) - Y(p))*(x - X(p))/(X(p+1) - X(p)) + Y(p);
}

#endif /* SCHNAR_MATHS_INTERPOLATE_INTERPOLATE1D_HPP_ */
