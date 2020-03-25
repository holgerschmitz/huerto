/*
 * interpolate2d.hpp
 *
 *  Created on: 23 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_MATHS_INTERPOLATE_INTERPOLATE2D_HPP_
#define HUERTO_MATHS_INTERPOLATE_INTERPOLATE2D_HPP_

#include "interpolate1d.hpp"


inline double linearInterpolate2d(const Grid1d &X, const Grid1d &Y, const Grid2d &T, double x, double y)
{

  int ix, ixp, iy, iyp;
  int xl = X.getLo(0), xh = X.getHi(0);
  int yl = Y.getLo(0), yh = Y.getHi(0);

  if (x<=X(xl))
  {
    ix = ixp = xl;
  }
  else if (x>=X(xh))
  {
    ix = ixp = xh;
  }
  else
  {
    ix = findInsertIndex(X, x);
    ixp = ix + 1;
  }

  if (y<=Y(yl))
  {
    iy = iyp = yl;
  }
  else if (y>=Y(yh))
  {
    iy = iyp = yh;
  }
  else
  {
    iy = findInsertIndex(Y, y);
    iyp = iy + 1;
  }

  double Xl = X(ix);
  double Xh = X(ixp);
  double Yl = Y(iy);
  double Yh = Y(iyp);

  double Tll = T(ix,  iy);
  double Tlh = T(ix,  iyp);
  double Thl = T(ixp, iy);
  double Thh = T(ixp, iyp);

  double xA = ix != ixp ? (x - Xl)/(Xh - Xl) : 0.0;
  double yA = iy != iyp ? (y - Yl)/(Yh - Yl) : 0.0;

  double Tli = (Thl - Tll)*xA + Tll;
  double Thi = (Thh - Tlh)*xA + Tlh;
  double result = (Thi - Tli)*yA + Tli;
  return result;
}


#endif /* HUERTO_MATHS_INTERPOLATE_INTERPOLATE2D_HPP_ */
