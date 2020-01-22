/*
 * table_lookup.cpp
 *
 *  Created on: 25 Jul 2019
 *      Author: Holger Schmitz
 */

#include "table_lookup.hpp"
#include "../maths.hpp"

#include <boost/make_shared.hpp>

/******************************************************************************/
/*    TableBlock                                                              */
/******************************************************************************/

void TableBlock::initParameters(schnek::BlockParameters &parameters)
{
  FileTableDataReader::initParameters(parameters);
}

void TableBlock::registerData()
{
  self = this;
  addData("TABLE_"+getName(), self);
}

void TableBlock::preInit()
{
  DataReader::preInit();
}


/******************************************************************************/
/*    TableLookup                                                             */
/******************************************************************************/

inline std::pair<int, int> findIndex(const Grid1d &X, double x)
{
  int lo = X.getLo(0);
  int hi = X.getHi(0);

  if (x<=X(lo)) return std::pair<int, int>(lo, lo);
  if (x>=X(hi)) return std::pair<int, int>(hi, hi);

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
  return std::pair<int, int>(lo, hi);
}

inline double doInterpolate(const Grid1d &X, const Grid1d &Y, double x)
{
  if (x<=X(X.getLo(0))) return Y(X.getLo(0));
  if (x>=X(X.getHi(0))) return Y(X.getHi(0));
  std::pair<int, int> p = findIndex(X, x);
  return (Y(p.first) - Y(p.second))*(x - X(p.second))/(X(p.first) - X(p.second)) + Y(p.second);
}

void TableLookup::init(TableBlock &table, int xIndex, int yIndex, bool cumulative)
{
  xValues = &table.getValues(xIndex);
  yValues = &table.getValues(yIndex);
  if (cumulative)
  {
    initCumulative();
  }
}

void TableLookup::init(TableLookup &tableA, double weightA, TableLookup &tableB, double weightB, bool cumulative)
{
  xValues = tableA.xValues;

  int lo = xValues->getLo(0);
  int hi = xValues->getHi(0);

  pGrid1d combined = boost::make_shared<Grid1d>(lo, hi);
  ownedData.push_back(combined);
  yValues = combined.get();

  Grid1d &yA = *(tableA.yValues);
  Grid1d &yB = *(tableB.yValues);

  for (int i=lo+1; i<=hi; ++i)
  {
    (*yValues)(i) = weightA*yA(i) + weightB*yB(i);
  }

  if (cumulative)
  {
    initCumulative();
  }
}

double TableLookup::interpolate(double x) const
{
  return doInterpolate(*xValues, *yValues, x);
}


void TableLookup::initCumulative()
{
  Grid1d &X = *xValues;
  Grid1d &Y = *yValues;
  int lo = X.getLo(0);
  int hi = X.getHi(0);

  pGrid1d yCumulativePtr = boost::make_shared<Grid1d>(lo, hi);
  ownedData.push_back(yCumulativePtr);
  yCumulative = yCumulativePtr.get();

  Grid1d &Yc = *yCumulative;

  // sum up values
  Yc(lo) = 0;
  double sum = 0;
  for (int i=lo+1; i<=hi; ++i)
  {
    double delta = Y(i) * (X(i) - X(i-1));
    if (delta<=0)
    {
      std::cerr << "Negative values in probability distribution!" << std::endl;
      exit(-1);
    }
    sum += delta;
    Yc(i) = sum;
  }

  // normalise cumulative distribution
  for (int i=lo; i<=hi; ++i)
  {
    Yc(i) /= sum;
  }
}

double TableLookup::randomDist()
{
  return doInterpolate(*yCumulative, *xValues, random_unit_interval(rng));
}


/******************************************************************************/
/*    TableLookup2d                                                           */
/******************************************************************************/


inline double doInterpolate2d(const Grid1d &X, const Grid1d &Y, const Grid2d &T, double x, double y)
{
  std::pair<int, int> pX = findIndex(X, x);
  std::pair<int, int> pY = findIndex(Y, y);

  double Xl = X(pX.first);
  double Xh = X(pX.second);
  double Yl = Y(pY.first);
  double Yh = Y(pY.second);

  double Tll = T(pX.first, pY.first);
  double Tlh = T(pX.first, pY.second);
  double Thl = T(pX.second, pY.first);
  double Thh = T(pX.second, pY.second);

  double xA = pX.first != pX.second ? (x - Xl)/(Xh - Xl) : 0.0;
  double yA = pY.first != pY.second ? (y - Yl)/(Yh - Yl) : 0.0;

  double Tli = (Tlh - Tll)*xA + Tll;
  double Thi = (Thh - Thl)*xA + Thl;
  return (Thi - Tli)*yA + Tli;
}


void TableLookup2d::init(TableBlock &tableBlock)
{
  int Nx = tableBlock.getValues(0).getDims(0) - 1;
  int Ny = tableBlock.getColumns() - 1;

  xValues.resize(0, Nx - 1);
  yValues.resize(0, Ny - 1);
  table.resize(Index2d(0, 0), Index2d(Nx-1, Ny-1));

  Grid1d &firstCol = tableBlock.getValues(0);
  for (int i=0; i<Nx; ++i)
  {
    xValues[i] = firstCol[i+1];
  }

  for (int j=0; j<Ny; ++j)
  {
    yValues[j] = tableBlock.getValues(j+1)[0];
  }

  for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j)
    {
      table(i,j) = tableBlock.getValues(j+1)[i+1];
    }
}


double TableLookup2d::interpolate(double x, double y) const
{
  return doInterpolate2d(xValues, yValues, table, x, y);
}
