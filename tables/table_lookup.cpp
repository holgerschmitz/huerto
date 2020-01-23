/*
 * table_lookup.cpp
 *
 *  Created on: 25 Jul 2019
 *      Author: Holger Schmitz
 */

#include "table_lookup.hpp"
#include "../maths/interpolate/interpolate1d.hpp"
#include "../maths/interpolate/interpolate2d.hpp"
#include "../maths/random.hpp"

#include <boost/make_shared.hpp>

#include <exception>

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
  FileTableDataReader::preInit();
}


/******************************************************************************/
/*    TableLookup                                                             */
/******************************************************************************/

void TableLookup::init(TableBlock &table, int xIndex, int yIndex, bool cumulative)
{
  xValues = &table.getValues(xIndex);
  yValues = &table.getValues(yIndex);

  if (!checkStrictlyAscending(*xValues))
  {
    throw std::runtime_error("In table "+table.getName()+": x axis not strictly ascending");
  }

  if (cumulative)
  {
    initCumulative();
  }
}

void TableLookup::init(TableLookup &tableA, double weightA, TableLookup &tableB, double weightB, bool cumulative)
{
  xValues = tableA.xValues;

  if (!checkStrictlyAscending(*xValues))
  {
    throw std::runtime_error("In table lookup: x axis not strictly ascending");
  }

  int lo = xValues->getLo(0);
  int hi = xValues->getHi(0);

  yValues = new Grid1d(lo, hi);
  ownedData.push_back(std::unique_ptr<Grid1d>(yValues));

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
  return linearInterpolate(*xValues, *yValues, x);
}


void TableLookup::initCumulative()
{
  Grid1d &X = *xValues;
  Grid1d &Y = *yValues;
  int lo = X.getLo(0);
  int hi = X.getHi(0);

  yCumulative = new Grid1d(lo, hi);

  Grid1d &Yc = *yCumulative;

  // sum up values
  Yc(lo) = 0;
  double sum = 0;
  for (int i=lo+1; i<=hi; ++i)
  {
    double delta = Y(i) * (X(i) - X(i-1));
    if (delta<=0)
    {
      throw std::runtime_error("Negative values in probability distribution!");
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
  return linearInterpolate(*yCumulative, *xValues, random_unit_interval(rng));
}


/******************************************************************************/
/*    TableLookup2d                                                           */
/******************************************************************************/


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
  return linearInterpolate2d(xValues, yValues, table, x, y);
}
