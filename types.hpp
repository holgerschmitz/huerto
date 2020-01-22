/*
 * types.hpp
 *
 *  Created on: 10 Jan 2020
 *      Author: Holger Schmitz
 */


#ifndef SCHNAR_TYPES_H
#define SCHNAK_TYPES_H

#include <schnek/grid.hpp>

#ifndef SCHNAR_GRID_CHECKER
#ifdef SCHNAR_NDEBUG
#define SchnarGridChecker schnek::GridNoArgCheck
#else
#define SchnarGridChecker schnek::GridAssertCheck
#endif
#else
#define SchnarGridChecker SCHNAR_GRID_CHECKER
#endif

#ifdef SCHNAR_ONE_DIM
static const size_t DIMENSION = 1;
enum Direction {west, east};
#endif

#ifdef SCHNAR_TWO_DIM
static const size_t DIMENSION = 2;
enum Direction {west, east, south, north};
#endif

#ifdef SCHNAR_THREE_DIM
static const size_t DIMENSION = 3;
enum Direction {west, east, south, north, down, up};
#endif

typedef schnek::Array<int, DIMENSION> Index;
typedef schnek::Array<double, DIMENSION> Vector;
typedef schnek::Grid<double, DIMENSION, SchnarGridChecker> Grid;
typedef boost::shared_ptr<Grid> pGrid;
typedef schnek::Field<double, DIMENSION, SchnarGridChecker> Field;
typedef boost::shared_ptr<Field> pField;
typedef schnek::Range<int, DIMENSION> Range;
typedef schnek::Array<bool, DIMENSION> Stagger;

typedef schnek::Array<int, 1> Index1d;
typedef schnek::Array<double, 1> Vector1d;
typedef schnek::Grid<double, 1, SchnarGridChecker> Grid1d;
typedef boost::shared_ptr<Grid1d> pGrid1d;
typedef schnek::Field<double, 1, SchnarGridChecker> Field1d;
typedef boost::shared_ptr<Field1d> pField1d;
typedef schnek::Range<int, 1> Range1d;
typedef schnek::Array<bool, 1> Stagger1d;

typedef schnek::Array<int, 2> Index2d;
typedef schnek::Array<double, 2> Vector2d;
typedef schnek::Grid<double, 2, SchnarGridChecker> Grid2d;
typedef boost::shared_ptr<Grid2d> pGrid2d;
typedef schnek::Field<double, 2, SchnarGridChecker> Field2d;
typedef boost::shared_ptr<Field2d> pField2d;
typedef schnek::Range<int, 2> Range2d;
typedef schnek::Array<bool, 2> Stagger2d;

typedef schnek::Array<int, 3> Index3d;
typedef schnek::Array<double, 3> Vector3d;
typedef schnek::Grid<double, 3, SchnarGridChecker> Grid3d;
typedef boost::shared_ptr<Grid3d> pGrid3d;
typedef schnek::Field<double, 3, SchnarGridChecker> Field3d;
typedef boost::shared_ptr<Field3d> pField3d;
typedef schnek::Range<int, 3> Range3d;
typedef schnek::Array<bool, 3> Stagger3d;


#endif // SCHNAR_TYPES_H
