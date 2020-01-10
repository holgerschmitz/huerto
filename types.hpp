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
#endif

#ifdef SCHNAR_TWO_DIM
static const size_t DIMENSION = 2;
#endif

#ifdef SCHNAR_THREE_DIM
static const size_t DIMENSION = 3;
#endif


typedef schnek::Array<int, DIMENSION> Index;
typedef schnek::Array<double, DIMENSION> Vector;
typedef schnek::Grid<double, DIMENSION, SchnarGridChecker> Grid;
typedef boost::shared_ptr<Grid> pGrid;
typedef schnek::Field<double, DIMENSION, SchnarGridChecker> Field;
typedef boost::shared_ptr<Field> pField;
typedef schnek::Range<int, DIMENSION> Range;
typedef schnek::Array<bool, DIMENSION> Stagger;

enum Direction {north, south, west, east, up, down};

#endif // MPULSE_TYPES_H
