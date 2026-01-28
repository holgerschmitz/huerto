/*
 * types.hpp
 *
 *  Created on: 10 Jan 2020
 *      Author: Holger Schmitz
 */


#ifndef HUERTO_TYPES_H
#define HUERTO_TYPES_H

#include <schnek/grid/array.hpp>
#include <schnek/grid/grid.hpp>
#include <schnek/grid/field.hpp>
#include <schnek/grid/iteration/range-iteration.hpp>
#include <schnek/decomposition/mpi_cartesian_decomposition.hpp>

#include <memory>

#if !defined(HUERTO_ONE_DIM) && !defined(HUERTO_TWO_DIM) && !defined(HUERTO_THREE_DIM)
std::static_assert(false, "Must define Huerto Dimension before including types.hpp");
#endif

#ifndef HUERTO_GRID_CHECKER
#ifdef HUERTO_NDEBUG
#define HuertoGridChecker schnek::GridNoArgCheck
#else
#define HuertoGridChecker schnek::GridAssertCheck
#endif
#else
#define HuertoGridChecker HUERTO_GRID_CHECKER
#endif

#ifdef HUERTO_ONE_DIM
static const size_t DIMENSION = 1;
enum Direction {west, east};
#endif

#ifdef HUERTO_TWO_DIM
static const size_t DIMENSION = 2;
enum Direction {west, east, south, north};
#endif

#ifdef HUERTO_THREE_DIM
static const size_t DIMENSION = 3;
enum Direction {west, east, south, north, down, up};
#endif

typedef schnek::Array<ptrdiff_t, DIMENSION> Index;
typedef schnek::Array<size_t, DIMENSION> Size;
typedef schnek::Array<double, DIMENSION> Vector;
typedef schnek::Grid<double, DIMENSION, HuertoGridChecker> Grid;
typedef std::shared_ptr<Grid> pGrid;
typedef std::unique_ptr<Grid> uGrid;
typedef std::reference_wrapper<Grid> rGrid;
typedef schnek::Field<double, DIMENSION, HuertoGridChecker> Field;
typedef std::shared_ptr<Field> pField;
typedef std::unique_ptr<Field> uField;
typedef std::reference_wrapper<Field> rField;
typedef schnek::Range<ptrdiff_t, DIMENSION> Range;
typedef schnek::Range<double, DIMENSION> Domain;
typedef schnek::Array<bool, DIMENSION> Stagger;

typedef schnek::Field<double, DIMENSION-1, HuertoGridChecker> Surface;
typedef std::shared_ptr<Surface> pSurface;
typedef std::unique_ptr<Surface> uSurface;
typedef std::reference_wrapper<Surface> rSurface;

typedef schnek::Array<ptrdiff_t, 1> Index1d;
typedef schnek::Array<size_t, 1> Size1d;
typedef schnek::Array<double, 1> Vector1d;
typedef schnek::Grid<double, 1, HuertoGridChecker> Grid1d;
typedef std::shared_ptr<Grid1d> pGrid1d;
typedef std::unique_ptr<Grid1d> uGrid1d;
typedef std::reference_wrapper<Grid1d> rGrid1d;
typedef schnek::Field<double, 1, HuertoGridChecker> Field1d;
typedef std::shared_ptr<Field1d> pField1d;
typedef std::unique_ptr<Field1d> uField1d;
typedef std::reference_wrapper<Field1d> rField1d;
typedef schnek::Range<ptrdiff_t, 1> Range1d;
typedef schnek::Range<double, 1> Domain1d;
typedef schnek::Array<bool, 1> Stagger1d;

typedef schnek::Array<ptrdiff_t, 2> Index2d;
typedef schnek::Array<size_t, 2> Size2d;
typedef schnek::Array<double, 2> Vector2d;
typedef schnek::Grid<double, 2, HuertoGridChecker> Grid2d;
typedef std::shared_ptr<Grid2d> pGrid2d;
typedef std::unique_ptr<Grid2d> uGrid2d;
typedef std::reference_wrapper<Grid2d> rGrid2d;
typedef schnek::Field<double, 2, HuertoGridChecker> Field2d;
typedef std::shared_ptr<Field2d> pField2d;
typedef std::unique_ptr<Field2d> uField2d;
typedef std::reference_wrapper<Field2d> rField2d;
typedef schnek::Range<ptrdiff_t, 2> Range2d;
typedef schnek::Range<double, 2> Domain2d;
typedef schnek::Array<bool, 2> Stagger2d;

typedef schnek::Array<ptrdiff_t, 3> Index3d;
typedef schnek::Array<size_t, 3> Size3d;
typedef schnek::Array<double, 3> Vector3d;
typedef schnek::Grid<double, 3, HuertoGridChecker> Grid3d;
typedef std::shared_ptr<Grid3d> pGrid3d;
typedef std::unique_ptr<Grid3d> uGrid3d;
typedef std::reference_wrapper<Grid3d> rGrid3d;
typedef schnek::Field<double, 3, HuertoGridChecker> Field3d;
typedef std::shared_ptr<Field3d> pField3d;
typedef std::unique_ptr<Field3d> uField3d;
typedef std::reference_wrapper<Field3d> rField3d;
typedef schnek::Range<ptrdiff_t, 3> Range3d;
typedef schnek::Range<double, 3> Domain3d;
typedef schnek::Array<bool, 3> Stagger3d;

typedef schnek::RangeCIterationPolicy<DIMENSION> FieldIterator;
typedef schnek::RangeCIterationPolicy<1> Field1dIterator;
typedef schnek::RangeCIterationPolicy<2> Field2dIterator;
typedef schnek::RangeCIterationPolicy<3> Field3dIterator;

#ifndef HUERTO_DECOMPOSITION
#define HUERTO_DECOMPOSITION

typedef schnek::MpiCartesianDomainDecomposition<DIMENSION> HuertoDecomposition;

#endif

#define LOG_ERROR 0
#define LOG_WARN 1
#define LOG_DEBUG 2
#define LOG_TRACE 3


#endif // HUERTO_TYPES_H
