/*
 * table_lookup.hpp
 *
 *  Created on: 22 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_TABLES_TABLE_LOOKUP_HPP_
#define HUERTO_TABLES_TABLE_LOOKUP_HPP_

#include "../types.hpp"
#include "../io/table_data_source.hpp"

#include <boost/scoped_ptr.hpp>

/**
 * A block that reads table data from a text file and allows lookup tables
 * to be created from the data.
 */
class TableBlock : public schnek::Block, public FileTableDataReader {
  private:
    /**
     * A pointer to the TableBlock object
     */
    TableBlock *self;
  protected:

    /**
     * This will register the TableBlock in the Block data repository.
     * The lookup name is constructed by the Block name prefixed with "TABLE_"
     */
    void registerData();

    /**
     * Registers the parameters to be read from the setup file
     */
    void initParameters(schnek::BlockParameters &parameters);

    /**
     * Calls FileTableDataReader::preInit() to read the data from the file
     */
    void preInit();
};

typedef boost::shared_ptr<TableBlock> pTableBlock;

/**
 * A table lookup that can be constructed from a table block
 *
 * The table uses linear interpolation.
 */
class TableLookup
{
  private:
    /**
     * A vector containing the grids that have been dynamically created within
     * this object. Ensuring RAII
     */
    std::vector<rGrid1d> ownedData;
  protected:

    /**
     * The values on the x-axis in the lookup table
     */
    Grid1d *xValues;

    /**
     * The y-values to be looked up in the lookup table
     */
    Grid1d *yValues;

    /**
     * The cumulative distribution of y-values normalised to 1
     *
     * This grid is only generated if `cumulative` was set to `true` in one of
     * the init() functions
     */
    std::unique_ptr<Grid1d> yCumulative;
  protected:
    /**
     * Create and initialise the cumulative distribution
     */
    void initCumulative();
  public:

    /**
     * Create a new lookup table from a table block
     *
     * @param table   the table block from which to obtain the data
     * @param xIndex  the column index of the x-axis
     * @param yIndex  the column index of the y-values
     * @param cumulative  if set to `true` a cumulative distribution will be
     *                    generated
     */
    void init(TableBlock &table, int xIndex, int yIndex, bool cumulative);

    /**
     * Create a new lookup table by combining two TableLookup objects
     *
     * The y-values are calculated by the weigthed average of the y-values of
     * the two tables, table A and table B.
     *
     * Table A provides the x-axis and table B must cover at least the same
     * index range as table A.
     *
     * @param tableA   the lookup table A
     * @param weightA  the weight of the values in lookup table A
     * @param tableB   the lookup table B
     * @param weightB  the weight of the values in lookup table B
     * @param cumulative  if set to `true` a cumulative distribution will be
     *                    generated
     */
    void init(TableLookup &tableA, double weightA, TableLookup &tableB, double weightB, bool cumulative);

    /**
     * Return the interpolated y-value from an x-value
     *
     * If x falls between points defined in the table, linear interpolation is
     * used.
     *
     * If x falls outside the range of the x-axis, the corresponding limiting
     * y-values are returned.
     *
     * @param x  the x-value to look up
     * @return   the interpolated y-value
     */
    double interpolate(double x) const;

    /**
     * Return a random x-value taken from the probability distribution given
     * by the y-values.
     *
     * This function requires the cumulative distribution function. This is
     * created by setting `cumulative` to `true` when calling on of the init()
     * functions.
     */
    double randomDist();
};

/**
 * A 2d table lookup that can be constructed from a table block
 *
 * The table uses linear interpolation.
 */
class TableLookup2d
{
  protected:

    /**
     * The values on the x-axis in the lookup table
     */
    Grid1d xValues;

    /**
     * The values on the y-axis in the lookup table
     */
    Grid1d yValues;

    /**
     * The 2d grid of values to be looked up in the lookup table
     */
    Grid2d table;
  public:

    /**
     * Create a 2d lookup table from a table block
     *
     * The lookup table assumes the following layout of data in the table block.
     *
     * * The first column in the table contains the x-axis values.
     * * The first row in the table contains the y-axis values.
     * * The first number in the file is unused.
     * * The remainder of the data contains the table values
     *
     * @param tableBlock  the table block
     */
    void init(TableBlock &tableBlock);

    /**
     * Return the interpolated value from a pair of (x,y)-values
     *
     * If x and y falls between points defined in the table, linear interpolation is
     * used.
     *
     * If x or y falls outside the range of the axes, the corresponding limiting
     * table values are returned.
     *
     * @param x  the x-value to look up
     * @param y  the y-value to look up
     * @return   the interpolated table value
     */
    double interpolate(double x, double y) const;
};

#endif /* HUERTO_TABLES_TABLE_LOOKUP_HPP_ */
