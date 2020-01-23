/*
 * table_lookup.hpp
 *
 *  Created on: 22 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef SCHNAR_TABLES_TABLE_LOOKUP_HPP_
#define SCHNAR_TABLES_TABLE_LOOKUP_HPP_

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

class TableLookup
{
  private:
    TableLookup *self;
    std::vector<std::unique_ptr<Grid1d>> ownedData;
  protected:
    Grid1d *xValues;
    Grid1d *yValues;

    const std::unique_ptr<Grid1d> yCumulative;
  protected:
    void initCumulative();
  public:
    void init(TableBlock &table, int xIndex, int yIndex, bool cumulative);
    void init(TableLookup &tableA, double weightA, TableLookup &tableB, double weightB, bool cumulative);
    double interpolate(double x) const;
    double randomDist();
};

class TableLookup2d
{
  protected:
    Grid1d xValues;
    Grid1d yValues;
    Grid2d table;
  public:
    void init(TableBlock &tableBlock);
    double interpolate(double x, double y) const;
};

#endif /* SCHNAR_TABLES_TABLE_LOOKUP_HPP_ */
