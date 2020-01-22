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

class TableBlock : public schnek::Block, public FileTableDataReader {
  private:
    TableBlock *self;
  protected:
    void registerData();
    /**
     * Registers the parameters to be read from the setup file
     */
    void initParameters(schnek::BlockParameters &parameters);
    void preInit();
};

class TableLookup
{
  private:
    TableLookup *self;
    std::vector<pGrid1d> ownedData;
  protected:
    Grid1d *xValues;
    Grid1d *yValues;

    Grid1d *yCumulative;
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
