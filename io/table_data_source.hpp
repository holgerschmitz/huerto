/*
 * table_data_source.hpp
 *
 *  Created on: 22 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_IO_TABLE_DATA_SOURCE_HPP_
#define HUERTO_IO_TABLE_DATA_SOURCE_HPP_

#include "../types.hpp"

#include <schnek/variables/block.hpp>
#include <schnek/grid.hpp>

#include <string>
#include <vector>

/**
 * A source for tabular data
 *
 * This class is intended to be used within a Block's lifecycle
 */
class TableDataSource
{
  public:
    /**
     * C++ virtual destructor
     */
    virtual ~TableDataSource() {}

    /**
     * Init any parameters to be read from the setup file
     */
    virtual void initParameters(schnek::BlockParameters &) {};

    /**
     * This is where data should be initialised
     */
    virtual void preInit() {};

    /**
     * Get the value array of column n
     */
    virtual Grid1d &getValues(size_t n) const = 0;

    /**
     * Get the number of columns
     */
    virtual size_t getColumns() const = 0;
};

/**
 * Used to read tabular data from a file
 */
class FileTableDataReader : public TableDataSource
{
  protected:
    /**
     * File name to read the data from
     */
    std::string fileName;

    /**
     * Field separator inside an ASCII file
     *
     * A set of characters that are regarded as field seperators.
     */
    std::string seperator;

    /**
     * A vector of value arrays. Each vector element corresponds to a column in
     * the file
     */
    std::vector<pGrid1d> values;
  public:

    /**
     * Registers the parameters to be read from the setup file
     */
    void initParameters(schnek::BlockParameters &parameters) override;

    /**
     * Loads the data from the file and stores it in the values member
     */
    void preInit() override;

    /**
     * Get the value array of column n
     */
    Grid1d &getValues(size_t n) const override { return *(values[n]); }

    /**
     * Get the number of columns
     */
    size_t getColumns() const override { return values.size(); }
};



#endif /* SOLENA_IO_DATA_READER_HPP_ */
