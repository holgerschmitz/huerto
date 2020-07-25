/*
 * test_table_data_source.cpp
 *
 *  Created on: 22 Jan 2020
 *      Author: Holger Schmitz
 */

#include "../test_types.hpp"

#include "../fixtures/block_runner.hpp"
#include "../../constants.hpp"
#include "../../io/table_data_source.hpp"

#include <cmath>

class FileTableDataReaderBlock : public schnek::Block, public FileTableDataReader
{
    void initParameters(schnek::BlockParameters &parameters) override
    {
      schnek::Block::initParameters(parameters);
      FileTableDataReader::initParameters(parameters);
    }

    void preInit() override
    {
      schnek::Block::preInit();
      FileTableDataReader::preInit();
    }
};

typedef boost::shared_ptr<FileTableDataReaderBlock> pFileTableDataReaderBlock;

BOOST_AUTO_TEST_SUITE( maths )

BOOST_AUTO_TEST_SUITE( io )

BOOST_FIXTURE_TEST_CASE( readTable, SingleBlockRunner )
{
  pFileTableDataReaderBlock block = createBlock<FileTableDataReaderBlock>("TableRader", "tests/io/test_table_data_source_1.setup");

  BOOST_CHECK_EQUAL(block->getColumns(), 3);

  Grid1d &col0 = block->getValues(0);
  Grid1d &col1 = block->getValues(1);
  Grid1d &col2 = block->getValues(2);

  BOOST_CHECK_EQUAL(col0.getLo(0), 0);
  BOOST_CHECK_EQUAL(col0.getHi(0), 100);
  BOOST_CHECK_EQUAL(col1.getLo(0), 0);
  BOOST_CHECK_EQUAL(col1.getHi(0), 100);
  BOOST_CHECK_EQUAL(col2.getLo(0), 0);
  BOOST_CHECK_EQUAL(col2.getHi(0), 100);

  for (int i=0; i<=100; i++)
  {
    double x = TWO_PI*i/100.;
    BOOST_CHECK_CLOSE(col0(i)+1, x+1, 1e-3);
    BOOST_CHECK_CLOSE(col1(i)+2, sin(x)+2, 1e-3);
    BOOST_CHECK_CLOSE(col2(i)+2, cos(x)+2, 1e-3);
  }

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
