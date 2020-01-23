/*
 * test_table_lookup.cpp
 *
 *  Created on: 23 Jan 2020
 *      Author: Holger Schmitz
 */


/*
 * test_table_data_source.cpp
 *
 *  Created on: 22 Jan 2020
 *      Author: Holger Schmitz
 */

#include "../test_types.hpp"

#include "../fixtures/block_runner.hpp"
#include "../../constants.hpp"
#include "../../tables/table_lookup.hpp"

#include <cmath>

BOOST_AUTO_TEST_SUITE( maths )

BOOST_AUTO_TEST_SUITE( tables )

BOOST_FIXTURE_TEST_CASE( test_read_TableBlock, SingleBlockRunner )
{
  pTableBlock block = createBlock<TableBlock>("TableBlock", "tests/io/test_table_data_source_1.setup");

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

  TableBlock *retrievedBlock;
  BOOST_CHECK_NO_THROW(
      block->retrieveData("TABLE_TableBlock", retrievedBlock);
  );
  BOOST_CHECK(retrievedBlock);
}

BOOST_FIXTURE_TEST_CASE( test_TableLookup, SingleBlockRunner )
{
  pTableBlock block = createBlock<TableBlock>("TableBlock", "tests/io/test_table_data_source_1.setup");

  TableLookup lookup;
  lookup.init()
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()


