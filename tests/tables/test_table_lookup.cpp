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

#include <boost/progress.hpp>

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

BOOST_FIXTURE_TEST_CASE( test_TableLookup_interpolate, SingleBlockRunner )
{
  pTableBlock block = createBlock<TableBlock>("TableBlock", "tests/tables/test_table_lookup_1.setup");

  TableLookup lookup;
  lookup.init(*block, 0, 1, false);

  // test support points
  for (int i = 0; i <= 100; i++)
  {
    double x = TWO_PI * i / 100.0;
    BOOST_CHECK_CLOSE(lookup.interpolate(x) + 2.0, sin(x) + 2.0, 5e-3);
  }

  // test intervale between support points
  for (int i = 0; i <= 400; i++)
  {
    double x = TWO_PI * i / 400.0;
    BOOST_CHECK_CLOSE(lookup.interpolate(x) + 2.0, sin(x) + 2.0, 5e-2);
  }

}

BOOST_FIXTURE_TEST_CASE( test_TableLookup_unordered, SingleBlockRunner )
{
  pTableBlock block = createBlock<TableBlock>("TableBlock", "tests/tables/test_table_lookup_1.setup");

  TableLookup lookup;

  BOOST_CHECK_THROW(lookup.init(*block, 1, 0, false), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE( test_TableLookup_cumulative, SingleBlockRunner )
{
  int N = 10000000;
  int H = 50;

  pTableBlock block = createBlock<TableBlock>("TableBlock", "tests/tables/test_table_lookup_1.setup");

  TableLookup lookup;
  lookup.init(*block, 0, 2, true);

  Grid1d hist(Index1d(0), Index1d(H));
  hist = 0.0;
  boost::progress_display show_progress(N);

  for (int i = 0; i < N; i++)
  {
    double x = lookup.randomDist();
    int index = (int) (H * x / TWO_PI);
    hist(index) += 1.0;

    ++show_progress;
  }

  for (int i = 0; i < H; i++)
  {
    double xm = TWO_PI * (i) / (double) H;
    double xp = TWO_PI * (i + 1.0) / (double) H;
    double Im = 0.5 * xm - 0.25 * sin(2 * xm);
    double Ip = 0.5 * xp - 0.25 * sin(2 * xp);

    BOOST_CHECK_CLOSE(PI * hist(i) / N + 1.0, Ip - Im + 1.0, 5e-2);
  }
}


BOOST_FIXTURE_TEST_CASE( test_TableLookup_cumulative_negative, SingleBlockRunner )
{
  pTableBlock block = createBlock<TableBlock>("TableBlock", "tests/tables/test_table_lookup_1.setup");

  TableLookup lookup;
  BOOST_CHECK_THROW(lookup.init(*block, 0, 1, true), std::runtime_error);
}


BOOST_FIXTURE_TEST_CASE( test_TableLookup2d_interpolate, SingleBlockRunner )
{
  pTableBlock block = createBlock<TableBlock>("TableBlock", "tests/tables/test_table_lookup_2.setup");

  TableLookup2d lookup;
  lookup.init(*block);
  boost::progress_display show_progress(617);

  // test support points
  for (int i = 0; i <= 123; i++)
  {
    double x = i / 246.0;
    for (int j = 0; j <= 100; j++)
    {
      double y = j / 100.0;
      BOOST_CHECK_CLOSE(lookup.interpolate(x, y) + 2.0, sin(x) * cos(y) + 2.0, 5e-3);
    }
    ++show_progress;
  }

  // test intervale between support points
  for (int i = 0; i <= 492; i++)
  {
    double x = i / 984.0;
    for (int j = 0; j <= 400; j++)
    {
      double y = j / 400.0;
      BOOST_CHECK_CLOSE(lookup.interpolate(x, y) + 2.0, sin(x) * cos(y) + 2.0, 5e-2);
    }
    ++show_progress;
  }

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()


