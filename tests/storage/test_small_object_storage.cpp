/*
 * test_small_object_storage.cpp
 *
 * Created on: 03 Mar 2023
 * Author: Holger Schmitz
 * Email: holger@notjustphysics.com
 */

#include "../../storage/small_object_storage.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

struct TestType 
{
  double x, y, z;
};

struct TestSmallObjectStorageTraits {
    static const size_t chunkSize = 1000;
};

typedef SmallObjectStorage<TestType, TestSmallObjectStorageTraits> TestStorage;

boost::random::mt19937 rGen;

const size_t NRepeat = 100;
const size_t PartSize = 50001;

BOOST_AUTO_TEST_SUITE( storage )

BOOST_AUTO_TEST_SUITE( small_object_storage )

std::pair<double, size_t> addElements(TestStorage &storage)
{
  boost::random::uniform_int_distribution<> randSize(1,PartSize);
  boost::random::uniform_real_distribution<> randValue(-1.0,1.0);

  size_t numElements = randSize(rGen);

  double sum = 0.0;
  for (size_t n = 0; n<numElements; ++n)
  {
    TestType &t = storage.addElement();

    double x = randValue(rGen);
    double y = randValue(rGen);
    double z = randValue(rGen);
    t.x = x;
    t.y = y;
    t.z = z;
    sum += x + y + z;
  }
  return {sum, numElements};
}

std::pair<double, size_t> removeElements(TestStorage &storage)
{
  boost::random::uniform_real_distribution<> mc(0.0,1.0);
  auto it = storage.begin();
  size_t removed = 0;
  double sum = 0.0;

  while (it!=storage.end())
  {
    if (mc(rGen) < 0.5)
    {
      sum += it->x + it->y + it->z;
      ++removed;

      it = storage.removeElement(it);
    } else
      ++it;
  }

  return {sum, removed};
}

double updateElements(TestStorage &storage)
{
  boost::random::uniform_real_distribution<> randValue(-1.0,1.0);
  double sum = 0.0;

  for (TestType &t : storage)
  {
    double x = randValue(rGen);
    double y = randValue(rGen);
    double z = randValue(rGen);
    t.x = x;
    t.y = y;
    t.z = z;
    sum += x + y + z;
  }

  return sum;
}

double readElements(TestStorage &storage)
{
  double sum = 0.0;

  for (TestType &t : storage)
  {
    sum += t.x + t.y + t.z;
  }

  return sum;
}


BOOST_AUTO_TEST_CASE( add_items )
{
  for (size_t i=0; i<NRepeat; ++i)
  {
    TestStorage storage;
    auto result = addElements(storage);
    BOOST_CHECK_EQUAL(storage.getCount(), result.second);
    double sum = readElements(storage);
    BOOST_CHECK_CLOSE(sum, result.first, 1e-8);
  }
}

BOOST_AUTO_TEST_CASE( update_items )
{
  for (size_t i=0; i<NRepeat; ++i)
  {
    TestStorage storage;
    addElements(storage);
    double result = updateElements(storage);
    double sum = readElements(storage);
    BOOST_CHECK_CLOSE(sum, result, 1e-8);
  }
}


BOOST_AUTO_TEST_CASE( remove_items )
{
  for (size_t i=0; i<NRepeat; ++i)
  {
    TestStorage storage;
    auto before = addElements(storage);
    auto removed = removeElements(storage);
    double after = readElements(storage);
    BOOST_CHECK_CLOSE(before.first - removed.first, after, 1e-8);
    BOOST_CHECK_EQUAL(before.second - removed.second, storage.getCount());
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()