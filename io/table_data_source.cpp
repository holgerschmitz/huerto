/*
 * table_data_source.cpp
 *
 *  Created on: 22 Jan 2020
 *      Author: Holger Schmitz
 */

#include "table_data_source.hpp"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>

void FileTableDataReader::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addParameter("file", &fileName);
  parameters.addParameter("seperator", &seperator, std::string(" "));
}

void FileTableDataReader::preInit()
{
  std::ifstream input(fileName);
  if (!input)
  {
    std::cerr << "Could not find table file '" << fileName << "'!" << std::endl;
    exit(-1);
  }

  std::vector<std::shared_ptr<std::vector<double>>> readValues;

  size_t count = 0;
  std::string line;
  while (std::getline(input,line)) {

    std::vector<std::string> splitVec;
    boost::trim(line);
    boost::split(splitVec, line, boost::is_any_of(seperator), boost::algorithm::token_compress_on);

    while (splitVec.size() > readValues.size())
    {
      readValues.push_back(std::make_shared<std::vector<double>>(count, 0.0));
    }

    for (size_t i=0; i<splitVec.size(); ++i)
    {
      readValues[i]->push_back(boost::lexical_cast<double>(splitVec[i]));
    }

    for (size_t i=splitVec.size(); i<readValues.size(); ++i)
    {
      readValues[i]->push_back(0.0);
    }

    ++count;
  }

  for (size_t i=0; i<readValues.size(); ++i)
  {
    pGrid1d column(std::make_shared<Grid1d>(Grid1d::IndexType(count)));
    std::vector<double> &readColumn = *readValues[i];
    for (size_t row=0; row<count; ++row)
    {
      (*column)(row) = readColumn[row];
    }
    values.push_back(column);
  }
}

