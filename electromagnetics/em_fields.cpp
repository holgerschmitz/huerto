/*
 * em_fields.cpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */

#include "em_fields.hpp"
#include "../types.hpp"

#include <schnek/util/logger.hpp>

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include <string>
#include <iostream>

void EMFields::registerData()
{
  SCHNEK_TRACE_ENTER_FUNCTION(LOG_DEBUG)
  pEx = boost::make_shared<Field>();
  pEy = boost::make_shared<Field>();
  pEz = boost::make_shared<Field>();

  pBx = boost::make_shared<Field>();
  pBy = boost::make_shared<Field>();
  pBz = boost::make_shared<Field>();

  addData("Ex", pEx);
  addData("Ey", pEy);
  addData("Ez", pEz);

  addData("Bx", pBx);
  addData("By", pBy);
  addData("Bz", pBz);
}

// end of main
