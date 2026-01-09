/*
 * current.cpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "current.hpp"
#include "../util/fieldsum.hpp"

#include <memory>

void CurrentContainer::addCurrent(pCurrent current)
{
  current->init();
  if (current->isValid()) {
    this->currents.push_back(current);
  }
}

void CurrentContainer::addMagCurrent(pCurrent current)
{
  current->init();
  if (current->isValid()) {
    this->magCurrents.push_back(current);
  }
}

void CurrentContainer::sumCurrents()
{
  Jx = 0;
  Jy = 0;
  Jz = 0;

  for (pCurrent current: this->currents)
  {
    Grid jx = current->getJx();
    Grid jy = current->getJy();
    Grid jz = current->getJz();

    FieldSum<Field, Grid> sumJx{Jx, jx};
    FieldSum<Field, Grid> sumJy{Jy, jy};
    FieldSum<Field, Grid> sumJz{Jz, jz};

    FieldIterator::forEach(jx.getRange(), sumJx);
    FieldIterator::forEach(jy.getRange(), sumJy);
    FieldIterator::forEach(jz.getRange(), sumJz);
  }
}

void CurrentContainer::sumMagCurrents()
{
  Mx = 0;
  My = 0;
  Mz = 0;

  for (pCurrent current: this->magCurrents)
  {
    Grid jx = current->getJx();
    Grid jy = current->getJy();
    Grid jz = current->getJz();

    FieldSum<Field, Grid> sumJx{Mx, jx};
    FieldSum<Field, Grid> sumJy{My, jy};
    FieldSum<Field, Grid> sumJz{Mz, jz};

    FieldIterator::forEach(jx.getRange(), sumJx);
    FieldIterator::forEach(jy.getRange(), sumJy);
    FieldIterator::forEach(jz.getRange(), sumJz);
  }
}

void CurrentContainer::init(SimulationContext &context)
{

  schnek::DomainSubdivision<Field> &subdivision = context.getSubdivision();
  Index lowIn = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  Domain domainSize = subdivision.getInnerExtent(context.getSize());
  Jx.resize(lowIn, highIn, domainSize, exStaggerYee, 2);
  Jy.resize(lowIn, highIn, domainSize, eyStaggerYee, 2);
  Jz.resize(lowIn, highIn, domainSize, ezStaggerYee, 2);

  Mx.resize(lowIn, highIn, domainSize, bxStaggerYee, 2);
  My.resize(lowIn, highIn, domainSize, byStaggerYee, 2);
  Mz.resize(lowIn, highIn, domainSize, bzStaggerYee, 2);
}


//===============================================================
//==========  CurrentBlock
//===============================================================

void CurrentBlock::preInit()
{
  schnek::ChildBlock<CurrentBlock>::preInit();
  SimulationEntity::init(this);
}

