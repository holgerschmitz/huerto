/*
 * current.cpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "current.hpp"
#include "../util/field_util.hpp"

#include <memory>

void CurrentContainer::addCurrent(pCurrent current)
{
  current->init();
  this->currents.push_back(current);
}

void CurrentContainer::addMagCurrent(pCurrent current)
{
  current->init();
  this->magCurrents.push_back(current);
}

void CurrentContainer::sumCurrents()
{
  auto &decomposition = context->getDecomposition();

  setField<Field>(decomposition, {Jx}, 0.0);
  setField<Field>(decomposition, {Jy}, 0.0);
  setField<Field>(decomposition, {Jz}, 0.0);

  for (pCurrent current: this->currents)
  {
    addToField<Field, Grid>(decomposition, Jx, current->getJx());
    addToField<Field, Grid>(decomposition, Jy, current->getJy());
    addToField<Field, Grid>(decomposition, Jz, current->getJz());
  }
}

void CurrentContainer::sumMagCurrents()
{
  auto &decomposition = context->getDecomposition();

  setField<Field>(decomposition, {Mx}, 0.0);
  setField<Field>(decomposition, {My}, 0.0);
  setField<Field>(decomposition, {Mz}, 0.0);

  for (pCurrent current: this->magCurrents)
  {
    addToField<Field, Grid>(decomposition, Mx, current->getJx());
    addToField<Field, Grid>(decomposition, My, current->getJy());
    addToField<Field, Grid>(decomposition, Mz, current->getJz());
  }
}

void CurrentContainer::init(SimulationContext &context)
{
  this->context = &context;
  auto &decomposition = context.getDecomposition();
  Jx = decomposition.registerField(schnek::GridFactory<Field>{exStaggerYee, 2});
  Jy = decomposition.registerField(schnek::GridFactory<Field>{eyStaggerYee, 2});
  Jz = decomposition.registerField(schnek::GridFactory<Field>{ezStaggerYee, 2});

  Mx = decomposition.registerField(schnek::GridFactory<Field>{bxStaggerYee, 2});
  My = decomposition.registerField(schnek::GridFactory<Field>{byStaggerYee, 2});
  Mz = decomposition.registerField(schnek::GridFactory<Field>{bzStaggerYee, 2});
}


//===============================================================
//==========  CurrentBlock
//===============================================================

void CurrentBlock::preInit()
{
  schnek::ChildBlock<CurrentBlock>::preInit();
  SimulationEntity::init(this);
}

