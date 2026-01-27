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
  auto gridContext = decomposition.getGridContext({Jx, Jy, Jz});
  gridContext.forEach([&](Range& /* range */, Field &Jx, Field &Jy, Field &Jz) {
    Jx = 0;
    Jy = 0;
    Jz = 0;
  });

  for (pCurrent current: this->currents)
  {
    auto currentGridContext = decomposition.getGridContext({current->getJx(), current->getJy(), current->getJz(), Jx, Jy, Jz});

    gridContext.forEach([&](Range& range, Field &jx, Field &jy, Field &jz, Field &Jx, Field &Jy, Field &Jz) {
        FieldSum<Field, Grid> sumJx{Jx, jx};
        FieldSum<Field, Grid> sumJy{Jy, jy};
        FieldSum<Field, Grid> sumJz{Jz, jz};

        FieldIterator::forEach(range, sumJx);
        FieldIterator::forEach(range, sumJy);
        FieldIterator::forEach(range, sumJz);
    });
  }
}

void CurrentContainer::sumMagCurrents()
{
  auto &decomposition = context->getDecomposition();
  auto gridContext = decomposition.getGridContext({Mx, My, Mz});
  gridContext.forEach([&](Range& /* range */, Field &Jx, Field &Jy, Field &Jz) {
    Jx = 0;
    Jy = 0;
    Jz = 0;
  });

  for (pCurrent current: this->magCurrents)
  {
    auto currentGridContext = decomposition.getGridContext({current->getJx(), current->getJy(), current->getJz(), Mx, My, Mz});

    gridContext.forEach([&](Range& range, Field &jx, Field &jy, Field &jz, Field &Jx, Field &Jy, Field &Jz) {
        FieldSum<Field, Grid> sumJx{Jx, jx};
        FieldSum<Field, Grid> sumJy{Jy, jy};
        FieldSum<Field, Grid> sumJz{Jz, jz};

        FieldIterator::forEach(range, sumJx);
        FieldIterator::forEach(range, sumJy);
        FieldIterator::forEach(range, sumJz);
    });
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

