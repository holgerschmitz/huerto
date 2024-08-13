/*
 * current.cpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "current.hpp"

#include <memory>

void CurrentContainer::addCurrent(pCurrent current)
{
  current->init();
  if (current->isValid())
    this->currents.push_back(current);
}

void CurrentContainer::addMagCurrent(pCurrent current)
{
  current->init();
  if (current->isValid())
    this->magCurrents.push_back(current);
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

    Index low = jx.getLo();
    Index high = jx.getHi();

#ifdef HUERTO_ONE_DIM
    for (int i=low[0]; i<=high[0]; ++i) {
      Jx(i) += jx(i);
      Jx(i) += jy(i);
      Jx(i) += jz(i);
    }
#endif

#ifdef HUERTO_TWO_DIM
    for (int i=low[0]; i<=high[0]; ++i) {
      for (int j=low[1]; j<=high[1]; ++j) {
        Jx(i,j) += jx(i,j);
        Jy(i,j) += jy(i,j);
        Jz(i,j) += jz(i,j);
      }
    }
#endif

#ifdef HUERTO_THREE_DIM
    for (int i=low[0]; i<=high[0]; ++i) {
      for (int j=low[1]; j<=high[1]; ++j) {
        for (int k=low[2]; k<=high[2]; ++k) {
          Jx(i,j,k) += jx(i,j,k);
          Jy(i,j,k) += jy(i,j,k);
          Jz(i,j,k) += jz(i,j,k);
        }
      }
    }
#endif

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

    Index low = jx.getLo();
    Index high = jx.getHi();

#ifdef HUERTO_ONE_DIM
    for (int i=low[0]; i<=high[0]; ++i){
      Mx(i) += jx(i);
      My(i) += jy(i);
      Mz(i) += jz(i);
    }
#endif

#ifdef HUERTO_TWO_DIM
    for (int i=low[0]; i<=high[0]; ++i) {
      for (int j=low[1]; j<=high[1]; ++j) {
        Mx(i,j) += jx(i,j);
        My(i,j) += jy(i,j);
        Mz(i,j) += jz(i,j);
      }
    }
#endif

#ifdef HUERTO_THREE_DIM
    for (int i=low[0]; i<=high[0]; ++i) {
      for (int j=low[1]; j<=high[1]; ++j) {
        for (int k=low[2]; k<=high[2]; ++k) {
          Mx(i,j,k) += jx(i,j,k);
          My(i,j,k) += jy(i,j,k);
          Mz(i,j,k) += jz(i,j,k);
        }
      }
    }
#endif
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

