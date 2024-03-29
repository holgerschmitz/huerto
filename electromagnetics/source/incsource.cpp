/*
 * incsource.t
 *
 *  Created on: 30 Nov 2020
 *      Author: Holger Schmitz
 */

#include "incsource.hpp"

#include "../fieldsolver.hpp"

#include <vector>
#include <sstream>

//===============================================================
//==========  IncidentSource
//===============================================================

void IncidentSource::initCurrents(CurrentContainer &container) {
  if (needCurrent(east)) {
    container.addCurrent(makeECurrent(distance, east));
    container.addMagCurrent(makeHCurrent(distance, east));
  }

  if (needCurrent(west)) {
    container.addCurrent(makeECurrent(distance, west));
    container.addMagCurrent(makeHCurrent(distance, west));
  }

#ifndef HUERTO_ONE_DIM
  if (needCurrent(north)) {
    container.addCurrent(makeECurrent(distance, north));
    container.addMagCurrent(makeHCurrent(distance, north));
  }

  if (needCurrent(south)) {
    container.addCurrent(makeECurrent(distance, south));
    container.addMagCurrent(makeHCurrent(distance, south));
  }
#endif

#ifdef HUERTO_THREE_DIM
  if (needCurrent(up)) {
    container.addCurrent(makeECurrent(distance, up));
    container.addMagCurrent(makeHCurrent(distance, up));
  }

  if (needCurrent(down)) {
    container.addCurrent(makeECurrent(distance, down));
    container.addMagCurrent(makeHCurrent(distance, down));
  }
#endif
}

void IncidentSource::initParameters(schnek::BlockParameters &blockPars) {
  CurrentBlock::initParameters(blockPars);

  blockPars.addParameter("d", &this->distance, 15);
}

//===============================================================
//==========  IncidentSourceCurrent
//===============================================================

IncidentSourceCurrent::IncidentSourceCurrent(int distance, Direction dir, bool isH, SimulationContext &context)
  : distance(distance), dir(dir), isH(isH), context(context)
{
  lowOffset = 0;
  highOffset = 0;

  dx = context.getDx();
  dt = context.getDt();

  switch (dir)
  {
    case east:
    case west:  dim = 0;
                transverse1 = 1;
                transverse2 = 2;
                JT[0] = Jy;
                JT[1] = Jz;
                dN = dx[0];
                break;
#ifndef HUERTO_ONE_DIM
    case north:
    case south: dim = 1;
                transverse1 = 2;
                transverse2 = 0;
                JT[0] = Jz;
                JT[1] = Jx;
                dN = dx[1];
                break;
#endif
#ifdef HUERTO_THREE_DIM
    case up:
    case down:  dim = 2;
                transverse1 = 0;
                transverse2 = 1;
                JT[0] = Jx;
                JT[1] = Jy;
                dN = dx[2];
                break;
#endif
  }

#ifdef HUERTO_ONE_DIM
  reverse = (dir==east);
#endif

#ifdef HUERTO_TWO_DIM
  reverse = ( (dir==east) || (dir==north) );
#endif

#ifdef HUERTO_THREE_DIM
  reverse = ( (dir==east) || (dir==north) || (dir==up) );
#endif
}

