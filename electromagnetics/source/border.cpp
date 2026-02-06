/*
 * border.cpp
 *
 *  Created on: 30 Nov 2020
 *      Author: Holger Schmitz
 */

#include "border.hpp"

#include <schnek/util/logger.hpp>

#ifdef HUERTO_ONE_DIM
bool getBorderExtent(Direction dir,
                     int thickness,
                     int distance,
                     Index &blow,
                     Index &bhigh,
                     bool isH,
                     SimulationContext &context,
                     bool /* restricted */)
{
  int distanceLow = distance - (isH?1:0);

  auto &decomposition = context.getDecomposition();
  Range grange = decomposition.getGlobalRange();
  Index glow  = grange.getLo();
  Index ghigh = grange.getHi();

  switch (dir)
  {
    case west:
      bhigh[0] = glow[0]+thickness-1+distanceLow;
      blow[0] = glow[0]+distanceLow;

      break;
    case east:
      blow[0] = ghigh[0]-thickness+2-distance;
      bhigh[0] = ghigh[0]-distance+1;

      break;
  }

  return true;
}
#endif

#ifdef HUERTO_TWO_DIM
bool getBorderExtent(Direction dir,
                     int thickness,
                     int distance,
                     Index &blow,
                     Index &bhigh,
                     bool isH,
                     SimulationContext &context,
                     bool restricted)
{
  int distanceLow = distance - (isH?1:0);
  int distanceHigh = distance - 0; //(isH?0:1);

  auto &decomposition = context.getDecomposition();
  Range grange = decomposition.getGlobalRange();
  Index glow  = grange.getLo();
  Index ghigh = grange.getHi();

  blow[0] = glow[0];
  blow[1] = glow[1];

  bhigh[0] = ghigh[0];
  bhigh[1] = ghigh[1];

  Index coords(0,1);

//  int borderfit = 0;

  switch (dir)
  {
    case west:
    case east:
      coords = Index(0,1);
      break;
    case south:
    case north:
      coords = Index(1,0);
      break;
  }

  int normal = coords[0];
  int t1 = coords[1];

  switch (dir)
  {
    case west:
    case south:
      bhigh[normal] = glow[normal]+thickness-1+distanceLow;
      blow[normal] = glow[normal]+distanceLow;

      break;
    case east:
    case north:
      blow[normal] = ghigh[normal]-thickness+2-distance;
      bhigh[normal] = ghigh[normal]-distance+1;

      break;
  }

  if (restricted) {
    bhigh[t1] = ghigh[t1]-distanceHigh;
    blow[t1] = glow[t1]+distance;
  }

  return true;
}
#endif

#ifdef HUERTO_THREE_DIM
bool getBorderExtent(Direction dir,
                     int thickness,
                     int distance,
                     Index &blow,
                     Index &bhigh,
                     bool isH,
                     SimulationContext &context,
                     bool restricted)
{
  int distanceLow = distance - (isH?1:0);
  int distanceHigh = distance - 0; //(isH?0:1);

  auto &decomposition = context.getDecomposition();
  Range grange = decomposition.getGlobalRange();
  Index glow  = grange.getLo();
  Index ghigh = grange.getHi();

  blow[0] = glow[0];
  blow[1] = glow[1];
  blow[2] = glow[2];

  bhigh[0] = ghigh[0];
  bhigh[1] = ghigh[1];
  bhigh[2] = ghigh[2];

  Index coords(0,1,2);

//  int borderfit = 0;

  switch (dir)
  {
    case west:
    case east:
      coords = Index(0,1,2);
      break;
    case south:
    case north:
      coords = Index(1,2,0);
      break;
    case down:
    case up:
      coords = Index(2,0,1);
      break;
  }

  int normal = coords[0];
  int t1 = coords[1];
  int t2 = coords[2];

  switch (dir)
  {
    case west:
    case south:
    case down:
      bhigh[normal] = glow[normal]+thickness-1+distanceLow;
      blow[normal] = glow[normal]+distanceLow;

      break;
    case east:
    case north:
    case up:
      blow[normal] = ghigh[normal]-thickness+2-distance;
      bhigh[normal] = ghigh[normal]-distance+1;

      break;
  }

  if (restricted) {
    bhigh[t1] = ghigh[t1]-distanceHigh;
    blow[t1] = glow[t1]+distance;

    bhigh[t2] = ghigh[t2]-distanceHigh;
    blow[t2] = glow[t2]+distance;
  }

  return true;
}
#endif

