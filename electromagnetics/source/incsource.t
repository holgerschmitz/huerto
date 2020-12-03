/*
 * incsource.t
 *
 *  Created on: 30 Nov 2020
 *      Author: Holger Schmitz
 */

#include "border.hpp"

#include <boost/make_shared.hpp>

//===============================================================
//==========  IncidentSourceECurrent
//===============================================================

template<class SourceFunc>
IncidentSourceECurrent<SourceFunc>::IncidentSourceECurrent(int distance, Direction dir, SimulationContext &context)
  : IncidentSourceCurrent(distance, dir, false, context),
    SourceFunc(dir, context)
{}

template<class SourceFunc>
void IncidentSourceECurrent<SourceFunc>::init() {
  Index blow, bhigh;

  if (!getBorderExtent(IncidentSourceCurrent::dir, 1, distance, blow, bhigh, false, context)) {
    return;
  }

  pJx = std::make_shared<Grid>(blow, bhigh);
  pJy = std::make_shared<Grid>(blow, bhigh);
  pJz = std::make_shared<Grid>(blow, bhigh);

  pGrid allJ[3];
  allJ[0] = pJx;
  allJ[1] = pJy;
  allJ[2] = pJz;

  pJ[0] = allJ[IncidentSourceCurrent::transverse1];
  pJ[1] = allJ[IncidentSourceCurrent::transverse2];

  if ((pJx!=0) && (pJy!=0) && (pJz!=0))
     this->initSourceFunc(pJx, pJy, pJz);
}

template<class SourceFunc>
void IncidentSourceECurrent<SourceFunc>::stepSchemeInit(double dt)
{
  stepScheme(0.5*dt);
}

template<class SourceFunc>
void IncidentSourceECurrent<SourceFunc>::stepScheme(double dt)
{
  Grid &J0 = *pJ[0];
  Grid &J1 = *pJ[1];

  Index low  = J0.getLo();
  Index high = J0.getHi();

  Index ind, indn;

  double Time = context.getTime();
  this->setTime(Time);

#ifdef HUERTO_ONE_DIM
  int off[1] = {0};
#endif
#ifdef HUERTO_TWO_DIM
  int off[2] = {0, 0};
#endif
#ifdef HUERTO_THREE_DIM
  int off[3] = {0, 0, 0};
#endif

  off[IncidentSourceCurrent::dim] = reverse?0:-1;
  double factor = reverse?1:-1;

#ifdef HUERTO_ONE_DIM
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
  {
    int x = ind[0]+off[0];
    Vector H = this->getHField(x, Time);
    J0(ind[0]) = factor*H[IncidentSourceCurrent::transverse2]/dN;
    J1(ind[0]) = -factor*H[IncidentSourceCurrent::transverse1]/dN;
  }
#endif

#ifdef HUERTO_TWO_DIM
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
  {
    int x = ind[0]+off[0];
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
    {
      int y = ind[1]+off[1];
      Vector H = this->getHField(x, y, Time);
      J0(ind[0], ind[1]) = factor*H[IncidentSourceCurrent::transverse2]/dN;
      J1(ind[0], ind[1]) = -factor*H[IncidentSourceCurrent::transverse1]/dN;
    }
  }
#endif

#ifdef HUERTO_THREE_DIM
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
  {
    int x = ind[0]+off[0];
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
    {
      int y = ind[1]+off[1];
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2])
      {
        int z = ind[2]+off[2];
        Vector H = this->getHField(x, y, z, Time);
        J0(ind[0], ind[1], ind[2]) = factor*H[IncidentSourceCurrent::transverse2]/dN;
        J1(ind[0], ind[1], ind[2]) = -factor*H[IncidentSourceCurrent::transverse1]/dN;
      }
    }
  }
#endif
}


//===============================================================
//==========  IncidentSourceHCurrent
//===============================================================

template<class SourceFunc>
IncidentSourceHCurrent<SourceFunc>::IncidentSourceHCurrent(int distance, Direction dir, SimulationContext &context)
  : IncidentSourceCurrent(distance, dir, true, context),
    SourceFunc(dir, context)
{}

template<class SourceFunc>
void IncidentSourceHCurrent<SourceFunc>::init()
{
  Index blow, bhigh;

  if (!getBorderExtent(IncidentSourceCurrent::dir, 1, distance, blow, bhigh, true, context)) return;

  pJx = std::make_shared<Grid>(blow, bhigh);
  pJy = std::make_shared<Grid>(blow, bhigh);
  pJz = std::make_shared<Grid>(blow, bhigh);

  pGrid allJ[3];
  allJ[0] = pJx;
  allJ[1] = pJy;
  allJ[2] = pJz;

  pJ[0] = allJ[IncidentSourceCurrent::transverse1];
  pJ[1] = allJ[IncidentSourceCurrent::transverse2];

  this->initSourceFunc(pJx, pJy, pJz);
}

template<class SourceFunc>
void IncidentSourceHCurrent<SourceFunc>::stepSchemeInit(double dt)
{
  stepScheme(0.5*dt);
}

template<class SourceFunc>
void IncidentSourceHCurrent<SourceFunc>::stepScheme(double dt)
{
  Grid &J0 = *pJ[0];
  Grid &J1 = *pJ[1];

  Index low  = J0.getLo();
  Index high = J0.getHi();

  Index ind, indn;

  double Time = context.getTime();
  this->setTime(Time);

#ifdef HUERTO_ONE_DIM
  int off[1] = {0};
#endif
#ifdef HUERTO_TWO_DIM
  int off[2] = {0, 0};
#endif
#ifdef HUERTO_THREE_DIM
  int off[3] = {0, 0, 0};
#endif

  off[IncidentSourceCurrent::dim] = reverse?0:1;
  double factor = reverse?1:-1;

#ifdef HUERTO_THREE_DIM
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
  {
    int x = ind[0]+off[0];
    Vector E = this->getEField(x, Time);
    J0(ind[0]) = -factor*E[IncidentSourceCurrent::transverse2]/dN;
    J1(ind[0]) = factor*E[IncidentSourceCurrent::transverse1]/dN;
  }
#endif

#ifdef HUERTO_TWO_DIM
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
  {
    int x = ind[0]+off[0];
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
    {
      int y = ind[1]+off[1];
      Vector E = this->getEField(x, y, Time);
      J0(ind[0], ind[1]) = -factor*E[IncidentSourceCurrent::transverse2]/dN;
      J1(ind[0], ind[1]) = factor*E[IncidentSourceCurrent::transverse1]/dN;
    }
  }
#endif

#ifdef HUERTO_THREE_DIM
  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0])
  {
    int x = ind[0]+off[0];
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1])
    {
      int y = ind[1]+off[1];
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2])
      {
        int z = ind[2]+off[2];
        Vector E = this->getEField(x,y,z,Time);
        J0(ind[0], ind[1], ind[2]) = -factor*E[IncidentSourceCurrent::transverse2]/dN;
        J1(ind[0], ind[1], ind[2]) = factor*E[IncidentSourceCurrent::transverse1]/dN;
      }
    }
  }
#endif
}

//===============================================================
//==========  GenericIncidentSourceESource
//===============================================================

template<class FieldFunc>
GenericIncidentSourceESource<FieldFunc>::GenericIncidentSourceESource(
        Vector k,
        Vector3d H,
        const Vector &front,
        double eps,
        SimulationContext &context)
        : k(k), H(H), front(front), context(context)
{
  dt = context.getDt();
  om = clight*norm(k)/sqrt(eps);

  dx = context.getDx();
}

#ifdef HUERTO_ONE_DIM
template<class FieldFunc>
Vector GenericIncidentSourceESource<FieldFunc>::getHField(int i, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx[0] - front[0];

  double posx = k[0]*x - om*realtime;
  double posy = k[0]*(x+0.5*dx[0]) - om*realtime;
  double posz = k[0]*(x+0.5*dx[0]) - om*realtime;

  double hx = fieldFunc(posx);
  double hy = fieldFunc(posy);
  double hz = fieldFunc(posz);

  return Vector(hx, hy, hz);
}
#endif

#ifdef HUERTO_TWO_DIM
template<class FieldFunc>
Vector GenericIncidentSourceESource<FieldFunc>::getHField(int i, int j, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];

  double posx = k[0]*x + k[1]*(y+0.5*dx[1]) - om*realtime;
  double posy = k[0]*(x+0.5*dx[0]) + k[1]*y - om*realtime;
  double posz = k[0]*(x+0.5*dx[0]) + k[1]*(y+0.5*dx[1]) - om*realtime;

  double hx = fieldFunc(posx);
  double hy = fieldFunc(posy);
  double hz = fieldFunc(posz);

  return Vector(hx, hy, hz);
}
#endif

#ifdef HUERTO_THREE_DIM
template<class FieldFunc>
Vector GenericIncidentSourceESource<FieldFunc>::getHField(int i, int j, int l, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];
  double z = l*dx[2] - front[2];

  double posx = k[0]*x + k[1]*(y+0.5*dx[1]) + k[2]*(z+0.5*dx[2]) - om*realtime;
  double posy = k[0]*(x+0.5*dx[0]) + k[1]*y + k[2]*(z+0.5*dx[2]) - om*realtime;
  double posz = k[0]*(x+0.5*dx[0]) + k[1]*(y+0.5*dx[1]) + k[2]*z - om*realtime;

  double hx = fieldFunc(posx);
  double hy = fieldFunc(posy);
  double hz = fieldFunc(posz);

  return Vector(hx, hy, hz);
}
#endif


//===============================================================
//==========  GenericIncidentSourceHSource
//===============================================================

template<class FieldFunc>
GenericIncidentSourceHSource<FieldFunc>::GenericIncidentSourceHSource(
        Vector k,
        Vector3d E,
        const Vector &front,
        double eps,
        SimulationContext &context)
        : k(k), E(E), front(front), context(context)
{
  dt = context.getDt();
  om = clight*norm(k)/sqrt(eps);

  dx = context.getDx();
}


#ifdef HUERTO_ONE_DIM
template<class FieldFunc>
Vector GenericIncidentSourceHSource<FieldFunc>::getEField(int i, double time)
{
  double realtime = time;

  double x = i*dx[0] - front[0];

  double posx = k[0]*(x+0.5*dx[0]) - om*realtime;
  double posy = k[0]*x - om*realtime;
  double posz = k[0]*x - om*realtime;

  double ex = fieldFunc(posx, ramp, E[0]);
  double ey = fieldFunc(posy, ramp, E[1]);
  double ez = fieldFunc(posz, ramp, E[2]);

  return Vector(ex, ey, ez);
}
#endif

#ifdef HUERTO_TWO_DIM
template<class FieldFunc>
Vector GenericIncidentSourceHSource<FieldFunc>::getEField(int i, int j, double time)
{
  double realtime = time;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];

  double posx = k[0]*(x+0.5*dx[0]) + k[1]*y - om*realtime;
  double posy = k[0]*x + k[1]*(y+0.5*dx[1]) - om*realtime;
  double posz = k[0]*x + k[1]*y - om*realtime;

  double ex = fieldFunc(posx, ramp, E[0]);
  double ey = fieldFunc(posy, ramp, E[1]);
  double ez = fieldFunc(posz, ramp, E[2]);

  return Vector(ex, ey, ez);
}
#endif

#ifdef HUERTO_THREE_DIM
template<class FieldFunc>
Vector GenericIncidentSourceHSource<FieldFunc>::getEField(int i, int j, int l, double time)
{
  double realtime = time;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];
  double z = l*dx[2] - front[2];

  double posx = k[0]*(x+0.5*dx[0]) + k[1]*y + k[2]*z - om*realtime;
  double posy = k[0]*x + k[1]*(y+0.5*dx[1]) + k[2]*z - om*realtime;
  double posz = k[0]*x + k[1]*y + k[2]*(z+0.5*dx[2]) - om*realtime;

  double ex = fieldFunc(posx);
  double ey = fieldFunc(posy);
  double ez = fieldFunc(posz);

  return Vector(ex, ey, ez);
}
#endif
