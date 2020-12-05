#include "plane_wave.hpp"

#include "../../constants.hpp"
#include "../../maths/vector/vector.hpp"

#include <cmath>

//===============================================================
//==========  Plane Wave
//===============================================================

inline double applyPlaneGaussField(double pos, double width, double F) {
  double r = pos/width;
  return F*exp(-r*r)*sin(pos);
}

inline Vector3d kCrossB(const Vector &k, const Vector3d &H) {
#ifdef HUERTO_ONE_DIM
  return Vector3d(0, -k[0]*H[2], k[0]*H[1]);
#endif

#ifdef HUERTO_TWO_DIM
  return Vector3d(
    k[1]*H[2],
    - k[0]*H[2],
    k[0]*H[1] - k[1]*H[0]
  );
#endif

#ifdef HUERTO_THREE_DIM
  return cross(k, H);
#endif
}


void PlaneWaveFieldFunc::setParam(double ramp) {
  this->ramp = ramp;
}

inline double PlaneWaveFieldFunc::fieldFunc(double pos, double F) {
  if (pos>0) return 0.0;
  double f = F*sin(pos);
  return (pos > -ramp) ? -pos/ramp*f : f;
}

pCurrent PlaneWaveSource::makeECurrent(int distance, Direction dir)
{
  typedef GenericIncidentSourceESource<PlaneWaveFieldFunc> PlaneWaveSourceEFunc;
  typedef IncidentSourceECurrent<PlaneWaveSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setGenericParam(k, H, origin, eps);
  cur->setParam(ramp);
  return pCurrent(cur);
}

pCurrent PlaneWaveSource::makeHCurrent(int distance, Direction dir)
{
  Vector3d E = kCrossB(k, H);

  double bmag = norm(H);
  double factor = -bmag/norm(E);

  E *= clight*factor/sqrt(eps);

  typedef GenericIncidentSourceHSource<PlaneWaveFieldFunc> PlaneWaveSourceHFunc;
  typedef IncidentSourceHCurrent<PlaneWaveSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setGenericParam(k, E, origin, eps);
  cur->setParam(ramp);
  return pCurrent(cur);
}

void PlaneWaveSource::initParameters(schnek::BlockParameters &blockPars)
{
  IncidentSource::initParameters(blockPars);

  blockPars.addArrayParameter("k", this->k, 0.0);
  blockPars.addArrayParameter("H", this->H, 0.0);

  blockPars.addParameter("ramp", &this->ramp, 0.5);
  blockPars.addParameter("eps", &this->eps, 1.0);
  blockPars.addArrayParameter("origin", this->origin, Vector(0));
}

//===============================================================
//==========  Plane Gaussian Wave Packet
//===============================================================

void PlaneGaussFieldFunc::setParam(double width) {
  this->width = width;
}

inline double PlaneGaussFieldFunc::fieldFunc(double pos, double F) {
  double r = pos/width;
  return F*exp(-r*r)*sin(pos);
}


pCurrent PlaneGaussSource::makeECurrent(int distance, Direction dir)
{
  typedef GenericIncidentSourceESource<PlaneGaussFieldFunc> PlaneGaussSourceEFunc;
  typedef IncidentSourceECurrent<PlaneGaussSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setGenericParam(k, H, origin, eps);
  cur->setParam(width);
  return pCurrent(cur);
}

pCurrent PlaneGaussSource::makeHCurrent(int distance, Direction dir)
{
  Vector3d E = kCrossB(k, H);

  double bmag = norm(H);
  double factor = -bmag/norm(E);

  E *= clight*factor/sqrt(eps);

  typedef GenericIncidentSourceHSource<PlaneGaussFieldFunc> PlaneGaussSourceHFunc;
  typedef IncidentSourceHCurrent<PlaneGaussSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setGenericParam(k, E, origin, eps);
  cur->setParam(width);
  return pCurrent(cur);
}

void PlaneGaussSource::initParameters(schnek::BlockParameters &blockPars)
{
  IncidentSource::initParameters(blockPars);

  blockPars.addArrayParameter("k", this->k, 0.0);
  blockPars.addArrayParameter("H", this->H, 0.0);

  blockPars.addParameter("width", &this->width, 10.0);

  blockPars.addParameter("eps", &this->eps, 1.0);
  blockPars.addArrayParameter("origin", this->origin, Vector(0));
}

