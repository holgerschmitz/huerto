#include "plane_wave.hpp"

#include "../../constants.hpp"
#include "../../maths/vector/vector.hpp"

#include <cmath>

//===============================================================
//==========  Plane Wave
//===============================================================

inline double applyPlaneWaveField(double pos, double ramp, double F) {
  if (pos>0) return 0.0;
  double f = F*sin(pos);
  return (pos > -ramp) ? -pos/ramp*f : f;
}

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

pCurrent PlaneWaveSource::makeECurrent(int distance, Direction dir)
{
  Vector3d E = kCrossB(k, H);

  double bmag = norm(H);
  double factor = -bmag/norm(E);

  E *= clight*factor/sqrt(eps);

  typedef GenericIncidentSourceHSource<PlaneWaveFieldFunc> PlaneWaveSourceEFunc
  typedef IncidentSourceECurrent<PlaneWaveSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setParam(k, E, H, ramp, eps, front);
  return pCurrent(cur);
}

pCurrent PlaneWaveSource::makeHCurrent(int distance, Direction dir)
{
  Vector3d E = kCrossB(k, H);

  double bmag = bmag = norm(H);
  double factor = -bmag/norm(E);

  E *= clight*factor/sqrt(eps);

  typedef IncidentSourceHCurrent<PlaneWaveSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setParam(k, E, H, ramp, eps, front);
  return pCurrent(cur);
}

void PlaneWaveSource::initParameters(schnek::BlockParameters &blockPars)
{
  IncidentSource::initParameters(blockPars);

  blockPars.addArrayParameter("k", this->k, 0.0);
  blockPars.addArrayParameter("H", this->H, 0.0);

  blockPars.addParameter("ramp", &this->ramp, 0.5);
  blockPars.addParameter("eps", &this->eps, 1.0);
  blockPars.addArrayParameter("front", this->front, Vector(0,0,0));
}


PlaneWaveSourceEFunc::PlaneWaveSourceEFunc(Direction dir, bool isH, SimulationContext &context)
  : dir(dir), isH(isH), context(context)
{}

void PlaneWaveSourceEFunc::setParam(Vector k, Vector E, Vector H, double ramp, double eps, const Vector &front)
{
  this->k = k;
  this->E = E;
  this->H = H;
  this->ramp = ramp;
  this->eps = eps;
  this->front = front;
  dt = context.getDt();
  om = clight*norm(k)/sqrt(eps);

  dx = context.getDx();
}

#ifdef HUERTO_ONE_DIM
Vector PlaneWaveSourceEFunc::getHField(int i, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx[0] - front[0];

  double posx = k[0]*x - om*realtime;
  double posy = k[0]*(x+0.5*dx[0]) - om*realtime;
  double posz = k[0]*(x+0.5*dx[0]) - om*realtime;

  double hx = applyPlaneWaveField(posx, ramp, H[0]);
  double hy = applyPlaneWaveField(posy, ramp, H[1]);
  double hz = applyPlaneWaveField(posz, ramp, H[2]);

  return Vector(hx, hy, hz);
}
#endif

#ifdef HUERTO_TWO_DIM
Vector PlaneWaveSourceEFunc::getHField(int i, int j, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];

  double posx = k[0]*x + k[1]*(y+0.5*dx[1]) - om*realtime;
  double posy = k[0]*(x+0.5*dx[0]) + k[1]*y - om*realtime;
  double posz = k[0]*(x+0.5*dx[0]) + k[1]*(y+0.5*dx[1]) - om*realtime;

  double hx = applyPlaneWaveField(posx, ramp, H[0]);
  double hy = applyPlaneWaveField(posy, ramp, H[1]);
  double hz = applyPlaneWaveField(posz, ramp, H[2]);

  return Vector(hx, hy, hz);
}
#endif

#ifdef HUERTO_THREE_DIM
Vector PlaneWaveSourceEFunc::getHField(int i, int j, int l, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];
  double z = l*dx[2] - front[2];

  double posx = k[0]*x + k[1]*(y+0.5*dx[1]) + k[2]*(z+0.5*dx[2]) - om*realtime;
  double posy = k[0]*(x+0.5*dx[0]) + k[1]*y + k[2]*(z+0.5*dx[2]) - om*realtime;
  double posz = k[0]*(x+0.5*dx[0]) + k[1]*(y+0.5*dx[1]) + k[2]*z - om*realtime;

  double hx = applyPlaneWaveField(posx, ramp, H[0]);
  double hy = applyPlaneWaveField(posy, ramp, H[1]);
  double hz = applyPlaneWaveField(posz, ramp, H[2]);

  return Vector(hx, hy, hz);
}
#endif


PlaneWaveSourceHFunc::PlaneWaveSourceHFunc(Direction dir, bool isH, SimulationContext &context)
  : dir(dir), isH(isH), context(context)
{}

void PlaneWaveSourceHFunc::setParam(Vector k, Vector E, Vector H, double ramp, double eps, const Vector &front)
{
  this->k = k;
  this->E = E;
  this->H = H;
  this->ramp = ramp;
  this->eps = eps;
  this->front = front;
  dt = context.getDt();
  om = clight*sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])/sqrt(eps);

  dx = context.getDx();
}

#ifdef HUERTO_ONE_DIM
Vector PlaneWaveSourceHFunc::getEField(int i, double time)
{
  double realtime = time;

  double x = i*dx[0] - front[0];

  double posx = k[0]*(x+0.5*dx[0]) - om*realtime;
  double posy = k[0]*x - om*realtime;
  double posz = k[0]*x - om*realtime;

  double ex = applyPlaneWaveField(posx, ramp, E[0]);
  double ey = applyPlaneWaveField(posy, ramp, E[1]);
  double ez = applyPlaneWaveField(posz, ramp, E[2]);

  return Vector(ex, ey, ez);
}
#endif

#ifdef HUERTO_TWO_DIM
Vector PlaneWaveSourceHFunc::getEField(int i, int j, double time)
{
  double realtime = time;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];

  double posx = k[0]*(x+0.5*dx[0]) + k[1]*y - om*realtime;
  double posy = k[0]*x + k[1]*(y+0.5*dx[1]) - om*realtime;
  double posz = k[0]*x + k[1]*y - om*realtime;

  double ex = applyPlaneWaveField(posx, ramp, E[0]);
  double ey = applyPlaneWaveField(posy, ramp, E[1]);
  double ez = applyPlaneWaveField(posz, ramp, E[2]);

  return Vector(ex, ey, ez);
}
#endif

#ifdef HUERTO_THREE_DIM
Vector PlaneWaveSourceHFunc::getEField(int i, int j, int l, double time)
{
  double realtime = time;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];
  double z = l*dx[2] - front[2];

  double posx = k[0]*(x+0.5*dx[0]) + k[1]*y + k[2]*z - om*realtime;
  double posy = k[0]*x + k[1]*(y+0.5*dx[1]) + k[2]*z - om*realtime;
  double posz = k[0]*x + k[1]*y + k[2]*(z+0.5*dx[2]) - om*realtime;

  double ex = applyPlaneWaveField(posx, ramp, E[0]);
  double ey = applyPlaneWaveField(posy, ramp, E[1]);
  double ez = applyPlaneWaveField(posz, ramp, E[2]);

  return Vector(ex, ey, ez);
}
#endif

//===============================================================
//==========  Plane Gaussian Wave Packet
//===============================================================

pCurrent PlaneGaussSource::makeECurrent(int distance, Direction dir)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);

  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);

  E[0] *= clight*factor/sqrt(eps);
  E[1] *= clight*factor/sqrt(eps);
  E[2] *= clight*factor/sqrt(eps);

  typedef IncidentSourceECurrent<PlaneGaussSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setParam(k, E, H, width, eps, front);
  return pCurrent(cur);
}

pCurrent PlaneGaussSource::makeHCurrent(int distance, Direction dir)
{
  Vector k(kx,ky,kz);
  Vector H(Hx,Hy,Hz);

  Vector E(ky*Hz-kz*Hy, kz*Hx-kx*Hz, kx*Hy-ky*Hx);

  double bmag = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
  double factor = -bmag/sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);

  E[0] *= clight*factor/sqrt(eps);
  E[1] *= clight*factor/sqrt(eps);
  E[2] *= clight*factor/sqrt(eps);

  typedef IncidentSourceHCurrent<PlaneGaussSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setParam(k, E, H, width, eps, front);
  return pCurrent(cur);
}

void PlaneGaussSource::initParameters(schnek::BlockParameters &blockPars)
{
  IncidentSource::initParameters(blockPars);

  blockPars.addParameter("kx", &this->kx, 0.0);
  blockPars.addParameter("ky", &this->ky, 0.0);
  blockPars.addParameter("kz", &this->kz, 1.0);

  blockPars.addParameter("Hx", &this->Hx, 0.0);
  blockPars.addParameter("Hy", &this->Hy, 0.0);
  blockPars.addParameter("Hz", &this->Hz, 0.0);

  blockPars.addParameter("width", &this->width, 10.0);

  blockPars.addParameter("eps", &this->eps, 1.0);
  blockPars.addArrayParameter("front", this->front, Vector(0,0,0));
}


PlaneGaussSourceEFunc::PlaneGaussSourceEFunc(Direction dir, bool isH, SimulationContext &context)
  : dir(dir), isH(isH), context(context)
{}

void PlaneGaussSourceEFunc::setParam(Vector k, Vector E, Vector H, double width, double eps, const Vector &front)
{
  this->k = k;
  this->E = E;
  this->H = H;
  double kn = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
  this->width = width*kn;
  this->eps = eps;
  this->front = front;
  dt = context.getDt();
  om = clight*kn/sqrt(eps);

  dx = context.getDx();
}


Vector PlaneGaussSourceEFunc::getHField(int i, int j, int l, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];
  double z = l*dx[2] - front[2];

  double posx = k[0]*x + k[1]*(y+0.5*dx[1]) + k[2]*(z+0.5*dx[2]) - om*realtime;
  double posy = k[0]*(x+0.5*dx[0]) + k[1]*y + k[2]*(z+0.5*dx[2]) - om*realtime;
  double posz = k[0]*(x+0.5*dx[0]) + k[1]*(y+0.5*dx[1]) + k[2]*z - om*realtime;

  double hx = applyPlaneGaussField(posx, width, H[0]);
  double hy = applyPlaneGaussField(posy, width, H[1]);
  double hz = applyPlaneGaussField(posz, width, H[2]);

  return Vector(hx, hy, hz);
}


PlaneGaussSourceHFunc::PlaneGaussSourceHFunc(Direction dir, bool isH, SimulationContext &context)
  : dir(dir), isH(isH), context(context)
{}

void PlaneGaussSourceHFunc::setParam(Vector k, Vector E, Vector H, double width, double eps, const Vector &front)
{
  this->k = k;
  this->E = E;
  this->H = H;
  double kn = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
  this->width = width*kn;
  this->eps = eps;
  this->front = front;
  dt = context.getDt();
  om = clight*kn/sqrt(eps);

  dx = context.getDx();
}


Vector PlaneGaussSourceHFunc::getEField(int i, int j, int l, double time)
{
  double realtime = time;

  double x = i*dx[0] - front[0];
  double y = j*dx[1] - front[1];
  double z = l*dx[2] - front[2];

  double posx = k[0]*(x+0.5*dx[0]) + k[1]*y + k[2]*z - om*realtime;
  double posy = k[0]*x + k[1]*(y+0.5*dx[1]) + k[2]*z - om*realtime;
  double posz = k[0]*x + k[1]*y + k[2]*(z+0.5*dx[2]) - om*realtime;

  double ex = applyPlaneGaussField(posx, width, E[0]);
  double ey = applyPlaneGaussField(posy, width, E[1]);
  double ez = applyPlaneGaussField(posz, width, E[2]);

  return Vector(ex, ey, ez);
}

