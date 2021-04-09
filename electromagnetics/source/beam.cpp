/*
 *  beam.cpp
 *
 *  Created on: 23 Dec 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#include "beam.hpp"

#include "../../maths/vector/vector.hpp"


#ifndef HUERTO_ONE_DIM

pCurrent GaussBeamSource::makeECurrent(int distance, Direction dir) {
  typedef IncidentSourceECurrent<GaussBeamSourceEFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  std::cout << "GaussBeamSource::makeECurrent " << dir << " (" << H << ")" << std::endl;
  cur->setParam(k, origin, H, waist, rise, offset, eps, circ, superGaussian);
  return pCurrent(cur);
}

pCurrent GaussBeamSource::makeHCurrent(int distance, Direction dir) {
#ifdef HUERTO_TWO_DIM
  Vector3d k3(k[0], k[1], 0.0);
#endif

#ifdef HUERTO_THREE_DIM
  Vector3d k3(k);
#endif

  Vector3d E = cross(k3, H);

  double bmag = norm(H);
  double factor = -clight*bmag/norm(E);

  E *= factor/sqrt(eps);

  std::cout << "GaussBeamSource::makeHCurrent " << dir << " (" << H << ") (" << k3 << ") (" << E << ")" << std::endl;

  typedef IncidentSourceHCurrent<GaussBeamSourceHFunc> CurrentType;
  CurrentType *cur = new CurrentType(distance, dir, getContext());
  cur->setParam(k, origin, E, waist, rise, offset, eps, circ, superGaussian);
  return pCurrent(cur);
}

void GaussBeamSource::initParameters(schnek::BlockParameters &blockPars) {
  IncidentSource::initParameters(blockPars);

  blockPars.addArrayParameter("k", this->k, 0.0);
  blockPars.addArrayParameter("origin", this->origin, 0.0);

  blockPars.addArrayParameter("H", this->H, 0.0);

  blockPars.addParameter("waist", &this->waist, 10.0);
  blockPars.addParameter("rise", &this->rise, 10.0);
  blockPars.addParameter("offset", &this->offset, 40.0);
  blockPars.addParameter("eps", &this->eps, 1.0);
  blockPars.addParameter("circ", &this->circ, 0.0);
  blockPars.addParameter("superGaussian", &this->superGaussian, 1);
}


GaussBeamSourceEFunc::GaussBeamSourceEFunc(Direction dir, SimulationContext &context)
  : dir(dir), context(context)
{}

void GaussBeamSourceEFunc::setParam(Vector k,
                                    Vector origin,
                                    Vector3d H,
                                    double waist,
                                    double rise,
                                    double offset,
                                    double eps,
                                    double circ,
                                    int superGaussian)
{
  this->k = k;
  kn = norm(k);

#ifdef HUERTO_TWO_DIM
  Vector3d k3(k[0], k[1], 0.0);
  kperp = Vector(k[1]/kn, -k[0]/kn);
#endif

#ifdef HUERTO_THREE_DIM
  Vector3d k3(k);
  kperpA = cross(k3, H);
  kperpA /= norm(kperpA);
  kperpB = cross(k, kperpA);
  kperpB /= norm(kperpB);
#endif

  this->origin = origin;

  this->H = H / mu_0;
  Vector3d kxH = cross(k3, H);
  this->Hp =  circ * norm(this->H) * kxH / norm(kxH);

  std::cout << "GaussBeamSourceEFunc::setParam " << dir << " (" << this->H << ") (" << kxH << ") " << circ << " (" << this->Hp << ")" << std::endl;

  this->waist = waist;
  this->rise = kn*rise;
  this->offset = kn*offset;
  this->eps = eps;
  this->superGaussian = superGaussian;

  dt = context.getDt();
  om = clight*kn/sqrt(eps);
  zr = kn*kn*waist*waist/2.0;

  dx = context.getDx();
}


#ifdef HUERTO_TWO_DIM
Vector3d GaussBeamSourceEFunc::getHField(int i, int j, double time)
{
  // The current time in s
  double realtime = time - 0.5*dt;

  // The current position relative to the origin in m
  double x = i*dx[0] - origin[0];
  double y = j*dx[1] - origin[1];

  // The normalised position along the beam axis
  Vector3d zA(k[0]*x + k[1]*(y+0.5*dx[1]),
              k[0]*(x+0.5*dx[0]) + k[1]*y,
              k[0]*(x+0.5*dx[0]) + k[1]*(y+0.5*dx[1]));

  // The position perpendicular to the beam axis in physical units
  Vector3d pA(kperp[0]*x + kperp[1]*(y+0.5*dx[1]),
              kperp[0]*(x+0.5*dx[0]) + kperp[1]*y,
              kperp[0]*(x+0.5*dx[0]) + kperp[1]*(y+0.5*dx[1]));

  Vector3d h;
  for (int d=0; d<3; ++d) {
    // z is normalised to the wavelength because k is the wavevector in 1/m
    double z = zA[d];

    // r is NOT normalised because kperp has been normalised to length 1
    double r = pA[d];

    // The current beam width away from the focal point in m
    double w = waist*sqrt(1 + z*z/(zr*zr));

    double Rinv = z / (z*z + zr*zr);
    double ph = 0.5*kn*kn*r*r*Rinv;

    z -= om*realtime;

    double zenv = (z + offset)/rise;
    double env = (zenv>0) ? exp(-zenv*zenv) : 1.0;

    double amp = sqrt(waist/w)*env*exp( -std::pow( r/w, 2*superGaussian) );

    h[d] = amp*(H[d]*sin(z + ph) + Hp[d]*cos(z + ph));

//    if (j == 1000 && i>1000) {
//      std::cout << "getHField " << i << ":" << d << " " << z << " " << r << " " << w << " | "
//              << Rinv << " "  << ph << " "  << zenv << " "  << env << " | "  << amp << " " << h[d] << std::endl;
//    }
  }

  return h;
}
#endif

#ifdef HUERTO_THREE_DIM
Vector3d GaussBeamSourceEFunc::getHField(int i, int j, int l, double time)
{
  double realtime = time - 0.5*dt;

  double x = i*dx[0] - origin[0];
  double y = j*dx[1] - origin[1];
  double z = l*dx[2] - origin[2];

  Vector3d zA(k[0]*x + k[1]*(y+0.5*dx[1]) + k[2]*(z+0.5*dx[2]),
              k[0]*(x+0.5*dx[0]) + k[1]*y + k[2]*(z+0.5*dx[2]),
              k[0]*(x+0.5*dx[0]) + k[1]*(y+0.5*dx[1]) + k[2]*z);

  Vector3d pA(kperpA[0]*x + kperpA[1]*(y+0.5*dx[1]) + kperpA[2]*(z+0.5*dx[2]),
              kperpA[0]*(x+0.5*dx[0]) + kperpA[1]*y + kperpA[2]*(z+0.5*dx[2]),
              kperpA[0]*(x+0.5*dx[0]) + kperpA[1]*(y+0.5*dx[1]) + kperpA[2]*z);

  Vector3d pB(kperpB[0]*x + kperpB[1]*(y+0.5*dx[1]) + kperpB[2]*(z+0.5*dx[2]),
              kperpB[0]*(x+0.5*dx[0]) + kperpB[1]*y + kperpB[2]*(z+0.5*dx[2]),
              kperpB[0]*(x+0.5*dx[0]) + kperpB[1]*(y+0.5*dx[1]) + kperpB[2]*z);

  Vector3d h;
  for (int d=0; d<3; ++d)
  {
    // z is normalised to the wavelength because k is the wavevector in 1/m
    double z = zA[d];
    // r is NOT normalised because kperp has been normalised to length 1
    double rA = pA[d];
    double rB = pB[d];
    double r2 = rA*rA + rB*rB;
    double r = sqrt(r2);
    double w = waist*sqrt(1 + z*z/(zr*zr));
    double Rinv = z / (z*z + zr*zr);
    double ph = 0.5*kn*kn*r*r*Rinv;

    z -= om*realtime;

    double zenv = (z + offset)/rise;
    double env = (zenv>0) ? exp(-zenv*zenv) : 1.0;

    double amp = sqrt(waist/w)*env*exp( -std::pow( r/w, 2*superGaussian) );

    h[d] = amp*(H[d]*sin(z + ph) + Hp[d]*cos(z + ph));
  }

  return h;
}
#endif


GaussBeamSourceHFunc::GaussBeamSourceHFunc(Direction dir, SimulationContext &context)
  : dir(dir), context(context)
{}

void GaussBeamSourceHFunc::setParam(Vector k,
                                    Vector origin,
                                    Vector3d E,
                                    double waist,
                                    double rise,
                                    double offset,
                                    double eps,
                                    double circ,
                                    int superGaussian)
{
  this->k = k;
  kn = norm(k);

#ifdef HUERTO_TWO_DIM
  Vector3d k3(k[0], k[1], 0.0);
  kperp = Vector(k[1]/kn, -k[0]/kn);
#endif

#ifdef HUERTO_THREE_DIM
  Vector3d k3(k);
  kperpA = cross(k3, E);
  kperpA /= norm(kperpA);
  kperpB = cross(k3, kperpA);
  kperpB /= norm(kperpB);
#endif

  this->origin = origin;

  this->E = E;
  Vector3d kxE = cross(k3, E);
  this->Ep =  circ * norm(E) * kxE/norm(kxE);

  std::cout << "GaussBeamSourceHFunc::setParam " << dir << " (" << this->E << ") (" << kxE << ") (" << this->Ep << ")" << std::endl;

  this->waist = waist;
  this->rise = kn*rise;
  this->offset = kn*offset;
  this->eps = eps;
  this->superGaussian = superGaussian;

  om = clight*kn/sqrt(eps);
  zr = kn*kn*waist*waist/2.0;

  dx = context.getDx();
}

#ifdef HUERTO_TWO_DIM
Vector3d GaussBeamSourceHFunc::getEField(int i, int j, double time) {
  double realtime = time;

  double x = i*dx[0] - origin[0];
  double y = j*dx[1] - origin[1];

  Vector3d zA(k[0]*(x+0.5*dx[0]) + k[1]*y,
              k[0]*x + k[1]*(y+0.5*dx[1]),
              k[0]*x + k[1]*y);

  Vector3d pA(kperp[0]*(x+0.5*dx[0]) + kperp[1]*y,
              kperp[0]*x + kperp[1]*(y+0.5*dx[1]),
              kperp[0]*x + kperp[1]*y);

  Vector3d e;
  for (int d=0; d<3; ++d) {
    // z is normalised to the wavelength because k is the wavevector in 1/m
    double z = zA[d];

    // r is NOT normalised because kperp has been normalised to length 1
    double r = pA[d];
    double w = waist*sqrt(1 + z*z/(zr*zr));
    double Rinv = z / (z*z + zr*zr);
    double ph = 0.5*kn*kn*r*r*Rinv;

    z -= om*realtime;

    double zenv = (z + offset)/rise;
    double env = (zenv>0) ? exp(-zenv*zenv) : 1.0;

    double amp = sqrt(waist/w)*env*exp( - std::pow( r/w, 2*superGaussian));
    e[d] = amp*(E[d]*sin(z + ph) + Ep[d]*cos(z + ph));

//    if (j == 1000 && i>1000) {
//      std::cout << "getEField " << i << ":" << d << " " << z << " " << r << " " << w << " | "
//              << Rinv << " "  << ph << " "  << zenv << " "  << env << " | "  << amp << " " << e[d] << std::endl;
//    }
  }

  return e;
}
#endif

#ifdef HUERTO_THREE_DIM
Vector3d GaussBeamSourceHFunc::getEField(int i, int j, int l, double time) {
  double realtime = time;

  double x = i*dx[0] - origin[0];
  double y = j*dx[1] - origin[1];
  double z = l*dx[2] - origin[2];

  Vector3d zA(k[0]*(x+0.5*dx[0]) + k[1]*y + k[2]*z,
              k[0]*x + k[1]*(y+0.5*dx[1]) + k[2]*z,
              k[0]*x + k[1]*y + k[2]*(z+0.5*dx[2]));

  Vector3d pA(kperpA[0]*(x+0.5*dx[0]) + kperpA[1]*y + kperpA[2]*z,
              kperpA[0]*x + kperpA[1]*(y+0.5*dx[1]) + kperpA[2]*z,
              kperpA[0]*x + kperpA[1]*y + kperpA[2]*(z+0.5*dx[2]));

  Vector3d pB(kperpB[0]*(x+0.5*dx[0]) + kperpB[1]*y + kperpB[2]*z,
              kperpB[0]*x + kperpB[1]*(y+0.5*dx[1]) + kperpB[2]*z,
              kperpB[0]*x + kperpB[1]*y + kperpB[2]*(z+0.5*dx[2]));

  Vector3d e;
  for (int d=0; d<3; ++d) {
    // z is normalised to the wavelength because k is the wavevector in 1/m
    double z = zA[d];
    // r is NOT normalised because kperp has been normalised to length 1
    double rA = pA[d];
    double rB = pB[d];
    double r2 = rA*rA + rB*rB;
    double r = sqrt(r2);
    double w = waist*sqrt(1 + z*z/(zr*zr));
    double Rinv = z / (z*z + zr*zr);
    double ph = 0.5*kn*kn*r*r*Rinv;

    z -= om*realtime;

    double zenv = (z + offset)/rise;
    double env = (zenv>0) ? exp(-zenv*zenv) : 1.0;

    double amp = sqrt(waist/w)*env*exp( - std::pow( r/w, 2*superGaussian));
    e[d] = amp*(E[d]*sin(z + ph) + Ep[d]*cos(z + ph));
  }

  return e;
}
#endif


#endif // not HUERTO_ONE_DIM

