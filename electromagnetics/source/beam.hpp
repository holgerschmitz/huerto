/*
 *  beam.hpp
 *
 *  Created on: 23 Dec 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_ELECTROMAGNETICS_SOURCE_BEAM_HPP_
#define HUERTO_ELECTROMAGNETICS_SOURCE_BEAM_HPP_

#include "incsource.hpp"

//===============================================================
//==========  Gaussian Beam Source
//===============================================================

#ifndef HUERTO_ONE_DIM

class GaussBeamSource : public IncidentSource
{
  public:
    ~GaussBeamSource() {}
  protected:
    pCurrent makeECurrent(int distance, Direction dir);
    pCurrent makeHCurrent(int distance, Direction dir);

    void initParameters(schnek::BlockParameters &blockPars);

    Vector k;
    Vector origin;
    Vector3d H;
    double waist;
    double rise;
    double offset;
    double eps;

    // degree to which the beam is circularly polarised
    // 0 means linear polarisation
    // 1 means circular polarisation
    double circ;
    int superGaussian;

};

class GaussBeamSourceEFunc
{
  public:
    GaussBeamSourceEFunc(Direction dir, SimulationContext &context);
    void setParam(Vector k,
                  Vector origin,
                  Vector3d H,
                  double waist,
                  double rise,
                  double offset,
                  double eps,
                  double circ,
                  int superGaussian);

#ifdef HUERTO_TWO_DIM
    Vector3d getHField(int i, int j, double time);
#endif

#ifdef HUERTO_THREE_DIM
    Vector3d getHField(int i, int j, int k, double time);
#endif

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(double) {}

  private:
    Vector k;
#ifdef HUERTO_TWO_DIM
    Vector kperp;
#endif

#ifdef HUERTO_THREE_DIM
    Vector kperpA;
    Vector kperpB;
#endif
    double kn;
    Vector origin;
    Vector3d H;
    Vector3d Hp;

    double dt;
    double om;
    double zr;
    double waist;
    double rise;
    double offset;
    double eps;
    int superGaussian;

    Vector dx;

    Direction dir;
    SimulationContext context;
};

class GaussBeamSourceHFunc
{
  public:
    GaussBeamSourceHFunc(Direction dir, SimulationContext &context);
    void setParam(Vector k,
                  Vector origin,
                  Vector3d E,
                  double waist,
                  double rise,
                  double offset,
                  double eps,
                  double circ,
                  int superGaussian);

#ifdef HUERTO_TWO_DIM
    Vector3d getEField(int i, int j, double time);
#endif

#ifdef HUERTO_THREE_DIM
    Vector3d getEField(int i, int j, int k, double time);
#endif

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(double) {}

  private:
    Vector k;
#ifdef HUERTO_TWO_DIM
    Vector kperp;
#endif

#ifdef HUERTO_THREE_DIM
    Vector kperpA;
    Vector kperpB;
#endif
    double kn;
    Vector origin;
    Vector3d E;
    Vector3d Ep;

    double dt;
    double om;
    double zr;
    double waist;
    double rise;
    double offset;
    double eps;
    int superGaussian;

    Vector dx;

    Direction dir;
    SimulationContext context;
};

#endif // not HUERTO_ONE_DIM

#endif /* HUERTO_ELECTROMAGNETICS_SOURCE_BEAM_HPP_ */
