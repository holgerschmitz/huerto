#ifndef HUERTO_EM_SOURCE_PLANE_WAVE_H
#define HUERTO_EM_SOURCE_PLANE_WAVE_H

#include "incsource.hpp"

//===============================================================
//==========  Plane Wave
//===============================================================

/**
 * An infinite plane wave source
 */
class PlaneWaveSource : public IncidentSource
{
  public:
    ~PlaneWaveSource() {}
  protected:
    pCurrent makeECurrent(int distance, Direction dir);
    pCurrent makeHCurrent(int distance, Direction dir);

    void initParameters(schnek::BlockParameters &blockPars);

    Vector k;
    Vector3d H;
    double ramp;
    double eps;
    Vector origin;
};

class PlaneWaveFieldFunc
{
  public:
    PlaneWaveFieldFunc() {}
    void setParam(double ramp);
    void initSourceFunc(Grid&, Grid&, Grid&) {}
    void setTime(int) {}

  protected:
    double fieldFunc(double pos, double F);
  private:
    double ramp;
};

//===============================================================
//==========  Plane Gauss Packet Source
//===============================================================

class PlaneGaussSource : public IncidentSource
{
  public:
    ~PlaneGaussSource() {}
  protected:
    pCurrent makeECurrent(int distance, Direction dir);
    pCurrent makeHCurrent(int distance, Direction dir);

    void initParameters(schnek::BlockParameters &blockPars);

    Vector k;
    Vector3d H;

    double length;
    double eps;

    Vector origin;
};


class PlaneGaussFieldFunc
{
  public:
    PlaneGaussFieldFunc() {}
    void setParam(double width);
    void initSourceFunc(Grid&, Grid&, Grid&) {}
    void setTime(int) {}

  protected:
    double fieldFunc(double pos, double F);
  private:
    double length;
};


#endif // HUERTO_EM_SOURCE_PLANE_WAVE_H
