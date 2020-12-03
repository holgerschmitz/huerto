#ifndef HUERTO_EM_SOURCE_PLANE_WAVE_H
#define HUERTO_EM_SOURCE_PLANE_WAVE_H

#include "incsource.hpp"

//===============================================================
//==========  Plane Wave
//===============================================================

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
    Vector front;
};

class PlaneWaveFieldFunc
{
  public:
    PlaneWaveFieldFunc();
    void setParam(double ramp);
    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(int) {}

  protected:
    double fieldFunc(double pos);
  private:
    double ramp;
};

class PlaneWaveSourceHFunc
{
  public:
    PlaneWaveSourceHFunc(Direction dir, SimulationContext &context);
    void setParam(Vector k, Vector E, Vector H, double ramp, double eps, const Vector &front);

    Vector getEField(int i, int j, int k, double time);

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(int) {}

  private:
    Vector k;
    Vector3d E;
    Vector3d H;

    double dt;
    double om;
    double ramp;
    double eps;
    Vector front;

    Vector dx;

    Direction dir;
    SimulationContext &context;
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

    double width;
    double eps;

    Vector front;
};

class PlaneGaussSourceEFunc
{
  public:
    PlaneGaussSourceEFunc(Direction dir, bool isH, SimulationContext &context);
    void setParam(Vector k, Vector E, Vector H, double width, double eps, const Vector &front);

    Vector getHField(int i, int j, int k, double time);

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(int) {}

  private:
    Vector k;
    Vector3d E;
    Vector3d H;

    double dt;
    double om;
    double width;
    double eps;
    Vector front;

    Vector dx;

    Direction dir;
    bool isH;
    SimulationContext &context;
};

class PlaneGaussSourceHFunc
{
  public:
    PlaneGaussSourceHFunc(Direction dir, bool isH, SimulationContext &context);
    void setParam(Vector k, Vector E, Vector H, double width, double eps, const Vector &front);

    Vector getEField(int i, int j, int k, double time);

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(int) {}

  private:
    Vector k;
    Vector3d E;
    Vector3d H;

    double dt;
    double om;
    double width;
    double eps;
    Vector front;

    Vector dx;

    Direction dir;
    bool isH;
    SimulationContext &context;
};


#endif // HUERTO_EM_SOURCE_PLANE_WAVE_H
