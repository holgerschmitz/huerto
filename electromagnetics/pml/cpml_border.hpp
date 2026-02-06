#ifndef MPULSE_CPML_BORDER_H
#define MPULSE_CPML_BORDER_H

#include "../../types.hpp"
#include "../current.hpp"

class CPMLBorder : public CurrentBlock
{
  public:
#ifdef HUERTO_ONE_DIM
    typedef schnek::GridRegistration GridLineRegistration;
#else
    /// 1D grid containing the coefficients
    typedef schnek::ProjectedGridRegistration<1, DIMENSION> GridLineRegistration;
#endif
    void initCurrents(CurrentContainer &container);
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void init();
  private:

    void initCoefficients();

    int thickness;
    double kappaMax;
    double aMax;
    double sigmaMax;
    double eps;
};

class CPMLBorderCurrent : public Current
{
  public:
    CPMLBorderCurrent(int thickness,
                      Direction dir,
                      bool isH,
                      double kappaMax,
                      double aMax,
                      double sigmaMax,
                      double eps,
                      CurrentBlock &borderBlock);
  protected:
    bool reverse;
    int thickness;

    size_t dim;
    int transverse1, transverse2;

    Direction dir;
    bool isH;
    int lowOffset;
    int highOffset;

    int zerolayer;

    double kappaMax;
    double aMax;
    double sigmaMax;
    double eps;

    Range borderRange;

    /// 1D grid containing the coefficients
    CPMLBorder::GridLineRegistration bCoeff;
    CPMLBorder::GridLineRegistration cCoeff;

    CurrentBlock &borderBlock;

    void makeCoeff();
};

class CPMLBorderECurrent : public CPMLBorderCurrent
{
  public:
    CPMLBorderECurrent(int thickness,
                       Direction dir,
                       double kappaMax,
                       double aMax,
                       double sigmaMax,
                       double eps,
                       CurrentBlock &borderBlock);

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:

    // The B-field components
    schnek::GridRegistration B[3];
    
    // The components of the Psi grid
    schnek::GridRegistration Psi[2];
    double dx;
};

class CPMLBorderHCurrent : public CPMLBorderCurrent
{
  public:
    CPMLBorderHCurrent(int thickness,
                       Direction dir,
                       double kappaMax,
                       double aMax,
                       double sigmaMax,
                       double eps,
                       CurrentBlock &borderBlock);

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:
    // The E-field components
    schnek::GridRegistration E[3];
    // The components of the Psi grid
    schnek::GridRegistration Psi[2];
    double dx;
};


#endif
