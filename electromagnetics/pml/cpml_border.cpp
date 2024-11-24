#include "cpml_border.hpp"

#include "../../constants.hpp"
#include "../fieldsolver.hpp"
#include "../source/border.hpp"

#include <schnek/tools/literature.hpp>

#include <memory>
#include <vector>

//===============================================================
//==========  CPMLBorder
//===============================================================

void CPMLBorder::initParameters(schnek::BlockParameters &blockPars) {
  CurrentBlock::initParameters(blockPars);

  blockPars.addParameter("d", &thickness, 8);
  blockPars.addParameter("kappaMax", &kappaMax, 15.0);
  blockPars.addParameter("aMax", &aMax, 0.0);
  blockPars.addParameter("sigmaMax", &sigmaMax, 1.0);
  blockPars.addParameter("eps", &eps, 1.0);
}


void CPMLBorder::init() {
  schnek::ChildBlock<CurrentBlock>::init();

  schnek::LiteratureArticle Roden2000("Roden2000", "Roden, J. A. and Gedney, S. D.",
      "An efficient fdtd implementation of the cfs-pml for arbitrary media",
      "Microwave and Optical Technology Letters", "2000", "27", "334--339");

  schnek::LiteratureManager::instance().addReference(
      "Implementation of the Convolution Perfectly Matched Layer (CPML) boundary condition", Roden2000);
}

void CPMLBorder::initCurrents(CurrentContainer &container)
{
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, east, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, west, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
#ifndef HUERTO_ONE_DIM
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, north, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, south, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
#endif
#ifdef HUERTO_THREE_DIM
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
  container.addCurrent(
    std::make_shared<CPMLBorderECurrent>(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
#endif

  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, east, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, west, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
#ifndef HUERTO_ONE_DIM
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, north, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, south, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
#endif
#ifdef HUERTO_THREE_DIM
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, up, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
  container.addMagCurrent(
    std::make_shared<CPMLBorderHCurrent>(thickness, down, this->kappaMax, this->aMax, this->sigmaMax, eps, boost::ref(*this))
  );
#endif

  initCoefficients();

}

void CPMLBorder::initCoefficients()
{
  schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();
  // initialize Kappas here

  Range gdomain = subdivision.getGlobalDomain();
  Index glow  = gdomain.getLo();
  Index ghigh = gdomain.getHi();

  Index low  = subdivision.getInnerLo();
  Index high = subdivision.getInnerHi();

  std::vector<Grid1d> KappaEdk(DIMENSION);
  std::vector<Grid1d> KappaHdk(DIMENSION);

  std::vector<Grid1d> CpmlSigmaE(DIMENSION);
  std::vector<Grid1d> CpmlSigmaH(DIMENSION);

  retrieveData("KappaEdx", KappaEdk[0]);
#ifndef HUERTO_ONE_DIM
  retrieveData("KappaEdy", KappaEdk[1]);
#endif
#ifdef HUERTO_THREE_DIM
  retrieveData("KappaEdz", KappaEdk[2]);
#endif

  retrieveData("KappaHdx", KappaHdk[0]);
#ifndef HUERTO_ONE_DIM
  retrieveData("KappaHdy", KappaHdk[1]);
#endif
#ifdef HUERTO_THREE_DIM
  retrieveData("KappaHdz", KappaHdk[2]);
#endif

  for (size_t dim = 0; dim<DIMENSION; ++dim)
  {
    std::cerr << "Dim " << dim << std::endl;

    KappaEdk[dim] = 1.0;
    KappaHdk[dim] = 1.0;


    Index blow, bhigh;
    Direction dir;

    switch (dim)
    {
      case 0: dir = west; break;
#ifndef HUERTO_ONE_DIM
      case 1: dir = south; break;
#endif
#ifdef HUERTO_THREE_DIM
      case 2: dir = down; break;
#endif
    }

    if (getBorderExtent(dir, thickness, 1, blow, bhigh, false, getContext(), false)) {
      int lowk  = blow[dim];
      int highk = bhigh[dim];
      int kLimit = highk-lowk + 1;

      for (int k=0; k<kLimit; ++k) {
        double x = 1 - double(k)/double(thickness);
        double x3 = x*x*x;

        KappaEdk[dim](lowk+k) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    if (getBorderExtent(dir, thickness, 1, blow, bhigh, true, getContext(), false)) {
      int lowk  = blow[dim];
      int highk = bhigh[dim];
      int kLimit = highk-lowk + 1;

      for (int k=0; k<kLimit; ++k) {
        double x = 1 - (double(k) - 0.5)/double(thickness);
        double x3 = x*x*x;

        KappaHdk[dim](lowk+k) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    switch (dim)
    {
      case 0: dir = east; break;
#ifndef HUERTO_ONE_DIM
      case 1: dir = north; break;
#endif
#ifdef HUERTO_THREE_DIM
      case 2: dir = up; break;
#endif
    }

    if (getBorderExtent(dir, thickness, 1, blow, bhigh, false, getContext(), false)) {
      int lowk  = blow[dim];
      int highk = bhigh[dim];
      int kLimit = highk-lowk + 1;

      for (int k=0; k<kLimit; ++k) {
        double x = 1 - double(k)/double(thickness);
        double x3 = x*x*x;

        KappaEdk[dim](highk-k) = 1 + (this->kappaMax - 1)*x3;
      }
    }

    if (getBorderExtent(dir, thickness, 1, blow, bhigh, true, getContext(), false)) {
      int lowk  = blow[dim];
      int highk = bhigh[dim];
      int kLimit = highk-lowk + 1;

      for (int k=0; k<kLimit; ++k) {
        double x = 1 - (double(k) - 0.5)/double(thickness);
        double x3 = x*x*x;

        KappaHdk[dim](highk-k) = 1 + (this->kappaMax - 1)*x3;
      }
    }

//
//
//    if (low[dim]<glow[dim]+thickness)
//    {
//
//      for (int i=0; i<=thickness; ++i)
//      {
//        double x  = 1 - double(i)/double(thickness);
//        double x3 = x*x*x;
//        (*pKappaEdk[dim])(low[dim]+i) = 1 + (this->kappaMax - 1)*x3;
//      }
//      for (int i=0; i<thickness; ++i)
//      {
//        double x  = 1 - (double(i)+0.5)/double(thickness);
//        double x3 = x*x*x;
//        (*pKappaHdk[dim])(low[dim]+i) = 1 + (this->kappaMax - 1)*x3;
//      }
//    }
//
//    if (high[dim]>ghigh[dim]-thickness)
//    {
//      for (int i=0; i<=thickness; ++i)
//      {
//        double x  = 1 - double(i)/double(thickness);
//        double x3 = x*x*x;
//        (*pKappaHdk[dim])(high[dim]-i) = 1 + (this->kappaMax - 1)*x3;
//      }
//      for (int i=0; i<thickness; ++i)
//
//      {
//        double x  = 1 - (double(i)+0.5)/double(thickness);
//        double x3 = x*x*x;
//        (*pKappaEdk[dim])(high[dim]-i) = 1 + (this->kappaMax - 1)*x3;
//      }
//    }
  }
}

//===============================================================
//==========  CPMLBorderCurrent
//===============================================================


CPMLBorderCurrent::CPMLBorderCurrent(int thickness, Direction dir, bool isH,
                                     double kappaMax, double aMax, double sigmaMax, double eps,
                                     CurrentBlock &borderBlock)
  : reverse(false), thickness(thickness), dir(dir), isH(isH), lowOffset(0), highOffset(0), zerolayer(0),
    kappaMax(kappaMax), aMax(aMax), sigmaMax(sigmaMax), eps(eps), borderBlock(borderBlock)
{
  switch (dir)
  {
    case east:
    case west:  dim = 0;
                transverse1 = 1;
                transverse2 = 2;
                break;
#ifndef HUERTO_ONE_DIM
    case north:
    case south: dim = 1;
                transverse1 = 0;
                transverse2 = 2;
                break;
#endif
#ifdef HUERTO_THREE_DIM
    case up:
    case down:  dim = 2;
                transverse1 = 0;
                transverse2 = 1;
                break;
#endif
  }

//  double eta = sqrt(mu_0/eps_0);
  sigmaMax = 10 * sigmaMax * clight / borderBlock.getContext().getDx()[dim];
}

void CPMLBorderCurrent::makeCoeff()
{
  Index low  = Jx.getLo();
  Index high = Jx.getHi();

  switch (dir)
  {
    case east:
#ifndef HUERTO_ONE_DIM
    case north:
#endif
#ifdef HUERTO_THREE_DIM
    case up:
#endif
      reverse = false;
      break;
    case west:
#ifndef HUERTO_ONE_DIM
    case south:
#endif
#ifdef HUERTO_THREE_DIM
    case down:
#endif
      reverse = true;
      break;
  }


  int lowk  = low[dim];
  int highk = high[dim];

  bCoeff.resize(lowk, highk);
  cCoeff.resize(lowk, highk);

  double offset = 0.0;
  lowOffset = 1;

  if (isH)
  {
    offset = 0.5;
    lowOffset = 0;
  }


  int kLimit = highk-lowk + 1;

  for (int k=0; k<kLimit; ++k) {
    double x = 1 - (double(k)-offset)/double(thickness);
    double x3 = x*x*x;

    int pos = reverse ? (lowk+k) : (highk-k);

    double sigma = x3*sigmaMax;
    double kappa = 1 + (kappaMax - 1)*x3;
    double a = aMax*(1-x);

    double b = exp(-(sigma/kappa + a));
    double c = sigma*(b-1)/(kappa*(sigma+kappa*a));

    bCoeff(pos) = b;
    cCoeff(pos) = c;
  }

}

//===============================================================
//==========  CPMLBorderECurrent
//===============================================================


CPMLBorderECurrent::CPMLBorderECurrent(int thickness,
                                       Direction dir,
                                       double kappaMax,
                                       double aMax,
                                       double sigmaMax,
                                       double eps,
                                       CurrentBlock &borderBlock)
  : CPMLBorderCurrent(thickness,dir,false,kappaMax,aMax,sigmaMax,eps,borderBlock),
    dx(0)
{}

void CPMLBorderECurrent::init()
{
  Index blow, bhigh;

  if (!getBorderExtent(dir, thickness, 1, blow, bhigh, false, borderBlock.getContext(), false)) return;

  Jx.resize(blow, bhigh);
  Jy.resize(blow, bhigh);
  Jz.resize(blow, bhigh);

  Jx = 0.0;
  Jy = 0.0;
  Jz = 0.0;

  switch (dir)
  {
    case east:
    case west:
      Psi[0] = Jy;
      Psi[1] = Jz;
      borderBlock.retrieveData("By", B[0]);
      borderBlock.retrieveData("Bz", B[1]);
      borderBlock.retrieveData("Bx", B[2]);

      dx = borderBlock.getContext().getDx()[0];
      break;
#ifndef HUERTO_ONE_DIM
    case north:
    case south:
      Psi[0] = Jz;
      Psi[1] = Jx;
      borderBlock.retrieveData("Bz", B[0]);
      borderBlock.retrieveData("Bx", B[1]);
      borderBlock.retrieveData("By", B[2]);

      dx = borderBlock.getContext().getDx()[1];
      break;
#endif
#ifdef HUERTO_THREE_DIM
    case up:
    case down:
      Psi[0] = Jx;
      Psi[1] = Jy;
      borderBlock.retrieveData("Bx", B[0]);
      borderBlock.retrieveData("By", B[1]);
      borderBlock.retrieveData("Bz", B[2]);

      dx = borderBlock.getContext().getDx()[2];
      break;
#endif
  }

  makeCoeff();
}

void CPMLBorderECurrent::stepSchemeInit(double /* dt */)
{}

void CPMLBorderECurrent::stepScheme(double /* dt */)
{
  Index low  = Psi[0].getLo();
  Index high = Psi[0].getHi();

  Grid &Psi0 = Psi[0];
  Grid &Psi1 = Psi[1];
  Field &B0 = B[0];
  Field &B1 = B[1];

  Index ind, indn;

  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0]) {
#ifndef HUERTO_ONE_DIM
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1]) {
#endif
#ifdef HUERTO_THREE_DIM
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2]) {
#endif

        int j = ind[dim];
        Index indm(ind);
        --indm[dim];

        Psi0[ind]
          = bCoeff(j)*Psi0[ind]
            + cCoeff(j)*(B1[ind]-B1[indm])/(mu_0*dx);
        Psi1[ind]
          = bCoeff(j)*Psi1[ind]
            - cCoeff(j)*(B0[ind]-B0[indm])/(mu_0*dx);
#ifdef HUERTO_THREE_DIM
      }
#endif
#ifndef HUERTO_ONE_DIM
    }
#endif
  }
}


//===============================================================
//==========  CPMLBorderHCurrent
//===============================================================

CPMLBorderHCurrent::CPMLBorderHCurrent(int thickness,
                                       Direction dir,
                                       double kappaMax,
                                       double aMax,
                                       double sigmaMax,
                                       double eps,
                                       CurrentBlock &borderBlock)
  : CPMLBorderCurrent(thickness,dir,true,kappaMax,aMax,sigmaMax,eps,borderBlock),
    dx(0)
{}

void CPMLBorderHCurrent::init()
{
  Index blow, bhigh;

  if (!getBorderExtent(dir, thickness, 1, blow, bhigh, true, borderBlock.getContext(), false)) return;

  Jx.resize(blow, bhigh);
  Jy.resize(blow, bhigh);
  Jz.resize(blow, bhigh);

  Jx = 0.0;
  Jy = 0.0;
  Jz = 0.0;

  switch (dir)
  {
    case east:
    case west:
      Psi[0] = Jy;
      Psi[1] = Jz;
      borderBlock.retrieveData("Ey", E[0]);
      borderBlock.retrieveData("Ez", E[1]);
      borderBlock.retrieveData("Ex", E[2]);

      dx = borderBlock.getContext().getDx()[0];
      break;
#ifndef HUERTO_ONE_DIM
    case north:
    case south:
      Psi[0] = Jz;
      Psi[1] = Jx;
      borderBlock.retrieveData("Ez", E[0]);
      borderBlock.retrieveData("Ex", E[1]);
      borderBlock.retrieveData("Ey", E[2]);

      dx = borderBlock.getContext().getDx()[1];
      break;
#endif
#ifdef HUERTO_THREE_DIM
    case up:
    case down:
      Psi[0] = Jx;
      Psi[1] = Jy;
      borderBlock.retrieveData("Ex", E[0]);
      borderBlock.retrieveData("Ey", E[1]);
      borderBlock.retrieveData("Ez", E[2]);

      dx = borderBlock.getContext().getDx()[2];
      break;
#endif
  }

  makeCoeff();

}

void CPMLBorderHCurrent::stepSchemeInit(double dt)
{
  stepScheme(0.5*dt);
}

void CPMLBorderHCurrent::stepScheme(double /* dt */)
{
  Index low  = Psi[0].getLo();
  Index high = Psi[0].getHi();

  Grid &Psi0 = Psi[0];
  Grid &Psi1 = Psi[1];
  Field &E0 = E[0];
  Field &E1 = E[1];

  Index ind, indn;

  for (ind[0]=low[0]; ind[0]<=high[0]; ++ind[0]) {
#ifndef HUERTO_ONE_DIM
    for (ind[1]=low[1]; ind[1]<=high[1]; ++ind[1]) {
#endif
#ifdef HUERTO_THREE_DIM
      for (ind[2]=low[2]; ind[2]<=high[2]; ++ind[2]) {
#endif
        int j = ind[dim];
        Index indp(ind);
        ++indp[dim];

        Psi0[ind]
          = bCoeff(j)*Psi0[ind]
            + cCoeff(j)*(E1[indp]-E1[ind])/dx;

        Psi1[ind]
          = bCoeff(j)*Psi1[ind]
            - cCoeff(j)*(E0[indp]-E0[ind])/dx;
#ifdef HUERTO_THREE_DIM
      }
#endif
#ifndef HUERTO_ONE_DIM
    }
#endif
  }
}


