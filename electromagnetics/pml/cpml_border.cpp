#include "cpml_border.hpp"

#include "../../constants.hpp"
#include "../../util/field_util.hpp"
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
  auto &decomposition = getContext().getDecomposition();
  // initialize Kappas here

  Range gdomain = decomposition.getGlobalRange();
  Index glow  = gdomain.getLo();
  Index ghigh = gdomain.getHi();

  std::vector<CPMLBorder::GridLineRegistration> KappaEdk(DIMENSION);
  std::vector<CPMLBorder::GridLineRegistration> KappaHdk(DIMENSION);

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
    ptrdiff_t lo_loE = glow[dim] + 1;
    ptrdiff_t lo_hiE = glow[dim] + thickness;

    ptrdiff_t lo_loH = glow[dim];
    ptrdiff_t lo_hiH = glow[dim] + thickness - 1;

    ptrdiff_t hi_lo = ghigh[dim] - thickness + 1;
    ptrdiff_t hi_hi = ghigh[dim];

#ifdef HUERTO_ONE_DIM
    auto gridContext = decomposition.getGridContext({KappaEdk[dim], KappaHdk[dim]});
#else
    auto gridContext = decomposition.getProjectedGridContext({KappaEdk[dim], KappaHdk[dim]});
#endif
    gridContext.forEach([&](Range1d range, Grid1d &kappaEdk, Grid1d &kappaHdk) {
      Field1dIterator::forEach(range, [=, &kappaEdk, &kappaHdk](Index1d posIndex){
        ptrdiff_t pos = posIndex[0];

        kappaEdk(pos) = 1.0;
        kappaHdk(pos) = 1.0;

        double x, x3;

        if (pos >= lo_loE && pos <= lo_hiE) {
          ptrdiff_t k = pos - lo_loE;
          x = 1 - double(k)/double(thickness);
          x3 = x*x*x;
          kappaEdk(pos) = 1 + (kappaMax - 1)*x3;
        }  

        if (pos >= lo_loH && pos <= lo_hiH) {
          ptrdiff_t k = pos - lo_loH;
          x = 1 - (double(k) - 0.5)/double(thickness);
          x3 = x*x*x;
          kappaHdk(pos) = 1 + (kappaMax - 1)*x3;
        }

        if (pos >= hi_lo && pos <= hi_hi) {
          ptrdiff_t k = hi_hi - pos;
          x = 1 - double(k)/double(thickness);
          x3 = x*x*x;
          kappaEdk(pos) = 1 + (kappaMax - 1)*x3;
          
          double x = 1 - (double(k) - 0.5)/double(thickness);
          double x3 = x*x*x;
          kappaHdk(pos) = 1 + (kappaMax - 1)*x3;
        }
      });
    });
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
  int lowk  = borderRange.getLo(dim);
  int highk = borderRange.getHi(dim);

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

  auto &decomposition = borderBlock.getContext().getDecomposition();
#ifdef HUERTO_ONE_DIM
  bCoeff = decomposition.registerField(schnek::GridFactory<Grid1d>{});
  cCoeff = decomposition.registerField(schnek::GridFactory<Grid1d>{});
#else
  bCoeff = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {dim});
  cCoeff = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {dim});
#endif

  double offset = 0.0;
  lowOffset = 1;

  if (isH)
  {
    offset = 0.5;
    lowOffset = 0;
  }

#ifdef HUERTO_ONE_DIM
  auto gridContext = decomposition.getGridContext({bCoeff, cCoeff});
#else
  auto gridContext = decomposition.getProjectedGridContext({bCoeff, cCoeff});
#endif
  gridContext.forEach([=](Range1d range, Grid1d &b_coeff, Grid1d &c_coeff) {
    Field1dIterator::forEach(range, [=, &b_coeff, &c_coeff](Index1d posIndex){
      ptrdiff_t pos = posIndex[0];
      double b, c;
      if (pos >= lowk && pos <= highk) {
        int k = reverse ? pos - lowk : highk - pos;
        double x = 1 - (double(k)-offset)/double(thickness);
        double x3 = x*x*x;

        double sigma = x3*sigmaMax;
        double kappa = 1 + (kappaMax - 1)*x3;
        double a = aMax*(1-x);

        b = exp(-(sigma/kappa + a));
        c = sigma*(b-1)/(kappa*(sigma+kappa*a));
      } else {
        b = 0.0;
        c = 0.0;
      }

      b_coeff(pos) = b;
      c_coeff(pos) = c;
    });
  });
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
  getBorderExtent(dir, thickness, 1, blow, bhigh, false, borderBlock.getContext(), false);

  auto &decomposition = borderBlock.getContext().getDecomposition();
  borderRange = Range{blow, bhigh};
  Jx = decomposition.registerField(schnek::GridFactory<Grid>{}, borderRange);
  Jy = decomposition.registerField(schnek::GridFactory<Grid>{}, borderRange);
  Jz = decomposition.registerField(schnek::GridFactory<Grid>{}, borderRange);

  setField<Grid>(decomposition, Jx, 0.0);
  setField<Grid>(decomposition, Jy, 0.0);
  setField<Grid>(decomposition, Jz, 0.0);

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
  auto gridContext = borderBlock.getContext().getDecomposition().getGridContext({Psi[0], Psi[1], B[0], B[1], bCoeff, cCoeff});
  gridContext.forEach([&](Range range, Grid &psi0, Grid &psi1, Field &B0, Field &B1, Grid1d &b_coeff, Grid1d &c_coeff) {
    FieldIterator::forEach(range, [=, &psi0, &psi1](Index ind) {
        int j = ind[dim];
        Index indm{ind};
        --indm[dim];

        psi0[ind]
          = b_coeff(j)*psi0[ind]
            + c_coeff(j)*(B1[ind]-B1[indm])/(mu_0*dx);
        psi1[ind]
          = b_coeff(j)*psi1[ind]
            - c_coeff(j)*(B0[ind]-B0[indm])/(mu_0*dx);
    });
  });
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
  getBorderExtent(dir, thickness, 1, blow, bhigh, true, borderBlock.getContext(), false);

  auto &decomposition = borderBlock.getContext().getDecomposition();
  borderRange = Range{blow, bhigh};
  Jx = decomposition.registerField(schnek::GridFactory<Grid>{}, borderRange);
  Jy = decomposition.registerField(schnek::GridFactory<Grid>{}, borderRange);
  Jz = decomposition.registerField(schnek::GridFactory<Grid>{}, borderRange);
  
  setField<Grid>(decomposition, Jx, 0.0);
  setField<Grid>(decomposition, Jy, 0.0);
  setField<Grid>(decomposition, Jz, 0.0);

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
  auto gridContext = borderBlock.getContext().getDecomposition().getGridContext({Psi[0], Psi[1], E[0], E[1], bCoeff, cCoeff});
  gridContext.forEach([&](Range range, Grid &psi0, Grid &psi1, Field &E0, Field &E1, Grid1d &b_coeff, Grid1d &c_coeff) {
    FieldIterator::forEach(range, [=, &psi0, &psi1](Index ind) {
        int j = ind[dim];
        Index indp(ind);
        ++indp[dim];

        psi0[ind]
          = b_coeff(j)*psi0[ind]
            + c_coeff(j)*(E1[indp]-E1[ind])/dx;

        psi1[ind]
          = b_coeff(j)*psi1[ind]
            - c_coeff(j)*(E0[indp]-E0[ind])/dx;
    });
  });
}
