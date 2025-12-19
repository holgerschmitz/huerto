/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "../../constants.hpp"

#include <schnek/grid.hpp>
#include <schnek/tools/literature.hpp>

#include "fdtd_plain.hpp"

#include <memory>

void FDTD_Plain::registerData()
{
  self = this;
  addData("FDTD_Plain", self);

#ifdef HUERTO_ONE_DIM
  addData("KappaEdx", KappaEdx);
  addData("KappaHdx", KappaHdx);
#endif

#ifdef HUERTO_TWO_DIM
  addData("KappaEdx", KappaEdx);
  addData("KappaEdy", KappaEdy);

  addData("KappaHdx", KappaHdx);
  addData("KappaHdy", KappaHdy);
#endif

#ifdef HUERTO_THREE_DIM
  addData("KappaEdx", KappaEdx);
  addData("KappaEdy", KappaEdy);
  addData("KappaEdz", KappaEdz);

  addData("KappaHdx", KappaHdx);
  addData("KappaHdy", KappaHdy);
  addData("KappaHdz", KappaHdz);
#endif
}

void FDTD_Plain::init()
{
  SimulationEntity::init(this);

  schnek::DomainSubdivision<Field> &subdivision = getContext().getSubdivision();
  Index low  = subdivision.getLo();
  Index high = subdivision.getHi();

#ifdef HUERTO_ONE_DIM
  KappaEdx.resize(schnek::Array<ptrdiff_t, 1>(low[0]), schnek::Array<ptrdiff_t, 1>(high[0]));
  KappaHdx.resize(schnek::Array<ptrdiff_t, 1>(low[0]), schnek::Array<ptrdiff_t, 1>(high[0]));
  KappaEdx = 1.0;
  KappaHdx = 1.0;
#endif

#ifdef HUERTO_TWO_DIM
  KappaEdx.resize(schnek::Array<ptrdiff_t, 1>(low[0]), schnek::Array<ptrdiff_t, 1>(high[0]));
  KappaEdy.resize(schnek::Array<ptrdiff_t, 1>(low[1]), schnek::Array<ptrdiff_t, 1>(high[1]));
  KappaEdx = 1.0;
  KappaEdy = 1.0;

  KappaHdx.resize(schnek::Array<ptrdiff_t, 1>(low[0]), schnek::Array<ptrdiff_t, 1>(high[0]));
  KappaHdy.resize(schnek::Array<ptrdiff_t, 1>(low[1]), schnek::Array<ptrdiff_t, 1>(high[1]));
  KappaHdx = 1.0;
  KappaHdy = 1.0;
#endif

#ifdef HUERTO_THREE_DIM
  KappaEdx.resize(schnek::Array<ptrdiff_t, 1>(low[0]), schnek::Array<ptrdiff_t, 1>(high[0]));
  KappaEdy.resize(schnek::Array<ptrdiff_t, 1>(low[1]), schnek::Array<ptrdiff_t, 1>(high[1]));
  KappaEdz.resize(schnek::Array<ptrdiff_t, 1>(low[2]), schnek::Array<ptrdiff_t, 1>(high[2]));
  KappaEdx = 1.0;
  KappaEdy = 1.0;
  KappaEdz = 1.0;

  KappaHdx.resize(schnek::Array<ptrdiff_t, 1>(low[0]), schnek::Array<ptrdiff_t, 1>(high[0]));
  KappaHdy.resize(schnek::Array<ptrdiff_t, 1>(low[1]), schnek::Array<ptrdiff_t, 1>(high[1]));
  KappaHdz.resize(schnek::Array<ptrdiff_t, 1>(low[2]), schnek::Array<ptrdiff_t, 1>(high[2]));
  KappaHdx = 1.0;
  KappaHdy = 1.0;
  KappaHdz = 1.0;
#endif

  retrieveData("Ex", Ex);
  retrieveData("Ey", Ey);
  retrieveData("Ez", Ez);

  retrieveData("Bx", Bx);
  retrieveData("By", By);
  retrieveData("Bz", Bz);

  for (pCurrentBlock current: schnek::BlockContainer<CurrentBlock>::childBlocks())
  {
    current->initCurrents(*this);
  }

  CurrentContainer::init(getContext());

  schnek::LiteratureArticle Yee1966("Yee1966", "Yee, K",
      "Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media.",
      "IEEE Transactions on Antennas and Propagation", "1966", "AP-14", "302--307");

  schnek::LiteratureManager::instance().addReference(
      "Integration of electrodynamic fields uses the Finite Difference Time Domain method.",
      Yee1966);
}

void FDTD_Plain::stepSchemeInit(double dt)
{
  for (pCurrent current: this->magCurrents)
  {
    current->stepSchemeInit(dt);
  }

  stepB(0.5*dt);

  for (pCurrent current: this->currents)
  {
    current->stepSchemeInit(dt);
  }
}

void FDTD_Plain::stepScheme(double dt)
{
  for (pCurrent current: this->currents)
  {
    current->stepScheme(dt);
  }

  stepE(dt);


  for (pCurrent current: this->magCurrents)
  {
    current->stepScheme(dt);
  }

  stepB(dt);
}

#ifdef HUERTO_ONE_DIM
void FDTD_Plain::stepE(double dt, 
                       Field* Ex_, 
                       Field* Ey_, 
                       Field* Ez_, 
                       Field* /* Bx_ */, 
                       Field* By_,
                       Field* Bz_, 
                       Grid1d* KappaEdx_)
{
  Field &Ex = Ex_ != NULL ? *Ex_ : this->Ex;
  Field &Ey = Ey_ != NULL ? *Ey_ : this->Ey;
  Field &Ez = Ez_ != NULL ? *Ez_ : this->Ez;

  Field &By = By_ != NULL ? *By_ : this->By;
  Field &Bz = Bz_ != NULL ? *Bz_ : this->Bz;

  Grid1d &rKappaEdx= KappaEdx_ != NULL ? *KappaEdx_ : this->KappaEdx;


  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
  {
    double jx = this->Jx(i);
    double jy = this->Jy(i);
    double jz = this->Jz(i);

    double kappaEdx = rKappaEdx(i)*dx[0];

    Ex(i) = Ex(i) - dt*jx/eps_0;

    Ey(i) = Ey(i)
      + dt*(
          clight2*(
          - (Bz(i) - Bz(i-1))/kappaEdx
          )
        - jy/eps_0
      );

    Ez(i) = Ez(i)
      + dt*(
          clight2*(
            (By(i) - By(i-1))/kappaEdx
          )
        - jz/eps_0
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_Plain::stepB(double dt)
{
  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i) {
    double jy = this->My(i);
    double jz = this->Mz(i);

    double kappaHdx = KappaHdx(i)*dx[0];

    By(i) = By(i)
      + dt*(
          (Ez(i+1) - Ez(i))/kappaHdx
        + jy
      );

    Bz(i) = Bz(i)
      + dt*(
        - (Ey(i+1) - Ey(i))/kappaHdx
        + jz
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif

#ifdef HUERTO_TWO_DIM
void FDTD_Plain::stepE(double dt, 
                       Field *Ex_, 
                       Field *Ey_, 
                       Field *Ez_, 
                       Field *Bx_, 
                       Field *By_,
                       Field *Bz_, 
                       Grid1d *KappaEdx_, 
                       Grid1d *KappaEdy_) {
  Field &Ex = Ex_ != NULL ? *Ex_ : this->Ex;
  Field &Ey = Ey_ != NULL ? *Ey_ : this->Ey;
  Field &Ez = Ez_ != NULL ? *Ez_ : this->Ez;

  Field &Bx = Bx_ != NULL ? *Bx_ : this->Bx;
  Field &By = By_ != NULL ? *By_ : this->By;
  Field &Bz = Bz_ != NULL ? *Bz_ : this->Bz;

  Grid1d &rKappaEdx= KappaEdx_ != NULL ? *KappaEdx_ : this->KappaEdx;
  Grid1d &rKappaEdy= KappaEdy_ != NULL ? *KappaEdy_ : this->KappaEdy;

  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
  {
    double jx = this->Jx(i,j);
    double jy = this->Jy(i,j);
    double jz = this->Jz(i,j);

    double kappaEdx = rKappaEdx(i)*dx[0];
    double kappaEdy = rKappaEdy(j)*dx[1];

    Ex(i,j) = Ex(i,j)
      + dt*(
          clight2*(
            (Bz(i,j) - Bz(i,j-1))/kappaEdy
          )
        - jx/eps_0
      );

    Ey(i,j) = Ey(i,j)
      + dt*(
          clight2*(
          - (Bz(i,j) - Bz(i-1,j))/kappaEdx
          )
        - jy/eps_0
      );

    Ez(i,j) = Ez(i,j)
      + dt*(
          clight2*(
            (By(i,j) - By(i-1,j))/kappaEdx
          - (Bx(i,j) - Bx(i,j-1))/kappaEdy
          )
        - jz/eps_0
      );

  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_Plain::stepB(double dt)
{
  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
  {
    double jx = this->Mx(i,j);
    double jy = this->My(i,j);
    double jz = this->Mz(i,j);

    double kappaHdx = KappaHdx(i)*dx[0];
    double kappaHdy = KappaHdy(j)*dx[1];

    Bx(i,j) = Bx(i,j)
      + dt*(
        - (Ez(i,j+1) - Ez(i,j))/kappaHdy
        + jx
      );

    By(i,j) = By(i,j)
      + dt*(
          (Ez(i+1,j) - Ez(i,j))/kappaHdx
        + jy
      );

    Bz(i,j) = Bz(i,j)
      + dt*(
          (Ex(i,j+1) - Ex(i,j))/kappaHdy
        - (Ey(i+1,j) - Ey(i,j))/kappaHdx
        + jz
      );

  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif

#ifdef HUERTO_THREE_DIM
void FDTD_Plain::stepE(double dt, 
                       Field *Ex_, 
                       Field *Ey_, 
                       Field *Ez_, 
                       Field *Bx_, 
                       Field *By_,
                       Field *Bz_, 
                       Grid1d *KappaEdx_, 
                       Grid1d *KappaEdy_, 
                       Grid1d *KappaEdz_)
{
  Field &Ex = Ex_ != NULL ? *Ex_ : this->Ex;
  Field &Ey = Ey_ != NULL ? *Ey_ : this->Ey;
  Field &Ez = Ez_ != NULL ? *Ez_ : this->Ez;

  Field &Bx = Bx_ != NULL ? *Bx_ : this->Bx;
  Field &By = By_ != NULL ? *By_ : this->By;
  Field &Bz = Bz_ != NULL ? *Bz_ : this->Bz;

  Grid1d &rKappaEdx= KappaEdx_ != NULL ? *KappaEdx_ : this->KappaEdx;
  Grid1d &rKappaEdy= KappaEdy_ != NULL ? *KappaEdy_ : this->KappaEdy;
  Grid1d &rKappaEdz= KappaEdz_ != NULL ? *KappaEdz_ : this->KappaEdz;



  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double jx = this->Jx(i,j,k);
    double jy = this->Jy(i,j,k);
    double jz = this->Jz(i,j,k);

    double kappaEdx = rKappaEdx(i)*dx[0];
    double kappaEdy = rKappaEdy(j)*dx[1];
    double kappaEdz = rKappaEdz(k)*dx[2];

    Ex(i,j,k) = Ex(i,j,k)
      + dt*(
          clight2*(
            (Bz(i,j,k) - Bz(i,j-1,k))/kappaEdy
          - (By(i,j,k) - By(i,j,k-1))/kappaEdz
          )
        - jx/eps_0
      );

    Ey(i,j,k) = Ey(i,j,k)
      + dt*(
          clight2*(
            (Bx(i,j,k) - Bx(i,j,k-1))/kappaEdz
          - (Bz(i,j,k) - Bz(i-1,j,k))/kappaEdx
          )
        - jy/eps_0
      );

    Ez(i,j,k) = Ez(i,j,k)
      + dt*(
          clight2*(
            (By(i,j,k) - By(i-1,j,k))/kappaEdx
          - (Bx(i,j,k) - Bx(i,j-1,k))/kappaEdy
          )
        - jz/eps_0
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}


void FDTD_Plain::stepB(double dt)
{
  Index low = Ex.getInnerLo();
  Index high = Ex.getInnerHi();

  Vector dx = getContext().getDx();

  sumMagCurrents();

  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
  {
    double jx = this->Mx(i,j,k);
    double jy = this->My(i,j,k);
    double jz = this->Mz(i,j,k);

    double kappaHdx = KappaHdx(i)*dx[0];
    double kappaHdy = KappaHdy(j)*dx[1];
    double kappaHdz = KappaHdz(k)*dx[2];

    Bx(i,j,k) = Bx(i,j,k)
      + dt*(
          (Ey(i,j,k+1) - Ey(i,j,k))/kappaHdz
        - (Ez(i,j+1,k) - Ez(i,j,k))/kappaHdy
        + jx
      );

    By(i,j,k) = By(i,j,k)
      + dt*(
          (Ez(i+1,j,k) - Ez(i,j,k))/kappaHdx
        - (Ex(i,j,k+1) - Ex(i,j,k))/kappaHdz
        + jy
      );

    Bz(i,j,k) = Bz(i,j,k)
      + dt*(
          (Ex(i,j+1,k) - Ex(i,j,k))/kappaHdy
        - (Ey(i+1,j,k) - Ey(i,j,k))/kappaHdx
        + jz
      );
  }

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}
#endif
