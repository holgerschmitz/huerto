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


struct FDTD_FieldContainer {
    Vector dx;
    double dt;
    Field Ex, Ey, Ez;
    Field Bx, By, Bz;
    Field Jx, Jy, Jz;
    Grid1d KappaDx;

#ifndef HUERTO_ONE_DIM
    Grid1d KappaDy;
#endif

#ifdef HUERTO_THREE_DIM
    Grid1d KappaDz;
#endif
};

#ifdef HUERTO_ONE_DIM
struct FDTD_StepE : public FDTD_FieldContainer {
  SCHNEK_INLINE void operator()(Index pos) {
    auto i = pos[0];
    double kappaEdx = KappaDx(i) * dx[0];

    Ex(i) += - dt * Jx(i) / eps_0;

    Ey(i) += dt * (
          clight2 * (
          - (Bz(i) - Bz(i-1)) / kappaEdx
          )
        - Jy(i) / eps_0
      );

    Ez(i) += dt * (
          clight2 * (
            (By(i) - By(i-1)) / kappaEdx
          )
        - Jz(i) / eps_0
      );
  }
};

struct FDTD_StepB : public FDTD_FieldContainer {
  SCHNEK_INLINE void operator()(Index pos) {
    auto i = pos[0];
    double kappaHdx = KappaDx(i) * dx[0];

    By(i) += dt * (
          (Ez(i+1) - Ez(i)) / kappaHdx
        + Jy(i)
      );

    Bz(i) += dt * (
        - (Ey(i+1) - Ey(i)) / kappaHdx
        + Jz(i)
      );
  }
};
#endif

#ifdef HUERTO_TWO_DIM
struct FDTD_StepE : public FDTD_FieldContainer {
  SCHNEK_INLINE void operator()(Index pos) {
    auto [i, j] = pos;
    double kappaEdx = KappaDx(i)*dx[0];
    double kappaEdy = KappaDy(j)*dx[1];

    Ex(i, j) += dt*(
          clight2*(
            (Bz(i, j) - Bz(i, j-1))/kappaEdy
          )
        - Jx(i, j) / eps_0
      );

    Ey(i, j) += dt*(
          clight2*(
          - (Bz(i, j) - Bz(i-1, j))/kappaEdx
          )
        - Jy(i, j) / eps_0
      );

    Ez(i, j) += dt*(
          clight2*(
            (By(i, j) - By(i-1, j))/kappaEdx
          - (Bx(i, j) - Bx(i, j-1))/kappaEdy
          )
        - Jz(i, j) / eps_0
      );
  }
};

struct FDTD_StepB : public FDTD_FieldContainer {
  SCHNEK_INLINE void operator()(Index pos) {
    auto [i, j] = pos;
    double kappaHdx = KappaDx(i)*dx[0];
    double kappaHdy = KappaDy(j)*dx[1];

    Bx(i, j) += dt*(
        - (Ez(i, j+1) - Ez(i, j))/kappaHdy
        + Jx(i, j)
      );

    By(i, j) += dt*(
          (Ez(i+1, j) - Ez(i, j))/kappaHdx
        + Jy(i, j)
      );

    Bz(i, j) += dt*(
          (Ex(i, j+1) - Ex(i, j))/kappaHdy
        - (Ey(i+1, j) - Ey(i, j))/kappaHdx
        + Jz(i, j)
      );
  }
};
#endif

#ifdef HUERTO_THREE_DIM
struct FDTD_StepE : public FDTD_FieldContainer {
  SCHNEK_INLINE void operator()(Index pos) {
    auto [i, j, k] = pos;
    double kappaEdx = KappaDx(i)*dx[0];
    double kappaEdy = KappaDy(j)*dx[1];
    double kappaEdz = KappaDz(k)*dx[2];

    Ex(i, j, k) += dt * (
          clight2 * (
            (Bz(i, j, k) - Bz(i, j-1, k)) / kappaEdy
          - (By(i, j, k) - By(i, j, k-1)) / kappaEdz
          )
        - Jx(i, j, k) / eps_0
      );

    Ey(i, j, k) += dt * (
          clight2 * (
            (Bx(i, j, k) - Bx(i, j, k-1)) / kappaEdz
          - (Bz(i, j, k) - Bz(i-1, j, k)) / kappaEdx
          )
        - Jy(i, j, k) / eps_0
      );

    Ez(i, j, k) += dt * (
          clight2 * (
            (By(i, j, k) - By(i-1, j, k)) / kappaEdx
          - (Bx(i, j, k) - Bx(i, j-1, k)) / kappaEdy
          )
        - Jz(i, j, k) / eps_0
      );
  }
};

struct FDTD_StepB : public FDTD_FieldContainer {
  SCHNEK_INLINE void operator()(Index pos) {
    auto [i, j, k] = pos;
    double kappaHdx = KappaDx(i)*dx[0];
    double kappaHdy = KappaDy(j)*dx[1];
    double kappaHdz = KappaDz(k)*dx[2];

    Bx(i, j, k) += dt * (
          (Ey(i, j, k+1) - Ey(i, j, k)) / kappaHdz
        - (Ez(i, j+1, k) - Ez(i, j, k)) / kappaHdy
        + Jx(i, j, k)
      );

    By(i, j, k) += dt * (
          (Ez(i+1, j, k) - Ez(i, j, k)) / kappaHdx
        - (Ex(i, j, k+1) - Ex(i, j, k)) / kappaHdz
        + Jy(i, j, k)
      );

    Bz(i, j, k) += dt*(
          (Ex(i, j+1, k) - Ex(i, j, k)) / kappaHdy
        - (Ey(i+1, j, k) - Ey(i, j, k)) / kappaHdx
        + Jz(i, j, k)
      );
  }
};
#endif


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

void FDTD_Plain::stepSchemeInit(double dt) {
  for (pCurrent current: this->magCurrents) {
    current->stepSchemeInit(dt);
  }

  stepB(0.5*dt);

  for (pCurrent current: this->currents) {
    current->stepSchemeInit(dt);
  }
}

void FDTD_Plain::stepScheme(double dt) {
  for (pCurrent current: this->currents) {
    current->stepScheme(dt);
  }

  stepE(dt);

  for (pCurrent current: this->magCurrents) {
    current->stepScheme(dt);
  }

  stepB(dt);
}

void FDTD_Plain::stepE(double dt)
{
  sumCurrents();
#ifdef HUERTO_ONE_DIM
  FDTD_StepE stepEFunc{getContext().getDx(), dt, Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, KappaEdx};
#endif
#ifdef HUERTO_TWO_DIM
  FDTD_StepE stepEFunc{getContext().getDx(), dt, Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, KappaEdx, KappaEdy};
#endif
#ifdef HUERTO_THREE_DIM
  FDTD_StepE stepEFunc{getContext().getDx(), dt, Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, KappaEdx, KappaEdy, KappaEdz};
#endif
  FieldIterator::forEach(Ex.getInnerRange(), stepEFunc);

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Ex);
  sub.exchange(Ey);
  sub.exchange(Ez);
}

void FDTD_Plain::stepB(double dt)
{
  sumMagCurrents();
#ifdef HUERTO_ONE_DIM
  FDTD_StepB stepBFunc{getContext().getDx(), dt, Ex, Ey, Ez, Bx, By, Bz, Mx, My, Mz, KappaHdx};
#endif
#ifdef HUERTO_TWO_DIM
  FDTD_StepB stepBFunc{getContext().getDx(), dt, Ex, Ey, Ez, Bx, By, Bz, Mx, My, Mz, KappaHdx, KappaHdy};
#endif
#ifdef HUERTO_THREE_DIM
  FDTD_StepB stepBFunc{getContext().getDx(), dt, Ex, Ey, Ez, Bx, By, Bz, Mx, My, Mz, KappaHdx, KappaHdy, KappaHdz};
#endif
  FieldIterator::forEach(Ex.getInnerRange(), stepBFunc);

  schnek::DomainSubdivision<Field> &sub = getContext().getSubdivision();

  sub.exchange(Bx);
  sub.exchange(By);
  sub.exchange(Bz);
}

