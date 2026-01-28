/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#include "../../constants.hpp"
#include "../../util/field_util.hpp"

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

  auto &decomposition = getContext().getDecomposition();

#ifdef HUERTO_ONE_DIM
  KappaEdx = decomposition.registerField(schnek::GridFactory<Grid1d>{});
  KappaHdx = decomposition.registerField(schnek::GridFactory<Grid1d>{});
  setField<Grid1d>(decomposition, KappaEdx, 1.0);
  setField<Grid1d>(decomposition, KappaHdx, 1.0);
#endif

#ifdef HUERTO_TWO_DIM
  KappaEdx = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {0});
  KappaEdy = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {1});
  setField1D<Grid1d>(decomposition, KappaEdx, 1.0, 0U);
  setField1D<Grid1d>(decomposition, KappaEdy, 1.0, 1U);

  KappaHdx = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {0});
  KappaHdy = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {1});
  setField1D<Grid1d>(decomposition, KappaHdx, 1.0, 0U);
  setField1D<Grid1d>(decomposition, KappaHdy, 1.0, 1U);
#endif

#ifdef HUERTO_THREE_DIM
  KappaEdx = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {0});
  KappaEdy = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {1});
  KappaEdz = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {2});
  setField1D<Grid1d>(decomposition, KappaEdx, 1.0, 0U);
  setField1D<Grid1d>(decomposition, KappaEdy, 1.0, 1U);
  setField1D<Grid1d>(decomposition, KappaEdz, 1.0, 2U);

  KappaHdx = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {0});
  KappaHdy = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {1});
  KappaHdz = decomposition.registerFieldProjection(schnek::GridFactory<Grid1d>{}, {2});
  setField1D<Grid1d>(decomposition, KappaHdx, 1.0, 0U);
  setField1D<Grid1d>(decomposition, KappaHdy, 1.0, 1U);
  setField1D<Grid1d>(decomposition, KappaHdz, 1.0, 2U);
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
  auto &decomposition = getContext().getDecomposition();

#ifdef HUERTO_ONE_DIM
  auto gridContext = decomposition.getGridContext({Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, KappaEdx});
  gridContext.forEach([&](
    Range range, 
    Field &ex, Field &ey, Field &ez, 
    Field &bx, Field &by, Field &bz,
    Field &jx, Field &jy, Field &jz,
    Grid1d &kappaEdx) {
      FDTD_StepE stepEFunc{getContext().getDx(), dt, ex, ey, ez, bx, by, bz, jx, jy, jz, kappaEdx};
      FieldIterator::forEach(range, stepEFunc);
  });
#endif
#ifdef HUERTO_TWO_DIM
  auto gridContext = decomposition.getGridContext({Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, KappaEdx, KappaEdy});
  gridContext.forEach([&](
    Range range, 
    Field &ex, Field &ey, Field &ez, 
    Field &bx, Field &by, Field &bz,
    Field &jx, Field &jy, Field &jz,
    Grid1d &kappaEdx, Grid1d &kappaEdy) {
      FDTD_StepE stepEFunc{getContext().getDx(), dt, ex, ey, ez, bx, by, bz, jx, jy, jz, kappaEdx, kappaEdy};
      FieldIterator::forEach(range, stepEFunc);
  });
#endif
#ifdef HUERTO_THREE_DIM
  auto gridContext = decomposition.getGridContext({Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, KappaEdx, KappaEdy, KappaEdz});
  gridContext.forEach([&](
    Range range, 
    Field &ex, Field &ey, Field &ez, 
    Field &bx, Field &by, Field &bz,
    Field &jx, Field &jy, Field &jz,
    Grid1d &kappaEdx, Grid1d &kappaEdy, Grid1d &kappaEdz) {
      FDTD_StepE stepEFunc{getContext().getDx(), dt, ex, ey, ez, bx, by, bz, jx, jy, jz, kappaEdx, kappaEdy, kappaEdz};
      FieldIterator::forEach(range, stepEFunc);
  });
#endif

  decomposition.exchange({Ex, Ey, Ez});
}

void FDTD_Plain::stepB(double dt)
{
  sumMagCurrents();
  auto &decomposition = getContext().getDecomposition();

#ifdef HUERTO_ONE_DIM
  auto gridContext = decomposition.getGridContext({Ex, Ey, Ez, Bx, By, Bz, Mx, My, Mz, KappaHdx});
  gridContext.forEach([&](
    Range range, 
    Field &ex, Field &ey, Field &ez, 
    Field &bx, Field &by, Field &bz,
    Field &mx, Field &my, Field &mz,
    Grid1d &kappaHdx) {
      FDTD_StepB stepBFunc{getContext().getDx(), dt, ex, ey, ez, bx, by, bz, mx, my, mz, kappaHdx};
      FieldIterator::forEach(range, stepBFunc);
  });
#endif
#ifdef HUERTO_TWO_DIM
  auto gridContext = decomposition.getGridContext({Ex, Ey, Ez, Bx, By, Bz, Mx, My, Mz, KappaHdx, KappaHdy});
  gridContext.forEach([&](
    Range range, 
    Field &ex, Field &ey, Field &ez, 
    Field &bx, Field &by, Field &bz,
    Field &mx, Field &my, Field &mz,
    Grid1d &kappaHdx, Grid1d &kappaHdy) {
      FDTD_StepB stepBFunc{getContext().getDx(), dt, ex, ey, ez, bx, by, bz, mx, my, mz, kappaHdx, kappaHdy};
      FieldIterator::forEach(range, stepBFunc);
  });
#endif
#ifdef HUERTO_THREE_DIM
  auto gridContext = decomposition.getGridContext({Ex, Ey, Ez, Bx, By, Bz, Mx, My, Mz, KappaHdx, KappaHdy, KappaHdz});
  gridContext.forEach([&](
    Range range, 
    Field &ex, Field &ey, Field &ez, 
    Field &bx, Field &by, Field &bz,
    Field &mx, Field &my, Field &mz,
    Grid1d &kappaHdx, Grid1d &kappaHdy, Grid1d &kappaHdz) {
      FDTD_StepB stepBFunc{getContext().getDx(), dt, ex, ey, ez, bx, by, bz, mx, my, mz, kappaHdx, kappaHdy, kappaHdz};
      FieldIterator::forEach(range, stepBFunc);
  });
#endif

  decomposition.exchange({Bx, By, Bz});
}

