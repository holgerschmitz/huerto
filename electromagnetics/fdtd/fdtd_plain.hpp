/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */


#ifndef HUERTO_FDTD_PLAIN_H
#define HUERTO_FDTD_PLAIN_H

#include "../fieldsolver.hpp"
#include "../current.hpp"

#include "../../simulation/simulation_context.hpp"

#include "../../types.hpp"


/**
 * Plain FDTD field solver
 *
 * The solver allows any number of electric and "magnetic" currents as well as
 */
class FDTD_Plain : public FieldSolver,
                   public CurrentContainer,
                   public schnek::BlockContainer<CurrentBlock>
{
  private:

    /**
     * References to the local grid of the electric field components,
     * \f$\mathbf{E}\f$
     */
    pField pEx, pEy, pEz;

    /**
     * References to the local grid of the magnetic field components,
     * \f$\mathbf{E}\f$
     */
    pField pBx, pBy, pBz;

#ifdef HUERTO_ONE_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaEdx;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaHdx;
#endif

#ifdef HUERTO_TWO_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaEdx, pKappaEdy;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaHdx, pKappaHdy;
#endif

#ifdef HUERTO_THREE_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaEdx, pKappaEdy, pKappaEdz;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    pGrid1d pKappaHdx, pKappaHdy, pKappaHdz;
#endif

  public:

    /**
     * Registers helper grids for sharing
     */
    void registerData();

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);

  private:
    void stepD(double dt);
    void stepB(double dt);

};

#endif
