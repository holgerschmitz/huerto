/*
 * fdtd_plain.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 * 
 * This file was initially created as part of the MPulse project before being added to Huerto
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
    FDTD_Plain *self;
    
    /**
     * The local grid of the electric field components,
     * \f$\mathbf{E}\f$
     */
    Field Ex, Ey, Ez;

    /**
     * The local grid of the magnetic field components,
     * \f$\mathbf{E}\f$
     */
    Field Bx, By, Bz;

#ifdef HUERTO_ONE_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    Grid1d KappaEdx;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    Grid1d KappaHdx;
#endif

#ifdef HUERTO_TWO_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    Grid1d KappaEdx, KappaEdy;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    Grid1d KappaHdx, KappaHdy;
#endif

#ifdef HUERTO_THREE_DIM
    /**
     * Stretch factors for the grid cell size of the electric fields; used by CPML
     * schemes.
     */
    Grid1d KappaEdx, KappaEdy, KappaEdz;

    /**
     * Stretch factors for the grid cell size of the magnetic fields; used by CPML
     * schemes.
     */
    Grid1d KappaHdx, KappaHdy, KappaHdz;
#endif

  public:

    /**
     * Registers helper grids for sharing
     */
    void registerData();

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);

#ifdef HUERTO_ONE_DIM
    void stepE(double dt, 
               Field *Ex = NULL, 
               Field *Ey = NULL, 
               Field *Ez = NULL, 
               Field *Bx = NULL, 
               Field *By = NULL, 
               Field *Bz = NULL, 
               Grid1d *KappaEdx = NULL);
#endif

#ifdef HUERTO_TWO_DIM
    void stepE(double dt, 
               Field *Ex = NULL, 
               Field *Ey = NULL, 
               Field *Ez = NULL, 
               Field *Bx = NULL, 
               Field *By = NULL, 
               Field *Bz = NULL, 
               Grid1d *KappaEdx = NULL, 
               Grid1d *KappaEdy = NULL);
#endif

#ifdef HUERTO_THREE_DIM
    void stepE(double dt, 
               Field *Ex = NULL, 
               Field *Ey = NULL, 
               Field *Ez = NULL, 
               Field *Bx = NULL, 
               Field *By = NULL, 
               Field *Bz = NULL, 
               Grid1d *KappaEdx = NULL, 
               Grid1d *KappaEdy = NULL, 
               Grid1d *KappaEdz = NULL);
#endif

    void stepB(double dt);
};

#endif
