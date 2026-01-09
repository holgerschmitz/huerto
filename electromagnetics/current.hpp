/*
 * current.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: holger
 */


#ifndef HUERTO_EM_CURRENT_H
#define HUERTO_EM_CURRENT_H

#include "em_fields.hpp"

#include "../simulation/simulation_context.hpp"

#include "../types.hpp"

#include <schnek/variables/blockcontainer.hpp>

class Current;
typedef std::shared_ptr<Current> pCurrent;

/**
 * A container for currents that distinguishes between electric and magnetic
 * currents. The currents are added up from all the Current references added
 * to the container.
 *
 * The magnetic currents are current-like fields that enter the Faraday equation
 * much like the electric current enters the Ampere equation. These can be used
 * to provide absorbing boundary layers. They are not real physical currents but
 * mathematical constructs.
 * 
 * Typically an implementation of an electromagnetic FieldSolver will also inherit 
 * from CurrentContainer. This allows the field solver to access the summed up 
 * currents and include them in the field equations.
 */
class CurrentContainer
{
  protected:
    /**
     * A list of pointers to currents
     */
    typedef std::list<pCurrent> CurrentList;

    /**
     * The list of electric currents
     */
    CurrentList currents;

    /**
     * The list of magnetic currents
     */
    CurrentList magCurrents;

    /**
     * Sum of all electric currents
     */
    Field Jx, Jy, Jz;

    /**
     * Sum of all magnetic currents
     */
    Field Mx, My, Mz;

    /**
     * Add all electric currents and store the result in #pJx, #pJy, #pJz
     *
     * Each current can be defined on a separate sub-domain of the local grid
     * domain
     */
    void sumCurrents();

    /**
     * Add all magnetic currents and store the result in #pMx, #pMy, #pMz
     *
     * Each current can be defined on a separate sub-domain of the local grid
     * domain
     */
    void sumMagCurrents();

  public:

    /**
     * Add an electric current to the list
     */
    void addCurrent(pCurrent current);

    /**
     * Add a magnetic current to the list
     */
    void addMagCurrent(pCurrent current);

    /**
     * Initialises the sum of the currents to the local grid domain
     */
    void init(SimulationContext &context);
};

/**
 * An abstract block type that calculates currents
 * 
 * The current block typically will be a child block of a CurrentContainer instance. 
 * The parent block will call #initCurrents function to add currents to the container.
 */
class CurrentBlock :
        public schnek::ChildBlock<CurrentBlock>,
        public SimulationEntity
{
  public:
    virtual ~CurrentBlock() {}

    /**
     * Initialise any currents with the current container
     *
     * This usually involves adding the currents to the container.
     *
     * This should be called by the parent block during the `Block.init` phase
     */
    virtual void initCurrents(CurrentContainer &container) = 0;

    /**
     * Initialises the current block
     */
    void preInit();
};

/**
 * A shared pointer to a current block
 */
typedef std::shared_ptr<CurrentBlock> pCurrentBlock;

/**
 * An abstract base class for currents
 */
class Current
{
  protected:
    /**
     * The grids containing the components of the current
     */
    Grid Jx, Jy, Jz;
  public:
    /**
     * Virtual destructor
     */
    virtual ~Current() {}

    /**
     * Initialise the currents
     *
     * This method will be called during the `Block.init` phase.
     */
    virtual void init() = 0;

    /**
     * Initialise the simulation step.
     *
     * This function will be called once during execution of the simulation.
     * It is not intended to perform any household tasks. Instead is should be used
     * to perform numerical computations that are part of the main scheme but
     * need only be computed before the first time step.
     */
    virtual void stepSchemeInit(double dt) = 0;

    /**
     * Perform a simulation time step
     */
    virtual void stepScheme(double dt) = 0;

    /**
     * Check if the current is valid
     *
     * The current is valid if all the components of the current have been
     * initialised
     */
    bool isValid() const { return (Jx.getSize() > 0) && (Jy.getSize() > 0) && (Jz.getSize() > 0); }

    /**
     * Get the x-component of the current
     */
    const Grid &getJx() { return Jx; }

    /**
     * Get the y-component of the current
     */
    const Grid &getJy() { return Jy; }

    /**
     * Get the z-component of the current
     */
    const Grid &getJz() { return Jz; }
};


#endif
