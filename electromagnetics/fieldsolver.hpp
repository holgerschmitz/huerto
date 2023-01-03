/*
 * fieldsolver.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_EM_FIELDSOLVER_H
#define HUERTO_EM_FIELDSOLVER_H

#include "../simulation/simulation_context.hpp"

#include <schnek/variables/blockcontainer.hpp>

/**
 * An interface for field solvers
 *
 * Field solvers form the basis of the EM field simulation.
 *
 * Field solvers should not create their own EM field grids. This is created in
 * the #EMFields class. They should instead acquire their copy through the
 * `retrieveData` method.
 */
class FieldSolver : public schnek::ChildBlock<FieldSolver>,
                    public SimulationEntity {
  public:

    /**
     * Initialise the simulation step.
     *
     * This function will be called once during execution of the simulation.
     * It is not intended to perform any household tasks. Instead is should be used
     * to perform numerical computations that are part of the main scheme but
     * need only be computed before the first time step.
     *
     * For the FDTD scheme this will perform the first half time step of the magnetic
     * field to initialise the correct time staggering of the fields.
     */
    virtual void stepSchemeInit(double dt) = 0;

    /**
     * Perform a simulation time step
     */
    virtual void stepScheme(double dt) = 0;

    /**
     * Initialise the field solver
     */
    void init()
    {
      SimulationEntity::init(this);
    }
};

typedef std::shared_ptr<FieldSolver> pFieldSolver;

#endif
