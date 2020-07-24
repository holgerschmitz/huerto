/*
 * simulation_context.hpp
 *
 *  Created on: 25 Mar 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_SIMULATION_SIMULATION_CONTEXT_HPP_
#define HUERTO_SIMULATION_SIMULATION_CONTEXT_HPP_

#include "../types.hpp"

#include <schnek/variables/block.hpp>

/**
 * A simulation context representing a regular grid domain simulated through time
 */
class SimulationContext
{
  protected:
    /**
     * The global grid size, \f$N\f$
     *
     * This quantity is specified through the `N` parameter in the setup file
     */
    Index gridSize;

    /**
     * The physical size of the global grid, \f$L\f$
     *
     * This quantity is specified through the `L` parameter in the setup file
     */
    Vector size;

    /**
     * The grid spacing
     *
     * Calculated as \f$dx_i = L_i / N_i\f$
     */
    Vector dx;

    /**
     * The time step, \f$dt\f$, of the simulation
     *
     * The time step is calculated as
     *
     * \f[
     * dt = f_{\mathrm{CLF}} \frac{\min_i(dx_i)}{c}
     * \f]
     */
    double dt;

    /**
     * The total physical simulation time
     *
     * This quantity is specified through the `tMax` parameter in the setup file
     */
    double tMax;

    /**
     * The current physical simulation time
     */
    double time;

    /**
     * The position vector that is exposed in the setup file and can be used for initialisation
     */
    Vector x;

    /**
     * The parser parameters associated with the position `x`
     */
    schnek::Array<schnek::pParameter, DIMENSION> x_parameters;

    /**
     * The Cartesian MPI subdivision of the grids containing the electromagnetic field
     */
    std::shared_ptr<schnek::DomainSubdivision<Field>> subdivision;

  public:
    /**
     * Get the global maximum grid index, #globalMax
     */
    Index getGridSize() { return gridSize; }

    /**
     * Get the grid spacing #dx
     */
    Vector getDx() { return dx; }

    /**
     * Get the time step, #dt
     */
    double getDt() { return dt; }

    /**
     * Get the physical size of the global grid, \f$L\f$
     */
    Vector getSize() { return size; }

    /**
     * Get the current simulation time
     */
    double getTime() { return time; }

    /**
     * Get a reference to the independent position vector
     */
    Vector &getX() { return x; }

    /**
     * The parser parameters associated with the position vector
     */
    schnek::Array<schnek::pParameter, DIMENSION> &getXParameter() { return x_parameters; }

    /**
     * Get the grid subdivision
     */
    schnek::DomainSubdivision<Field> &getSubdivision() { return *subdivision; };

    /**
     * Initialise the global parameters exposed in the setup file
     *
     * This function needs to be called explicitly by the simulation that is inheriting from the
     * simulation context
     */
    void initParameters(schnek::BlockParameters &blockPars) {
      x_parameters = blockPars.addArrayParameter("", x, schnek::BlockParameters::readonly);
    }

    void init() {

    }
};

/**
 * A simulation entity that has access to a simulation context
 */
class SimulationEntity
{
  private:
    SimulationContext *context;

  protected:
    void init(schnek::Block *entity)
    {
      context = NULL;

      while (context == NULL && entity != NULL)
      {
        context = dynamic_cast<SimulationContext*>(entity);
        entity = entity->getParent().get();
      }
    }

    SimulationContext &getContext() { return *context; }
};



#endif /* HUERTO_SIMULATION_SIMULATION_CONTEXT_HPP_ */
