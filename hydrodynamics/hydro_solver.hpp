/*
 * hydro_solver.hpp
 *
 *  Created on: 24 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_HYDRODYNAMICS_HYDRO_SOLVER_HPP_
#define HUERTO_HYDRODYNAMICS_HYDRO_SOLVER_HPP_

#include "../boundary/boundary.hpp"
#include "../simulation/simulation_context.hpp"

#include <schnek/variables/blockcontainer.hpp>

/**
 * Base class for hydro solvers
 */
class HydroSolver :
        public schnek::ChildBlock<HydroSolver>,
        public SimulationEntity
{
  public:
    typedef std::shared_ptr<HydroSolver> pSolver;

    virtual ~HydroSolver() {}
    void init() override { SimulationEntity::init(this); }

    virtual double maxDt() = 0;
    virtual void timeStep(double dt) = 0;
};

typedef std::shared_ptr<HydroSolver> pHydroSolver;

#endif /* HUERTO_HYDRODYNAMICS_HYDRO_SOLVER_HPP_ */
