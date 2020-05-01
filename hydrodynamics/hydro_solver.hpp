/*
 * hydro_solver.hpp
 *
 *  Created on: 24 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_HYDRODYNAMICS_HYDRO_SOLVER_HPP_
#define HUERTO_HYDRODYNAMICS_HYDRO_SOLVER_HPP_

#include "../boundary/boundary.hpp"

#include <schnek/variables/blockcontainer.hpp>

template<class Field, size_t dimension>
class HydroSolver :
        public schnek::ChildBlock<HydroSolver>,
        public schnek::BlockContainer<BoundaryCondition<Field, dimension>>
{
  public:
    typedef boost::shared_ptr<HydroSolver> pSolver;

    virtual ~HydroSolver() {}
    virtual double maxDt() = 0;
    virtual void timeStep(double dt) = 0;
};



#endif /* HUERTO_HYDRODYNAMICS_HYDRO_SOLVER_HPP_ */
