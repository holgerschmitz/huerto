/*
 * knp_scheme.hpp
 *
 *  Created on: 16 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_MATHS_INTEGRATE_HYPERBOLIC_KNP_SCHEME_HPP_
#define HUERTO_MATHS_INTEGRATE_HYPERBOLIC_KNP_SCHEME_HPP_

#include "../../../types.hpp"

#include <schnek/grid/array.hpp>

/**
 * Kurganov, Noelle, Petrova solver for conservative equations
 *
 * The solver defines rhs() that calculates the fluxes and returns the right hand side of the time
 * advance equation \f$\partial_t M = \nabla F\f$. These can then be used in a time stepping algorithm
 * such as Runge Kutta
 *
 * The `rank` template argument defines the dimensional rank of the simulation domain.
 *
 * The `dim` template argument defines the number of fields, i.e. the dimensionality of the conservation equation.
 *
 * The solver does not define the fluxes but is templated with a Model.
 *
 * Model must implement
 * * `flux_function` a function that calculates
 * * `calc_internal_vars` should calculate any additional variables from the local fluid values,
 *   for example pressure, temperature, etc
 * * `flow_speed` calculate the fluid flow speed
 * * `sound_speed` calculate the fluid sound speed, i.e. the speed of the fastest travelling wave
 */
template<int rank, int dim, template<int, int> class Model>
class KurganovNoellePetrova : public Model<rank, dim>
{
  public:
    typedef schnek::Field<double, rank, SchnarGridChecker> Field;
    typedef std::reference_wrapper<Field> rField;

    typedef schnek::Array<double, dim> FluidValues;
    typedef schnek::Array<int, rank> Index;
    typedef schnek::Array<double, Model<rank, dim>::internalDim> InternalVars;

    schnek::Array<rField, dim> fields;
  private:
    void setField(int d, Field &field);

    double van_leer(double u, double up, double um);
    void reconstruct(size_t direction, const Index &pos, int dir, FluidValues &u);
    void flux(size_t direction, const Index &pos, FluidValues& flux);
    void rhs(Index p, FluidValues &dudt);
    void minmax_local_speed(size_t direction,
                            const FluidValues &uW,
                            const FluidValues &uE,
                            const InternalVars &pW,
                            const InternalVars &pE,
                            double &ap,
                            double &am);

    void operator()(Index p, FluidValues &dudt) { rhs(p, dudt); }
};

#include "knp_scheme.t"

#endif /* HUERTO_MATHS_INTEGRATE_HYPERBOLIC_KNP_SCHEME_HPP_ */
