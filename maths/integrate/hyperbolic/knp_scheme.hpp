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

template<int rank, int dimension, int internalDimension>
struct KurganovNoellePetrovaTypes
{
    static const int dim = dimension;
    static const int internalDim = internalDimension;
    typedef schnek::Field<double, rank, SchnarGridChecker> Field;

    typedef schnek::Array<double, dimension> FluidValues;
    typedef schnek::Array<int, rank> Index;
    typedef schnek::Array<double, internalDimension> InternalVars;
};


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
template<int rank, template<int> class Model>
class KurganovNoellePetrova : public Model<rank>
{
  public:
    typedef KurganovNoellePetrovaTypes<rank, Model<rank>::dim, Model<rank>::internalDim> KNP;
    static const int dim = KNP::dim;
    static const int internalDim = KNP::internalDim;
    typedef typename KNP::Field Field;
    typedef typename KNP::FluidValues FluidValues;
    typedef typename KNP::InternalVars InternalVars;

    typedef schnek::Array<int, rank> Index;
    typedef std::shared_ptr<Field> pField;

  private:
    schnek::Array<Field*, dim> fields;
  public:
    void setField(int d, Field &field);

    double van_leer(double u, double up, double um) const;
    void reconstruct(size_t direction, const Index &pos, int dir, FluidValues &u) const;
    void flux(size_t direction, const Index &pos, FluidValues& flux) const;
    void rhs(Index p, FluidValues &dudt) const;
    void minmax_local_speed(size_t direction,
                            const FluidValues &uW,
                            const FluidValues &uE,
                            const InternalVars &pW,
                            const InternalVars &pE,
                            double &ap,
                            double &am) const;

    void operator()(Index p, FluidValues &dudt) const { rhs(p, dudt); }
};

#include "knp_scheme.t"

#endif /* HUERTO_MATHS_INTEGRATE_HYPERBOLIC_KNP_SCHEME_HPP_ */
