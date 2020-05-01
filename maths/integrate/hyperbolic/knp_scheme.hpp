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
 * Model must implement
 * * `flux_function`
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
    void reconstruct(size_t dim, const Index &pos, int dir, FluidValues &u);
    void flux(size_t dim, const Index &pos, FluidValues& flux);
    void rhs(Index p, FluidValues &dudt);
    void minmax_local_speed(size_t dim,
                            const FluidValues &uW,
                            const FluidValues &uE,
                            const InternalVars &pW,
                            const InternalVars &pE,
                            double &ap,
                            double &am);

};

#include "knp_scheme.t"

#endif /* HUERTO_MATHS_INTEGRATE_HYPERBOLIC_KNP_SCHEME_HPP_ */
