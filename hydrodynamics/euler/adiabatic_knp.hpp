/*
 * adiabatic_knp.hpp
 *
 *  Created on: 16 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_HYDRODYNAMICS_EULER_ADIABATIC_KNP_HPP_
#define HUERTO_HYDRODYNAMICS_EULER_ADIABATIC_KNP_HPP_

#include "../hydro_solver.hpp"

#include "../../types.hpp"
#include "../../maths/integrate/hyperbolic/knp_scheme.hpp"

template<int rank>
class AdiabaticKnpModel
{
  public:
    static const int dim = rank + 2;
    static const int internalDim = 1;
    typedef KurganovNoellePetrova<rank, dim, AdiabaticKnpModel> KNP;
    typedef typename KNP::Field Field;
    typedef typename KNP::rField rField;
    typedef typename KNP::FluidValues FluidValues;
    typedef typename KNP::InternalVars InternalVars;

  private:
    double adiabaticGamma;

    static const int C_RHO = 0;
    static const int C_E   = 1;
    static const int C_M[]  = {2, 3, 4};
  protected:
    double flow_speed(size_t direction, const FluidValues &u, const InternalVars &p);
    double sound_speed(const FluidValues &u, const InternalVars &p);
    void calc_internal_vars(const FluidValues &u, InternalVars &p);
    void flux_function(size_t direction, const FluidValues &u, const InternalVars &p, FluidValues &f);
  public:
    void setAdiabaticGamma(double adiabaticGamma);
};


template<int rank>
class AdiabaticKnp : public HydroSolver<typename AdiabaticKnpModel<rank>::Field, AdiabaticKnpModel<rank>::dim>
{
  public:
    static const int dim = AdiabaticKnpModel<rank>::dim;
  private:
    KurganovNoellePetrova<rank, dim, AdiabaticKnpModel> scheme;
    FieldRungeKutta4<rank, dim> integrator;

    double adiabaticGamma;

    pField Rho;
    schnek::Array<pField, DIMENSION> M;
    pField E;
  public:
    /**
     * Initialise the parameters available through the setup file
     */
    void initParameters(schnek::BlockParameters &parameters) override;

    void init() override;
    void postInit() override;
    void timeStep(double dt) override;
    double maxDt() override;
};

#include "adiabatic_knp.t"

#endif /* HUERTO_HYDRODYNAMICS_EULER_ADIABATIC_KNP_HPP_ */
