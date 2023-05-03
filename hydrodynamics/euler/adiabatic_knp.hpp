/*
 * adiabatic_knp.hpp
 *
 *  Created on: 16 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_HYDRODYNAMICS_ADIABATIC_KNP_HPP_
#define HUERTO_HYDRODYNAMICS_ADIABATIC_KNP_HPP_

#include "../hydro_solver.hpp"

#include "../../types.hpp"
#include "../../maths/integrate/hyperbolic/knp_scheme.hpp"
#include "../../maths/integrate/runge_kutta.hpp"

/**
 * The model for the Euler equations
 */
template<int rank>
class AdiabaticKnpModel
{
  public:
    typedef KurganovNoellePetrovaTypes<rank, rank + 1, 1> KNP;
    static const int dim = KNP::dim;
    static const int internalDim = KNP::internalDim;
    typedef typename KNP::Field Field;
    typedef typename KNP::FluidValues FluidValues;
    typedef typename KNP::InternalVars InternalVars;

    static const int C_RHO = 0;
    static const int C_M[];
  private:
    double adiabaticGamma;

    /**
     * The proportionality factor in the adiabatic equation.
     * 
     * This is really
     * \f$\frac{p_0}{\rho_0^\gamma}\f$
     */
    double p0;
    schnek::Array<double, rank> dx;
  protected:
    double flow_speed(size_t direction, const FluidValues &u, const InternalVars &p) const;
    void flux_function(size_t direction, const FluidValues &u, const InternalVars &p, FluidValues &f) const;

    const schnek::Array<double, rank> &getDx() const { return dx; }
  public:
    double sound_speed(const FluidValues &u, const InternalVars &p) const;
    void calc_internal_vars(const FluidValues &u, InternalVars &p) const;
    void setParameters(double adiabaticGamma, double p0, const schnek::Array<double, rank> &dx);
};


template<int rank>
class AdiabaticKnp : 
        public HydroSolver,
        public schnek::BlockContainer<BoundaryCondition<
          typename AdiabaticKnpModel<rank>::Field, 
          AdiabaticKnpModel<rank>::dim
        >>
{
  public:
    static const int dim = AdiabaticKnpModel<rank>::dim;
    typedef typename AdiabaticKnpModel<rank>::Field Field;
    typedef typename AdiabaticKnpModel<rank>::FluidValues FluidValues;
    typedef typename AdiabaticKnpModel<rank>::InternalVars InternalVars;
  private:
    typedef HydroSolver Super;

    KurganovNoellePetrova<rank, AdiabaticKnpModel> scheme;
    FieldRungeKuttaHeun<rank, dim> integrator;
    BoundaryApplicator<Field, dim> boundary;

    double adiabaticGamma;
    double p0;

    Field Rho;
    schnek::Array<Field, rank> M;

    schnek::Array<double, rank> dx;
  public:
    /**
     * Initialise the parameters available through the setup file
     */
    void initParameters(schnek::BlockParameters &parameters) override;

    void init() override;
    double maxDt() override;
    void timeStep(double dt) override;
};

#include "adiabatic_knp.t"

#endif /* HUERTO_HYDRODYNAMICS_ADIABATIC_KNP_HPP_ */
