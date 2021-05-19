/*
 * adiabatic_knp.hpp
 *
 *  Created on: 16 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_HYDRODYNAMICS_EULER_KNP_HPP_
#define HUERTO_HYDRODYNAMICS_EULER_KNP_HPP_

#include "../hydro_solver.hpp"

#include "../../types.hpp"
#include "../../maths/integrate/hyperbolic/knp_scheme.hpp"
#include "../../maths/integrate/runge_kutta.hpp"

/**
 * The model for the Euler equations
 */
template<int rank>
class EulerKnpModel
{
  public:
    typedef KurganovNoellePetrovaTypes<rank, rank + 2, 1> KNP;
    static const int dim = KNP::dim;
    static const int internalDim = KNP::internalDim;
    typedef typename KNP::Field Field;
    typedef typename KNP::FluidValues FluidValues;
    typedef typename KNP::InternalVars InternalVars;

    static const int C_RHO = 0;
    static const int C_E   = 1;
    static const int C_M[];
  private:
    double adiabaticGamma;
    schnek::Array<double, rank> dx;
  protected:
    double flow_speed(size_t direction, const FluidValues &u, const InternalVars &p) const;
    void flux_function(size_t direction, const FluidValues &u, const InternalVars &p, FluidValues &f) const;

    const schnek::Array<double, rank> &getDx() const { return dx; }
  public:
    double sound_speed(const FluidValues &u, const InternalVars &p) const;
    void calc_internal_vars(const FluidValues &u, InternalVars &p) const;
    void setParameters(double adiabaticGamma, const schnek::Array<double, rank> &dx);
};


template<int rank>
class EulerKnp : 
        public HydroSolver,
        public schnek::BlockContainer<BoundaryCondition<
          typename EulerKnpModel<rank>::Field, 
          EulerKnpModel<rank>::dim
        >>
{
  public:
    static const int dim = EulerKnpModel<rank>::dim;
    typedef typename EulerKnpModel<rank>::Field Field;
    typedef std::shared_ptr<Field> pField;
    typedef typename EulerKnpModel<rank>::FluidValues FluidValues;
    typedef typename EulerKnpModel<rank>::InternalVars InternalVars;
  private:
    typedef HydroSolver Super;

    KurganovNoellePetrova<rank, EulerKnpModel> scheme;
    FieldRungeKutta4<rank, dim> integrator;
    BoundaryApplicator<Field, dim> boundary;

    double adiabaticGamma;

    pField Rho;
    schnek::Array<pField, rank> M;
    pField E;

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

#include "euler_knp.t"

#endif /* HUERTO_HYDRODYNAMICS_EULER_KNP_HPP_ */
