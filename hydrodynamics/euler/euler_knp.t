/*
 * adiabatic_knp.t
 *
 *  Created on: 29 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#include "../../constants.hpp"

#include <schnek/tools/literature.hpp>

template<int rank>
inline double EulerKnpModel<rank>::flow_speed(size_t direction, const FluidValues &u, const InternalVars &p) const
{
  return u[C_M[direction]] / u[C_RHO];
}

template<int rank>
inline double EulerKnpModel<rank>::sound_speed(const FluidValues &u, const InternalVars &p) const
{
  return (p[0]>0.0)?(0.5*sqrt(4.0*adiabaticGamma*p[0]/u[C_RHO])):0.0;
}

template<int rank>
inline void EulerKnpModel<rank>::calc_internal_vars(const FluidValues &u, InternalVars &p) const
{
  double sqrU = 0.0;
  for (size_t i=0; i<rank; ++i)
  {
    sqrU += u[C_M[i]]*u[C_M[i]];
  }

  double internal_energy = std::max(0.0, u[C_E] - 0.5*sqrU/u[C_RHO]);

  p[0] = (adiabaticGamma-1.0)*internal_energy;
}

template<int rank>
void EulerKnpModel<rank>::flux_function(size_t direction,
                                            const FluidValues &u,
                                            const InternalVars &p,
                                            FluidValues &f) const
{
  double rho = u[C_RHO];
  double mdir = u[C_M[direction]];
  double engy = u[C_E];

  f[C_RHO]   = mdir;
  f[C_E]     = (engy + p[0])*mdir/rho;
  for (size_t i=0; i<rank; ++i)
  {
    f[C_M[i]] = mdir*u[C_M[i]]/rho;
  }
  f[C_M[direction]] += p[0];

}

template<int rank>
void EulerKnpModel<rank>::setParameters(double adiabaticGamma, const schnek::Array<double, rank> &dx)
{
  this->adiabaticGamma = adiabaticGamma;
  this->dx = dx;
}

template<int rank>
void EulerKnp<rank>::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addParameter("gamma", &adiabaticGamma, 1.4);
}


template<int rank>
void EulerKnp<rank>::init()
{
  Super::init();
  boundary.setSubdivision(this->getContext().getSubdivision());

  this->retrieveData("Rho", Rho);
  scheme.setField(EulerKnpModel<rank>::C_RHO, Rho);
  integrator.setField(EulerKnpModel<rank>::C_RHO, Rho);
  boundary.setField(EulerKnpModel<rank>::C_RHO, &Rho);

  this->retrieveData("E", E);
  scheme.setField(EulerKnpModel<rank>::C_E, E);
  integrator.setField(EulerKnpModel<rank>::C_E, E);
  boundary.setField(EulerKnpModel<rank>::C_E, &E);


  for (size_t i=0; i<rank; ++i)
  {
    this->retrieveData(indexToCoord(i, "M"), M[i]);
    scheme.setField(EulerKnpModel<rank>::C_M[i], M[i]);
    integrator.setField(EulerKnpModel<rank>::C_M[i], M[i]);
    boundary.setField(EulerKnpModel<rank>::C_M[i], &M[i]);
  }

  auto boundaries = schnek::BlockContainer<BoundaryCondition<Field, dim> >::childBlocks();
  boundary.addBoundaries(boundaries.begin(), boundaries.end());

  dx = this->getContext().getDx();
  scheme.setParameters(adiabaticGamma, dx);

  schnek::LiteratureArticle Kurganov2001("Kurganov2001", "A. Kurganov and S. Noelle and G. Petrova",
      "Semidiscrete central-upwind schemes for hyperbolic conservation laws and Hamilton--Jacobi equations",
      "SIAM J. Sci. Comput.", "2001", "23", "707");

  schnek::LiteratureManager::instance().addReference(
      "Semidiscrete central-upwind scheme for hyperbolic conservation laws", Kurganov2001);
}

template<int rank>
double EulerKnp<rank>::maxDt()
{
  schnek::DomainSubdivision<Field> &subdivision = this->getContext().getSubdivision();

  Index lo = Rho.getInnerLo();
  Index hi = Rho.getInnerHi();
  Range range(lo, hi);
  Range::iterator range_end = range.end();

  FluidValues u;
  double max_speed = 0.0;

  double min_dx = dx[0];

  for (size_t i=1; i<rank; i++)
  {
    min_dx = std::min(min_dx, dx[i]);
  }

  for (Range::iterator it = range.begin();
       it != range_end;
       ++it)
  {
    const Index &p = *it;
    u[EulerKnpModel<rank>::C_RHO] = Rho[p];

    double maxU = 0.0;
    for (size_t i=0; i<rank; ++i)
    {
      u[EulerKnpModel<rank>::C_M[i]]    = (*M[i])[p];
      maxU = std::max(maxU, fabs(u[EulerKnpModel<rank>::C_M[i]]));
    }
    u[EulerKnpModel<rank>::C_E] = E[p];

    InternalVars pressure = -1.0;
    scheme.calc_internal_vars(u, pressure);

    double v_max = maxU/u[EulerKnpModel<rank>::C_RHO];

    max_speed = std::max(max_speed,scheme.sound_speed(u, pressure) + v_max);
  }

  max_speed = subdivision.maxReduce(max_speed);
  return min_dx/max_speed;
}

template<int rank>
void EulerKnp<rank>::timeStep(double dt)
{
  integrator.integrateStep(dt, scheme, boundary);
}

