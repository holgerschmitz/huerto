/*
 * adiabatic_knp.t
 *
 *  Created on: 29 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#include "adiabatic_knp.hpp"

#include "../../constants.hpp"

#include <schnek/tools/literature.hpp>

template<int rank>
inline double AdiabaticKnpModel<rank>::flow_speed(size_t direction, const FluidValues &u, const InternalVars &p) const
{
  return u[C_M[direction]] / u[C_RHO];
}

template<int rank>
inline double AdiabaticKnpModel<rank>::sound_speed(const FluidValues &u, const InternalVars &p) const
{
  return (p[0]>0.0)?(0.5*sqrt(4.0*adiabaticGamma*p[0]/u[C_RHO])):0.0;
}

template<int rank>
inline void AdiabaticKnpModel<rank>::calc_internal_vars(const FluidValues &u, InternalVars &p) const
{
  p[0] = p0*pow(u[C_RHO] / rho0, adiabaticGamma);
}

template<int rank>
void AdiabaticKnpModel<rank>::flux_function(size_t direction,
                                            const FluidValues &u,
                                            const InternalVars &p,
                                            FluidValues &f) const
{
  double rho = u[C_RHO];
  double mdir = u[C_M[direction]];

  f[C_RHO]   = mdir;
  for (size_t i=0; i<rank; ++i)
  {
    f[C_M[i]] = mdir*u[C_M[i]]/rho;
  }
  f[C_M[direction]] += p[0];

}

template<int rank>
void AdiabaticKnpModel<rank>::setParameters(double adiabaticGamma, double rho0, const schnek::Array<double, rank> &dx)
{
  this->adiabaticGamma = adiabaticGamma;
  this->rho0 = rho0;
  this->dx = dx;
}

template<int rank>
void AdiabaticKnp<rank>::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addParameter("gamma", &adiabaticGamma, 1.4);
  parameters.addParameter("rho0", &rho0, 1);
}


template<int rank>
void AdiabaticKnp<rank>::init()
{
  Super::init();
  boundary.setSubdivision(this->getContext().getSubdivision());

  this->retrieveData("Rho", Rho);
  scheme.setField(AdiabaticKnpModel<rank>::C_RHO, *Rho);
  integrator.setField(AdiabaticKnpModel<rank>::C_RHO, *Rho);
  boundary.setField(AdiabaticKnpModel<rank>::C_RHO, &(*Rho));


  for (size_t i=0; i<rank; ++i)
  {
    this->retrieveData(indexToCoord(i, "M"), M[i]);
    scheme.setField(AdiabaticKnpModel<rank>::C_M[i], *M[i]);
    integrator.setField(AdiabaticKnpModel<rank>::C_M[i], *M[i]);
    boundary.setField(AdiabaticKnpModel<rank>::C_M[i], &(*M[i]));
  }

  auto boundaries = schnek::BlockContainer<BoundaryCondition<Field, dim> >::childBlocks();
  boundary.addBoundaries(boundaries.begin(), boundaries.end());

  dx = this->getContext().getDx();
  scheme.setParameters(adiabaticGamma, rho0, dx);

  schnek::LiteratureArticle Kurganov2001("Kurganov2001", "A. Kurganov and S. Noelle and G. Petrova",
      "Semidiscrete central-upwind schemes for hyperbolic conservation laws and Hamilton--Jacobi equations",
      "SIAM J. Sci. Comput.", "2001", "23", "707");

  schnek::LiteratureManager::instance().addReference(
      "Semidiscrete central-upwind scheme for hyperbolic conservation laws", Kurganov2001);
}

template<int rank>
double AdiabaticKnp<rank>::maxDt()
{
  schnek::DomainSubdivision<Field> &subdivision = this->getContext().getSubdivision();

  Field &Rho = *(this->Rho);
  schnek::Array<Field*, rank> M;
  for (size_t i=0; i<rank; ++i)
  {
    M[i] = &(*this->M[i]);
  }

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
    u[AdiabaticKnpModel<rank>::C_RHO]    = Rho[p];

    double maxU = 0.0;
    for (size_t i=0; i<rank; ++i)
    {
      u[AdiabaticKnpModel<rank>::C_M[i]]    = (*M[i])[p];
      maxU = std::max(maxU, fabs(u[AdiabaticKnpModel<rank>::C_M[i]]));
    }

    InternalVars pressure = -1.0;
    scheme.calc_internal_vars(u, pressure);

    double v_max = maxU/u[AdiabaticKnpModel<rank>::C_RHO];

    max_speed = std::max(max_speed,scheme.sound_speed(u, pressure) + v_max);
  }

  max_speed = subdivision.maxReduce(max_speed);
  return min_dx/max_speed;
}

template<int rank>
void AdiabaticKnp<rank>::timeStep(double dt)
{
  integrator.integrateStep(dt, scheme, boundary);
}

