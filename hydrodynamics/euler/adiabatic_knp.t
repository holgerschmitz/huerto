/*
 * adiabatic_knp.t
 *
 *  Created on: 29 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

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
  double sqrU = 0.0;
  for (size_t i=0; i<dim; ++i)
  {
    sqrU += u[C_M[i]]*u[C_M[i]];
  }

  double internal_energy = std::max(0.0, u[C_E] - 0.5*sqrU/u[C_RHO]);

  p[0] = (adiabaticGamma-1.0)*internal_energy;
}

template<int rank>
void AdiabaticKnpModel<rank>::flux_function(size_t direction,
                                            const FluidValues &u,
                                            const InternalVars &p,
                                            FluidValues &f) const
{
  double rho = u[C_RHO];
  double mdir = u[C_M[direction]];
  double engy = u[C_E];

  f[C_RHO]   = mdir;
  f[C_E]     = (engy + p[0])*mdir/rho;
  for (size_t i=0; i<dim; ++i)
  {
    f[C_M[i]] = mdir*u[C_M[i]]/rho;
  }
  f[C_M[direction]] += p[0];

}

template<int rank>
void AdiabaticKnpModel<rank>::setParameters(double adiabaticGamma, const schnek::Array<double, rank> &dx)
{
  this->adiabaticGamma = adiabaticGamma;
  this->dx = dx;
}

template<int rank>
inline double AdiabaticKnpModel<rank>::speed_cf(double rho, double p)
{
 return (p>0.0)?(0.5*sqrt(4.0*adiabaticGamma*p/rho)):0.0;
}



template<int rank>
void AdiabaticKnp<rank>::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addParameter("gamma", &adiabaticGamma, 1.4);
}


template<int rank>
void AdiabaticKnp<rank>::init()
{
  Super::init();

  this->retrieveData("Rho", Rho);
  scheme.setField(AdiabaticKnpModel<rank>::C_RHO, *Rho);
  integrator.setField(AdiabaticKnpModel<rank>::C_RHO, *Rho);
  boundary.setField(AdiabaticKnpModel<rank>::C_RHO, &(*Rho));

  this->retrieveData("E", E);
  scheme.setField(AdiabaticKnpModel<rank>::C_E, *E);
  integrator.setField(AdiabaticKnpModel<rank>::C_E, *E);
  boundary.setField(AdiabaticKnpModel<rank>::C_E, &(*E));


  for (size_t i=0; i<DIMENSION; ++i)
  {
    this->retrieveData(indexToCoord(i, "M"), M[i]);
    scheme.setField(AdiabaticKnpModel<rank>::C_M[i], *M[i]);
    integrator.setField(AdiabaticKnpModel<rank>::C_M[i], *M[i]);
    boundary.setField(AdiabaticKnpModel<rank>::C_M[i], &(*M[i]));
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
double AdiabaticKnp<rank>::maxDt()
{
  schnek::DomainSubdivision<Field> &subdivision = this->getContext().getSubdivision();

  Field &Rho = *(this->Rho);
  Field &E = *(this->E);
  schnek::Array<Field*, DIMENSION> M;
  for (size_t i=0; i<DIMENSION; ++i)
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

  for (size_t i=1; i<DIMENSION; i++)
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
    for (size_t i=0; i<DIMENSION; ++i)
    {
      u[AdiabaticKnpModel<rank>::C_M[i]]    = (*M[i])[p];
      maxU = std::max(maxU, fabs(u[AdiabaticKnpModel<rank>::C_M[i]]));
    }
    u[AdiabaticKnpModel<rank>::C_E] = E[p];

    InternalVars pressure;
    scheme.calc_internal_vars(u, pressure);

    double v_max = maxU/u[AdiabaticKnpModel<rank>::C_RHO];

    max_speed = std::max(max_speed,scheme.speed_cf(u[AdiabaticKnpModel<rank>::C_RHO], pressure[0])+v_max);
  }

  max_speed = subdivision.maxReduce(max_speed);
  return min_dx/max_speed;
}

template<int rank>
void AdiabaticKnp<rank>::timeStep(double dt)
{
  integrator.integrateStep(dt, scheme, boundary);
}

