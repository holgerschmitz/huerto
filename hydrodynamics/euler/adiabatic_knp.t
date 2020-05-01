/*
 * adiabatic_knp.t
 *
 *  Created on: 29 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */


template<int rank>
inline double AdiabaticKnpModel<rank>::flow_speed(size_t direction, const FluidValues &u, const InternalVars &p)
{
  return u[C_M[direction]] / u[C_RHO];
}

template<int rank>
inline double AdiabaticKnpModel<rank>::sound_speed(const FluidValues &u, const InternalVars &p)
{
  return (p[0]>0.0)?(0.5*sqrt(4.0*adiabaticGamma*p/u[C_RHO])):0.0;
}

template<int rank>
inline void AdiabaticKnpModel<rank>::calc_internal_vars(const FluidValues &u, InternalVars &p)
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
                                            FluidValues &f)
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
void AdiabaticKnp<rank>::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addParameter("gamma", &adiabaticGamma, 1.4);
}
