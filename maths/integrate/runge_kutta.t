/*
 * runge_kutta.t
 *
 *  Created on: 27 Mar 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */


#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <cmath>

template<int rank, int dim>
void FieldRungeKutta4<rank, dim>::setField(int d, Field &field)
{
  fields[d] = field;
  fields_tmp[d] = Field(field);
}


template<int rank, int dim>
template<typename RHS, typename BC>
void FieldRungeKutta4<rank, dim>::integrateStep(double dt, RHS rhs, BC boundary)
{
  typename Field::IndexType lo = fields[0].getInnerLo();
  typename Field::IndexType hi = fields[0].getInnerHi();
  typename Field::RangeType range(lo, hi);

  schnek::Array<double, dim> dudt;
  typename Field::IndexType p;

  // First step
  BOOST_FOREACH(p, range)
  {
    rhs(p, dudt);

    for (int d=0; d<dim; ++d)
    {
      fields_tmp[d][p] = fields[d][p] + dt*dudt[d];
    }
  }

  // Swap starred fields and the unstarred fields
  for (int d=0; d<dim; ++d)
  {
    Field &f = *fields[d];
    Field &f_tmp = *fields_tmp[d];
    BOOST_FOREACH(p, range)
    {
      std::swap(f[p], f_tmp[p]);
    }
  }

  boundary();

  // Second step
  BOOST_FOREACH(p, range)
  {
    rhs(p, dudt);

    for (int d=0; d<dim; ++d)
    {
      fields_tmp[d][p] =
              0.5*(fields[d][p] + fields_tmp[d][p])
              + dt*dudt[d];
    }
  }

  // Copy starred fields back into the unstarred fields
  for (int d=0; d<dim; ++d)
  {
    Field &f = fields[d];
    Field &f_tmp = fields_tmp[d];
    BOOST_FOREACH(p, range)
    {
      f[p] = f_tmp[p];
    }
  }

  boundary();

}
