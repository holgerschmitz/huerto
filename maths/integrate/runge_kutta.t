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
void FieldRungeKutta4<rank, dim>::setField(int pos, pField field)
{
  fields[pos] = field;
  fields_tmp[pos] = boost::make_shared<Field>(field);

  fields_fast[pos] = *(fields[pos]);
  fields_tmp_fast[pos] = *(fields_tmp[pos]);
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
    rhs(p, dudt, 0.5);

    for (int d=0; d<dim; ++d)
    {
      (*fields_tmp_fast[d])[p] = (*fields_fast[d])[p] + dt*dudt[d];
    }
  }

  // Swap starred fields and the unstarred fields
  for (int d=0; d<dim; ++d)
  {
    Field &f = *fields_fast[d];
    Field &f_tmp = *fields_tmp_fast[d];
    BOOST_FOREACH(p, range)
    {
      std::swap(f[p], f_tmp[p]);
    }
  }

  boundary();

  // Second step
  BOOST_FOREACH(p, range)
  {
    rhs(p, dudt, 0.0);

    for (int d=0; d<dim; ++d)
    {
      (*fields_tmp_fast[d])[p] =
              0.5*((*fields_fast[d])[p] + (*fields_tmp_fast[d])[p])
              + dt*dudt[d];
    }
  }

  // Copy starred fields back into the unstarred fields
  for (int d=0; d<dim; ++d)
  {
    Field &f = *fields_fast[d];
    Field &f_tmp = *fields_tmp_fast[d];
    BOOST_FOREACH(p, range)
    {
      f[p] = f_tmp[p];
    }
  }

  boundary();

}
