/*
 * runge_kutta.t
 *
 *  Created on: 27 Mar 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */


#include <memory>
#include <cmath>

template<size_t rank, size_t dim>
void FieldRungeKuttaHeun<rank, dim>::setField(size_t d, Field &field)
{
  fields[d] = &field;
  fields_tmp[d] = std::unique_ptr<Field>(new Field(field));
}


template<size_t rank, size_t dim>
template<typename RHS, typename BC>
void FieldRungeKuttaHeun<rank, dim>::integrateStep(double dt, const RHS &rhs, BC boundary)
{
  typename Field::IndexType lo = fields[0]->getInnerLo();
  typename Field::IndexType hi = fields[0]->getInnerHi();
  typename schnek::Range<int, rank> range(lo, hi);

  schnek::Array<double, dim> dudt;

  // First step
  for(auto p: range)
  {
    rhs(p, dudt, 0.0);

    for (size_t d=0; d<dim; ++d)
    {
      (*fields_tmp[d])[p] = (*fields[d])[p] + dt*dudt[d];
    }
  }

  // Swap starred fields and the unstarred fields
  for (size_t d=0; d<dim; ++d)
  {
    Field &f = *fields[d];
    Field &f_tmp = *fields_tmp[d];
    for(auto p: range)
    {
      std::swap(f[p], f_tmp[p]);
    }
  }

  boundary();

  // Second step
  for(auto p: range)
  {
    rhs(p, dudt, dt);

    for (size_t d=0; d<dim; ++d)
    {
      (*fields_tmp[d])[p] =
              0.5*((*fields[d])[p] + (*fields_tmp[d])[p]
              + dt*dudt[d]);
    }
  }

  // Copy starred fields back into the unstarred fields
  for (size_t d=0; d<dim; ++d)
  {
    Field &f = *fields[d];
    Field &f_tmp = *fields_tmp[d];
    for(auto p: range)
    {
      f[p] = f_tmp[p];
    }
  }

  boundary();
}

template<size_t rank, size_t dim>
template<typename RHS, typename BC, typename STEPPER>
void FieldRungeKuttaHeun<rank, dim>::integrateStep(double dt, const RHS &rhs, BC boundary, STEPPER &stepper)
{
  typename Field::IndexType lo = fields[0]->getInnerLo();
  typename Field::IndexType hi = fields[0]->getInnerHi();
  typename schnek::Range<int, rank> range(lo, hi);

  schnek::Array<double, dim> dudt;

  // First step
  stepper.step(0);
  for(auto p: range)
  {
    rhs(p, dudt, 0.0);

    for (size_t d=0; d<dim; ++d)
    {
      (*fields_tmp[d])[p] = (*fields[d])[p] + dt*dudt[d];
    }
  }

  // Swap starred fields and the unstarred fields
  for (size_t d=0; d<dim; ++d)
  {
    Field &f = *fields[d];
    Field &f_tmp = *fields_tmp[d];
    for(auto p: range)
    {
      std::swap(f[p], f_tmp[p]);
    }
  }

  boundary();

  // Second step
  stepper.step(1);
  for(auto p: range)
  {
    rhs(p, dudt, dt);

    for (size_t d=0; d<dim; ++d)
    {
      (*fields_tmp[d])[p] =
              0.5*((*fields[d])[p] + (*fields_tmp[d])[p]
              + dt*dudt[d]);
    }
  }

  // Copy starred fields back into the unstarred fields
  for (size_t d=0; d<dim; ++d)
  {
    Field &f = *fields[d];
    Field &f_tmp = *fields_tmp[d];
    for(auto p: range)
    {
      f[p] = f_tmp[p];
    }
  }

  boundary();
}