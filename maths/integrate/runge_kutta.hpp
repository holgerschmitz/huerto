/*
 * runge_kutta.hpp
 *
 *  Created on: 27 Mar 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_MATHS_INTEGRATE_RUNGE_KUTTA_HPP_
#define HUERTO_MATHS_INTEGRATE_RUNGE_KUTTA_HPP_

#include "../../types.hpp"

#include <schnek/grid/array.hpp>
#include <schnek/grid/domainsubdivision.hpp>

template<int rank, int dim>
class FieldRungeKutta4
{
  public:
    typedef schnek::Field<double, rank, HuertoGridChecker> Field;
    typedef std::shared_ptr<Field> pField;
    typedef std::reference_wrapper<Field> rField;
  private:
    schnek::Array<rField, dim> fields;
    schnek::Array<rField, dim> fields_tmp;
  public:
    void setField(int d, Field &field);

    template<typename RHS, typename BC>
    void integrateStep(double dt, RHS rhs, BC boundary);
};

#include "runge_kutta.t"



#endif /* HUERTO_MATHS_INTEGRATE_RUNGE_KUTTA_HPP_ */
