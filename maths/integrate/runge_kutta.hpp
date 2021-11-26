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
  private:
    schnek::Array<Field*, dim> fields;
    schnek::Array<std::unique_ptr<Field>, dim> fields_tmp;
  public:
    void setField(int d, Field &field);

    template<typename RHS, typename BC>
    void integrateStep(double dt, const RHS &rhs, BC boundary);

    template<typename RHS, typename BC, typename STEPPER>
    void integrateStep(double dt, const RHS &rhs, BC boundary, STEPPER &stepper);
};

#include "runge_kutta.t"



#endif /* HUERTO_MATHS_INTEGRATE_RUNGE_KUTTA_HPP_ */
