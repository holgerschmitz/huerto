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
    typedef schnek::Field<double, rank, SchnarGridChecker> Field;
    typedef boost::shared_ptr<Field> pField;
  private:
    schnek::Array<pField, dim> fields;
    schnek::Array<pField, dim> fields_tmp;

    schnek::Array<Field*, dim> fields_fast;
    schnek::Array<Field*, dim> fields_tmp_fast;
  public:
    void setField(int pos, pField field);

    template<typename RHS, typename BC>
    void integrateStep(double dt, RHS rhs, BC boundary);
};

#include "runge_kutta.t"



#endif /* HUERTO_MATHS_INTEGRATE_RUNGE_KUTTA_HPP_ */
