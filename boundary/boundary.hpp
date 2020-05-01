/*
 * boundary.hpp
 *
 *  Created on: 16 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_BOUNDARY_BOUNDARY_HPP_
#define HUERTO_BOUNDARY_BOUNDARY_HPP_

#include "../simulation/simulation_context.hpp"

#include <schnek/variables/block.hpp>
#include <schnek/variables/blockcontainer.hpp>

#include <cstdlib>

template<class Field>
class ZeroNeumannBoundary
{
  public:
    void applyLo(size_t dim, Field &f);
    void applyHi(size_t dim, Field &f);
};

template<class Field>
class ZeroDirichletBoundary
{
  public:
    void applyLo(size_t dim, Field &f);
    void applyHi(size_t dim, Field &f);
};

template<class Field, size_t dimension>
class BoundaryCondition : public schnek::ChildBlock<BoundaryCondition>
{
  private:
    Index applyLo;
    Index applyHi;

    SimulationContext &context;
  public:
    typedef Field* pField;

    BoundaryCondition(SimulationContext &context);
    virtual ~BoundaryCondition();

    void initParameters(schnek::BlockParameters &blockPars) override;
    void apply(schnek::Array<pField, dimension> fields);

    virtual void applyLoDim(schnek::Array<pField, dimension> fields) = 0;
    virtual void applyHiDim(schnek::Array<pField, dimension> fields) = 0;
};

template<class Field, size_t dimension>
class ZeroNeumannBoundaryBlock : public BoundaryCondition<Field, dimension>
{
  private:
    ZeroNeumannBoundary<Field> boundary;
  public:
    void applyLoDim(schnek::Array<pField, dimension> fields) override;
    void applyHiDim(schnek::Array<pField, dimension> fields) override;
};

#include "boundary.t"

#endif /* HUERTO_BOUNDARY_BOUNDARY_HPP_ */
