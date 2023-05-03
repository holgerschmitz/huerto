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
/**
 * A base class for boundary conditions
 *
 * Field: the field type to apply boundary conditions to
 * dimension: the number of fields
 */
template<class Field, size_t dimension>
class BoundaryCondition : public schnek::ChildBlock<BoundaryCondition<Field, dimension> >,
                          public SimulationEntity
{
  private:
    Index applyLo;
    Index applyHi;
  public:

    virtual ~BoundaryCondition() {};

    void initParameters(schnek::BlockParameters &blockPars) override;
    void init() override;

    void apply(schnek::Array<Field, dimension> &fields);

    virtual void applyLoDim(int dim, schnek::Array<Field, dimension> &fields) = 0;
    virtual void applyHiDim(int dim, schnek::Array<Field, dimension> &fields) = 0;
};

template<class Field, size_t dimension>
class ZeroNeumannBoundaryBlock : public BoundaryCondition<Field, dimension>
{
  private:
    ZeroNeumannBoundary<Field> boundary;
  public:
    void applyLoDim(int dim, schnek::Array<Field, dimension> &fields) override;
    void applyHiDim(int dim, schnek::Array<Field, dimension> &fields) override;
};

template<class Field, size_t dimension>
class BoundaryApplicator
{
  private:
    schnek::Array<Field, dimension> fields;
    std::list<std::shared_ptr<BoundaryCondition<Field, dimension>>> boundaryConditions;
    schnek::DomainSubdivision<Field> *subdivision;
  public:
    BoundaryApplicator() : subdivision(NULL) {}
    template<class iterator>
    void addBoundaries(iterator start, iterator end);
    void setField(int dim, Field &f);
    void setSubdivision(schnek::DomainSubdivision<Field> &subdivision);
    void operator()();
};

#include "boundary.t"

#endif /* HUERTO_BOUNDARY_BOUNDARY_HPP_ */
