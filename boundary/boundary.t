/*
 * boundary.t
 *
 *  Created on: 16 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */



template<class Field>
void ZeroNeumannBoundary<Field>::applyLo(size_t dim, Field& f)
{
  Index lo_inner = f.getInnerLo();

  Range range(f.getLo(), f.getHi());

  int innerLo = f.getInnerLo()[dim];
  range.getHi()[dim] =  innerLo - 1;

  const Range::iterator itEnd = range.end();

  for (Range::iterator it = range.begin();
       it != itEnd;
       ++it)
  {
    Index src = *it;
    src[dim] = innerLo;
    f[*it] = f[src];
  }
}


template<class Field>
void ZeroNeumannBoundary<Field>::applyHi(size_t dim, Field& f)
{
  Index lo_inner = f.getInnerLo();

  Range range(f.getLo(), f.getHi());

  int innerHi = f.getInnerHi()[dim];
  range.getLo()[dim] =  innerHi + 1;

  const Range::iterator itEnd = range.end();

  for (Range::iterator it = range.begin();
       it != itEnd;
       ++it)
  {
    Index src = *it;
    src[dim] = innerHi;
    f[*it] = f[src];
  }
}


template<class Field>
void ZeroDirichletBoundary<Field>::applyLo(size_t dim, Field& f)
{
  Index lo_inner = f.getInnerLo();

  Range range(f.getLo(), f.getHi());

  int innerLo = f.getInnerLo()[dim];
  range.getHi()[dim] =  innerLo - 1;

  const Range::iterator itEnd = range.end();

  for (Range::iterator it = range.begin();
       it != itEnd;
       ++it)
  {
    f[*it] = 0.0;
  }
}

template<class Field>
void ZeroDirichletBoundary<Field>::applyHi(size_t dim, Field& f)
{
  Index lo_inner = f.getInnerLo();

  Range range(f.getLo(), f.getHi());

  int innerHi = f.getInnerHi()[dim];
  range.getLo()[dim] =  innerHi + 1;

  const Range::iterator itEnd = range.end();

  for (Range::iterator it = range.begin();
       it != itEnd;
       ++it)
  {
    f[*it] = 0.0;
  }
}

template<class Field, size_t dimension>
void BoundaryCondition<Field, dimension>::initParameters(schnek::BlockParameters &blockPars)
{
  blockPars.addArrayParameter("low_", applyLo, Index(0));
  blockPars.addArrayParameter("high_", applyHi, Index(0));
}


template<class Field, size_t dimension>
void BoundaryCondition<Field, dimension>::init() {
    schnek::ChildBlock<BoundaryCondition<Field, dimension> >::init();
    SimulationEntity::init(this);
}

template<class Field, size_t dimension>
void BoundaryCondition<Field, dimension>::apply(schnek::Array<pField, dimension> &fields)
{
  schnek::DomainSubdivision<Field> &subdivision = this->getContext().getSubdivision();

  for (size_t i=0; i<DIMENSION; i++)
  {
    if (bool(applyLo[i]) && subdivision.isBoundLo(i)) applyLoDim(i, fields);
    if (bool(applyHi[i]) && subdivision.isBoundHi(i)) applyHiDim(i, fields);
  }
}

template<class Field, size_t dimension>
void ZeroNeumannBoundaryBlock<Field, dimension>::applyLoDim(int dim, schnek::Array<pField, dimension> &fields)
{
  for (size_t i=0; i<dimension; i++)
  {
    boundary.applyLo(dim, *fields[i]);
  }
}

template<class Field, size_t dimension>
void ZeroNeumannBoundaryBlock<Field, dimension>::applyHiDim(int dim, schnek::Array<pField, dimension> &fields)
{
  for (size_t i=0; i<dimension; i++)
  {
    boundary.applyHi(dim, *fields[i]);
  }
}


template<class Field, size_t dimension>
template<class iterator>
void BoundaryApplicator<Field, dimension>::addBoundaries(iterator start, iterator end) {
  boundaryConditions.insert(boundaryConditions.end(), start, end);
}

template<class Field, size_t dimension>
void BoundaryApplicator<Field, dimension>::setField(int dim, pField f)
{
  fields[dim] = f;
}

template<class Field, size_t dimension>
void BoundaryApplicator<Field, dimension>::setSubdivision(schnek::DomainSubdivision<Field> &sub)
{
  subdivision = &sub;
}

template<class Field, size_t dimension>
void BoundaryApplicator<Field, dimension>::operator()()
{
  if (subdivision != NULL)
  {
    for (size_t i=0; i<dimension; i++)
    {
      subdivision->exchange(*fields[i]);
    }
  }

  for (auto boundary : boundaryConditions)
  {
    boundary->apply(fields);
  }
}
