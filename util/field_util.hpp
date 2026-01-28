/*
 * field_util.hpp
 *
 *  Created on: 09 Jan 2026
 *      Author: Holger Schmitz
 */
#include <schnek/macros.hpp>

#include "../types.hpp"

template<typename FieldType, typename IndexType = Index>
struct SetField {
  FieldType field;
  typename FieldType::value_type val;
  SCHNEK_INLINE void operator()(typename FieldType::IndexType pos) {
    field[pos] = val;
  }
};

template<typename FieldType, typename Decomposition>
void setField(Decomposition &decomposition, schnek::GridRegistration fieldReg, typename FieldType::value_type val) {
    auto gridContext = decomposition.getGridContext({fieldReg});
    gridContext.forEach([&](Range range, FieldType &field) {
        SetField<FieldType> set{field, val};
        FieldIterator::forEach(range, set);
    });
}

template<typename FieldType, typename Decomposition, size_t rank>
void setField1D(Decomposition &decomposition, schnek::ProjectedGridRegistration<1, rank> fieldReg, typename FieldType::value_type val, size_t dim) {
    auto gridContext = decomposition.getGridContext({fieldReg});
    gridContext.forEach([&](Range range, FieldType &field) {
        Range1d range1d{Index1d{range.getLo(dim)}, Index1d{range.getHi(dim)}};
        SetField<FieldType> set{field, val};
        Field1dIterator::forEach(range1d, set);
    });
}

template<typename AccType, typename FieldType, typename IndexType = Index>
struct AddToField {
  AccType accumulator;
  FieldType field;
  SCHNEK_INLINE void operator()(typename FieldType::IndexType pos) {
    accumulator[pos] += field[pos];
  }
};

template<typename AccType, typename FieldType, typename Decomposition>
void addToField(Decomposition &decomposition, schnek::GridRegistration accReg, schnek::GridRegistration fieldReg) {
    auto gridContext = decomposition.getGridContext({fieldReg, accReg});
    gridContext.forEach([&](Range range, FieldType &field, AccType &acc) {
        AddToField<AccType, FieldType> sum{acc, field};
        FieldIterator::forEach(range, sum);
    });
}
