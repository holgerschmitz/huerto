#include "../types.hpp"

#include <memory>

template<typename T>
struct InitialsedField {
    std::shared_ptr<schnek::Field<T, DIMENSION, SchnarGridChecker>> field;
    schnek::pParameter parameter;
    T value;
};

