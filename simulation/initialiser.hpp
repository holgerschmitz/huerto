
#ifndef HUERTO_SIMULATION_INITIALISER_CONTEXT_HPP_
#define HUERTO_SIMULATION_INITIALISER_CONTEXT_HPP_

#include "../types.hpp"

#include <schnek/variables/block.hpp>

#include <memory>

/**
 * @brief A named field that can be initialised through the setup file
 * 
 * @tparam T The value type in the field, e.g `double`
 * @tparam D The dimension of the field (default DIMENSION)
 */
template<typename T, int D = DIMENSION>
struct InitialsedField {
    std::shared_ptr<schnek::Field<T, D, HuertoGridChecker>> field;
    schnek::pParameter parameter;
    T value;
};

/**
 * @brief A data reference to named data in a different block
 * 
 * @tparam T The value type of the data, e.g. `Field`
 * @tparam P The pointer type of the data (default T*)
 */
template<typename T, typename P = T*>
class DataReference {
  protected:
    std::string refName;
    P data;
  public:
    void initParameter(schnek::BlockParameters &blockPars, std::string name);
    void init(schnek::Block &block);
    bool isSet() { return refName != ""; }
    T& operator*() { return *data; }
    P operator->() { return data; }
};

template<typename T, typename P>
void DataReference<T, P>::initParameter(schnek::BlockParameters &blockPars, std::string name) {
    blockPars.addParameter(name, &refName, std::string(""));
}

template<typename T, typename P>
void DataReference<T, P>::init(schnek::Block &block) {
    if (isSet()) {
      block.retrieveData(refName, data);
    }
}

#endif // HUERTO_SIMULATION_INITIALISER_CONTEXT_HPP_