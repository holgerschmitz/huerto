#include "../types.hpp"

#include <schnek/variables/block.hpp>

#include <memory>

template<typename T>
struct InitialsedField {
    std::shared_ptr<schnek::Field<T, DIMENSION, HuertoGridChecker>> field;
    schnek::pParameter parameter;
    T value;
};

template<typename T, typename P = T*>
class DataReference {
  protected:
    std::string refName;
    P data;
  public:
    void initParameter(schnek::BlockParameters &blockPars, std::string name);
    void init(schnek::Block &block);
    T& operator*() { return *data; }
    P operator->() { return data; }
};

template<typename T, typename P>
void DataReference<T, P>::initParameter(schnek::BlockParameters &blockPars, std::string name) {
    blockPars.addParameter(name, &refName);
}

template<typename T, typename P>
void DataReference<T, P>::init(schnek::Block &block) {
    block.retrieveData(refName, data);
}