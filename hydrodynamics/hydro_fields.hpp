/*
 * hydro_fields.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef VELLAMO_HYDRO_FIELDS_H
#define VELLAMO_HYDRO_FIELDS_H

#include "../types.hpp"
#include "../simulation/simulation_context.hpp"
#include "../simulation/initialiser.hpp"

#include <schnek/variables/blockcontainer.hpp>

class Vellamo;

/** @brief A container block for hydrodynamic fields
 *
 *  Multiple sets of fields can be defined
 */
class HydroFields :
        public schnek::ChildBlock<HydroFields>,
        public SimulationEntity
{
  private:
    friend class Vellamo;
    /// The fluid mass density
    InitialisedField<double> Rho;
    /// The fluid momentum
    schnek::Array<InitialisedField<double>, DIMENSION> M;
    /// The internal energy of the fluid
    InitialisedField<double> E;

    /**
     * Fill the field values from the expressions provided in the setup file
     */
    void fillValues();
  public:

    /** @brief Constructor taking an optional parent block
     *
     * Constructor is compatible with the schnek::Block constructor
     */
    HydroFields(schnek::pBlock parent = schnek::pBlock()) : schnek::ChildBlock<HydroFields>(parent)
    {}

    /** @brief Register the hydrodynamic fields
     *
     * The hydrodynamic fields are created and registered with the Block storage.
     */
    void registerData();

    /**
     * Initialise the parameters available through the setup file
     */
    void initParameters(schnek::BlockParameters &parameters);

    /**
     * Initialise the simulation data
     */
    void init();
};

//inline void diagField(Field &f, std::string msg) {
//  std::cout << "DIAG " << msg << std::endl;
//  for (int i=f.getLo(0); i<=f.getHi(0); ++i)
//    for (int j=f.getLo(1); j<=f.getHi(1); ++j) {
//      std::cout << i << " " << j << " " << f(i,j) << std::endl;
//    }
//}

#endif
