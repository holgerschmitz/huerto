/*
 * em_fields.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef HUERTO_EM_FIELDS_H
#define HUERTO_EM_FIELDS_H

#include "../types.hpp"
#include "../simulation/simulation_context.hpp"
#include "../simulation/initialiser.hpp"

#include <schnek/variables/blockcontainer.hpp>

class MPulse;

/**
 * @brief A container block for electromagnetic fields
 *
 * Multiple sets of fields can be defined
 */
class EMFields :
        public schnek::ChildBlock<EMFields>,
        public SimulationEntity
{
  private:
    /// The electric field in V/m
    schnek::Array<InitialsedField<double>, 3> E;

    /// The magnetic field in Tesla
    schnek::Array<InitialsedField<double>, 3> B;

    /**
     * Fill the field values from the expressions provided in the setup file
     */
    void fillValues();
  public:

    /**
     * @brief Constructor taking an optional parent block
     *
     * Constructor is compatible with the schnek::Block constructor
     *
     */
    EMFields(schnek::pBlock parent = schnek::pBlock()) : schnek::ChildBlock<EMFields>(parent)
    {}

    /**
     * @brief Register the electromagnetic fields
     *
     * The electromagnetic fields are created and registered with the Block storage.
     */
    void registerData();

    /**
     * Initialise the parameters available through the setup file
     */
    void initParameters(schnek::BlockParameters &parameters);

    /**
     * Make preInit lifecycle hook public
     */
    void preInit();

    /**
     * Initialise the simulation data
     */
    void init();
};


#ifdef HUERTO_ONE_DIM
static const Stagger exStaggerYee(true );
static const Stagger eyStaggerYee(false);
static const Stagger ezStaggerYee(false);
static const Stagger bxStaggerYee(false);
static const Stagger byStaggerYee(true );
static const Stagger bzStaggerYee(true );
#endif

#ifdef HUERTO_TWO_DIM
static const Stagger exStaggerYee(true,  false);
static const Stagger eyStaggerYee(false, true );
static const Stagger ezStaggerYee(false, false);
static const Stagger bxStaggerYee(false, true );
static const Stagger byStaggerYee(true,  false);
static const Stagger bzStaggerYee(true,  true );
#endif

#ifdef HUERTO_THREE_DIM
static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true );
static const Stagger bxStaggerYee(false, true,  true );
static const Stagger byStaggerYee(true,  false, true );
static const Stagger bzStaggerYee(true,  true,  false);
#endif


#endif
