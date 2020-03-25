/*
 * em_fields.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef HUERTO_EM_FIELDS_H
#define HUERTO_EM_FIELDS_H

#include "../../types.hpp"

#include <schnek/variables/blockcontainer.hpp>

class MPulse;

/**
 * @brief A container block for electromagnetic fields
 *
 * Multiple sets of fields can be defined
 */
class EMFields : public schnek::ChildBlock<EMFields>
{
  private:
    friend class MPulse;
    /// The x-component of the electric field in V/m
    pField pEx;
    /// The y-component of the electric field in V/m
    pField pEy;
    /// The z-component of the electric field in V/m
    pField pEz;
    /// The x-component of the magnetic field in Tesla
    pField pBx;
    /// The y-component of the magnetic field in Tesla
    pField pBy;
    /// The z-component of the magnetic field in Tesla
    pField pBz;
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
