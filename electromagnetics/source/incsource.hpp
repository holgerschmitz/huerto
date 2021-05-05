/*
 * incsource.hpp
 *
 *  Created on: 30 Nov 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_EM_SOURCE_INCSOURCE_H
#define HUERTO_EM_SOURCE_INCSOURCE_H

#include "../../types.hpp"
#include "../current.hpp"

class IncidentSourceCurrent;

//===============================================================
//==========  Base Classes
//===============================================================

/**
 * An abstract base class for sources on the simulation boundary
 *
 * Sources for time-dependent electromagnetic fields are represented by
 * transverse currents on the boundary.
 */
class IncidentSource : public CurrentBlock {
  public:

    /**
     * Initialise any currents
     *
     * This will check all 6 faces of the boundary box. If #needCurrent indicates
     * that currents are needed on that boundary, #makeECurrent and #makeHCurrent
     * are used to create the transverse currents.
     *
     * The current created in this way will usually extend #IncidentSourceCurrent
     */
    void initCurrents(CurrentContainer &container);

  protected:

    /**
     * Create an #IncidentSourceECurrent for a single boundary at a given distance
     */
    virtual pCurrent makeECurrent(int distance, Direction dir) = 0;

    /**
     * Create an #IncidentSourceBCurrent for a single boundary at a given distance
     */
    virtual pCurrent makeHCurrent(int distance, Direction dir) = 0;

    /**
     * This can be implemented to indicate that some boundaries do not need currents
     */
    virtual bool needsCurrent(Direction dir) { return true; }

    /**
     * Initialise the setup parameters
     */
    void initParameters(schnek::BlockParameters &blockPars);
  private:

    /**
     * The distance of the currents from the outer boundary of the simulation box
     */
    int distance;
};

/**
 * An incident source current
 *
 * Incident source currents are intended to be used with #IncidentSource current
 * blocks. They define a transverse current on a boundary plane.
 */
class IncidentSourceCurrent : public Current {
  public:

    /**
     * Initialise the incident source current
     *
     * @param distance_  the distance of the injection plane from the simulation
     *                   border
     * @param dir_       the direction of the border from which the wave should
     *                   be injected
     * @param isH_       a flag indicating whether the current is a magnetic
     *                   current.
     */
    IncidentSourceCurrent(int distance, Direction dir, bool isH, SimulationContext &context);

  protected:
    /**
     * A flag indicating whether the wave should propagate in the negative
     * coordinate direction
     */
    bool reverse;

    /**
     * The distance of the injection plane from the simulation
     * border
     */
    int distance;

    /**
     * The propagation axis
     */
    int dim;

    /**
     * The transverse axes
     */
    int transverse1, transverse2;

    /**
     * The direction of the border from which the wave should be injected
     */
    Direction dir;

    /**
     * A flag indicating whether the current is a magnetic current.
     */
    bool isH;
    int lowOffset;
    int highOffset;

    /**
     * The grid spacing
     */
    Vector dx;

    /**
     * The time step
     */
    double dt;

    /**
     * References to the transverse components of the current from #Current.pJx,
     * #Current.pJy, and #Current.pJz
     */
    pGrid pJ[2];

    /**
     * The grid spacing in the direction normal to the injection plane
     */
    double dN;

    /**
     * The simulation context
     */
    SimulationContext &context;
};

/**
 * Specialisation of #IncidentSourceCurrent for electric currents
 *
 * This function is templated with a source function that should provide the
 * value of the **magnetic field**. `SourceFunc` must expose a method
 * `getHField` to obtain this field and a method `initSourceFunc` to initialise
 * the source function.
 */
template<class SourceFunc>
class IncidentSourceECurrent : public IncidentSourceCurrent, public SourceFunc {
  public:

    /**
     * Create a new IncidentSourceECurrent instance
     *
     * @param distance_  the distance of the injection plane from the simulation
     *                   border
     * @param dir_       the direction of the border from which the wave should
     *                   be injected
     */
    IncidentSourceECurrent(int distance, Direction dir, SimulationContext &context);

    void init();
    void stepSchemeInit(double dt);
    void stepScheme(double dt);
};

/**
 * Specialisation of #IncidentSourceCurrent for magnetic currents
 *
 * This function is templated with a source function that should provide the
 * value of the **electric field**. `SourceFunc` must expose a method
 * `getEField` to obtain this field and a method `initSourceFunc` to initialise
 * the source function.
 */
template<class SourceFunc>
class IncidentSourceHCurrent : public IncidentSourceCurrent, public SourceFunc {
  public:

    /**
     * Create a new IncidentSourceBCurrent instance
     *
     * @param distance_  the distance of the injection plane from the simulation
     *                   border
     * @param dir_       the direction of the border from which the wave should
     *                   be injected
     */
    IncidentSourceHCurrent(int distance, Direction dir, SimulationContext &context);

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);
};

template<class FieldFunc>
class GenericIncidentSourceESource : public FieldFunc {
  public:
    GenericIncidentSourceESource(Direction dir, SimulationContext &context);
    void setGenericParam(Vector k, Vector3d H, const Vector &origin, double eps);
#ifdef HUERTO_ONE_DIM
    Vector3d getHField(int i, double time);
#endif

#ifdef HUERTO_TWO_DIM
    Vector3d getHField(int i, int j, double time);
#endif

#ifdef HUERTO_THREE_DIM
    Vector3d getHField(int i, int j, int k, double time);
#endif

  private:
    Vector k;
    Vector3d H;

    double dt;
    double om;
    Vector origin;
    Vector dx;
    SimulationContext &context;
};

template<class FieldFunc>
class GenericIncidentSourceHSource : public FieldFunc {
  public:
    GenericIncidentSourceHSource(Direction dir, SimulationContext &context);
    void setGenericParam(Vector k, Vector3d E, const Vector &origin, double eps);
#ifdef HUERTO_ONE_DIM
    Vector3d getEField(int i, double time);
#endif

#ifdef HUERTO_TWO_DIM
    Vector3d getEField(int i, int j, double time);
#endif

#ifdef HUERTO_THREE_DIM
    Vector3d getEField(int i, int j, int k, double time);
#endif

  private:
    Vector k;
    Vector3d E;

    double dt;
    double om;
    Vector origin;
    Vector dx;
    SimulationContext &context;
};

#include "incsource.t"


#endif
