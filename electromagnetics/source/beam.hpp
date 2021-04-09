/*
 *  beam.hpp
 *
 *  Created on: 23 Dec 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_ELECTROMAGNETICS_SOURCE_BEAM_HPP_
#define HUERTO_ELECTROMAGNETICS_SOURCE_BEAM_HPP_

#include "incsource.hpp"

//===============================================================
//==========  Gaussian Beam Source
//===============================================================

#ifndef HUERTO_ONE_DIM

class GaussBeamSource : public IncidentSource
{
  public:
    ~GaussBeamSource() {}
  protected:
    pCurrent makeECurrent(int distance, Direction dir) override;
    pCurrent makeHCurrent(int distance, Direction dir) override;

    void initParameters(schnek::BlockParameters &blockPars) override;

    bool needsCurrent(Direction dir) override;

    /// The wavevector in 1/m
    Vector k;

    /// The position of the focal point in m
    Vector origin;

    /// The maximum magnetic field amplitude in physical units [Tesla]
    Vector3d H;

    /// The beam waist width at the focal point in m
    double waist;

    /// The length of the initial amplitude rise of the beam in m
    double rise;

    /**
     * The offset of the point of full amplitude of the beam along the beam's
     * axis with respoect to the origin of the beam. In physical units [m].
     */
    double offset;

    /// The relative permittivity (defaults to 1)
    double eps;

    /** 
     * degree to which the beam is circularly polarised
     * 0 means linear polarisation
     * 1 means circular polarisation
     */
    double circ;

    /// Multiplaction coefficient for the exponent in the transverse Gaussian profile (default 1)
    int superGaussian;

};

class GaussBeamSourceEFunc
{
  public:
    GaussBeamSourceEFunc(Direction dir, SimulationContext &context);
    void setParam(Vector k,
                  Vector origin,
                  Vector3d H,
                  double waist,
                  double rise,
                  double offset,
                  double eps,
                  double circ,
                  int superGaussian);

#ifdef HUERTO_TWO_DIM
    Vector3d getHField(int i, int j, double time);
#endif

#ifdef HUERTO_THREE_DIM
    Vector3d getHField(int i, int j, int k, double time);
#endif

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(double) {}

  private:
    /// The wavevector in 1/m
    Vector k;

#ifdef HUERTO_TWO_DIM
    /// A unit vector perpendicular to the wavevector
    Vector kperp;
#endif

#ifdef HUERTO_THREE_DIM
    /// A unit vector perpendicular to the wavevector and H
    Vector kperpA;
    /// A unit vector perpendicular to the wavevector and kperpA
    Vector kperpB;
#endif

    /// The wavenumber (norm of the wavevector) in 1/m
    double kn;

    /// The position of the focal point in m
    Vector origin;

    /// The maximum magnetic field amplitude normalised in Tesla / mu_0
    Vector3d H;

    /**
     * The maximum magnetic field amplitude in the perpendicular direction normalised in Tesla / mu_0.
     * 
     * Creates a wave polarised in the perpendicular direction with pi/2 phase shift to create elliptically
     * polarised waves.
     * 
     * In Tesla / mu_0
     */
    Vector3d Hp;

    /// The simulation time step in s
    double dt;

    /// The frequency of the beam in 1/s
    double om;

    /// The normalised Rayleigh length
    double zr;

    /// The beam waist width at the focal point in m
    double waist;

    /// The length of the initial amplitude rise of the beam normalised to the inverse wavenumber
    double rise;

    /**
     * The offset of the point of full amplitude of the beam along the beam's
     * axis with respoect to the origin of the beam. Normalised to the inverse wavelength.
     **/
    double offset;

    /// The relative permittivity (defaults to 1)
    double eps;

    /// Multiplaction coefficient for the exponent in the transverse Gaussian profile (default 1)
    int superGaussian;

    /// The grid spacing in m
    Vector dx;

    Direction dir;
    SimulationContext context;
};

class GaussBeamSourceHFunc
{
  public:
    GaussBeamSourceHFunc(Direction dir, SimulationContext &context);
    void setParam(Vector k,
                  Vector origin,
                  Vector3d E,
                  double waist,
                  double rise,
                  double offset,
                  double eps,
                  double circ,
                  int superGaussian);

#ifdef HUERTO_TWO_DIM
    Vector3d getEField(int i, int j, double time);
#endif

#ifdef HUERTO_THREE_DIM
    Vector3d getEField(int i, int j, int k, double time);
#endif

    void initSourceFunc(pGrid, pGrid, pGrid) {}
    void setTime(double) {}

  private:
    /// The wavevector in 1/m
    Vector k;

#ifdef HUERTO_TWO_DIM
    /// A unit vector perpendicular to the wavevector
    Vector kperp;
#endif

#ifdef HUERTO_THREE_DIM
    /// A unit vector perpendicular to the wavevector and E
    Vector kperpA;
    /// A unit vector perpendicular to the wavevector and kperpA
    Vector kperpB;
#endif

    /// The wavenumber (norm of the wavevector) in 1/m
    double kn;

    /// The position of the focal point in m
    Vector origin;
    Vector3d E;
    Vector3d Ep;

    /// The frequency of the beam in 1/s
    double om;

    /// The normalised Rayleigh length
    double zr;

    /// The beam waist width at the focal point in m
    double waist;

    /// The length of the initial amplitude rise of the beam normalised to the inverse wavenumber
    double rise;

    /**
     * The offset of the point of full amplitude of the beam along the beam's
     * axis with respoect to the origin of the beam. Normalised to the inverse wavelength.
     **/
    double offset;

    /// The relative permittivity (defaults to 1)
    double eps;

    /// Multiplaction coefficient for the exponent in the transverse Gaussian profile (default 1)
    int superGaussian;

    /// The grid spacing in m
    Vector dx;

    Direction dir;
    SimulationContext context;
};

#endif // not HUERTO_ONE_DIM

#endif /* HUERTO_ELECTROMAGNETICS_SOURCE_BEAM_HPP_ */
