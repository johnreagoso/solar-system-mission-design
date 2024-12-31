////  =====================================================================
////  These data are to be used as described in the related document
////  titled "Keplerian Elements for Approximate Positions of the
////  Major Planets" by E.M. Standish (JPL/Caltech) available from
////  the JPL Solar System Dynamics web site (http://ssd.jpl.nasa.gov/).
////  =====================================================================
////  Keplerian elements and their rates, with respect to the Mean Ecliptic
////  and Mean Equinox of 2000 (EME2000), valid for the time-interval 1800 AD - 2050 AD
//
////  This .cpp-file will compute the specific planetry ephemeris data for a
////  specific Julian Date
//
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

struct PlanetParameters {

    double a_mean_rtn;
    double planet_mu_rtn;
    double radius_planet_rtn;
    double a0rtn, a1rtn, e0rtn, e1rtn, i0rtn, i1rtn, L0rtn, L1rtn, w_bar0rtn, w_bar1rtn, omega0rtn, omega1rtn;
    double vMinAtmAltitude; //km

} stPlanetSelectParam;

// Function Prototype ====================================================================================================
// =======================================================================================================================

struct PlanetParameters PlanetSelectOrbitParameters(int vPlanetSelect)
{
    double earth_mu = 3.9865e5; 	// km^3/sec^2
    double earth_radius = 6378;		// km  (Earth equatorial radius)

    double planet_mu, radius_planet, a_mean, a0, a1, e0, e, e1, i0, i1, L0, L, L1, w_bar0, w_bar, w_bar1, omega, omega0, omega1;
    double vPlanetMinAtmAlt; //km

    double vAU = 1.495978e8;  //km 

    switch (vPlanetSelect)
    {
    case 3:		// Earth
        a_mean = 1.0;
        planet_mu = 1 * earth_mu;
        radius_planet = 1 * earth_radius;
        a0 = 1.00000261;
        a1 = 0.00000562;
        e0 = 0.01671123;
        e1 = -0.00004392;
        i0 = -0.00001531;
        i1 = -0.01294668;
        L0 = 100.46457166;
        L1 = 35999.37344981;
        w_bar0 = 102.93768193;
        w_bar1 = 0.32327364;
        omega0 = 0.00;
        omega1 = 0.00;
        vPlanetMinAtmAlt = 6578; //km

        stPlanetSelectParam.a_mean_rtn = a_mean * vAU;

        stPlanetSelectParam.planet_mu_rtn = planet_mu;

        stPlanetSelectParam.radius_planet_rtn = radius_planet;

        stPlanetSelectParam.a0rtn = a0;
        stPlanetSelectParam.a1rtn = a1;
        stPlanetSelectParam.e0rtn = e0;
        stPlanetSelectParam.e1rtn = e1;
        stPlanetSelectParam.i0rtn = i0;
        stPlanetSelectParam.i1rtn = i1;
        stPlanetSelectParam.L0rtn = L0;
        stPlanetSelectParam.L1rtn = L1;
        stPlanetSelectParam.w_bar0rtn = w_bar0;
        stPlanetSelectParam.w_bar1rtn = w_bar1;
        stPlanetSelectParam.omega0rtn = omega0;
        stPlanetSelectParam.omega1rtn = omega1;
        stPlanetSelectParam.vMinAtmAltitude = vPlanetMinAtmAlt;
        break;

    case 2:		//Venus

        a_mean = 0.72333566;
        planet_mu = 0.8149 * earth_mu;
        radius_planet = 0.949 * earth_radius;
        a0 = 0.72333566;
        a1 = 0.00000390;
        e0 = 0.00677672;
        e1 = -0.00004107;
        i0 = 3.39467605;
        i1 = -0.00078890;
        L0 = 181.97909950;
        L1 = 58517.81538729;
        w_bar0 = 131.60246718;
        w_bar1 = 0.00268329;
        omega0 = 76.67984255;
        omega1 = -0.27769418;
        vPlanetMinAtmAlt = radius_planet; //km

        stPlanetSelectParam.a_mean_rtn = a_mean * vAU;

        stPlanetSelectParam.planet_mu_rtn = planet_mu;

        stPlanetSelectParam.radius_planet_rtn = radius_planet;

        stPlanetSelectParam.a0rtn = a0;
        stPlanetSelectParam.a1rtn = a1;
        stPlanetSelectParam.e0rtn = e0;
        stPlanetSelectParam.e1rtn = e1;
        stPlanetSelectParam.i0rtn = i0;
        stPlanetSelectParam.i1rtn = i1;
        stPlanetSelectParam.L0rtn = L0;
        stPlanetSelectParam.L1rtn = L1;
        stPlanetSelectParam.w_bar0rtn = w_bar0;
        stPlanetSelectParam.w_bar1rtn = w_bar1;
        stPlanetSelectParam.omega0rtn = omega0;
        stPlanetSelectParam.omega1rtn = omega1;
        stPlanetSelectParam.vMinAtmAltitude = vPlanetMinAtmAlt;

        break;

    case 4:		//Mars

        a_mean = 1.52371034;
        planet_mu = 0.1074 * earth_mu;
        radius_planet = 0.532 * earth_radius;
        a0 = 1.52371034;
        a1 = 0.00001847;
        e0 = 0.09339410;
        e1 = 0.00007882;
        i0 = 1.84969142;
        i1 = -0.00813131;
        L0 = -4.553343205;
        L1 = 19140.30268499;
        w_bar0 = -23.94362959;
        w_bar1 = 0.44441088;
        omega0 = 49.55953891;
        omega1 = -0.29257343;
        vPlanetMinAtmAlt = radius_planet + 100.0; //km

        stPlanetSelectParam.a_mean_rtn = a_mean * vAU;

        stPlanetSelectParam.planet_mu_rtn = planet_mu;

        stPlanetSelectParam.radius_planet_rtn = radius_planet;

        stPlanetSelectParam.a0rtn = a0;
        stPlanetSelectParam.a1rtn = a1;
        stPlanetSelectParam.e0rtn = e0;
        stPlanetSelectParam.e1rtn = e1;
        stPlanetSelectParam.i0rtn = i0;
        stPlanetSelectParam.i1rtn = i1;
        stPlanetSelectParam.L0rtn = L0;
        stPlanetSelectParam.L1rtn = L1;
        stPlanetSelectParam.w_bar0rtn = w_bar0;
        stPlanetSelectParam.w_bar1rtn = w_bar1;
        stPlanetSelectParam.omega0rtn = omega0;
        stPlanetSelectParam.omega1rtn = omega1;
        stPlanetSelectParam.vMinAtmAltitude = vPlanetMinAtmAlt;

        break;
    }

    //stPlanetSelectParam.a_mean_rtn = a_mean * vAU;

    //stPlanetSelectParam.planet_mu_rtn = planet_mu;

    //stPlanetSelectParam.radius_planet_rtn = radius_planet;

    //stPlanetSelectParam.a0rtn = a0;
    //stPlanetSelectParam.a1rtn = a1;
    //stPlanetSelectParam.e0rtn = e0;
    //stPlanetSelectParam.e1rtn = e1;
    //stPlanetSelectParam.i0rtn = i0;
    //stPlanetSelectParam.i1rtn = i1;
    //stPlanetSelectParam.L0rtn = L0;
    //stPlanetSelectParam.L1rtn = L1;
    //stPlanetSelectParam.w_bar0rtn = w_bar0;
    //stPlanetSelectParam.w_bar1rtn = w_bar1;
    //stPlanetSelectParam.omega0rtn = omega0;
    //stPlanetSelectParam.omega1rtn = omega1;
    //stPlanetSelectParam.vMinAtmAltitude = vPlanetMinAtmAlt;

    return stPlanetSelectParam;
}