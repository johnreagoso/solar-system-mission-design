//  ===============================================================================================================================
/*	PROCEDURE NAME: SSMD_GetPlanetPosVel.cpp

	Purpose- return planetary heliocentric position and velocity at a specified epoch (mean ecliptic J2000)

	===============================================================================================================================
	Input:
	1. User-selection input:					int vPlanetSelect (Planet no., i.e. Earth-3, Mercury-1
												double vJDdate (Julian Day epoch)

	2. Output:									(structure) PlanetPosVel stPlanetPosVel, which contains:
													a. 	double JDdate		(epoch)
													b.  double aPosArray[3] (Position vector in Cartesian)
													c.  double aVelArray[3] (Velocity vector in Cartesian)
	===============================================================================================================================
	External Refs:
	NAME										PURPOSE

	AngleMinimize()								Ensures that the angle returned is less than 2pi
	SSMD_KeplerRad								Computes eccentirc anomaly via Mean anomaly and eccentricity input

	===============================================================================================================================
	Development History
	Name					Date				Description of Change
	Reagoso, J

	These data are to be used as described in the related document titled "Keplerian Elements for Approximate Positions of the
	Major Planets" by E.M. Standish (JPL/Caltech) available from the JPL Solar System Dynamics web site (http://ssd.jpl.nasa.gov/).
	If interested, please refer to http://articles.adsabs.harvard.edu//full/1990A%26A...233..252S/0000270.000.html

	================================================================================================================================
	Keplerian elements and their rates, with respect to the mean ecliptic and equinox of J2000, valid for the time-interval
	1800 AD - 2050 AD. These initial condition and coefficients represent the best Least Squares fit of thousand of over 50000 types
	of planetary observations/measurements. (DE200 has also been the basis for the calculation of Astronomical Almanac planetary tables
	since 1984. While not as accurate as DE405, it is suitable for our mission design needs. */

	//  This m-file will compute the specific planetry ephemeris data for a specific Julian Date.

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

struct PlanetPosVelForJD {
	//struct PlanetPosVel{
	double JDdate;
	double aPosArray[3];
	double aVelArray[3];
} stPlanetPosVel;

struct PlanetParameters
{
	double a_mean_rtn;
	double planet_mu_rtn;
	double radius_planet_rtn;
	double a0rtn, a1rtn, e0rtn, e1rtn, i0rtn, i1rtn, L0rtn, L1rtn, w_bar0rtn, w_bar1rtn, omega0rtn, omega1rtn;

} *stPlanetSelectParamType;

// Function Prototype ====================================================================================================
double Kepler(double* variable1, double* variable2);
double AngleMinimize(double* variable3);

// =======================================================================================================================
double const pi = 3.1415926535897932;	double const vMuSun = 1.32714e11;// km^3/ sec^2 
double AU = 1.495978e8;  //km                      
double vRadConv = pi / 180;
double argmt_periapsis; double vJDdate;

//double planet_mu, radius_planet, a, a0, a1, e0, e, e1, i, i0, i1, L0,L, L1, M, w_bar0, w_bar, w_bar1, omega, omega0, omega1;
double E, h, nu;
double xPrime, yPrime, zPrime, xEcl, yEcl, zEcl;
double f, g, gDot;
double vxEcl, vyEcl, vzEcl;

double vJulianCentEpoch;
struct PlanetPosVelForJD SSMD_GetPlanetPosVel(int vPlanetSelect, double vJDdate, struct PlanetParameters* stPlanetSelectParam)
	//struct PlanetPosVel SSMD_GetPlanetPosVel(int vPlanetSelect, double vJDdate, struct PlanetParameters *stPlanetSelectParam)
{
	double planet_mu, radius_planet, a, a_mean_in, a_mean, a0, a1, e0, e, e1, i, i0, i1, L0, L, L1, M, w_bar0, w_bar, w_bar1, omega, omega0, omega1;

	double earth_radius = 6378.12;  // km
	double earth_mu = 3.986e5;  // km^3/ sec^2

	vJulianCentEpoch = ((vJDdate - 2451545.0) / 36525);

	a_mean_in = stPlanetSelectParam->a_mean_rtn;

	a0 = stPlanetSelectParam->a0rtn;		a1 = stPlanetSelectParam->a1rtn;

	e0 = stPlanetSelectParam->e0rtn;		e1 = stPlanetSelectParam->e1rtn;

	i0 = stPlanetSelectParam->i0rtn;		i1 = stPlanetSelectParam->i1rtn;

	omega0 = stPlanetSelectParam->omega0rtn;	omega1 = stPlanetSelectParam->omega1rtn;

	w_bar0 = stPlanetSelectParam->w_bar0rtn;	w_bar1 = stPlanetSelectParam->w_bar1rtn;

	L0 = stPlanetSelectParam->L0rtn;		L1 = stPlanetSelectParam->L1rtn;

	a_mean = a_mean_in * AU;	a = a0 * AU + a1 * vJulianCentEpoch * AU;

	e = e0 + e1 * vJulianCentEpoch;

	//  JPL planetary ephemeris data are provided in degrees/degrees per century..
	//  convert to radians as well:

	i = vRadConv * i0 + vRadConv * i1 * vJulianCentEpoch;

	omega = vRadConv * omega0 + vRadConv * omega1 * (vJulianCentEpoch);        omega = AngleMinimize(&omega);

	w_bar = vRadConv * w_bar0 + vRadConv * w_bar1 * vJulianCentEpoch;         w_bar = AngleMinimize(&w_bar);

	argmt_periapsis = w_bar - omega;										argmt_periapsis = AngleMinimize(&argmt_periapsis);

	L = vRadConv * L0 + vRadConv * L1 * vJulianCentEpoch;						L = AngleMinimize(&L);

	M = (L - w_bar);														M = AngleMinimize(&M);

	// Compute planetary heliocentric vectors =================================================================================================

	E = Kepler(&M, &e);
	E = AngleMinimize(&E);

	nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)); //nu in radians

	nu = AngleMinimize(&nu);

	h = sqrt(vMuSun * a * (1 - pow(e, 2)));

	xPrime = a * (cos(E) - e);
	yPrime = a * sqrt(1 - pow(e, 2)) * sin(E);
	zPrime = 0.00;


	xEcl = (cos(argmt_periapsis) * cos(omega) - sin(argmt_periapsis) * sin(omega) * cos(i)) * xPrime + (-1 * sin(argmt_periapsis) * cos(omega) - cos(argmt_periapsis) * sin(omega) * cos(i)) * yPrime;
	yEcl = (cos(argmt_periapsis) * sin(omega) + sin(argmt_periapsis) * cos(omega) * cos(i)) * xPrime + (-1 * sin(argmt_periapsis) * sin(omega) + cos(argmt_periapsis) * cos(omega) * cos(i)) * yPrime;
	zEcl = sin(argmt_periapsis) * sin(i) * xPrime + cos(argmt_periapsis) * sin(i) * yPrime;


	vxEcl = (-1 * vMuSun / h) * cos(omega) * sin(argmt_periapsis + nu) - (vMuSun / h) * e * cos(omega) * sin(argmt_periapsis) + (-vMuSun / h) * sin(omega) * cos(i) * cos(argmt_periapsis + nu) + (-vMuSun / h) * sin(omega) * cos(i) * e * cos(argmt_periapsis);
	vyEcl = (-1 * vMuSun / h) * sin(omega) * sin(argmt_periapsis + nu) + (-1 * vMuSun / h) * e * sin(omega) * sin(argmt_periapsis) + (vMuSun / h) * cos(omega) * cos(i) * cos(argmt_periapsis + nu) + (vMuSun / h) * cos(omega) * cos(i) * e * cos(argmt_periapsis);
	vzEcl = (vMuSun / h) * cos(argmt_periapsis + nu) * sin(i) + (vMuSun / h) * e * cos(argmt_periapsis) * sin(i);


	stPlanetPosVel.aPosArray[0] = xEcl;		stPlanetPosVel.aPosArray[1] = yEcl;		stPlanetPosVel.aPosArray[2] = zEcl;
	stPlanetPosVel.aVelArray[0] = vxEcl;	stPlanetPosVel.aVelArray[1] = vyEcl;	stPlanetPosVel.aVelArray[2] = vzEcl;

	stPlanetPosVel.JDdate = vJDdate;

	return stPlanetPosVel;
}

double AngleMinimize(double* vAngle)
{
	double vNoAngleRev, vWholeMultiple, Nangle;
	double vAngleOutput;

	if (abs(*vAngle) > (2 * pi))
	{

		vNoAngleRev = abs(*vAngle) / (2 * pi);

		vWholeMultiple = floor(vNoAngleRev);

		Nangle = vNoAngleRev - vWholeMultiple;

		if (*vAngle < 0) {
			Nangle = Nangle * -1;
		}
		else {
			Nangle = Nangle * 1;
		}

		vAngleOutput = Nangle * 2 * pi;
	}
	else
	{
		vAngleOutput = *vAngle * 1;
	}

	return vAngleOutput;
}