#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

struct PlanetPosVel {
	double vMJDno;
	double aPosVector1[3];
	double aVelVector[3];

} aPos1, aPos2;

struct	LambertData {

	double aVelVector1[3];
	double aVelVector2[3];
	double vIterationsReq;

} stLambertReturn;

double const pi = 3.1415926535897932;
double const vMaxIterations = 1000;

// Prototype Functions ===================================================================================================================================================================
double vec_magntd(double aArray[3]);
double Fz(double variable, double variable2, double variable3, double variable4, double variable5, double variable6);
double y(double variable, double variable2, double variable3, double variable4);
double dFz(double variable, double variable2, double variable3, double variable4);
double StumpfS(double variable);
double StumpfC(double variable);
// =======================================================================================================================================================================================

int const vMax_iterations = 100;
double const vTol = 1.00e-7;
double const vNmax = 600;

//struct LambertData SSMD_LambertsSolverUnivVariables(double aPosVector1[3], double aPosVector2[3], double vThetaSwept, double vTOF, double vMu)
struct LambertData SSMD_LambertsSolverUnivVariables(struct PlanetPosVel* aPos1, struct PlanetPosVel* aPos2, double* vThetaSweptInput, double* vTOFInput)
{
	double vThetaSwept = *vThetaSweptInput;
	double vTOF = *vTOFInput;
	double vMu = 1.32714e11;

	int ii = 0;
	int kk = 0;
	double vRatio = 1.00;

	double vR1;	double vR2;	double vA;	double z;	double vFz; double dFzDenom;

	double aVelVector1comp[3];
	double aVelVector2comp[3];

	double f;	double g;	double gDot;
	double aPosVector1[3], aPosVector2[3];

	aPosVector1[0] = aPos1->aPosVector1[0];
	aPosVector1[1] = aPos1->aPosVector1[1];
	aPosVector1[2] = aPos1->aPosVector1[2];

	aPosVector2[0] = aPos2->aPosVector1[0];
	aPosVector2[1] = aPos2->aPosVector1[1];
	aPosVector2[2] = aPos2->aPosVector1[2];

	vR1 = vec_magntd(aPosVector1);
	vR2 = vec_magntd(aPosVector2);

	vA = sin(vThetaSwept) * sqrt((vR1 * vR2) / (1 - cos(vThetaSwept)));
	z = 0.010;

	vFz = Fz(vR1, vR2, vA, z, vTOF, vMu);
	// ======================================================================
	// Provides for initial guess to assist with Newton Raphson iteration====

	while (vFz < 0.00)
	{
		z = z + 0.10;
		vFz = Fz(vR1, vR2, vA, z, vTOF, vMu);
		ii = ii + 1;
	}

	// Newton-Raphson Solver to determine (z)=================================
	// =======================================================================
	while ((abs(vRatio) > vTol) && (kk <= vNmax))
	{

		dFzDenom = 1 / dFz(vR1, vR2, vA, z);

		vRatio = Fz(vR1, vR2, vA, z, vTOF, vMu) * dFzDenom;

		z = z - vRatio;
		//cout << z  << endl; 

		kk = kk + 1;
	}
	double vIterationsReq = kk - 1;
	//cout << kk  << endl; 
// ============================================================================
/* f_g Series to determine velocity vectors ===================================
   at position 1 and 2, respectively. ========================================= */

	double vR1denom = 1 / vR1;
	double vR2denom = 1 / vR2;

	f = 1 - y(vR1, vR2, vA, z) * vR1denom;

	g = vA * sqrt(y(vR1, vR2, vA, z) / vMu);

	gDot = 1 - y(vR1, vR2, vA, z) * vR2denom;

	double vGdenom = 1 / g;

	aVelVector1comp[0] = aPosVector2[0] * vGdenom - f * aPosVector1[0] * vGdenom;
	aVelVector1comp[1] = aPosVector2[1] * vGdenom - f * aPosVector1[1] * vGdenom;
	aVelVector1comp[2] = aPosVector2[2] * vGdenom - f * aPosVector1[2] * vGdenom;

	aVelVector2comp[0] = gDot * aPosVector2[0] * vGdenom - aPosVector1[0] * vGdenom;
	aVelVector2comp[1] = gDot * aPosVector2[1] * vGdenom - aPosVector1[1] * vGdenom;
	aVelVector2comp[2] = gDot * aPosVector2[2] * vGdenom - aPosVector1[2] * vGdenom;

	stLambertReturn.aVelVector1[0] = aVelVector1comp[0];
	stLambertReturn.aVelVector1[1] = aVelVector1comp[1];
	stLambertReturn.aVelVector1[2] = aVelVector1comp[2];

	stLambertReturn.aVelVector2[0] = aVelVector2comp[0];
	stLambertReturn.aVelVector2[1] = aVelVector2comp[1];
	stLambertReturn.aVelVector2[2] = aVelVector2comp[2];

	stLambertReturn.vIterationsReq = vIterationsReq;

	return stLambertReturn;
}

double Fz(double vR1, double vR2, double vA, double z, double vTOF, double vMu)
{
	double dummy2;
	dummy2 = pow((y(vR1, vR2, vA, z) / StumpfC(z)), 1.5) * StumpfS(z) + vA * sqrt(y(vR1, vR2, vA, z)) - sqrt(vMu) * (vTOF);

	return dummy2;
}

double y(double vR1, double vR2, double vA, double z)
{
	double dummy3;
	double vCdiv;

	vCdiv = 1 / sqrt(StumpfC(z));

	dummy3 = vR1 + vR2 + vA * z * StumpfS(z) * vCdiv - vA * vCdiv;

	return dummy3;
}

double dFz(double vR1, double vR2, double vA, double z)
{

	double dummy6;
	double vAdiv8 = vAinput / 8;

	if (zInput == 0)
	{
		dummy6 = 0.025 * pow(2, 0.5) * pow(y(vR1, vR2, vA, 0), 1.5) + vAdiv8 * (sqrt(y(vR1, vR2, vA, 0)) + vA * sqrt(1 / 2 / y(vR1, vR2, vA, 0)));
	}
	else
	{
		dummy6 = pow((y(vR1, vR2, vA, z) / StumpfC(z)), 1.5) * (1 / 2 / (z) * (StumpfC(z) - 3 * StumpfS(z) / 2 / StumpfC(z)) + 3 * pow(StumpfS(z), 2) / 4 / StumpfC(z)) +
			vAdiv8 * (3 * StumpfS(z) / StumpfC(z) * sqrt(y(vR1, vR2, vA, z)) + vA * sqrt(StumpfC(z) / y(vR1, vR2, vA, z)));
	}
	return dummy6;
}

double StumpfS(double z)
{
	double dummy4;
	double vZlowDenom;
	double vZhighDenom;

	if (z > 0) {

		vZlowDenom = 1 / pow(sqrt(z), 3);

		dummy4 = vZlowDenom * sqrt(z) - vZlowDenom * sin(sqrt(z));

	}
	else if (z < 0) {

		vZhighDenom = 1 / pow(sqrt(-z), 3);

		dummy4 = vZhighDenom * sinh(sqrt(-z)) - vZhighDenom * sqrt(-z);

	}
	else {

		dummy4 = 0.166666670;
	}
	return dummy4;
}

double StumpfC(double z)
{
	double dummy5;
	double vZdenom;

	if (z > 0) {

		vZdenom = 1 / (z);
		dummy5 = vZdenom - cos(sqrt(z)) * vZdenom;
	}
	else if (z < 0)
	{
		vZdenom = 1 / (z);
		dummy5 = -1 * cosh(sqrt(-z)) * vZdenom + vZdenom;
	}

	else {
		dummy5 = 0.50;
	}

	return dummy5;
}










