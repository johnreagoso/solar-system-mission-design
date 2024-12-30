#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
//#include "SpiceUsr.h"

using namespace std;

struct PlanetPosVelForJD {
	double vJulianDay;
	double aPosVector[3];
	double aVelVector[3];
	struct stPlPostVelforJD;

} stEarthPosVelForJD[9000], stVenusPosVelForJD[9000], stMarsPosVelForJD[9000];

//struct PorkChop
//{
//	double vJDdepart;
//	double vJDarrive;
//	double vDeltaV;
//	double vC3;
//
//} stPorkChopData1[1000000]; 

struct LambertData
{
	double aVelVector1[3];
	double aVelVector2[3];
	int vIterationsReq;

} LambertReturn1, LambertReturn2, LambertReturn3;

struct CrossProdCompute
{
	double aArray[3];

} stSpecAngMomVector, stNodeVector, CrossProdOut;

struct OrbitElementsCompute
{
	double hVector[3];
	double h;
	double e;
	double i;
	double RAAN;
	double ArgPer;
	double TrueAnom;
	double TrueLong;

} stTrajElementRtn1, stTrajElementRtn2, stTrajElementRtn3, stTrajPrototype;

struct Flyby1Data {

	double vTotalTOF;
	double vTOF1;
	double vTOF2;

	double vC3;
	double deltaV;

	double vThetaSwept;
	double vDepartJD;
	double vArrive1JD;
	double vArrive2JD;

	double vPlanetFlyby;
	double vFlybyAlt;
	double vVinfDeltaMeasured;
	double vTurnAngleRad;

}stFlyby1Data[1000000]; //, stFlyby2Data[1000000], stFlyby3Data[1000000];

struct PlanetParameters
{
	double a_mean_rtn;
	double planet_mu_rtn;
	double radius_planet_rtn;
	double a0rtn, a1rtn, e0rtn, e1rtn, i0rtn, i1rtn, L0rtn, L1rtn, w_bar0rtn, w_bar1rtn, omega0rtn, omega1rtn;
	double vMinAtmAltitude; //km

} stPlanetSelect1Param, stPlanetSelect2Param, stPlanetSelect3Param;

struct PlanetPosVel {
	double vMJDlook;
	double aPosVector[3];
	double aVelVector[3];

}stPlanetPosVel1, stPlanetPosVel2, stPlanetPosVel3, stPlanetPosVel4;

struct PorkChop {

	double vMJDdepart;
	double vMJDarrive;
	double vTOF;
	double vDeltaV;
	double vC3;
	double velInfVector_vx;
	double velInfVector_vy;
	double velInfVector_vz;
	double velInfVector_vx2;
	double velInfVector_vy2;
	double velInfVector_vz2;

}stPorkChopData1[1000000];  

double const pi = 3.1415926535897932;
double const vMuEarth = 3.9865e5;
double const vMuSun = 1.32714e11;
double earth_mu = 3.9865e5;
double earth_radius = 6378;
double radius_planet = 6378;
double const vMaxIterations = 10000;

// ==================================================================================================================================================
/* Function Prototype(s) declared: ================================================================================================================*/
struct OrbitElementsCompute state_elements(struct PlanetPosVel* worker3, struct LambertData* worker4, double variable5);

struct PlanetParameters PlanetSelectOrbitParameters(int vPlanetSelect);

struct LambertData SSMD_LambertsSolverUnivVariables(struct PlanetPosVel* worker1, struct PlanetPosVel* worker2, double* vVariable2, double* vVariable3);

struct OrbitElementsCompute state_elements(struct PlanetPosVel* worker3, double worker1, double worker2, double worker6, double variable5);

double SSMD_JulianDateCompute(int vYearInput, int vMonthInput, int vDayInput);

// alternate method:
struct PlanetPosVelForJD SSMD_GetPlanetPosVel(int vPlanetSelect, double vVariable2, struct PlanetParameters* struct1);

double dot_prod(double vector1[3], double vector2[3]);

string JulianDay2Greg(double variable1);

double vec_magntd(double vector1[3]);

struct CrossProdCompute cross_prod(double vector1[3], double vector2[3]);

double AngleSwept(struct PlanetPosVel* strct1, struct PlanetPosVel* strct2);
// ==================================================================================================================================================
// ==================================================================================================================================================

int const vMax_iterations = 1000;
double const vTflight = 100;
double const vTol = 1.00e-8;
double const vNmax = 1000;

/*
**************************************************************************************
**************************************************************************************
					************************************
					*								   *
					*	 SSMD MISSION USER INPUTS	   *
					*								   *
					************************************
*/
int vPlanetSelect1 = 4;
int vPlanetSelect2 = 3; 
int vPlanetSelect3 = 5;  // Assumes only (1) gravity-assist for this code version... 

double vDeltaVmax = 10;
//******************************************************************
int vFlybyNo = 0;
//******************************************************************

int vYearInput = 2028;			int vMonthInput = 1;		int vDayInput = 1;
int vYearInputLast = 2032;		int vMonthInputLast = 12;	int vDayInputLast = 31;
//int vYearInputLast = 2014;		int vMonthInputLast = 1;	int vDayInputLast = 31;

double vArc1Min = 90;				double vArc1Max = 500.0; //double vArc1Max = 365*1;

double vArc2Min = 365*2;			double vArc2Max = 365*10;
double vArc3Min = 365*2;			double vArc3Max = 365*10;


/***************************************************************************************
**************************************************************************************/

double vTotTrajLook = (vArc1Max - vArc1Min) * (vArc2Max - vArc2Min);
double vDeltaV;
double vDayCountIndex;
double vDayCountBase;

string sGregDateOutputDep;
string sGregDateOutputArr1;
string sGregDateOutputArr2;

int main()
{
	clock_t startTime = clock();

	double aPosVector1[3], aPosVector2[3], aPosVector3[3];;
	double aVelVector1[3], aVelVector2[3], aVelVector3[3];
	double vTemp, vReturn, vJulianDay;

	double vDepartMJDstart, vDepartMJDlast, vArrival1_MJDlast, vArrival2_MJDlast;
	double iiDepartMJD, vArrive1MJD, vArrive2MJD, vArrive3MJD;

	int ii, jj, kk, mm, vTextNum, vDepartMJDspan;

	int	vFlyby3Count = 0;

	int jjCount = 0;	int kkCount = 0;
	int	vDepartCount = 0;	int vArc1no = 0;

	int vArc2no;			int vFlyby1Count = 0;	int vFlyby2Count = 0;
	int vArc1noOut = 0;		int vArc2noOut = 0;

	double vTOF1, vTOF2, vTOF3;
	double vTOF1secs, vTOF2secs, vTOF3secs, vAngleSweptRad1, vAngleSweptRad2, vAngleSweptRad3;

	double vSMA1, vEnergy1, vSMA2, vSMA3, vEnergy2, vEnergy3;

	double vSCarc1DepInvVecMagntd;

	double vC3Mgntd;					double vVscP;
	double vSCVelEarthOrbit;   			double vDeltaV;

	double aSCarc1DepInfVelVector[3];	double aSCarc2ArrInfVelVector[3];

	double aPreFlyby1VinfVel[3];		double vPreFlyby1VinfVelMag;
	double aPostFlyby1VinfVel[3];		double vPostFlyby1VinfVelMag;

	double aPreFlyby2VinfVel[3];		double vPreFlyby2VinfVelMag;
	double aPostFlyby2VinfVel[3];		double vPostFlyby2VinfVelMag;

	double vDeltaVinf1;					double vDeltaVinf2;
	double vDeltaVinfAvg;				double vDeltaVinfFlyby1, vDeltaVinfFlyby2;

	double vSCTurnAngleRad;		double vEccenTurnRad;	double vSCTurnCompute2;
	double vFlybyVinfAvgMag;	double vMuPlanetDenom; 		double vMaxSCTurn;
	double vRpFlybyDistActual;	double vSCTurnCompute;		double vTotalFlight1;

	double vTotalFlight2;
	double vCross;

	vJulianDay = SSMD_JulianDateCompute(vYearInput, vMonthInput, vDayInput);

	vDepartMJDstart = SSMD_JulianDateCompute(vYearInput, vMonthInput, vDayInput);
	vDepartMJDlast = SSMD_JulianDateCompute(vYearInputLast, vMonthInputLast, vDayInputLast);

	vDepartMJDspan = (vDepartMJDlast - vDepartMJDstart) + 1;

	double vMissionSpanDays = (vDepartMJDlast + vArc1Max + vArc2Max) - vDepartMJDstart;
	double vJDindexStart;
	double vJDindexEnd;
	int vJDindex;
	int	vDepartCounter = 0;
	double vDepartMJD;

	string output;
	string output1;
	string output2;
	string output3;

	stPlanetSelect1Param = PlanetSelectOrbitParameters(vPlanetSelect1);
	stPlanetSelect2Param = PlanetSelectOrbitParameters(vPlanetSelect2);
	stPlanetSelect3Param = PlanetSelectOrbitParameters(vPlanetSelect3);

	// here we populate the 10 years worth of planetary position and velocity vectors at once: 
	for (vJDindex = 0; vJDindex < vMissionSpanDays; vJDindex++)
	{
		vDepartMJD = vJDindex + vDepartMJDstart;
		stEarthPosVelForJD[vJDindex] = SSMD_GetPlanetPosVel(3, vDepartMJD, &stPlanetSelect1Param);
		stMarsPosVelForJD[vJDindex] = SSMD_GetPlanetPosVel(4, vDepartMJD, &stPlanetSelect2Param);
		stVenusPosVelForJD[vJDindex] = SSMD_GetPlanetPosVel(2, vDepartMJD, &stPlanetSelect2Param);
	}

	stPlanetSelect1Param = PlanetSelectOrbitParameters(vPlanetSelect1);
	stPlanetSelect2Param = PlanetSelectOrbitParameters(vPlanetSelect2);
	stPlanetSelect3Param = PlanetSelectOrbitParameters(vPlanetSelect3);

	double vDepartMJDarray[5000]; // vDepartMJD; // this array-size must be modified as a user-input as well. 

	for (ii = 0; ii < vDepartMJDspan; ii++)  // this for-loop populates the MJD array of departured dates, for individual indexed MJDs to be pulled during ops. 
	{
		vDepartMJDarray[ii] = vDepartMJDstart + ii;
	}

	//======================================================================================================================================================================
	// Main Body ===========================================================================================================================================================
	for (vDepartCount = 0; vDepartCount < vDepartMJDspan; vDepartCount = vDepartCount + 2)
	{
		vDepartCounter = vDepartCounter + 1;
		vDepartMJD = vDepartMJDarray[vDepartCount];
		string output;
		output = JulianDay2Greg(vJulianDay);

		stPlanetPosVel1.aPosVector[0] = stEarthPosVelForJD[vDepartCount].aPosVector[0];
		stPlanetPosVel1.aPosVector[1] = stEarthPosVelForJD[vDepartCount].aPosVector[1];
		stPlanetPosVel1.aPosVector[2] = stEarthPosVelForJD[vDepartCount].aPosVector[2];

		stPlanetPosVel1.aVelVector[0] = stEarthPosVelForJD[vDepartCount].aVelVector[0];
		stPlanetPosVel1.aVelVector[1] = stEarthPosVelForJD[vDepartCount].aVelVector[1];
		stPlanetPosVel1.aVelVector[2] = stEarthPosVelForJD[vDepartCount].aVelVector[2];

		stPlanetPosVel1.vMJDlook = vDepartMJD;

		//stPlanetPosVel1 = SSMD_GetPlanetPosVel(vPlanetSelect1, vDepartMJD, &stPlanetSelect1Param); // Here we will extract the Departure planet's heliocentric position and velocity..

		cout << "Total number of Departs are:" << endl;
		cout << vDepartCounter << endl;

		//for (jj = vArc1Min; jj <= vArc1Max; jj++)
		for (jj = vArc1Min; jj <= vArc1Max; jj = jj + 2)
		{
			vTOF1 = jj;

			vTOF1secs = vTOF1*24*60*60;

			vArrive1MJD = vDepartMJD + jj;

			//stPlanetPosVel2.aPosVector[0] = stVenusPosVelForJD[vDepartCount+jj].aPosVector[0];
			//stPlanetPosVel2.aPosVector[1] = stVenusPosVelForJD[vDepartCount+jj].aPosVector[1];
			//stPlanetPosVel2.aPosVector[2] = stVenusPosVelForJD[vDepartCount+jj].aPosVector[2];

			//stPlanetPosVel2.aVelVector[0] = stVenusPosVelForJD[vDepartCount+jj].aVelVector[0];
			//stPlanetPosVel2.aVelVector[1] = stVenusPosVelForJD[vDepartCount+jj].aVelVector[1];
			//stPlanetPosVel2.aVelVector[2] = stVenusPosVelForJD[vDepartCount+jj].aVelVector[2];


			stPlanetPosVel2.aPosVector[0] = stMarsPosVelForJD[vDepartCount + jj].aPosVector[0];
			stPlanetPosVel2.aPosVector[1] = stMarsPosVelForJD[vDepartCount + jj].aPosVector[1];
			stPlanetPosVel2.aPosVector[2] = stMarsPosVelForJD[vDepartCount + jj].aPosVector[2];

			stPlanetPosVel2.aVelVector[0] = stMarsPosVelForJD[vDepartCount + jj].aVelVector[0];
			stPlanetPosVel2.aVelVector[1] = stMarsPosVelForJD[vDepartCount + jj].aVelVector[1];
			stPlanetPosVel2.aVelVector[2] = stMarsPosVelForJD[vDepartCount + jj].aVelVector[2];

			stPlanetPosVel2.vMJDlook = vDepartMJD + vTOF1;

			sGregDateOutputDep = JulianDay2Greg(vDepartMJD);

			vAngleSweptRad1 = AngleSwept(&stPlanetPosVel1, &stPlanetPosVel2);

			LambertReturn1 = SSMD_LambertsSolverUnivVariables(&stPlanetPosVel1, &stPlanetPosVel2, &vAngleSweptRad1, &vTOF1secs);

			vSMA1 = pow(2 / vec_magntd(stPlanetPosVel1.aPosVector) / vMuSun - pow(vec_magntd(LambertReturn1.aVelVector1), 2) / vMuSun, -1);

			vEnergy1 = -1 *vMuSun / 2* vSMA1;

			stTrajElementRtn1 = state_elements(&stPlanetPosVel1, &LambertReturn1, vMuSun);

			aSCarc1DepInfVelVector[0] = LambertReturn1.aVelVector1[0] - stPlanetPosVel1.aVelVector[0];
			aSCarc1DepInfVelVector[1] = LambertReturn1.aVelVector1[1] - stPlanetPosVel1.aVelVector[1];
			aSCarc1DepInfVelVector[2] = LambertReturn1.aVelVector1[2] - stPlanetPosVel1.aVelVector[2];

			aSCarc2ArrInfVelVector[0] = LambertReturn1.aVelVector2[0] - stPlanetPosVel2.aVelVector[0];
			aSCarc2ArrInfVelVector[1] = LambertReturn1.aVelVector2[1] - stPlanetPosVel2.aVelVector[1];
			aSCarc2ArrInfVelVector[2] = LambertReturn1.aVelVector2[2] - stPlanetPosVel2.aVelVector[2];

			vSCarc1DepInvVecMagntd = vec_magntd(aSCarc1DepInfVelVector);

			vC3Mgntd = pow(vec_magntd(aSCarc1DepInfVelVector), 2);

			// vVscP = sqrt(vC3Mgntd + 2 * vMuEarth / 6578);
			// vSCVelEarthOrbit = sqrt(vMuEarth / 6578);

			vVscP = sqrt(vC3Mgntd + 2*42828.37 / 3739.2);
			vSCVelEarthOrbit = sqrt(42828.37 / 3739.2);

			vDeltaV = vVscP - vSCVelEarthOrbit;

//		//	if (vSCarc1DepInvVecMagntd <= 9.00)
//			{
//				stPorkChopData1[vArc1noOut].vMJDdepart	= vDepartMJD;
//				stPorkChopData1[vArc1noOut].vMJDarrive	= vArrive1MJD;
//				stPorkChopData1[vArc1noOut].vC3			= vC3Mgntd;
//				stPorkChopData1[vArc1noOut].vDeltaV		= vDeltaV;
//				stPorkChopData1[vArc1noOut].vTOF = vTOF1;
//				stPorkChopData1[vArc1noOut].velInfVector_vx = aSCarc1DepInfVelVector[0];
//				stPorkChopData1[vArc1noOut].velInfVector_vy = aSCarc1DepInfVelVector[1];
//				stPorkChopData1[vArc1noOut].velInfVector_vz = aSCarc1DepInfVelVector[2];
////
//		//		vArc1noOut = vArc1noOut + 1;
//		//	}

			stPorkChopData1[vArc1noOut].vMJDdepart = vDepartMJD;
			stPorkChopData1[vArc1noOut].vMJDarrive = vArrive1MJD;
			stPorkChopData1[vArc1noOut].vC3 = vC3Mgntd;
			stPorkChopData1[vArc1noOut].vDeltaV = vDeltaV;
			stPorkChopData1[vArc1noOut].vTOF = vTOF1;
			stPorkChopData1[vArc1noOut].velInfVector_vx = aSCarc1DepInfVelVector[0];
			stPorkChopData1[vArc1noOut].velInfVector_vy = aSCarc1DepInfVelVector[1];
			stPorkChopData1[vArc1noOut].velInfVector_vz = aSCarc1DepInfVelVector[2];
			
			stPorkChopData1[vArc1noOut].velInfVector_vx2 = aSCarc2ArrInfVelVector[0];
			stPorkChopData1[vArc1noOut].velInfVector_vy2 = aSCarc2ArrInfVelVector[1];
			stPorkChopData1[vArc1noOut].velInfVector_vz2 = aSCarc2ArrInfVelVector[2];

			vArc1noOut = vArc1noOut + 1;




			// compute the deltaV and C3 components.. then determine if we should look at a flyby opportunity.. 

			//if (vDeltaV <= vDeltaVmax && vFlybyNo > 0)

/************************************************************************************************************************************
 ***********************************************************************************************************************************
								**********************************
								*								 *
								*	      Flyby-1 Loop			 *
								*								 *
								**********************************/
			if (vSCarc1DepInvVecMagntd <= 12.0 && vFlybyNo > 0)
			//if (vDeltaV <= 6.00 && vFlybyNo > 0)
			{
				//for (kk = vArc2Min; kk <= vArc2Max; kk++)
				for (kk = vArc2Min; kk <= vArc2Max; kk = kk + 1)
				{
					vArrive2MJD = vArrive1MJD + kk;
					output3 = JulianDay2Greg(vArrive2MJD);

					vTOF2 = kk;
					vTOF2secs = vTOF2 * 24 * 60 * 60;

					/*stPlanetPosVel3.aPosVector[0] = stMarsPosVelForJD[vDepartCount+jj+kk].aPosVector[0];
					stPlanetPosVel3.aPosVector[1] = stMarsPosVelForJD[vDepartCount+jj+kk].aPosVector[1];
					stPlanetPosVel3.aPosVector[2] = stMarsPosVelForJD[vDepartCount+jj+kk].aPosVector[2];

					stPlanetPosVel3.aVelVector[0] = stMarsPosVelForJD[vDepartCount+jj+kk].aVelVector[0];
					stPlane0tPosVel3.aVelVector[1] = stMarsPosVelForJD[vDepartCount+jj+kk].aVelVector[1];
					stPlanetPosVel3.aVelVector[2] = stMarsPosVelForJD[vDepartCount+jj+kk].aVelVector[2];*/

					stPlanetPosVel3.aPosVector[0] = stEarthPosVelForJD[vDepartCount + jj + kk].aPosVector[0];
					stPlanetPosVel3.aPosVector[1] = stEarthPosVelForJD[vDepartCount + jj + kk].aPosVector[1];
					stPlanetPosVel3.aPosVector[2] = stEarthPosVelForJD[vDepartCount + jj + kk].aPosVector[2];

					stPlanetPosVel3.aVelVector[0] = stEarthPosVelForJD[vDepartCount + jj + kk].aVelVector[0];
					stPlanetPosVel3.aVelVector[1] = stEarthPosVelForJD[vDepartCount + jj + kk].aVelVector[1];
					stPlanetPosVel3.aVelVector[2] = stEarthPosVelForJD[vDepartCount + jj + kk].aVelVector[2];

					stPlanetPosVel3.vMJDlook = vDepartMJD + vTOF1 + vTOF2;

					//stPlanetPosVel3	= SSMD_GetPlanetPosVel(vPlanetSelect3, vArrive2MJD, &stPlanetSelect3Param); // Here we will extract the Arrival planet's heliocentric position and velocity based on the vDepart date plus jj:

					vAngleSweptRad2 = AngleSwept(&stPlanetPosVel2, &stPlanetPosVel3);

					LambertReturn2 = SSMD_LambertsSolverUnivVariables(&stPlanetPosVel2, &stPlanetPosVel3, &vAngleSweptRad2, &vTOF2secs);

					aPreFlyby1VinfVel[0] = LambertReturn1.aVelVector2[0] - stPlanetPosVel2.aVelVector[0];
					aPreFlyby1VinfVel[1] = LambertReturn1.aVelVector2[1] - stPlanetPosVel2.aVelVector[1];
					aPreFlyby1VinfVel[2] = LambertReturn1.aVelVector2[2] - stPlanetPosVel2.aVelVector[2];

					aPostFlyby1VinfVel[0] = LambertReturn2.aVelVector1[0] - stPlanetPosVel2.aVelVector[0];
					aPostFlyby1VinfVel[1] = LambertReturn2.aVelVector1[1] - stPlanetPosVel2.aVelVector[1];
					aPostFlyby1VinfVel[2] = LambertReturn2.aVelVector1[2] - stPlanetPosVel2.aVelVector[2];

					vPreFlyby1VinfVelMag = vec_magntd(aPreFlyby1VinfVel);
					vPostFlyby1VinfVelMag = vec_magntd(aPostFlyby1VinfVel);

					vDeltaVinf1 = fabs(1 - fabs(vPreFlyby1VinfVelMag / vPostFlyby1VinfVelMag));
					vDeltaVinf2 = fabs(1 - fabs(vPostFlyby1VinfVelMag / vPreFlyby1VinfVelMag));

					vDeltaVinfAvg = 0.5 * (vDeltaVinf1)+0.5 * (vDeltaVinf2);

					//			    vDeltaVinf		=  1 - abs(norm(aSCPostFlyby_InfVelVector)/norm(aSCPreFlyby_InfVelVector));
					//				vDeltaVinfAvg	= 0.5*(abs(vDeltaVinf) + abs(vDeltaVinf2));

					vDeltaVinfFlyby1 = abs(1 - abs(vPostFlyby1VinfVelMag / vPreFlyby1VinfVelMag));   // computes the percentage (in decimals) of the difference in magnitudes between the pre/post flyby inf-velocity vectors.

					vCross = aPreFlyby1VinfVel[0] * aPostFlyby1VinfVel[1] - aPreFlyby1VinfVel[1] * aPostFlyby1VinfVel[0];

					vSCTurnAngleRad = acos(dot_prod(aPreFlyby1VinfVel, aPostFlyby1VinfVel) / (vec_magntd(aPreFlyby1VinfVel) * vec_magntd(aPostFlyby1VinfVel)));

					vEccenTurnRad = 1 / sin(0.5 * vSCTurnAngleRad);

					vFlybyVinfAvgMag = 0.5 * vec_magntd(aPreFlyby1VinfVel) + 0.5 * vec_magntd(aPostFlyby1VinfVel);
					vMuPlanetDenom = 1 / stPlanetSelect2Param.planet_mu_rtn;					//1/vMuPlanet;

					vMaxSCTurn = 2 * asin(1 / (1 + stPlanetSelect2Param.vMinAtmAltitude * pow(vFlybyVinfAvgMag, 2) * vMuPlanetDenom));

					vRpFlybyDistActual = (vEccenTurnRad - 1) * stPlanetSelect2Param.planet_mu_rtn / pow(vFlybyVinfAvgMag, 2);

					vSCTurnCompute = 2 * asin(1 / (1 + vRpFlybyDistActual * pow(vFlybyVinfAvgMag, 2) * vMuPlanetDenom));

					vSCTurnCompute2 = pi - 2 * acos(1 / (1 + vRpFlybyDistActual * pow(vFlybyVinfAvgMag, 2) * vMuPlanetDenom));

					//vSCTurnCompute		= 2*asin(1/(1 + vRpFlybyDistActual* pow(vFlybyVinfAvgMag,2)*vMuPlanetDenom));

					vTotalFlight1 = vTOF1 + vTOF2;

					if (vDeltaVinfFlyby1 <= 0.025 && vRpFlybyDistActual >= stPlanetSelect2Param.vMinAtmAltitude)
					{
						stFlyby1Data[vFlyby1Count].vDepartJD = vDepartMJD;
						stFlyby1Data[vFlyby1Count].vArrive1JD = vArrive1MJD;
						stFlyby1Data[vFlyby1Count].vArrive2JD = vArrive2MJD;
						stFlyby1Data[vFlyby1Count].vC3 = vC3Mgntd;
						
						stFlyby1Data[vFlyby1Count].deltaV = vDeltaV;
						stFlyby1Data[vFlyby1Count].vTOF1  = vTOF1;
						stFlyby1Data[vFlyby1Count].vTOF2  = vTOF2;

						stFlyby1Data[vFlyby1Count].vTotalTOF = vTotalFlight1;
						stFlyby1Data[vFlyby1Count].vFlybyAlt = vRpFlybyDistActual;
						stFlyby1Data[vFlyby1Count].vTurnAngleRad = vSCTurnCompute;
						stFlyby1Data[vFlyby1Count].vVinfDeltaMeasured = vDeltaVinfFlyby1;

						vFlyby1Count = vFlyby1Count + 1;
						//break;
					}
				}

				vArc1no = vArc1no + 1;
			}
		}
		cout << sGregDateOutputDep.c_str() << endl;
	}
	// End Main Body ===========================================================================================================================================================

	clock_t endTime = clock();

	cout << "Total number of executed Lamberts runs are:" << endl;
	cout << vDepartCounter * (vArc1Max - vArc1Min) << endl;
	cout << "Total script execution speed for this run is:" << endl;
	cout << double(clock() - startTime) / 1000 << endl;

	ofstream PorkChopFile;
	PorkChopFile.open ("PorkChop_Mars_Earth_Jan2028_2_Dec2032.txt");
	PorkChopFile << "Depart Date (Greg)" << setw(30)  <<   "Arrive-1 Date (Greg)" <<  setw(30)  <<   "Depart Date (JD)" <<  setw(30)  <<" Arrive-1 Date (JD)" << setw(30) << " C3 Magntd (km^2/sec^2)" << setw(20)  << "DeltaV"  << setw(30) << "TOF"  << setw(30) << "VinfVel_vx" << setw(30) << "VinfVel_vy" << setw(30) << "VinfVel_vz" << "Arrival VinfVel_vx" << setw(30) << "Arrival VinfVel_vy" << setw(30) << "Arrival VinfVel_vz" << endl;
	PorkChopFile << setprecision(12);

		for (kk = 0; kk< vArc1noOut ; kk++)
		{
		sGregDateOutputDep  = JulianDay2Greg(stPorkChopData1[kk].vMJDdepart);
		sGregDateOutputArr1 = JulianDay2Greg(stPorkChopData1[kk].vMJDarrive);
		
		PorkChopFile <<  sGregDateOutputDep.c_str() <<  setw(30)      << sGregDateOutputArr1.c_str() <<  setw(30) << stPorkChopData1[kk].vMJDdepart  << setw(30)  << stPorkChopData1[kk].vMJDarrive << setw(30) << stPorkChopData1[kk].vC3 <<  setw(30)  << stPorkChopData1[kk].vDeltaV  <<setw(30) << stPorkChopData1[kk].vTOF << setw(30) << stPorkChopData1[kk].velInfVector_vx << setw(30) << stPorkChopData1[kk].velInfVector_vy << setw(30) << stPorkChopData1[kk].velInfVector_vz << setw(30) <<
			stPorkChopData1[kk].velInfVector_vx2 << setw(30) << stPorkChopData1[kk].velInfVector_vy2 << setw(30) << stPorkChopData1[kk].velInfVector_vz2 << setw(30) << kk <<endl;
		}

	PorkChopFile.close();
 	return 0;

//******************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************
//                 Flyby Output .CSV Export
//******************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************
	ofstream myfile;
	myfile.open("mars_earth_2030_2032.txt");
	myfile << "Depart Date (Greg)" << setw(30) << "Arrive-1 Date (Greg)" << setw(30) << "Arrive-1 Date (Greg)" << setw(30) << "Depart Date (JD)" << setw(30) << " Arrive-1 Date (JD)" << setw(30) << " Arrive-2 Date (JD)" << setw(30) << " C3 Magntd (km^2/sec^2)" << setw(30) << "  delta-v (km/sec) "  << setw(30) << " TOF1 (days) " << setw(30) << "TOF2 (days) " << setw(30) << "Total TOF" << setw(30) << " Turn Angle (rad)" << setw(30) << "Vinf Vel Vector Delta Measured" << setw(30) << "Flyby Altitude" << endl;
	myfile << setprecision(12);

	for (kk = 0; kk < vFlyby1Count; kk++)
	{
		sGregDateOutputDep = JulianDay2Greg(stFlyby1Data[kk].vDepartJD);
		sGregDateOutputArr1 = JulianDay2Greg(stFlyby1Data[kk].vArrive1JD);
		sGregDateOutputArr2 = JulianDay2Greg(stFlyby1Data[kk].vArrive2JD);

		myfile << sGregDateOutputDep.c_str() << setw(30) << sGregDateOutputArr1.c_str() << setw(30) << sGregDateOutputArr2.c_str() << setw(30) << stFlyby1Data[kk].vDepartJD << setw(30) << stFlyby1Data[kk].vArrive1JD << setw(30) << stFlyby1Data[kk].vArrive2JD << setw(30) << stFlyby1Data[kk].vC3 << setw(30) << stFlyby1Data[kk].deltaV << setw(30) << stFlyby1Data[kk].vTOF1 << setw(30) << stFlyby1Data[kk].vTOF2 << setw(30) << stFlyby1Data[kk].vTotalTOF << setw(30) << stFlyby1Data[kk].vTurnAngleRad << setw(30) << stFlyby1Data[kk].vVinfDeltaMeasured << setw(30) << stFlyby1Data[kk].vFlybyAlt << kk << endl;
	}

	myfile.close();

	//******************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************
	//******************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************	
	return 0;
}

struct OrbitElementsCompute state_elements(struct PlanetPosVel* aPos1state, struct LambertData* aVel1state, double vMu)
{
	double const pi = 3.1415926535897932;
	double aPosVector[3];
	double aVelVector[3];
	double vRangeMagntd;
	double vVelMagntd;
	double aSpecAngMomVector[3];
	double vSpecAngMomMagtnd;
	double vIncl;
	double aEccenVector[3];
	double e;
	double aNodeVector[3];
	double eVector[3];
	double vArgPer;		double vTrueAnom;	double vRAAN;	double vTrueLong;
	double aZvector[3] = { 0,0,1 };

	aPosVector[0] = aPos1state->aPosVector[0];	aPosVector[1] = aPos1state->aPosVector[1];	aPosVector[2] = aPos1state->aPosVector[2];
	aVelVector[0] = aVel1state->aVelVector1[0];	aVelVector[1] = aVel1state->aVelVector1[1];	aVelVector[2] = aVel1state->aVelVector1[2];

	vRangeMagntd = vec_magntd(aPosVector);
	vVelMagntd = vec_magntd(aVelVector);

	stSpecAngMomVector = cross_prod(aPosVector, aVelVector);
	vSpecAngMomMagtnd = vec_magntd(stSpecAngMomVector.aArray);

	stNodeVector = cross_prod(aZvector, stSpecAngMomVector.aArray);

	eVector[0] = (1 / vMu) * ((pow(vVelMagntd, 2) - vMu / vRangeMagntd) * aPosVector[0] - dot_prod(aPosVector, aVelVector) * aVelVector[0]);
	eVector[1] = (1 / vMu) * ((pow(vVelMagntd, 2) - vMu / vRangeMagntd) * aPosVector[1] - dot_prod(aPosVector, aVelVector) * aVelVector[1]);
	eVector[2] = (1 / vMu) * ((pow(vVelMagntd, 2) - vMu / vRangeMagntd) * aPosVector[2] - dot_prod(aPosVector, aVelVector) * aVelVector[2]);

	e = vec_magntd(eVector);

	if (stNodeVector.aArray[1] < 0.0) {
		vRAAN = 2 * pi - acos(stNodeVector.aArray[0] / vec_magntd(stNodeVector.aArray));
	}
	else {
		vRAAN = acos(stNodeVector.aArray[0] / vec_magntd(stNodeVector.aArray));
	}

	if (eVector[2] < 0.00) {
		vArgPer = vRAAN = 2 * pi - acos(dot_prod(stNodeVector.aArray, eVector) / (vec_magntd(stNodeVector.aArray) * vec_magntd(eVector)));
	}
	else {
		vArgPer = vRAAN = acos(dot_prod(stNodeVector.aArray, eVector) / (vec_magntd(stNodeVector.aArray) * vec_magntd(eVector)));
	}

	if (dot_prod(aPosVector, aVelVector) < 0.00) {
		vTrueAnom = 2 * pi - acos(dot_prod(eVector, aPosVector) / e * vRangeMagntd);
	}
	else {
		vTrueAnom = acos(dot_prod(eVector, aPosVector) / (e * vRangeMagntd));
	}

	vIncl = 10;

	vTrueLong = vRAAN + vArgPer + vTrueAnom;

	stTrajPrototype.hVector[0] = stSpecAngMomVector.aArray[0];
	stTrajPrototype.hVector[1] = stSpecAngMomVector.aArray[1];
	stTrajPrototype.hVector[2] = stSpecAngMomVector.aArray[2];

	stTrajPrototype.h = vSpecAngMomMagtnd;
	stTrajPrototype.e = e;
	stTrajPrototype.i = vIncl;
	stTrajPrototype.RAAN = vRAAN;
	stTrajPrototype.ArgPer = vArgPer;
	stTrajPrototype.TrueAnom = vTrueAnom;
	stTrajPrototype.TrueLong = vTrueLong;

	return stTrajPrototype;
}

double AngleSwept(struct PlanetPosVel* strct1, struct PlanetPosVel* strct2)
{
	double vThetaSwept;
	double vCross;
	double aPosVector1[3];
	double aPosVector2[3];
	double vPosVectorMagntd1;
	double vPosVectorMagntd2;

	aPosVector1[0] = strct1->aPosVector[0];
	aPosVector1[1] = strct1->aPosVector[1];
	aPosVector1[2] = strct1->aPosVector[2];

	aPosVector2[0] = strct2->aPosVector[0];
	aPosVector2[1] = strct2->aPosVector[1];
	aPosVector2[2] = strct2->aPosVector[2];

	vPosVectorMagntd1 = vec_magntd(aPosVector1);
	vPosVectorMagntd2 = vec_magntd(aPosVector2);

	// Here we compute the angle swept:
	vCross = aPosVector1[0] * aPosVector2[1] - aPosVector1[1] * aPosVector2[0];


	// Spacecraft angle swept (theta_swept) computation:
	if (vCross >= 0) {
		vThetaSwept = acos((dot_prod(aPosVector1, aPosVector2) / (vPosVectorMagntd1 * vPosVectorMagntd2)));
	}
	else {
		vThetaSwept = 2 * pi - acos((dot_prod(aPosVector1, aPosVector2) / (vPosVectorMagntd1 * vPosVectorMagntd2)));
	}

	return vThetaSwept;
}

