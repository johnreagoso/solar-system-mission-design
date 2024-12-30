//#include <iostream>
//#include <cmath>
//#include <vector>
//#include <fstream>
//#include <iomanip>
//
//using namespace std;
//
//struct StateCompute{
//	double hVector[3];
//	double h_out;
//	double i;
//	double e_out;
//	double h;
//	double RAAN;
//	double ArgPer;
//	double TrueAnom;
//	double TrueLong;
//
//} stStateElementReturn;
//
//struct PlanetPosVel{
//	double vJDdate;
//	double aPosArray[3];
//	double aVelArray[3];
//
//}aPos1state;
//
//struct LambertInput{
//	
//	double aVelArray1[3];
//	double aVelArray2[3];
//	int vVariableNotUsed;
//
//}aVel1state;
//
//double vec_magntd(double aPosVecInput1[3]);
//double cross_prod(double aPosVecInput2[3], double aVelVecInput1[3]);
//double dot_prod(double aPosVecInput3[3], double aVelVecInput2[3]);
//
////double StateCompute SSMD_state_elements(struct PlanetPosVel *aPos1state, struct LambertInput *aVel1state)	 
//double state_elements(struct PlanetPosVel *aPos1state, struct LambertInput *aVel1state)	 
//{
//double const pi = 3.1415926535897932 ;	
//double vMu = 1.3271e11;
//double aPosVector[3];
//double aVelVector[3];	
//
//	aPosVector[0] = aPos1state->aPosArray[0];
//	aPosVector[1] = aPos1state->aPosArray[1];
//	aPosVector[2] = aPos1state->aPosArray[2];
//		 
//	aVelVector[0] = aVel1state->aVelArray1[0];
//	aVelVector[1] = aVel1state->aVelArray1[1];
//	aVelVector[2] = aVel1state->aVelArray1[2];
//
//
//double vRangeMagntd;	
//double vVelMagntd;
//double aSpecAngMomVector[3];	
//double vSpecAngMomMagtnd;
//
//double aNodeVector[3];
//double vIncl;
//double aEccenVector[3];
//double e;
//double eVector[3];
//double aZvector[3] = {0,0,1};
//
//double vArgPer;		double vTrueAnom;	double vRAAN;	double vTrueLong;
//
//vRangeMagntd = vec_magntd(aPosVector);
//vVelMagntd   = vec_magntd(aVelVector);
//
//aSpecAngMomVector[3]	= cross_prod(aPosVector, aVelVector);
//vSpecAngMomMagtnd		= vec_magntd(aSpecAngMomVector);
//
//aNodeVector[3] = cross_prod(aZvector, aSpecAngMomVector);
//
//eVector[0] = (1/vMu)* ((pow(vVelMagntd,2) - vMu/vRangeMagntd)* aPosVector[0] - dot_prod(aPosVector, aVelVector)* aVelVector[0]);
//eVector[1] = (1/vMu)* ((pow(vVelMagntd,2) - vMu/vRangeMagntd)* aPosVector[1] - dot_prod(aPosVector, aVelVector)* aVelVector[1]);
//eVector[2] = (1/vMu)* ((pow(vVelMagntd,2) - vMu/vRangeMagntd)* aPosVector[2] - dot_prod(aPosVector, aVelVector)* aVelVector[2]);
//
//e = vec_magntd(eVector);
//
//if (aNodeVector[1] < 0.0)
//{
//	vRAAN = 2*pi - acos(aNodeVector[0]/vec_magntd(aNodeVector));
//}
//else 
//{
//	vRAAN = acos(aNodeVector[0]/vec_magntd(aNodeVector));
//}
//
//
//if (eVector[2] < 0.00)
//{
//	vArgPer = vRAAN = 2*pi - acos(dot_prod(aNodeVector, eVector)/ (vec_magntd(aNodeVector)*vec_magntd(eVector)));
//}
//else 
//{
//	vArgPer = vRAAN = acos(dot_prod(aNodeVector, eVector)/ (vec_magntd(aNodeVector)*vec_magntd(eVector)));
//}
//
//
//if (dot_prod(aPosVector, aVelVector) <0.00)
//{
//	vTrueAnom = 2*pi - acos(dot_prod(eVector, aPosVector)/e*vRangeMagntd);
//}
//else
//{
//	vTrueAnom = acos(dot_prod(eVector, aPosVector)/e*vRangeMagntd);
//}
//
//vTrueLong = vRAAN + vArgPer + vTrueAnom;
//
//return vTrueLong;
////
////stStateElementReturn.hVector[0] = aSpecAngMomVector[0];
////stStateElementReturn.hVector[1] = aSpecAngMomVector[1];
////stStateElementReturn.hVector[2] = aSpecAngMomVector[2];
////
////stStateElementReturn.h_out = vSpecAngMomMagtnd;
////stStateElementReturn.e_out = e;
////stStateElementReturn.i = vIncl;
////stStateElementReturn.RAAN = vRAAN;
////stStateElementReturn.ArgPer = vArgPer;
////stStateElementReturn.TrueAnom = vTrueAnom;
////stStateElementReturn.TrueLong = vTrueLong;
//
////return stStateElementReturn;
//
//}