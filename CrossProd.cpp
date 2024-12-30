#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

struct CrossProdCompute {
    double aArray[3];
} stCrossProd;

double aCrossVector[3];

struct CrossProdCompute cross_prod(double aVector1[3], double aVector2[3])
{
    aCrossVector[0] = aVector1[1] * aVector2[2] - aVector1[2] * aVector2[1];

    aCrossVector[1] = aVector1[2] * aVector2[0] - aVector1[0] * aVector2[2];

    aCrossVector[2] = aVector1[0] * aVector2[1] - aVector1[1] * aVector2[0];


    stCrossProd.aArray[0] = aCrossVector[0];
    stCrossProd.aArray[1] = aCrossVector[1];
    stCrossProd.aArray[2] = aCrossVector[2];

    return stCrossProd;
}