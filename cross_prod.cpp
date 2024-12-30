#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

double dot_prod(double vector1[3], double vector2[3])
{
    double dummy;
    dummy = vector1[1] * vector2[1] + vector1[2] * vector2[2] + vector1[3] * vector2[3];

    return dummy;
}

function[aVectorCrossOut] = Fcross(aVector1, aVector2)

aVectorCrossOut(1, 1) = aVector1(2) * aVector2(3) - aVector1(3) * aVector2(2);

aVectorCrossOut(2, 1) = aVector1(3) * aVector2(1) - aVector1(1) * aVector2(3);

aVectorCrossOut(3, 1) = aVector1(1) * aVector2(2) - aVector1(2) * aVector2(1);

return