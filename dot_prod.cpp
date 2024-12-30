#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

double dot_prod(double vector1[3], double vector2[3])
{
	double dummy;
	dummy = vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];

	return dummy;
}