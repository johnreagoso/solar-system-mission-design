#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

double vec_magntd(double vector[3])
{
	double dummy1;
	dummy1 = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);

	return dummy1;
}