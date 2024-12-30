#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

double Kepler(double* M_input, double* e_input)	//(double M, double e)- what this function requires for input.. 
{
	double delta_M, delta_E, E;
	double const pi = 3.1415926535897932;

	double e = *e_input;	double M = *M_input;

	if (e < 1)
	{
		E = M;
	}
	else
	{
		E = M + e * M;
	}

	do
	{
		delta_M = M - double(E - e * sin(E));  // "double" is called casting...

		delta_E = delta_M / double(1 - e * cos(E));

		E = E + delta_E;

	} while (abs(delta_E) > 0.0000001);


	if (E < 0.0)
	{
		E = E + 2 * pi;
	}
	else
	{
		E = E;
	}

	return E;
}
