#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

/* This c++ code will compute the Julian date for a specific Gregorian calendar date. This script will be
specifically used departure planet ephemeris data, and then used for the arrival planet(s) ephemeris
arrival date, and for specific iterations for the transfer orbit elements.							*/

double SSMD_JulianDateCompute(int vYearInput, int vMonthInput, int vDayInput)
{

    double vYear;
    double vMonth;
    double vA;
    double vB;
    double vJulianDay;

    if (vMonthInput > 2) {
        vYear = vYearInput;
    }
    else {
        vYear = vYearInput - 1;
    }

    if (vMonthInput > 2) {
        vMonth = vMonthInput;
    }
    else {
        vMonth = vMonthInput + 12;
    }

    vA = floor(vYear / 100.0);	vB = 2 + floor(vA / 4) - vA;

    vJulianDay = floor(365.25 * (vYear + 4716)) + floor(30.6001 * (vMonth + 1)) + vDayInput + vB - 1524.5;

    return vJulianDay;
}

