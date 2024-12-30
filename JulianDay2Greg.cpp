#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream> 

using namespace std;

double floor0(double variable1);

string JulianDay2Greg(double vJulianDay)	//(double M, double e)- what this function requires for input.. 
{
    double z;
    double w_top;

    double w;
    double x;
    double a;
    double b;
    double c;
    double c_top;
    double d;
    double e_top;
    double e;
    double f;

    int vGregYear;
    double vGregDay;

    string Greg_day;
    string sGregMonth;
    string Greg_year;
    //string sGregDate;

    int month;
    int month1;
    int month2;


    z = vJulianDay; //+ 0.50 ;

    w_top = z - 1867216.25;

    w = floor0(w_top / 36524.25); // 'round') ;

    x = floor0(w / 4);

    a = z + 1 + w - x;

    b = a + 1524;

    c_top = b - 122.1;

    c = floor0(c_top / 365.25);

    d = floor0(365.25 * c);

    e_top = b - d;

    e = floor0(e_top / 30.6001);

    f = floor0(30.6001 * e);

    vGregDay = b - d - f + 0.5;

    month1 = e - 1;
    month2 = e - 13;

    if (month1 < 1 || month1 > 12)
    {
        month = month2;
    }
    else
    {
        month = month1;
    }


    if (month == 1 || month == 2)
    {
        vGregYear = c - 4715;
    }
    else
    {
        vGregYear = c - 4716;
    }


    switch (month)
    {
    case 1:

        sGregMonth = "Jan";
        break;

    case 2:

        sGregMonth = "Feb";
        break;

    case 3:

        sGregMonth = "Mar";
        break;

    case 4:

        sGregMonth = "Apr";
        break;

    case 5:

        sGregMonth = "May";
        break;

    case 6:

        sGregMonth = "Jun";
        break;

    case 7:

        sGregMonth = "July";
        break;

    case 8:

        sGregMonth = "Aug";
        break;

    case 9:

        sGregMonth = "Sep";
        break;

    case 10:

        sGregMonth = "Oct";
        break;

    case 11:

        sGregMonth = "Nov";
        break;

    case 12:

        sGregMonth = "Dec";
        break;
    }

    std::ostringstream strs1;
    strs1 << vGregDay;
    std::string sGregDay = strs1.str();

    std::ostringstream strs2;
    strs2 << vGregYear;
    std::string sGregYear = strs2.str();

    std::stringstream ss;
    ss << sGregMonth << " " << sGregDay << " " << sGregYear;
    std::string sGregDate = ss.str();

    return sGregDate;
}

double floor0(double value)
{
    if (value < 0.0) {
        return ceil(value);
    }
    else
    {
        return floor(value);
    }
}