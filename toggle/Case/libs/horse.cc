//This library have the general purpose functions used in the programs getci.cc and feconf.cc
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>
#include <iostream>
#include "horse.h"

using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::string;

Horse::Horse()
{
}

Horse::Horse(int semilla)
{
    start_rng(semilla);
}

void Horse::start_rng(int semilla)
{
    r = gsl_rng_alloc(gsl_rng_default); // Media en 0
    gsl_rng_set(r, semilla);
}

double Horse::randreal()
{
  //returns random double in the [0, 1) interval
  double res = gsl_rng_uniform (r);
  return res;
}

void Horse::close_rng()
{
  gsl_rng_free (r);
}

void Horse::open_ifstream(ifstream& fe, string nomb)
{
	fe.open(nomb.c_str());
	if (!fe.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Horse::open_ifstream.\n";
		exit(1);
	}
}

void Horse::open_ofstream(ofstream& fs, string nomb)
{
	fs.open(nomb.c_str(), ios::app);
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Horse::open_ofstream.\n";
		exit(1);
	}
}

void Horse::open_ofstream_new(ofstream& fs, string nomb)
{
	fs.open(nomb.c_str());
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Horse::open_ofstream.\n";
		exit(1);
	}
}

string Horse::inttostring(int num)
{
	char buff[50];
	int j;
	j = sprintf(buff, "%d", num);
	if (j <0) {
		cout << "[Error]: Transformation of string to int was not possible using Horse::inttostring.\n";
		exit(1);
	}
	string res(buff);
	return res;
}
