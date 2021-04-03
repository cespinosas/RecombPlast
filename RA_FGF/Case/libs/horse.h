//g++ -Wall -c horse.cc
#ifndef HORSE_H
#define HORSE_H

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <iostream>

using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ifstream;
using std::string;
using std::ios;

class Horse 
{
public:
    gsl_rng * r;
    Horse();
    Horse(int semilla);
    void start_rng(int semilla);
    double randreal();
    void close_rng();
    void open_ifstream(ifstream& fe, string nomb);
    void open_ofstream(ofstream& fs, string nomb);
    void open_ofstream_new(ofstream& fs, string nomb);
    string inttostring(int num);
};

#endif
