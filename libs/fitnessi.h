#ifndef FITNESSI_H
#define FITNESSI_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "graphi.h"

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ifstream;
using std::istream;
using std::string;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;

class FitnessI
{
public:
	FitnessI();
  FitnessI(Alea& jacta);
  void start_rng(Alea& jacta);
  double strict(GraphI &red, int *goal);
  double strict(GraphI &red, int **goal, int lencyc);
  
  double distance(int *goal, int *atr, int tam);
  double distance(int *goal, int **atr, int numatr, int tam);
  double distance(int **goal, int nugs, int **atr, int numatr, int tam);
  
private:
  Alea est;
  Basics basic;
  double distance_aux(int **goal, int nugs, int **atr, int numatr, int tam);

};

#endif
