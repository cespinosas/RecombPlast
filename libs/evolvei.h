#ifndef EVOLVEI_H
#define EVOLVEI_H

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
#include "fitnessi.h"

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

class EvolveI
{
public:
	EvolveI();
  EvolveI(Alea& jacta);
  void start_rng(Alea& jacta);
  void set_params(double muratep, double wolp, int popsizep, int nugenp);
  void print_params(ostream& sal);
  void start_pop(GraphI &founder);
  void one_generation_sex();
  void assign_w(int cual, double cuanto);
  void calc_meanw();
  double return_meanw();
  double return_maxw();
  void clear();
  void clear_lop();
  double get_murate();
  double get_wol();
  void new_wol(double nwol);
  int get_popsize();
  int get_nugen();
  double get_w(int cual);
  
  
  GraphI *population;
  double mean_gen_dist, mean_gen_dist_thr, max_gen_dist, max_gen_dist_thr, *all_gdists, *all_gdists_thr;
  double mean_ph_dist, max_ph_dist, *all_pdists;
  
  int ****last_opt_phen; //
  int **perlop;
  int *nulop;
  
private:

  Alea est;
  Basics basic;
  double murate;
  double wol;
  int popsize;
  int nugen;
  bool yaparams;
  double meanw;
  double maxw;
  double *w;
  GraphI *nepop;
  bool yapop;
  
  bool yalop;//
  int vanop;//, generini;

};

#endif
