#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include <string>
#include "basics.h"
#include "graphi.h"
#include "fitnessi.h"
#include "evolvei.h"

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

EvolveI::EvolveI()
{
  yaparams = false;
  yapop = false;
  yalop = false;
}

EvolveI::EvolveI(Alea& jacta) //creates instance of class EvolveI and assigns rng to jacta
{
  start_rng(jacta);
  yaparams = false;
  yapop = false;
  yalop = false;
}

void EvolveI::start_rng(Alea& jacta) //assigns rng to jacta 
{
  est = jacta;
  basic.start_rng(est);
}

void EvolveI::set_params(double muratep, double wolp, int popsizep, int nugenp) //sets parameters for evolutionary simulations
{
  murate = muratep;
  wol = wolp;
  popsize = popsizep;
  nugen = nugenp;
  yaparams = true;
}


void EvolveI::print_params(ostream& sal) //prints parameters
{
  sal << "Integer weights\n";
  sal << "Number of genes: " << nugen << endl;
  sal << "Population size: " << popsize << endl;
  sal << "Mutation rate: " << murate << endl;
  sal << "Expected connectivity: " << wol << endl;
  sal << endl << endl;
}

void EvolveI::start_pop(GraphI &founder)
{
  if (!yaparams) {
    cout << "[Error]: Parameters had not been set when EvolveI::start_pop was called.\n";
    exit(1);
  }
  if (yapop) {
    cout << "[Error]: Population already constructed when EvolveI::start_pop was called.\n";
    exit(1);
  }
  population = new GraphI[popsize];
  nepop = new GraphI[popsize];
  w = new double[popsize];
  basic.fillv0(w, popsize);
  meanw =0;
  maxw = 0;
  int i;
  for (i=0; i<popsize; i++) {
    population[i].start_rng(est);
    nepop[i].start_rng(est);
    population[i].copy(founder);
    population[i].consider_mutation(murate, wol);
  }
  mean_gen_dist=0;
  mean_gen_dist_thr = 0;
  max_gen_dist = 0;
  max_gen_dist_thr = 0;
  int nucomb = (popsize*(popsize-1))/2;
  all_gdists = new double[nucomb];
  basic.fillv0(all_gdists, nucomb);
  all_gdists_thr = new double[nucomb];
  basic.fillv0(all_gdists_thr, nucomb);
  mean_ph_dist = 0;
  max_ph_dist = 0;
  all_pdists = new double[nucomb];
  yapop = true;
}

void EvolveI::one_generation_sex()
{
  int i,j,k;
  double frog, toad, tadp;
  for (i=0; i<popsize; i++) {
    frog = est.randreal()*meanw*popsize;
    toad = 0;
    for (j=0; j<popsize; j++) {
      toad = toad + w[j];
      if (toad >= frog)
        break;
    }
    frog = est.randreal()*meanw*popsize;
    tadp = 0;
    for (k=0; k< popsize; k++) {
      tadp = tadp + w[k];
      if (tadp >= frog)
        break;
    }
    nepop[i].mate(population[j], population[k]);
    nepop[i].consider_mutation(murate, wol);
  }
  for (i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
}

void EvolveI::assign_w(int cual, double cuanto)
{
  if ((cual >= popsize) || (cual < 0)) {
    cout << "[Error]: Individual does not exist. EvolveI::assign_w.\n";
    exit(1);
  }
  w[cual] = cuanto;
}

void EvolveI::calc_meanw()
{
  meanw = basic.get_mean(w, popsize);
  maxw = basic.find_max(w, popsize);
}

double EvolveI::return_meanw()
{
  return meanw;
}

double EvolveI::return_maxw()
{
  return maxw;
}

void EvolveI::clear()
{
  int i;
  if (yapop) {
    for (i=0; i<popsize; i++)
      population[i].clear();
    delete [] population;
    delete [] nepop;
    delete [] all_gdists_thr;
    delete [] all_gdists;
    delete [] w;
    delete [] all_pdists;
  }
  yapop = false;
  yaparams = false;
  if (yalop)
    clear_lop();
}

void EvolveI::clear_lop()
{
  int i,j,k;
  for (i=0; i < vanop; i++) {
    for (j=0; j < nulop[i]; j++) {
      for (k=0; k < perlop[i][j]; k++)
        delete [] last_opt_phen[i][j][k];
      delete [] last_opt_phen[i][j];
    }
    delete [] perlop[i];
    delete [] last_opt_phen[i];
  }
  delete [] nulop;
  delete [] perlop;
  delete [] last_opt_phen;
  vanop = 0;
  yalop = false;
}

double EvolveI::get_murate()
{
  return murate;
}

double EvolveI::get_wol()
{
  return wol;
}

void EvolveI::new_wol(double nwol) {
  wol = nwol;
}

int EvolveI::get_popsize()
{
  return popsize;
}

int EvolveI::get_nugen()
{
  return nugen;
}

double EvolveI::get_w(int cual) {
  return w[cual];
}


