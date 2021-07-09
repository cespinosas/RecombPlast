#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "fitnessi.h"
#include <string>
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

FitnessI::FitnessI()
{
}

FitnessI::FitnessI(Alea& jacta) //creates instance of class FitnessI and assigns rng to jacta
{
  start_rng(jacta);
}

void FitnessI::start_rng(Alea& jacta) //assigns rng to jacta
{
  est = jacta;
  basic.start_rng(est);
}

double FitnessI::strict(GraphI &red, int *goal) //strict fitness function
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::strict was called.\n";
    exit(1);
  }
  double w= 0;
  if (red.attractor_size() == 1) {
    int i;
    w=1;
    for (i=0; i< red.number_of_nodes(); i++)
      if (goal[i] != red.attractor_element(0,i)) {
        w = 0;
        break;
      }
  }
  return w;
}

double FitnessI::strict(GraphI &red, int **goal, int lencyc)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::strict was called.\n";
    exit(1);
  }
  double w= 0;
  if (red.attractor_size() == lencyc) {
    int ini,i,j;
    int **elatra;
    elatra = new int*[lencyc];
    for (i=0; i< lencyc; i++)
      elatra[i] = new int[red.number_of_nodes()];
    for (i=0; i< lencyc; i++)
      for (j=0; j< red.number_of_nodes(); j++)
        elatra[i][j] = red.attractor_element(i,j);
    ini = -1;
    for (i=0; i< lencyc; i++)
      if (basic.eqvec(goal[0], red.number_of_nodes(), elatra[i], red.number_of_nodes())) {
        ini = i;
        break;
      }
    if (ini >= 0) {
      w = 1;
      for (i=0; i < lencyc; i++)
        if (!basic.eqvec(goal[i], red.number_of_nodes(), elatra[(ini+i)%lencyc], red.number_of_nodes())) {
          w = 0;
          break;
        }
    }
    for (i=0; i< lencyc; i++)
      delete [] elatra[i];
    delete [] elatra;
  }
  return w;
}

double FitnessI::distance(int *goal, int *atr, int tam)
{
  int dif = 0,i;
  for (i=0; i<tam; i++)
    if (goal[i] != atr[i])
      dif++;
  return double(dif)/double(tam);
}

double FitnessI::distance(int *goal, int **atr, int numatr, int tam)
{
  int i,j;
  double dif = 0;
  for (i=0; i<numatr; i++)
    for (j=0; j < tam; j++)
      if (goal[j] != atr[i][j])
        dif = dif + (1.0/double(numatr));
  return dif/double(tam);
}

double FitnessI::distance(int **goal, int nugs, int **atr, int numatr, int tam)
{
  int i,j,k;
  double dis;
  if (nugs != numatr) {
    int mulng, mulna;
    if (nugs < numatr) {
      if ((numatr%nugs)==0) {
        mulna = 1;
        mulng = numatr/nugs;
      }
      else {
        mulng = numatr;
        mulna = nugs;
      }
    }
    else {
      if ((nugs%numatr)==0) {
        mulng = 1;
        mulna = nugs/numatr;
      }
      else {
        mulng = numatr;
        mulna = nugs;
      }
    }
    int **nugoal, **nuatr;
    nugoal = new int*[nugs*mulng];
    for (i= 0; i<(nugs*mulng); i++)
      nugoal[i] = new int[tam];
    nuatr = new int*[numatr*mulna];
    for (i=0; i <(numatr*mulna); i++)
      nuatr[i] = new int[tam];
    for (i=0; i<mulng; i++)
      for (j=0; j<nugs; j++)
        for (k=0; k < tam; k++)
          nugoal[(nugs*i)+j][k] = goal[j][k];
    for (i=0; i < mulna; i++)
      for (j=0; j<numatr; j++)
        for (k=0; k<tam;k++)
          nuatr[(numatr*i)+j][k] = atr[j][k];
    dis = distance_aux(nugoal, (nugs*mulng), nuatr, (numatr*mulna), tam);
    for (i= 0; i<(nugs*mulng); i++)
      delete [] nugoal[i];
    delete [] nugoal;
    for (i=0; i <(numatr*mulna); i++)
      delete [] nuatr[i];
    delete [] nuatr;
  }
  else
    dis = distance_aux(goal, nugs, atr, numatr, tam);
  return dis;
}

//private
double FitnessI::distance_aux(int **goal, int nugs, int **atr, int numatr, int tam)
{
  if (nugs != numatr) {
    cout << "[Error]: Wrong number of matrix rows in FitnessI::distance_aux.\n";
    exit(1);
  }
  int i,j,k;
  double *dists,d;
  dists = new double[nugs];
  basic.fillv0(dists, nugs);
  for (i=0; i<nugs; i++) {
    dists[i] =0;
    for (j=0; j<nugs; j++) {
      for (k=0;k<tam;k++)
        if (goal[j][k] != atr[(j+i)%nugs][k])
          dists[i] = dists[i] + (1.0/double(nugs));
    }
  }
  d = basic.find_min(dists, nugs);
  d = d/double(tam);
  delete [] dists;
  return d;
}
