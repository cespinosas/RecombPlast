#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphb.h"
#include "graphi.h"
#include "graphc.h"


using std::cout;
using std::cin;
using std::getline;
using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

int main(int argc, char **argv) {
  ifstream fe;
  ofstream fs;

  int seed, nuredes, nodes, edges, redesxarch, numarch, mame, mult, dcia, dcib, dab;
  seed = atoi(argv[1]);
  nuredes = atoi(argv[2]); //por job
  nodes = atoi(argv[3]);
  edges = atoi(argv[4]);
  redesxarch = atoi(argv[5]);
  mame = atoi(argv[6]);
  mult = atoi(argv[7]);
  dcia = atoi(argv[8]);
  dcib = atoi(argv[9]);
  dab = atoi(argv[10]);
  
  
  if ((nuredes%redesxarch)!=0) {
    cout << "[Error]: The number of networks is not a multiple integer of the number of networks per file.\n";
    exit(1);
  }
  numarch = nuredes/redesxarch;

  Alea jacta(seed);
  Basics bas(jacta);
  int i;
  int *ci, *phenA, *phenB;
  ci = new int[nodes];
  phenA = new int[nodes];
  phenB = new int[nodes];
  
  for (i=0; i < nodes; i++) {
    if ((i%2)==0)
      ci[i] = 1;
    else
      ci[i] = -1;
    if (i < (nodes-dcia))
      phenA[i] = ci[i];
    else
      phenA[i] = ci[i]*(-1);
    phenB[i] = ci[i];
  }
  if ((dcia+dcib)==dab) {
    for (i = (nodes-dcia)-1; i >= (nodes-dab); i--)
      phenB[i]*=(-1);
  } else {
    if ((dcia+dcib)>dab) {
      //overlap: ((dcia+dcib) -dab)/2
      for (i = (nodes - ((dcia+dcib+dab)/2)); i < (nodes + ((dcib-dab-dcia)/2)); i++) 
        phenB[i]*=(-1);
    } else {
      cout << "[Error]: Impossible phenotypes.\n";
      exit(1);
    }
  }


  if ((bas.difsinvecs(ci, nodes, phenA, nodes) != dcia) || (bas.difsinvecs(ci, nodes, phenB, nodes) != dcib) ||
      (bas.difsinvecs(phenB, nodes, phenA, nodes) != dab)) {
    cout << "[Error]: Distances do not match.\n";
    exit(1);
  }
  
  GraphI co(jacta);
  co.rwalk_2fp1t(nodes, edges, phenA, phenB, ci, nuredes, "networks/", mame, mult, true, numarch, ((seed-1)*numarch)+1);
  
  delete [] ci;
  delete [] phenB;
  delete [] phenA;
  jacta.close_rng();
  return 0;
}
