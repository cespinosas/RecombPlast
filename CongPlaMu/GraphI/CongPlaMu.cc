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
  ofstream fs;
  ifstream fe;
  int seed = 1;
  Alea jacta(seed);
  Basics bas(jacta);
  string orden;
  
  int nodes = atoi(argv[1]);
  int edges = atoi(argv[2]);
  int dcia = atoi(argv[3]);
  int dcib = atoi(argv[4]);
  int dab = atoi(argv[5]);
  int numarch = atoi(argv[6]);
  int nuredesperarch = atoi(argv[7]);
  int i, j,k;
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
  
  string tras;
  GraphI co(jacta);

  int totredes = nuredesperarch*numarch;
  
  //para plasticidad
  double *pn1B, *pn2B, *maB;
  pn1B = new double[totredes];
  pn2B = new double[totredes];
  maB = new double[totredes];
  bas.fillv0(pn1B, totredes);
  bas.fillv0(pn2B, totredes);
  bas.fillv0(maB, totredes);
  
  //para robustez
  double *pn1A, *pn2A, *maA;
  pn1A = new double[totredes];
  pn2A = new double[totredes];
  maA = new double[totredes];
  bas.fillv0(pn1A, totredes);
  bas.fillv0(pn2A, totredes);
  bas.fillv0(maA, totredes);
  
  int temA, temB;
  
  for (i=1; i <= numarch; i++) {
    bas.open_ifstream(fe, "networks/"+bas.inttostring(i)+".txt");
    for (j=0; j < nuredesperarch; j++) {
      co.get_dir_nw_from_file(nodes, fe, edges);
      fe >> tras;
      k = ((i-1)*nuredesperarch)+j;
      co.penetrance_n1_2p(ci, phenA, phenB, temA, temB);
      pn1A[k] = temA;///(nodes*1.0);
      pn1B[k] = temB;///(nodes*1.0);
      co.penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
      pn2A[k] = temA;//(2.0*(pn1A[k]+temA))/(nodes*(nodes+1.0));
      pn2B[k] = temB;//(2.0*(pn1B[k]+temB))/(nodes*(nodes+1.0));
      co.mutational_access_2p(ci, phenA, phenB, temA, temB);
      maA[k] = temA/((2.0*nodes*nodes)-edges);
      maB[k] = temB/((2.0*nodes*nodes)-edges);
      co.clear();
    }
    fe.close();
  }
  
  bas.open_ofstream(fs, "Results/readme.txt");
  fs << "B through plasticity (n1)\n";
  fs << "B through plasticity (n2)\n";
  fs << "B through mutations\n";
  fs << "A through plasticity (n1)\n";
  fs << "A through plasticity (n2)\n";
  fs << "A through mutations\n";
  fs.close();
  
  bas.open_ofstream(fs, "Results/CongPlaMu.txt");
  for (i=0; i<totredes; i++) {
    fs << pn1B[i] << "\t" << pn2B[i] << "\t" << maB[i] << "\t" << pn1A[i] << "\t" << pn2A[i] << "\t" << maA[i] << endl;
  }
  fs.close();
  
 
  
  delete [] ci;
  delete [] phenA;
  delete [] phenB;
  delete [] pn1B;
  delete [] pn2B;
  delete [] maB;
  delete [] pn1A;
  delete [] pn2A;
  delete [] maA;
  jacta.close_rng();
  return 0;
}


