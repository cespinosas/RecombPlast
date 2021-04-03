//This program takes the data from fami.m in order to separate the random networks that produce two GAPs with 0 offspring with native B and parents with 1 or more offspring with native B.

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"


using std::cout;
using std::cin;
using std::getline;
using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

int main(int argc, char **argv) {
  
  int seed = atoi(argv[1]);
  if (argc != 2) {
    cout << "[Error].\n";
    exit(1);
  }
  Alea jacta(seed);
  Basics bas(jacta);
  string eldir = "S"+bas.inttostring(seed)+"/";
  ofstream fs0, fs1, fst;
  ifstream fe;
  double gar;
  double pci, has;
  int cuecer = 0;
  int cueun = 0;
  
  bas.open_ofstream(fst, eldir+"has_pci.txt");
  bas.open_ofstream(fs0, eldir+"has_en_pci0.txt");
  bas.open_ofstream(fs1, eldir+"has_en_pci1.txt");
  bas.open_ifstream(fe, eldir+"pci_pmu_hci_hmu_has.txt");
  while (fe >> pci) {
    fe >> gar;
    fe >> gar;
    fe >> gar;
    fe >> has;
    if (has == 0) {
      cuecer++;
      fs0 << pci << endl;
      fst << 0 << "\t" << pci << endl;
    } else {
      cueun++;
      fs1 << pci << endl;
      fst << 1 << "\t" << pci << endl;
    }
  }
  fe.close();
  fs1.close();
  fs0.close();
  fst.close();
  bas.open_ofstream(fs0, eldir+"count_as.txt");
  fs0 << "Pairs with 0 offspring with native B: " << cuecer << endl;
  fs0 << "Pairs with 1 or more offspring with native B: " << cueun << endl;
  fs0.close();

  jacta.close_rng();
  return 0;
}


