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
  
  int nunws = atoi(argv[1]);
  int nwsperfile = atoi(argv[2]);
  if ((nunws%nwsperfile) != 0) {
    cout << "[Error]: The number of networks is not an integer multiple of the number of networks per file.\n";
    exit(1);
  }
  int nufiles = nunws/nwsperfile;
  int desde, para, i, j;
  ofstream fs;
  ifstream fe;
  int seed = 1;
  Alea jacta(seed);
  Basics bas(jacta);
  string orden;
  
  
  string losdir = "GraphI/pn1/12N/36e/3_3_6";
  
  
  string losdp = "GraphI/sinplas/12N/36e/3_3_6";
  
  bas.open_ifstream(fe, losdir+"/nexttowr.txt");
  fe >> para;
  fe.close();
  
  bas.open_ifstream(fe, "../Sample/"+losdp+"/grabthis.txt"); //dp
  fe >> desde;
  fe.close();
  for (j=0; j < nufiles; j++) {
    orden = "mv ../Sample/"+losdp+"/networks/"+bas.inttostring(desde+j)+".txt "+losdir+"/networks/Male/"+bas.inttostring(para+j)+".txt";//dp losdirmale
    system(orden.c_str());
  }
  for (j=0; j < nufiles; j++) {
    orden = "mv ../Sample/"+losdp+"/networks/"+bas.inttostring(desde+nufiles+j)+".txt "+losdir+"/networks/Wife1/"+bas.inttostring(para+j)+".txt"; //dp losdirw1
    system(orden.c_str());
  }
  bas.open_ofstream(fs, "../Sample/"+losdp+"/grabthis.txt"); //dp
  fs << desde+nufiles+nufiles;
  fs.close();
  
  bas.open_ifstream(fe, "../Sample/"+losdir+"/grabthis.txt"); //dir
  fe >> desde;
  fe.close();
  for (j=0; j < nufiles; j++) {
    orden = "mv ../Sample/"+losdir+"/networks/"+bas.inttostring(desde+j)+".txt "+losdir+"/networks/Wife2/"+bas.inttostring(para+j)+".txt"; //dir losdirw2
    system(orden.c_str());
  }
  bas.open_ofstream(fs, "../Sample/"+losdir+"/grabthis.txt"); //dir
  fs << desde+nufiles;
  fs.close();
  
  bas.open_ofstream(fs, losdir+"/nexttowr.txt");
  fs << para+nufiles;
  fs.close();
  
  
  jacta.close_rng();
  return 0;
}


