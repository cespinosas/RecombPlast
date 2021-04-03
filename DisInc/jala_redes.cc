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
  
  
  string losdir;
  losdir = "GraphI/pn1/12N/36e/3_3_6";
    bas.open_ifstream(fe, "../Sample/"+losdir+"/grabthis.txt");
    fe >> desde;
    fe.close();
    bas.open_ifstream(fe, losdir+"/nexttowr.txt");
    fe >> para;
    fe.close();
    for (j=0; j < nufiles; j++) {
      orden = "mv ../Sample/"+losdir+"/networks/"+bas.inttostring(desde+j)+".txt "+losdir+"/networks/Male/"+bas.inttostring(para+j)+".txt";
      system(orden.c_str());
    }
    bas.open_ofstream(fs, "../Sample/"+losdir+"/grabthis.txt");
    fs << desde+nufiles;
    fs.close();
    bas.open_ofstream(fs, losdir+"/nexttowr.txt");
    fs << para+nufiles;
    fs.close();
  
  

  jacta.close_rng();
  return 0;
}


