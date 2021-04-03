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
  
  
  string losdir[16];
  
  losdir[0] = "GraphC/pn1/12N/36e/3_3_6";
  
  losdir[1] = "GraphI/pfp/12N/36e/3_3_6";
  
  losdir[2] = "GraphI/pn1/12N/24e/3_3_6";
  
  losdir[3] = "GraphI/pn1/12N/36e/2_3_5";
  losdir[4] = "GraphI/pn1/12N/36e/3_2_5";
  losdir[5] = "GraphI/pn1/12N/36e/3_3_4";
  losdir[6] = "GraphI/pn1/12N/36e/3_3_6";
  losdir[7] = "GraphI/pn1/12N/36e/3_4_7";
  losdir[8] = "GraphI/pn1/12N/36e/4_3_7";
  losdir[9] = "GraphI/pn1/12N/36e/4_4_8";
  
  losdir[10] = "GraphI/pn1/12N/48e/3_3_6";
  
  losdir[11] = "GraphI/pn1/16N/48e/4_4_8";
  
  losdir[12] = "GraphI/pn1/8N/24e/2_2_4";
  
  losdir[13] = "GraphI/pn2/12N/36e/3_3_6";
  
  losdir[14] = "GraphI/sinplas/12N/36e/3_3_6";

  losdir[15] = "GraphI0/pn1/12N/36e/3_3_6";
  
  for (i = 0; i < 16; i++) {
    bas.open_ifstream(fe, "../Sample/"+losdir[i]+"/grabthis.txt");
    fe >> desde;
    fe.close();
    bas.open_ifstream(fe, losdir[i]+"/nexttowr.txt");
    fe >> para;
    fe.close();
    for (j=0; j < nufiles; j++) {
      orden = "mv ../Sample/"+losdir[i]+"/networks/"+bas.inttostring(desde+j)+".txt "+losdir[i]+"/networks/"+bas.inttostring(para+j)+".txt";
      system(orden.c_str());
    }
    bas.open_ofstream(fs, "../Sample/"+losdir[i]+"/grabthis.txt");
    fs << desde+nufiles;
    fs.close();
    bas.open_ofstream(fs, losdir[i]+"/nexttowr.txt");
    fs << para+nufiles;
    fs.close();
  }
  
  

  jacta.close_rng();
  return 0;
}


