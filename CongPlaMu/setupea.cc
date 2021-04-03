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
  int i;
  ofstream fs;
  int seed = 1;
  Alea jacta(seed);
  Basics bas(jacta);
  string orden;
  
  
  string losdir[16];
  losdir[0] = "GraphI/pfp/12N/36e/3_3_6/";
  losdir[1] = "GraphI/pn1/12N/24e/3_3_6/";
  losdir[2] = "GraphI/pn1/12N/36e/2_3_5/";
  losdir[3] = "GraphI/pn1/12N/36e/3_2_5/";
  losdir[4] = "GraphI/pn1/12N/36e/3_3_4/";
  losdir[5] = "GraphI/pn1/12N/36e/3_3_6/";
  losdir[6] = "GraphI/pn1/12N/36e/3_4_7/";
  losdir[7] = "GraphI/pn1/12N/36e/4_3_7/";
  losdir[8] = "GraphI/pn1/12N/36e/4_4_8/";
  losdir[9] = "GraphI/pn1/12N/48e/3_3_6/";
  
  losdir[10] = "GraphI/pn1/16N/48e/4_4_8/";
  
  losdir[11] = "GraphI/pn1/8N/24e/2_2_4/";
  
  losdir[12] = "GraphI/pn2/12N/36e/3_3_6/";
  
  losdir[13] = "GraphI/sinplas/12N/36e/3_3_6/";
  
  losdir[14] = "GraphC/pn1/12N/36e/3_3_6/";
  losdir[15] = "GraphI0/pn1/12N/36e/3_3_6/";

  
  int nunods[16], nuedges[16], dcia[16], dcib[16], dab[16];
  for (i=0; i<16; i++) {
    nunods[i] = 12;
    nuedges[i] = 36;
    dcia[i] = 3;
    dcib[i] = 3;
    dab[i] = 6;
  }
  nunods[10] = 16;
  nunods[11] = 8;
  
  nuedges[1] = 24;
  nuedges[9] = 48;
  nuedges[10] = 48;
  nuedges[11] = 24;
  dcia[2] = 2;
  dcia[7] = 4;
  dcia[8] = 4;
  dcia[10] = 4;
  dcia[11] = 2;
  
  dcib[3] = 2;
  dcib[6] = 4;
  dcib[8] = 4;
  dcib[10] = 4;
  dcib[11] = 2;
  
  dab[2] = 5;
  dab[3] = 5;
  dab[4] = 4;
  dab[6] = 7;
  dab[7] = 7;
  dab[8] = 8;
  dab[10] = 8;
  dab[11] = 4;
  
  for (i=0; i< 16; i++) {
    orden = "rm "+losdir[i]+"*.sh";
    system(orden.c_str());
    
    orden = "rm -r "+losdir[i]+"Figs/";
    system(orden.c_str());
    orden = "rm -r "+losdir[i]+"Results/";
    system(orden.c_str());
    orden = "rm "+losdir[i]+"networks/*.txt";
    system(orden.c_str());
    orden = "rm "+losdir[i]+"nexttowr.txt";
    system(orden.c_str());
    bas.open_ofstream(fs, losdir[i]+"nexttowr.txt");
    fs << 1 << endl;
    fs.close();
    
    bas.open_ofstream(fs, losdir[i]+"correle.sh");
    fs << "../../../../CongPlaMu " << nunods[i] << " " << nuedges[i] << " " << dcia[i] << " " << dcib[i] << " " <<  dab[i] << " " << 150 <<  " " << 20 << endl;
    fs.close();
    orden = "chmod +x "+losdir[i]+"correle.sh";
    system(orden.c_str());
    orden = "mkdir "+losdir[i]+"Results";
    system(orden.c_str());
  }
  
  
  
  jacta.close_rng();
  return 0;
}


