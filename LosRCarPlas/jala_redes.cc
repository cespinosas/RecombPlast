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
  
  //hacer directorio networkspn1 y sinplas
  
  int nunws = atoi(argv[1]);
  int nwsperfile = atoi(argv[2]);
  string cualdir = argv[3]; //incluir diagonal

  if ((nunws%nwsperfile) != 0) {
    cout << "[Error]: The number of networks is not an integer multiple of the number of networks per file.\n";
    exit(1);
  }
  int nodes = 12;
  int edges = 36;
  
  
  
  int nufiles = nunws/nwsperfile;
  int desde, para, i, j;
  ofstream fs;
  ifstream fe;
  int seed = 1;
  Alea jacta(seed);
  Basics bas(jacta);
  string orden;
  
  GraphI co(jacta);
  string entpn1 = "../SampleSP/";
  string entsinp = "../Sample/GraphI/sinplas/12N/36e/3_3_6/";
  
  
//  string losdir = "GraphI/pn1/12N/36e/3_3_6";
//  string sinplas = "GraphI/sinplas/12N/36e/3_3_6";
  
  para = 1;
  bas.run_command("mkdir "+cualdir+"/networks");
  bas.run_command("mkdir "+cualdir+"/networks/pn1");
  bas.run_command("mkdir "+cualdir+"/networks/sinplas");

//  bas.open_ifstream(fe, "nexttowr.txt");
//  fe >> para;
//  fe.close();
  
  bas.open_ifstream(fe, entpn1+"grabthis.txt");
  fe >> desde;
  fe.close();
  for (i=0; i < nufiles; i++) {
    bas.open_ifstream(fe, entpn1+"networks/"+bas.inttostring(desde+i)+".txt");
    for (j=0; j < nwsperfile; j++) {
      co.get_dir_nw_from_file(nodes, fe, edges);
      fe >> orden;
      co.export_nw(cualdir+"/networks/pn1/"+bas.inttostring((i*nwsperfile)+j+para)+".txt");
      co.clear();
    }
    fe.close();
    orden = "rm "+entpn1+"networks/"+bas.inttostring(desde+i)+".txt";
    system(orden.c_str());
  }
  bas.open_ofstream(fs, entpn1+"grabthis.txt");
  fs << desde+nufiles;
  fs.close();
  
  
  bas.open_ifstream(fe, entsinp+"grabthis.txt");
  fe >> desde;
  fe.close();
  for (i=0; i < nufiles; i++) {
    bas.open_ifstream(fe, entsinp+"networks/"+bas.inttostring(desde+i)+".txt");
    for (j=0; j < nwsperfile; j++) {
      co.get_dir_nw_from_file(nodes, fe, edges);
      fe >> orden;
      co.export_nw(cualdir+"/networks/sinplas/"+bas.inttostring((i*nwsperfile)+j+para)+".txt");
      co.clear();
    }
    fe.close();
    orden = "rm "+entsinp+"networks/"+bas.inttostring(desde+i)+".txt";
    system(orden.c_str());
  }
  bas.open_ofstream(fs, entsinp+"grabthis.txt");
  fs << desde+nufiles;
  fs.close();
  

  jacta.close_rng();
  return 0;
}


