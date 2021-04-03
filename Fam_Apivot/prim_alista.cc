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
  
  
  string losdir = "GraphI/pn1/12N/36e/3_3_6/";
  
  
  
  
  int nunods = 12, nuedges = 36, dcia=3, dcib=3, dab=6;
  
////  for (i=0; i< 15; i++) {
//    orden = "rm "+losdir+"*.sh";
//    system(orden.c_str());
//    orden = "rm "+losdir+"grabthis.txt";
//    system(orden.c_str());
    orden = "mkdir "+losdir+"Results";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Figs";
    system(orden.c_str());
    orden = "mkdir "+losdir+"networks/Male";
    system(orden.c_str());
    orden = "mkdir "+losdir+"networks/Wife1";
    system(orden.c_str());
    orden = "mkdir "+losdir+"networks/Wife2";
    system(orden.c_str());
    
    orden = "mkdir "+losdir+"Results/Male";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/Wife1";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/Wife2";
    system(orden.c_str());
    
    
    orden = "mkdir "+losdir+"Results/Midparent1";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/Midparent2";
    system(orden.c_str());
    
    orden = "mkdir "+losdir+"Results/Dist_dad_wife2";
    system(orden.c_str());
    
    orden = "mkdir "+losdir+"Results/Dist_dad_wife1";
    system(orden.c_str());
    ///
    orden = "mkdir "+losdir+"Results/offspring_f1";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/mean";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/median";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/max";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/min";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/sd";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/1q";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/3q";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/Anative";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/Bnative";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f1/Cnative";
    system(orden.c_str());
    ///
    orden = "mkdir "+losdir+"Results/offspring_f2";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/mean";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/median";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/max";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/min";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/sd";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/1q";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/3q";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/Anative";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/Bnative";
    system(orden.c_str());
    orden = "mkdir "+losdir+"Results/offspring_f2/Cnative";
    system(orden.c_str());
    
    bas.open_ofstream(fs, losdir+"nexttowr.txt");
    fs << 1 << endl;
    fs.close();
    
    bas.open_ofstream(fs, losdir + "run_Ap_" + bas.inttostring(nunods) + "_" + bas.inttostring(nuedges) + "_" + bas.inttostring(dcia) + "_" + bas.inttostring(dcib) + "_" + bas.inttostring(dab)+".sh");
    fs << "#!/bin/bash\n";
    fs << "#PBS -l walltime=98:00:00\n";
    fs << "#PBS -l nice=19\n";
    fs << "#PBS -l nodes=1:ppn=1\n";
    fs << "#PBS -q batch\n";
    fs << "cd $PBS_O_WORKDIR\n";
    fs << "../../../../fam_apiv ${PBS_ARRAYID} " << nunods << " " << nuedges << " " << dcia << " " << dcib << " " << dab << " " << "20\n";
    fs.close();
//  }
  
  
  
  jacta.close_rng();
  return 0;
}


