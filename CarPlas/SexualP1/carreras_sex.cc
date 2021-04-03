#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphi.h"
#include "fitnessi.h"
#include "evolvei.h"


using std::cout;
using std::cin;
using std::getline;
using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

int main(int argc, char **argv) {
  if ((argc != 10) && (argc != 16)) {
    cout << "[Error]; Wrong number of arguments in program. It should be ./asu seed nodes edges popsize murate generations pertrate selcoefpg gfromintsctn_one ...\n";
    exit(1);
  }
  ofstream fs;
  int seed = atoi(argv[1]);
  int nodes = atoi(argv[2]);
  int edges = atoi(argv[3]);
  int popsize = atoi(argv[4]);
  double murate = atof(argv[5]);
  int generations = atoi(argv[6]);
  double wol = double(edges)/double(nodes*nodes);
  double pertrate = atof(argv[7]);
  double selcoefpg = atof(argv[8]);
  int numint = atoi(argv[9]); //numint
  bool genfromint = false;
  if (numint == 1)
    genfromint = true;
  
  ifstream fe;
  int dcia, dcib, dab, dcic, dbc, dac;
  if (argc == 10) {
    dcia = 3;
    dcib = 3;
    dcic = 3;
    dab = 6;
    dbc = 6;
    dac = 6;
  } else {
    dcia = atoi(argv[10]);
    dcib = atoi(argv[11]);
    dcic = atoi(argv[12]);
    dab = atoi(argv[13]);
    dbc = atoi(argv[14]);
    dac = atoi(argv[15]);
  }
  
  Alea jacta(seed);
  Basics bas(jacta);
  EvolveI chw(jacta);
  FitnessI fifu(jacta);
  string orden;
  
  int i, j,k,l;
  double *defaultA, *defaultB, *defaultC, *defaultO, *wmax, *meanedges;
  
  defaultA = new double[generations];
  defaultB = new double[generations];
  defaultC = new double[generations];
  defaultO = new double[generations];
  wmax = new double[generations];
  meanedges = new double[generations];
  
  bas.fillvm1(defaultA, generations);
  bas.fillvm1(defaultB, generations);
  bas.fillvm1(defaultC, generations);
  bas.fillvm1(defaultO, generations);
  bas.fillvm1(wmax, generations);
  bas.fillvm1(meanedges, generations);
  int tpC = -1, tpB = -1, t50C = -1, t50B = -1, tdC = -1, tdB = -1;
  
  
  int *ci, *phenA, *phenB, *phenC;
  ci = new int[nodes];
  phenA = new int[nodes];
  phenB = new int[nodes];
  phenC = new int[nodes];
  for (i=0; i < nodes; i++) {
    if ((i%2)==0)
      ci[i] = 1;
    else
      ci[i] = -1;
    if (i < (nodes-dcia))
      phenA[i] = ci[i];
    else
      phenA[i] = ci[i]*(-1);
    if (i < dcic)
      phenC[i] = ci[i]*(-1);
    else
      phenC[i]  = ci[i];
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
      (bas.difsinvecs(phenB, nodes, phenA, nodes) != dab) || (bas.difsinvecs(ci, nodes, phenC, nodes) != dcic) || (bas.difsinvecs(phenC, nodes, phenB, nodes) != dbc) || (bas.difsinvecs(phenC, nodes, phenA, nodes) != dac)) {
    cout << "[Error]: Distances do not match.\n";
    exit(1);
  }
  
  string cuerda;
  
  if (genfromint)
    cuerda = "Results/pn1/";
  else
    cuerda = "Results/sinplas/";
  
  chw.set_params(murate, wol, popsize, nodes);
  bas.open_ofstream(fs, cuerda+"pars/"+bas.inttostring(seed)+"pars.txt");
  chw.print_params(fs);
  fs.close();
  
  GraphI adam(jacta);
  if (genfromint)
    cuerda = "networks/pn1/"+bas.inttostring(seed)+".txt";
  else
    cuerda = "networks/sinplas/"+bas.inttostring(seed)+".txt";
  
  bas.open_ifstream(fe, cuerda);
  adam.get_dir_nw_from_file(nodes, fe, edges);
  fe.close();
  
  int **matatr;
  
  chw.start_pop(adam);//nodes, edges, ics[0], goals[0][0]);
  double w, wb, wc, prope, disb, disc;
  int defa, defb, defc, sume;
  
  int *oci;
  oci = new int[nodes];
  for (i=0; i< nodes; i++)
    oci[i] = ci[i];
  
  for (i=0; i<generations; i++) {
    if ((i%100)==0)
      cout << i << endl;
    defa = 0;
    defb = 0;
    defc = 0;
    sume = 0;
    for (j=0; j < popsize; j++) {
      sume += chw.population[j].number_of_edges();
      //para contar
      for (k=0; k < nodes; k++)
        oci[k] = ci[k];
      chw.population[j].set_as_state(oci);
      chw.population[j].find_an_attractor();
      if (fifu.strict(chw.population[j], phenA) > 0) {
        defa++;
      }
      else {
        if (fifu.strict(chw.population[j], phenB) > 0) {
          defb++;
        }
        else {
          if (fifu.strict(chw.population[j], phenC) > 0) {
            defc++;
          }
        }
      }
      chw.population[j].clear_attractor();
      //paraevolu:
      for (k=0; k < nodes; k++) {
        prope = jacta.randreal();
        if (prope < pertrate)
          oci[k] = ci[k]*(-1);
        else
          oci[k] = ci[k];
      }
      chw.population[j].set_as_state(oci);
      chw.population[j].find_an_attractor();
      bas.create_array(matatr, chw.population[j].attractor_size(), nodes);
      for (k=0; k < chw.population[j].attractor_size(); k++)
        for (l=0; l< nodes; l++)
          matatr[k][l] = chw.population[j].attractor_element(k,l);
      disb = nodes*fifu.distance(phenB, matatr, chw.population[j].attractor_size(), nodes);
      disc = nodes*fifu.distance(phenC, matatr, chw.population[j].attractor_size(), nodes);
      for (k=0; k < chw.population[j].attractor_size(); k++)
        delete [] matatr[k];
      delete [] matatr; ////
      wb = pow((1 - selcoefpg), disb);
      wc = pow((1 - selcoefpg), disc);
      if (wb > wc)
        w = wb;
      else
        w = wc;
      chw.assign_w(j, w);
    }
    chw.calc_meanw();
    defaultC[i] = defc/(1.0*popsize);
    defaultB[i] = defb/(1.0*popsize);
    defaultA[i] = defa/(1.0*popsize);
    defaultO[i] = (popsize - defa - defb - defc)/(1.0*popsize);
    wmax[i] = chw.return_maxw();
    meanedges[i] = sume/(1.0*popsize);
    if ((defc > 0) && (tpC < 0)) {
      tpC = i;
      tdC = i;
    }
    if ((defc == 0) && (tpC >= 0))
      tdC = -1;
    if ((defc > 0) && (tpC >= 0))
      tdC = i;
    
    if ((defb > 0) && (tpB < 0)) {
      tpB = i;
      tdB = i;
    }
    if ((defb == 0) && (tpB >= 0))
      tdB = -1;
    if ((defb > 0) && (tpB >= 0) && (tdB < 0))
      tdB = i;
    if (defaultC[i] >= 0.50)
      t50C = i;
    
    if (defaultB[i] >= 0.50)
      t50B = i;
    
    if ((t50B > 0) || (t50C > 0))
      break;
    
    if ((i < (generations-1)) && (t50C < 0) && (t50B < 0))
      chw.one_generation_sex();
  }
  //checar
  
  int hasta = i;
  
  
  if (genfromint)
    cuerda = "Results/pn1/";
  else
    cuerda = "Results/sinplas/";
  
  
  
  if (hasta == generations) {
    bas.open_ofstream(fs, cuerda+"errors/"+bas.inttostring(seed)+"neverended.txt");
    fs << "never ended.\n";
    fs.close();
  }
  hasta++;
  
  bas.open_ofstream(fs, cuerda+"defaultA/"+bas.inttostring(seed)+"defaultA.txt"); //voyaquiiiiiiiii
  for (i=0; i< hasta; i++)
    fs << i << "\t" << defaultA[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"defaultB/"+bas.inttostring(seed)+"defaultB.txt");
  for (i=0; i< hasta; i++)
    fs << i << "\t" << defaultB[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"defaultC/"+bas.inttostring(seed)+"defaultC.txt");
  for (i=0; i< hasta; i++)
    fs << i << "\t" << defaultC[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"defaultO/"+bas.inttostring(seed)+"defaultO.txt");
  for (i=0; i< hasta; i++)
    fs << i << "\t" << defaultO[i] << endl;
  fs.close();
  
  
  
  bas.open_ofstream(fs, cuerda+"wmax/"+bas.inttostring(seed)+"wmax.txt");
  for (i=0; i< hasta; i++)
    fs << i << "\t" << wmax[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"meanedges/"+bas.inttostring(seed)+"meanedges.txt");
  for (i=0; i< hasta; i++)
    fs << i << "\t" << meanedges[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"tpC/"+bas.inttostring(seed)+"tpC.txt");
  fs << tpC << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"tdC/"+bas.inttostring(seed)+"tdC.txt");
  fs << tdC << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"t50C/"+bas.inttostring(seed)+"t50C.txt");
  fs << t50C << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"tpB/"+bas.inttostring(seed)+"tpB.txt");
  fs << tpB << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"tdB/"+bas.inttostring(seed)+"tdB.txt");
  fs << tdB << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"t50B/"+bas.inttostring(seed)+"t50B.txt");
  fs << t50B << endl;
  fs.close();
  
  //  int tpC = -1, tpB = -1, t50C = -1, t50B = -1, tdC = -1, tdB = -1;
  
  delete [] phenA;
  delete [] phenB;
  delete [] phenC;
  delete [] ci;
  delete [] oci;
  delete [] defaultA;
  delete [] defaultB;
  delete [] defaultC;
  delete [] defaultO;
  delete [] wmax;
  delete [] meanedges;
  
  jacta.close_rng();
  cout << 1 << endl;
  return 0;
}


