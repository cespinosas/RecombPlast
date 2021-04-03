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

void getmomAB(int indice, Alea& est, int tam, GraphI &pap, GraphI &cual, int *conin, int *fa, int *fb);
void getmomA(int indice, Alea& est, int tam, GraphI &pap, GraphI &cual, int *conin, int *fa);
int quitaponAB(Alea& est, int tam, GraphI &pap, GraphI &cual, int *conin, int *fa, int *fb, int disprev);
int quitaponA(Alea& est, int tam, GraphI &pap, GraphI &cual, int *conin, int *fa, int disprev);


int main(int argc, char **argv) {
  ofstream fs;
  ifstream fe;
  
  int seed = atoi(argv[1]);
  int nodes = atoi(argv[2]);
  int edges = atoi(argv[3]);
  int dcia = atoi(argv[4]);
  int dcib = atoi(argv[5]);
  int dab = atoi(argv[6]);
  int nuredes = atoi(argv[7]);
  Alea jacta(seed);
  Basics bas(jacta);
  string orden;
  
  int i, j,l,m, inmo;
  
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
  GraphI *gwife1;
  GraphI *gwife2;
  GraphI gdad(jacta);
  GraphI hijo(jacta);
  
  gwife1 = new GraphI[5];
  gwife2 = new GraphI[5];
  for (i=0; i < 5; i++) {
    gwife1[i].start_rng(jacta);
    gwife2[i].start_rng(jacta);
  }

  
  int numhijos = pow(2, nodes);
  
  double ***wife1, ***wife2, **dad, ***midparent1, ***midparent2;
  
  bas.create_array(wife1, 5, 6, nuredes);
  bas.create_array(wife2, 5, 6, nuredes);
  bas.create_array(midparent1, 5, 6, nuredes);
  bas.create_array(midparent2, 5, 6, nuredes);
  bas.create_array(dad, 6, nuredes);
  
  double **dist_pm1, **dist_pm2;
  bas.create_array(dist_pm1, 5, nuredes);
  bas.create_array(dist_pm2, 5, nuredes);
  
  
  double ***offspring_f2_mean, ***offspring_f2_median, ***offspring_f2_max, ***offspring_f2_min, ***offspring_f2_sd, ***offspring_f2_1q, ***offspring_f2_3q;
  bas.create_array(offspring_f2_mean, 5, 6, nuredes);
  bas.create_array(offspring_f2_median, 5, 6, nuredes);
  bas.create_array(offspring_f2_max, 5, 6, nuredes);
  bas.create_array(offspring_f2_min, 5, 6, nuredes);
  bas.create_array(offspring_f2_sd, 5, 6, nuredes);
  bas.create_array(offspring_f2_1q, 5, 6, nuredes);
  bas.create_array(offspring_f2_3q, 5, 6, nuredes);
  
  int **offspring_f2_Bnat, **offspring_f2_Anat, **offspring_f2_Cnat;
  bas.create_array(offspring_f2_Anat, 5, nuredes);
  bas.create_array(offspring_f2_Bnat, 5, nuredes);
  bas.create_array(offspring_f2_Cnat, 5, nuredes);
  
  double ***offspring_f1_mean, ***offspring_f1_median, ***offspring_f1_max, ***offspring_f1_min, ***offspring_f1_sd, ***offspring_f1_1q, ***offspring_f1_3q;
  bas.create_array(offspring_f1_mean, 5, 6, nuredes);
  bas.create_array(offspring_f1_median, 5, 6, nuredes);
  bas.create_array(offspring_f1_max, 5, 6, nuredes);
  bas.create_array(offspring_f1_min, 5, 6, nuredes);
  bas.create_array(offspring_f1_sd, 5, 6, nuredes);
  bas.create_array(offspring_f1_1q, 5, 6, nuredes);
  bas.create_array(offspring_f1_3q, 5, 6, nuredes);
  
  int **offspring_f1_Bnat, **offspring_f1_Anat, **offspring_f1_Cnat;
  bas.create_array(offspring_f1_Anat, 5, nuredes);
  bas.create_array(offspring_f1_Bnat, 5, nuredes);
  bas.create_array(offspring_f1_Cnat, 5, nuredes);
  
  double **hijos; //aqui
  bas.create_array(hijos, 6, numhijos);
  
  int temA, temB;
  ifstream fepap;//, femAB, femA;
  
  double *hijosor;
  hijosor = new double[numhijos];
  int *atratemp;
  atratemp = new int[nodes];
  
  bool *dema;
  dema = new bool[nodes];
  
  int cualhij;
  bas.open_ifstream(fepap, "networks/Male/"+bas.inttostring(seed)+".txt");
  for (j=0; j < nuredes; j++) {
    gdad.get_dir_nw_from_file(nodes, fepap, edges);
    fepap >> tras;
    
    gdad.penetrance_n1_2p(ci, phenA, phenB, temA, temB);
    dad[0][j] = temB;
    dad[3][j] = temA;
    gdad.penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
    dad[1][j] = temB;
    dad[4][j] = temA;
    gdad.mutational_access_2p(ci, phenA, phenB, temA, temB);
    dad[2][j] = temB/((2.0*nodes*nodes)-edges);
    dad[5][j] = temA/((2.0*nodes*nodes)-edges);
    
    for (inmo = 0; inmo < 5; inmo++) {
      getmomAB(inmo, jacta, nodes, gdad, gwife1[inmo], ci, phenA, phenB);
      getmomA(inmo, jacta, nodes, gdad, gwife2[inmo], ci, phenA);

      gwife2[inmo].penetrance_n1_2p(ci, phenA, phenB, temA, temB);
      wife2[inmo][0][j] = temB;
      wife2[inmo][3][j] = temA;
      gwife2[inmo].penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
      wife2[inmo][1][j] = temB;
      wife2[inmo][4][j] = temA;
      gwife2[inmo].mutational_access_2p(ci, phenA, phenB, temA, temB);
      wife2[inmo][2][j] = temB/((2.0*nodes*nodes)-edges);
      wife2[inmo][5][j] = temA/((2.0*nodes*nodes)-edges);
      dist_pm2[inmo][j] = gdad.distance_from_nw(gwife2[inmo]);
      
      gwife1[inmo].penetrance_n1_2p(ci, phenA, phenB, temA, temB);
      wife1[inmo][0][j] = temB;
      wife1[inmo][3][j] = temA;
      gwife1[inmo].penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
      wife1[inmo][1][j] = temB;
      wife1[inmo][4][j] = temA;
      gwife1[inmo].mutational_access_2p(ci, phenA, phenB, temA, temB);
      wife1[inmo][2][j] = temB/((2.0*nodes*nodes)-edges);
      wife1[inmo][5][j] = temA/((2.0*nodes*nodes)-edges);
      dist_pm1[inmo][j] = gdad.distance_from_nw(gwife1[inmo]);
      
      for (l=0; l < 6; l++) {
        midparent2[inmo][l][j] = (dad[l][j]+wife2[inmo][l][j])/2.0;
        midparent1[inmo][l][j] = (dad[l][j]+wife1[inmo][l][j])/2.0;
      }
      //acabe con padres
      bas.fillv0(dema, nodes);
      cualhij = 0;
      offspring_f2_Anat[inmo][j] = 0;
      offspring_f2_Bnat[inmo][j] = 0;
      offspring_f2_Cnat[inmo][j] = 0;
      //hijos maA
      while(cualhij < numhijos){ //}!bas.last_vector(dema, nodes)) {
        hijo.build_offspring(gwife2[inmo], gdad, dema);
        hijo.penetrance_n1_2p(ci, phenA, phenB, temA, temB);
        hijos[0][cualhij] = temB;
        hijos[3][cualhij] = temA;
        hijo.penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
        hijos[1][cualhij] = temB;
        hijos[4][cualhij] = temA;
        hijo.mutational_access_2p(ci, phenA, phenB, temA, temB);
        hijos[2][cualhij] = temB/((2.0*nodes*nodes)-hijo.number_of_edges());
        hijos[5][cualhij] = temA/((2.0*nodes*nodes)-hijo.number_of_edges());
        hijo.set_as_state(ci);
        hijo.find_an_attractor();
        if (hijo.attractor_size()==1) {
          for (m=0; m < nodes; m++)
          atratemp[m] = hijo.attractor_element(0, m);
          if (bas.eqvec(atratemp, nodes, phenA, nodes))
            offspring_f2_Anat[inmo][j]++;
          else {
            if (bas.eqvec(atratemp, nodes, phenB, nodes))
              offspring_f2_Bnat[inmo][j]++;
            else
              offspring_f2_Cnat[inmo][j]++;
          }
        } else
          offspring_f2_Cnat[inmo][j]++;
        hijo.clear();
        cualhij++;
        if (cualhij < numhijos)
          bas.next_in_vector(dema, nodes);
      }
      for (l=0; l < 6; l++) {
        offspring_f2_mean[inmo][l][j] = bas.get_mean(hijos[l], numhijos);
        offspring_f2_sd[inmo][l][j] = bas.get_pop_stddev(hijos[l], numhijos, offspring_f2_mean[inmo][l][j]);
        bas.sort(hijos[l], hijosor, numhijos);
        offspring_f2_min[inmo][l][j] = hijosor[0];
        offspring_f2_max[inmo][l][j] = hijosor[numhijos-1];
        offspring_f2_median[inmo][l][j] = bas.get_midpoint(hijosor, numhijos);
        offspring_f2_1q[inmo][l][j] = bas.get_q1(hijosor, numhijos);
        offspring_f2_3q[inmo][l][j] = bas.get_q3(hijosor, numhijos);
      }
      
      //hijos mAB
      bas.fillv0(dema, nodes);
      cualhij = 0;
      offspring_f1_Anat[inmo][j] = 0;
      offspring_f1_Bnat[inmo][j] = 0;
      offspring_f1_Cnat[inmo][j] = 0;
      while(cualhij < numhijos){ //!bas.last_vector(dema, nodes)) {
        hijo.build_offspring(gwife1[inmo], gdad, dema);
        hijo.penetrance_n1_2p(ci, phenA, phenB, temA, temB);
        hijos[0][cualhij] = temB;
        hijos[3][cualhij] = temA;
        hijo.penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
        hijos[1][cualhij] = temB;
        hijos[4][cualhij] = temA;
        hijo.mutational_access_2p(ci, phenA, phenB, temA, temB);
        hijos[2][cualhij] = temB/((2.0*nodes*nodes)-hijo.number_of_edges());
        hijos[5][cualhij] = temA/((2.0*nodes*nodes)-hijo.number_of_edges());
        hijo.set_as_state(ci);
        hijo.find_an_attractor();
        if (hijo.attractor_size()==1) {
          for (m=0; m < nodes; m++)
          atratemp[m] = hijo.attractor_element(0, m);
          if (bas.eqvec(atratemp, nodes, phenA, nodes))
            offspring_f1_Anat[inmo][j]++;
          else {
            if (bas.eqvec(atratemp, nodes, phenB, nodes))
              offspring_f1_Bnat[inmo][j]++;
            else
              offspring_f1_Cnat[inmo][j]++;
          }
        } else
          offspring_f1_Cnat[inmo][j]++;
        hijo.clear();
        cualhij++;
        if (cualhij < numhijos)
          bas.next_in_vector(dema, nodes);
      }
      for (l=0; l < 6; l++) {
        offspring_f1_mean[inmo][l][j] = bas.get_mean(hijos[l], numhijos);
        offspring_f1_sd[inmo][l][j] = bas.get_pop_stddev(hijos[l], numhijos, offspring_f1_mean[inmo][l][j]);
        bas.sort(hijos[l], hijosor, numhijos);
        offspring_f1_min[inmo][l][j] = hijosor[0];
        offspring_f1_max[inmo][l][j] = hijosor[numhijos-1];
        offspring_f1_median[inmo][l][j] = bas.get_midpoint(hijosor, numhijos);
        offspring_f1_1q[inmo][l][j] = bas.get_q1(hijosor, numhijos);
        offspring_f1_3q[inmo][l][j] = bas.get_q3(hijosor, numhijos);
      }
      
      gwife1[inmo].clear();
      gwife2[inmo].clear();
    }
    gdad.clear();
  }
  fepap.close();
//  femA.close();
//  femAB.close();
  
  bas.open_ofstream(fs, "Results/readme.txt");
  fs << "B through plasticity (n1)\n";
  fs << "B through plasticity (n2)\n";
  fs << "B through mutations\n";
  fs << "A through plasticity (n1)\n";
  fs << "A through plasticity (n2)\n";
  fs << "A through mutations\n";
  fs.close();
  
  bas.open_ofstream(fs, "Results/Male/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << dad[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  
  int dista[5] = {2,4,8,16,32};
  
  for (inmo=0; inmo < 5; inmo++) {
    bas.open_ofstream(fs, "Results/Wife2/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << wife2[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/Wife1/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << wife1[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    
    
    bas.open_ofstream(fs, "Results/Midparent2/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << midparent2[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/Midparent1/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << midparent1[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    
    bas.open_ofstream(fs, "Results/Dist_dad_wife2/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++)
    fs << dist_pm2[inmo][i] << endl;
    fs.close();
    bas.open_ofstream(fs, "Results/Dist_dad_wife1/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++)
    fs << dist_pm1[inmo][i] << endl;
    fs.close();
    
    
    bas.open_ofstream(fs, "Results/offspring_f2/mean/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f2_mean[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/median/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f2_median[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/max/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f2_max[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/min/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f2_min[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/sd/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f2_sd[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/1q/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f2_1q[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/3q/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f2_3q[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/Anative/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++)
    fs << offspring_f2_Anat[inmo][i] << endl;
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/Bnative/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++)
    fs << offspring_f2_Bnat[inmo][i] << endl;
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f2/Cnative/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++)
    fs << offspring_f2_Cnat[inmo][i] << endl;
    fs.close();//
    //
    //HIJOSab:
    bas.open_ofstream(fs, "Results/offspring_f1/mean/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f1_mean[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/median/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f1_median[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/max/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f1_max[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/min/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f1_min[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/sd/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f1_sd[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/1q/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f1_1q[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/3q/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++) {
      for (j=0; j < 6; j++)
      fs << offspring_f1_3q[inmo][j][i] << "\t";
      fs << endl;
    }
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/Anative/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++)
    fs << offspring_f1_Anat[inmo][i] << endl;
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/Bnative/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++)
    fs << offspring_f1_Bnat[inmo][i] << endl;
    fs.close();
    bas.open_ofstream(fs, "Results/offspring_f1/Cnative/"+bas.inttostring(seed)+"_d"+bas.inttostring(dista[inmo])+".txt");
    for (i= 0; i < nuredes; i++)
    fs << offspring_f1_Cnat[inmo][i] << endl;
    fs.close();
  }
  
  delete [] gwife1;
  delete [] gwife2;
  
  for (inmo=0; inmo < 5; inmo++) {
    for (i=0; i<6; i++) {
      delete [] wife1[inmo][i];
      delete [] wife2[inmo][i];
      delete [] midparent1[inmo][i];
      delete [] midparent2[inmo][i];
      delete [] offspring_f2_mean[inmo][i];
      delete [] offspring_f2_median[inmo][i];
      delete [] offspring_f2_max[inmo][i];
      delete [] offspring_f2_min[inmo][i];
      delete [] offspring_f2_sd[inmo][i];
      delete [] offspring_f2_1q[inmo][i];
      delete [] offspring_f2_3q[inmo][i];
      delete [] offspring_f1_mean[inmo][i];
      delete [] offspring_f1_median[inmo][i];
      delete [] offspring_f1_max[inmo][i];
      delete [] offspring_f1_min[inmo][i];
      delete [] offspring_f1_sd[inmo][i];
      delete [] offspring_f1_1q[inmo][i];
      delete [] offspring_f1_3q[inmo][i];
    }
    delete [] wife1[inmo];
    delete [] wife2[inmo];
    delete [] midparent1[inmo];
    delete [] midparent2[inmo];
    delete [] offspring_f2_mean[inmo];
    delete [] offspring_f2_median[inmo];
    delete [] offspring_f2_max[inmo];
    delete [] offspring_f2_min[inmo];
    delete [] offspring_f2_sd[inmo];
    delete [] offspring_f2_1q[inmo];
    delete [] offspring_f2_3q[inmo];
    delete [] offspring_f1_mean[inmo];
    delete [] offspring_f1_median[inmo];
    delete [] offspring_f1_max[inmo];
    delete [] offspring_f1_min[inmo];
    delete [] offspring_f1_sd[inmo];
    delete [] offspring_f1_1q[inmo];
    delete [] offspring_f1_3q[inmo];
  }
  
  
  for (i=0; i<6; i++) {
    delete [] dad[i];
    delete [] hijos[i];
  }

  
  delete [] dist_pm1;
  delete [] dist_pm2;
  delete [] offspring_f2_Bnat;
  delete [] offspring_f2_Anat;
  delete [] offspring_f2_Cnat;
  delete [] offspring_f1_Bnat;
  delete [] offspring_f1_Anat;
  delete [] offspring_f1_Cnat;
  delete [] hijosor;
  delete [] dema;
  delete [] wife1;
  delete [] wife2;
  delete [] dad;
  delete [] midparent1;
  delete [] midparent2;
  delete [] offspring_f2_mean;
  delete [] offspring_f2_median;
  delete [] offspring_f2_max;
  delete [] offspring_f2_min;
  delete [] offspring_f2_sd;
  delete [] offspring_f2_1q;
  delete [] offspring_f2_3q;
  delete [] offspring_f1_mean;
  delete [] offspring_f1_median;
  delete [] offspring_f1_max;
  delete [] offspring_f1_min;
  delete [] offspring_f1_sd;
  delete [] offspring_f1_1q;
  delete [] offspring_f1_3q;
  delete [] hijos;
  delete [] ci;
  delete [] phenA;
  delete [] phenB;
  
  jacta.close_rng();
  return 0;
}

void getmomAB(int indice, Alea& est, int tam, GraphI &pap, GraphI &cual, int *conin, int *fa, int *fb) {
  if (indice > 4) {
    cout << "[Error]: Index is not right.\n";
    exit(1);
  }
  int disprev = 0, tries = 0;
  int distancias[5] = {2, 4, 8, 16, 32};
  int ciclos[5] = {1, 2, 4, 8, 16};
  cual.copy(pap);
  int i, mien;
  for (i = 0; i < ciclos[indice]; i++) {
    mien = quitaponAB(est, tam, pap, cual, conin, fa, fb, disprev);
    if (mien < 0) {
      cual.clear();
      cual.copy(pap);
      i = -1;
      disprev = 0;
      tries++;
    }
    else
      disprev = mien;
    if (tries >= (4*tam*tam)) {
      cout << "[Error]: Random walk lead to nowhere.\n";
      exit(1);
    }
  }
  if (cual.distance_from_nw(pap) != distancias[indice]) {
    cout << "[Error]: Distance between parents is not right.\n";
    exit(1);
  }
}

void getmomA(int indice, Alea& est, int tam, GraphI &pap, GraphI &cual, int *conin, int *fa) {
  if (indice > 4) {
    cout << "[Error]: Index is not right.\n";
    exit(1);
  }
  int disprev = 0, tries = 0;
  int distancias[5] = {2, 4, 8, 16, 32};
  int ciclos[5] = {1, 2, 4, 8, 16};
  cual.copy(pap);
  int i, mien;
  for (i = 0; i < ciclos[indice]; i++) {
    mien = quitaponA(est, tam, pap, cual, conin, fa, disprev);
    if (mien < 0) {
      cual.clear();
      cual.copy(pap);
      i = -1;
      disprev = 0;
      tries++;
    }
    else
      disprev = mien;
    if (tries >= (4*tam*tam)) {
      cout << "[Error]: Random walk lead to nowhere.\n";
      exit(1);
    }
  }
  if (cual.distance_from_nw(pap) != distancias[indice]) {
    cout << "[Error]: Distance between parents is not right.\n";
    exit(1);
  }
}



int quitaponAB(Alea& est, int tam, GraphI &pap, GraphI &cual, int *conin, int *fa, int *fb, int disprev){
  int i, j, ii, jj, kk, ov, ll, dis, intentos = 0;
  bool pom, buena = true, buena2;
  
  int *modci;
  modci = new int[tam];
  for (i = 0; i < tam; i++) {
    modci[i] = conin[i];
  }
  
  do {
    do {
      ii = est.randint(0, tam);
      jj = est.randint(0, tam);
    } while (cual.weight(ii, jj) == 0);
    ov = cual.weight(ii, jj);
    cual.force_interaction(ii, jj, 0);
    do {
      kk = est.randint(0, tam);
      ll = est.randint(0, tam);
    } while (cual.weight(kk, ll) != 0);
    pom = est.toss();
    if (pom)
      cual.force_interaction(kk, ll, 1);
    else
      cual.force_interaction(kk, ll, -1);
    buena = true;
    dis = cual.distance_from_nw(pap);
    if (dis != (disprev+2)) {
      buena = false;
      intentos++; //esto!!!!
      cual.force_interaction(ii,jj, ov);//esto!!!!
      cual.force_interaction(kk, ll, 0);//esto!!!!
    }
    else {
      cual.set_as_state(conin);
      cual.find_an_attractor();
      if (cual.attractor_size() != 1) {
        buena = false;
      } else {
        for (i=0; i < tam; i++) {
          if (cual.attractor_element(0, i) != fa[i]) {
            buena = false;
            break;
          }
        }
      }
      cual.clear_attractor();
      if (buena) {
        buena2 = false;
        for (j=0; j<tam; j++) {
          modci[j] *= (-1);
          cual.set_as_state(modci);
          cual.find_an_attractor();
          if (cual.attractor_size() == 1) {
            for (i=0; i < tam; i++) {
              if (cual.attractor_element(0, i) != fb[i]) {
                break;
              }
            }
            if (i== tam) {
              buena2 = true;
              cual.clear_attractor(); //esto!!!!
              modci[j] *= (-1); //esto!!!!
              break;
            }
          }
          cual.clear_attractor();
          modci[j] *= (-1);
        }
      }
      intentos++;
      if ((!buena) || (!buena2)) {
        cual.force_interaction(ii,jj, ov);
        cual.force_interaction(kk, ll, 0);
      }
    }
  }while (((!buena) || (!buena2)) && (intentos < (2*tam*tam)));
  if ((!buena) || (!buena2))
    dis = -1;
  for (j = 0; j < tam; j++) { //esto!!
    if (modci[j] != conin[j]) {
      cout << "ahjijo\n";
      exit(1);
    }
  }
  return dis;
}

int quitaponA(Alea& est, int tam, GraphI &pap, GraphI &cual, int *conin, int *fa, int disprev) {
  int i, ii, jj, kk, ov, ll, dis, intentos = 0;
  bool pom, buena = true;
  
  do {
    do {
      ii = est.randint(0, tam);
      jj = est.randint(0, tam);
    } while (cual.weight(ii, jj) == 0);
    ov = cual.weight(ii, jj);
    cual.force_interaction(ii, jj, 0);
    do {
      kk = est.randint(0, tam);
      ll = est.randint(0, tam);
    } while (cual.weight(kk, ll) != 0);
    pom = est.toss();
    if (pom)
      cual.force_interaction(kk, ll, 1);
    else
      cual.force_interaction(kk, ll, -1);
    buena = true;
    dis = cual.distance_from_nw(pap);
    if (dis != (disprev+2)) {
      buena = false;
    } else {
      cual.set_as_state(conin);
      cual.find_an_attractor();
      if (cual.attractor_size() != 1) {
        buena = false;
      } else {
        for (i=0; i < tam; i++) {
          if (cual.attractor_element(0, i) != fa[i]) {
            buena = false;
            break;
          }
        }
      }
      cual.clear_attractor();
    }
    intentos++;
    if (!buena) {
      cual.force_interaction(ii,jj, ov);
      cual.force_interaction(kk, ll, 0);
    }
  }while ((!buena) && (intentos < (2*tam*tam)));
  if (!buena)
    dis = -1;
  return dis;
}
