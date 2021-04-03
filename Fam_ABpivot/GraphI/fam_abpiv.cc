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
  
  int i, j,l,m;
  
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
  GraphI gwife1(jacta);
  GraphI gwife2(jacta);
  GraphI gdad(jacta);
  GraphI hijo(jacta);
  
  
  int numhijos = pow(2, nodes);
  
  double **wife1, **wife2, **dad, **midparent1, **midparent2;
  
  bas.create_array(wife1, 6, nuredes);
  bas.create_array(wife2, 6, nuredes);
  bas.create_array(midparent1, 6, nuredes);
  bas.create_array(midparent2, 6, nuredes);
  bas.create_array(dad, 6, nuredes);
  
  double *dist_pm1, *dist_pm2;
  dist_pm2 = new double[nuredes];
  dist_pm1 = new double[nuredes];
  
  
  double **offspring_f2_mean, **offspring_f2_median, **offspring_f2_max, **offspring_f2_min, **offspring_f2_sd, **offspring_f2_1q, **offspring_f2_3q;
  bas.create_array(offspring_f2_mean, 6, nuredes);
  bas.create_array(offspring_f2_median, 6, nuredes);
  bas.create_array(offspring_f2_max, 6, nuredes);
  bas.create_array(offspring_f2_min, 6, nuredes);
  bas.create_array(offspring_f2_sd, 6, nuredes);
  bas.create_array(offspring_f2_1q, 6, nuredes);
  bas.create_array(offspring_f2_3q, 6, nuredes);
  
  int *offspring_f2_Bnat, *offspring_f2_Anat, *offspring_f2_Cnat;
  offspring_f2_Anat = new int[nuredes];
  offspring_f2_Bnat = new int[nuredes];
  offspring_f2_Cnat = new int[nuredes];
  
  
  double **offspring_f1_mean, **offspring_f1_median, **offspring_f1_max, **offspring_f1_min, **offspring_f1_sd, **offspring_f1_1q, **offspring_f1_3q;
  bas.create_array(offspring_f1_mean, 6, nuredes);
  bas.create_array(offspring_f1_median, 6, nuredes);
  bas.create_array(offspring_f1_max, 6, nuredes);
  bas.create_array(offspring_f1_min, 6, nuredes);
  bas.create_array(offspring_f1_sd, 6, nuredes);
  bas.create_array(offspring_f1_1q, 6, nuredes);
  bas.create_array(offspring_f1_3q, 6, nuredes);
  
  int *offspring_f1_Bnat, *offspring_f1_Anat, *offspring_f1_Cnat;
  offspring_f1_Anat = new int[nuredes];
  offspring_f1_Bnat = new int[nuredes];
  offspring_f1_Cnat = new int[nuredes];
  
  
  double **hijos;
  bas.create_array(hijos, 6, numhijos);
  
  int temA, temB;
  ifstream fepap, femAB, femA;
  
  double *hijosor;
  hijosor = new double[numhijos];
  int *atratemp;
  atratemp = new int[nodes];
  
  bool *dema;
  dema = new bool[nodes];
  
  int cualhij;
  bas.open_ifstream(fepap, "networks/Male/"+bas.inttostring(seed)+".txt");
  bas.open_ifstream(femAB, "networks/Wife1/"+bas.inttostring(seed)+".txt");
  bas.open_ifstream(femA, "networks/Wife2/"+bas.inttostring(seed)+".txt");
  for (j=0; j < nuredes; j++) {
    gdad.get_dir_nw_from_file(nodes, fepap, edges);
    fepap >> tras;
    gwife2.get_dir_nw_from_file(nodes, femA, edges);
    femA >> tras;
    gwife1.get_dir_nw_from_file(nodes, femAB, edges);
    femAB >> tras;
    
    gdad.penetrance_n1_2p(ci, phenA, phenB, temA, temB);
    dad[0][j] = temB;
    dad[3][j] = temA;
    gdad.penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
    dad[1][j] = temB;
    dad[4][j] = temA;
    gdad.mutational_access_2p(ci, phenA, phenB, temA, temB);
    dad[2][j] = temB/((2.0*nodes*nodes)-edges);
    dad[5][j] = temA/((2.0*nodes*nodes)-edges);
    
    gwife2.penetrance_n1_2p(ci, phenA, phenB, temA, temB);
    wife2[0][j] = temB;
    wife2[3][j] = temA;
    gwife2.penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
    wife2[1][j] = temB;
    wife2[4][j] = temA;
    gwife2.mutational_access_2p(ci, phenA, phenB, temA, temB);
    wife2[2][j] = temB/((2.0*nodes*nodes)-edges);
    wife2[5][j] = temA/((2.0*nodes*nodes)-edges);
    dist_pm2[j] = gdad.distance_from_nw(gwife2);
    
    gwife1.penetrance_n1_2p(ci, phenA, phenB, temA, temB);
    wife1[0][j] = temB;
    wife1[3][j] = temA;
    gwife1.penetrance_n2_ex_2p(ci, phenA, phenB, temA, temB);
    wife1[1][j] = temB;
    wife1[4][j] = temA;
    gwife1.mutational_access_2p(ci, phenA, phenB, temA, temB);
    wife1[2][j] = temB/((2.0*nodes*nodes)-edges);
    wife1[5][j] = temA/((2.0*nodes*nodes)-edges);
    dist_pm1[j] = gdad.distance_from_nw(gwife1);
    
    for (l=0; l < 6; l++) {
      midparent2[l][j] = (dad[l][j]+wife2[l][j])/2.0;
      midparent1[l][j] = (dad[l][j]+wife1[l][j])/2.0;
    }
    //acabe con padres
    bas.fillv0(dema, nodes);
    cualhij = 0;
    offspring_f2_Anat[j] = 0;
    offspring_f2_Bnat[j] = 0;
    offspring_f2_Cnat[j] = 0;
    //hijos maA
    while(cualhij < numhijos){//}!bas.last_vector(dema, nodes)) {
      hijo.build_offspring(gwife2, gdad, dema);
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
          offspring_f2_Anat[j]++;
        else {
          if (bas.eqvec(atratemp, nodes, phenB, nodes))
            offspring_f2_Bnat[j]++;
          else
            offspring_f2_Cnat[j]++;
        }
      } else
        offspring_f2_Cnat[j]++;
      hijo.clear();
      cualhij++;
      if (cualhij < numhijos)
        bas.next_in_vector(dema, nodes);
    }
    for (l=0; l < 6; l++) {
      offspring_f2_mean[l][j] = bas.get_mean(hijos[l], numhijos);
      offspring_f2_sd[l][j] = bas.get_pop_stddev(hijos[l], numhijos, offspring_f2_mean[l][j]);
      bas.sort(hijos[l], hijosor, numhijos);
      offspring_f2_min[l][j] = hijosor[0];
      offspring_f2_max[l][j] = hijosor[numhijos-1];
      offspring_f2_median[l][j] = bas.get_midpoint(hijosor, numhijos);
      offspring_f2_1q[l][j] = bas.get_q1(hijosor, numhijos);
      offspring_f2_3q[l][j] = bas.get_q3(hijosor, numhijos);
    }
    
    
    //hijos mAB
    bas.fillv0(dema, nodes);
    cualhij = 0;
    offspring_f1_Anat[j] = 0;
    offspring_f1_Bnat[j] = 0;
    offspring_f1_Cnat[j] = 0;
    while(cualhij < numhijos) {//}!bas.last_vector(dema, nodes)) {
      hijo.build_offspring(gwife1, gdad, dema);
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
          offspring_f1_Anat[j]++;
        else {
          if (bas.eqvec(atratemp, nodes, phenB, nodes))
            offspring_f1_Bnat[j]++;
          else
            offspring_f1_Cnat[j]++;
        }
      } else
        offspring_f1_Cnat[j]++;
      hijo.clear();
      cualhij++;
      if (cualhij < numhijos)
        bas.next_in_vector(dema, nodes);
    }
    for (l=0; l < 6; l++) {
      offspring_f1_mean[l][j] = bas.get_mean(hijos[l], numhijos);
      offspring_f1_sd[l][j] = bas.get_pop_stddev(hijos[l], numhijos, offspring_f1_mean[l][j]);
      bas.sort(hijos[l], hijosor, numhijos);
      offspring_f1_min[l][j] = hijosor[0];
      offspring_f1_max[l][j] = hijosor[numhijos-1];
      offspring_f1_median[l][j] = bas.get_midpoint(hijosor, numhijos);
      offspring_f1_1q[l][j] = bas.get_q1(hijosor, numhijos);
      offspring_f1_3q[l][j] = bas.get_q3(hijosor, numhijos);
    }
    gwife1.clear();
    gwife2.clear();
    gdad.clear();
    
  }
  fepap.close();
  femA.close();
  femAB.close();
  
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
  bas.open_ofstream(fs, "Results/Wife2/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << wife2[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/Wife1/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << wife1[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  
  
  bas.open_ofstream(fs, "Results/Midparent2/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << midparent2[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/Midparent1/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << midparent1[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  
  
  
  
  bas.open_ofstream(fs, "Results/Dist_dad_wife2/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++)
    fs << dist_pm2[i] << endl;
  fs.close();
  bas.open_ofstream(fs, "Results/Dist_dad_wife1/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++)
    fs << dist_pm1[i] << endl;
  fs.close();
  
  
  bas.open_ofstream(fs, "Results/offspring_f2/mean/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f2_mean[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/median/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f2_median[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/max/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f2_max[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/min/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f2_min[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/sd/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f2_sd[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/1q/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f2_1q[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/3q/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f2_3q[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/Anative/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++)
    fs << offspring_f2_Anat[i] << endl;
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/Bnative/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++)
    fs << offspring_f2_Bnat[i] << endl;
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f2/Cnative/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++)
    fs << offspring_f2_Cnat[i] << endl;
  fs.close();//
  //
  //HIJOSab:
  bas.open_ofstream(fs, "Results/offspring_f1/mean/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f1_mean[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/median/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f1_median[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/max/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f1_max[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/min/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f1_min[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/sd/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f1_sd[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/1q/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f1_1q[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/3q/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++) {
    for (j=0; j < 6; j++)
      fs << offspring_f1_3q[j][i] << "\t";
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/Anative/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++)
    fs << offspring_f1_Anat[i] << endl;
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/Bnative/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++)
    fs << offspring_f1_Bnat[i] << endl;
  fs.close();
  bas.open_ofstream(fs, "Results/offspring_f1/Cnative/"+bas.inttostring(seed)+".txt");
  for (i= 0; i < nuredes; i++)
    fs << offspring_f1_Cnat[i] << endl;
  fs.close();
    
  for (i=0; i<6; i++) {
    delete [] wife1[i];
    delete [] wife2[i];
    delete [] dad[i];
    delete [] midparent1[i];
    delete [] midparent2[i];
    delete [] offspring_f2_mean[i];
    delete [] offspring_f2_median[i];
    delete [] offspring_f2_max[i];
    delete [] offspring_f2_min[i];
    delete [] offspring_f2_sd[i];
    delete [] offspring_f2_1q[i];
    delete [] offspring_f2_3q[i];
    delete [] offspring_f1_mean[i];
    delete [] offspring_f1_median[i];
    delete [] offspring_f1_max[i];
    delete [] offspring_f1_min[i];
    delete [] offspring_f1_sd[i];
    delete [] offspring_f1_1q[i];
    delete [] offspring_f1_3q[i];
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


