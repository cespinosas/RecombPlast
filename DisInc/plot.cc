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

void impomat(ifstream& fli, string archi, double **mat, int colfil, int rofil);
void impovec(ifstream& fli, string archi, int *vec, int tam);
void impovec(ifstream& fli, string archi, double *vec, int tam);

void wilcoxon(ofstream& flusal, string archstat, int colx, int coly, char tgol, string descx, string descy);//primera que que segunda


Alea jacta(1);
Basics bas(jacta);

int main(int argc, char **argv) {
  int cuantos = atoi(argv[1]);
  int nodes = atoi(argv[2]);
  ifstream fe;
  ofstream fs;
  if ((cuantos%20)!=0) {
    cout << "[Error]: Wrong number of files.\n";
    exit(1);
  }
  int nufiles = cuantos/20;
  int i, j, k, l;
  //poner cada cosa en un solo archivo
  string losdir[27];
  losdir[0] = "Results/Male/";
  losdir[1] = "Results/Wife1/";
  losdir[2] = "Results/Wife2/";
  losdir[3] = "Results/Dist_dad_wife1/";
  losdir[4] = "Results/Dist_dad_wife2/";
  losdir[5] = "Results/offspring_f1/mean/";
  losdir[6] = "Results/offspring_f1/median/";
  losdir[7] = "Results/offspring_f1/max/";
  losdir[8] = "Results/offspring_f1/min/";
  losdir[9] = "Results/offspring_f1/sd/";
  losdir[10] = "Results/offspring_f1/1q/";
  losdir[11] = "Results/offspring_f1/3q/";
  losdir[12] = "Results/offspring_f1/Anative/";
  losdir[13] = "Results/offspring_f1/Bnative/";
  losdir[14] = "Results/offspring_f1/Cnative/";
  losdir[15] = "Results/offspring_f2/mean/";
  losdir[16] = "Results/offspring_f2/median/";
  losdir[17] = "Results/offspring_f2/max/";
  losdir[18] = "Results/offspring_f2/min/";
  losdir[19] = "Results/offspring_f2/sd/";
  losdir[20] = "Results/offspring_f2/1q/";
  losdir[21] = "Results/offspring_f2/3q/";
  losdir[22] = "Results/offspring_f2/Anative/";
  losdir[23] = "Results/offspring_f2/Bnative/";
  losdir[24] = "Results/offspring_f2/Cnative/";
  losdir[25] = "Results/Midparent1/";
  losdir[26] = "Results/Midparent2/";
  string orden;
  int inmo;
  int dista[5] ={2,4,8,16,32};
  inmo = 0;
  i = 3;
  for (inmo = 0; inmo < 5; inmo++) {
    for (i = 3; i < 27; i++) {
        orden = "cat "+losdir[i]+"1"+"_d"+bas.inttostring(dista[inmo])+".txt "+losdir[i]+"2"+"_d"+bas.inttostring(dista[inmo])+".txt > "+losdir[i]+"all"+"_d"+bas.inttostring(dista[inmo])+".txt";
        system(orden.c_str());
        for (j=3; j <= nufiles; j++) {
          orden = "cat "+losdir[i]+bas.inttostring(j)+"_d"+bas.inttostring(dista[inmo])+".txt >> "+losdir[i]+"all_d"+bas.inttostring(dista[inmo])+".txt";
          system(orden.c_str());
        }
    }
  }


  double ****off_f1, ****off_f2;
  off_f1 = new double***[5];
  off_f2 = new double***[5];

  for (inmo = 0; inmo < 5; inmo++) {
    //orden es: mean, median, max, min, sd, 1q, 3q,
    bas.create_array(off_f1[inmo], 2, 6, cuantos);
    bas.create_array(off_f2[inmo], 2, 6, cuantos);
  }
  int ***native_f1, ***native_f2;
  //orden es A,B,C
  bas.create_array(native_f1, 5,2, cuantos);
  bas.create_array(native_f2, 5,2, cuantos);
  
  
  
  for (inmo=0; inmo < 5; inmo++) {
    for (i=0; i< 2; i++) {//era 7
      impomat(fe, losdir[5+i]+"all_d"+bas.inttostring(dista[inmo])+".txt", off_f1[inmo][i], 6, cuantos);
      impomat(fe, losdir[15+i]+"all_d"+bas.inttostring(dista[inmo])+".txt", off_f2[inmo][i], 6, cuantos);
    }
    for (i=0; i<2; i++) { //era 2
      impovec(fe, losdir[12+i]+"all_d"+bas.inttostring(dista[inmo])+".txt", native_f1[inmo][i], cuantos);
      impovec(fe, losdir[22+i]+"all_d"+bas.inttostring(dista[inmo])+".txt", native_f2[inmo][i], cuantos);
    }
  }

    
//  parte para wilcoxon
    

  orden = "mkdir Figs";
  system(orden.c_str());
  orden = "mkdir Figs/mean";
  system(orden.c_str());
  orden = "mkdir Figs/native";
  system(orden.c_str());
  string quecosa[2] = {"mean", "median"};
  string prop[3] = {"n1", "n2", "mu"};
  string losfen[2] = {"B", "A"}; //para nat, alreves
  double losdiv[3] = {1, (nodes*(nodes-1))/2, 1};
  ofstream fs2, fsd;
  double *vecor;
  vecor = new double[cuantos];
  i=0;
  j=0;
      for (k = 0; k < 3; k++) { //prop
        bas.open_ofstream(fs, "Figs/"+quecosa[i]+"/desc"+losfen[j]+prop[k]+".txt");
        for (inmo=0; inmo < 5; inmo++) {
          fs << "*********Distance " << dista[inmo] << "******************\n";
          fs << "**Offspring 1 (ABxAB):\n";
          fs << "Mean: " << bas.get_mean(off_f1[inmo][i][(j*3)+k], cuantos)/losdiv[k] << endl;
          fs << "SD: " << bas.get_sample_stddev(off_f1[inmo][i][(j*3)+k], cuantos)/losdiv[k] << endl;
          bas.sort(off_f1[inmo][i][(j*3)+k], vecor, cuantos);
          fs << "Median: " << bas.get_midpoint(vecor, cuantos)/losdiv[k] << endl;
          fs << "Min: " << vecor[0]/losdiv[k] << endl;
          fs << "Q1: " << bas.get_q1(vecor, cuantos)/losdiv[k] << endl;
          fs << "Q3: " << bas.get_q3(vecor, cuantos)/losdiv[k] << endl;
          fs << "Max: " << vecor[cuantos-1]/losdiv[k] << endl;
          fs << endl;
          fs << "Offspring 2 (ABxA):\n";
          fs << "Mean: " << bas.get_mean(off_f2[inmo][i][(j*3)+k], cuantos)/losdiv[k] << endl;
          fs << "SD: " << bas.get_sample_stddev(off_f2[inmo][i][(j*3)+k], cuantos)/losdiv[k] << endl;
          bas.sort(off_f2[inmo][i][(j*3)+k], vecor, cuantos);
          fs << "Median: " << bas.get_midpoint(vecor, cuantos)/losdiv[k] << endl;
          fs << "Min: " << vecor[0]/losdiv[k] << endl;
          fs << "Q1: " << bas.get_q1(vecor, cuantos)/losdiv[k] << endl;
          fs << "Q3: " << bas.get_q3(vecor, cuantos)/losdiv[k] << endl;
          fs << "Max: " << vecor[cuantos-1]/losdiv[k] << endl;
        }
        fs.close();
      }
  ofstream fst;
  int *vecori; //native_f1[inmo][j]
  vecori = new int[cuantos];
  j = 0;
    bas.open_ofstream(fs, "Figs/native/desc"+losfen[j]+"nat.txt");
    for (inmo=0; inmo < 5; inmo++) {
      fs << "*********Distance " << dista[inmo] << "******************\n";
      fs << "**native 1 (ABxAB):\n";
      fs << "Mean: " << bas.get_mean(native_f1[inmo][1-j], cuantos) << endl;
      fs << "SD: " << bas.get_sample_stddev(native_f1[inmo][1-j], cuantos) << endl;
      bas.sort(native_f1[inmo][1-j], vecori, cuantos);
      fs << "Median: " << bas.get_midpoint(vecori, cuantos) << endl;
      fs << "Min: " << vecori[0] << endl;
      fs << "Q1: " << bas.get_q1(vecori, cuantos) << endl;
      fs << "Q3: " << bas.get_q3(vecori, cuantos) << endl;
      fs << "Max: " << vecori[cuantos-1] << endl;
      fs << endl;
      fs << "native 2 (ABxA):\n";
      fs << "Mean: " << bas.get_mean(native_f2[inmo][1-j], cuantos) << endl;
      fs << "SD: " << bas.get_sample_stddev(native_f2[inmo][1-j], cuantos) << endl;
      bas.sort(native_f2[inmo][1-j], vecori, cuantos);
      fs << "Median: " << bas.get_midpoint(vecori, cuantos) << endl;
      fs << "Min: " << vecori[0]/losdiv[k] << endl;
      fs << "Q1: " << bas.get_q1(vecori, cuantos) << endl;
      fs << "Q3: " << bas.get_q3(vecori, cuantos) << endl;
      fs << "Max: " << vecori[cuantos-1] << endl;
    }
    fs.close();

  string palab[3] = {"plasticity", "plasticity", "mutation"};
  i=0;
    bas.open_ofstream(fs2, "Figs/"+quecosa[i]+"/paR.sh");
    fs2 << "library(ggplot2)\n";
  j = 0;
      for (k = 0; k < 3; k++) { //prop
        for (inmo=0; inmo < 5; inmo++) {
          bas.open_ofstream(fs, "Figs/"+quecosa[i]+"/d"+bas.inttostring(dista[inmo])+"_"+losfen[j]+prop[k]+".txt");
          for (l=0; l < cuantos; l++) {
            fs << off_f1[inmo][i][(j*3)+k][l] << "\t" << off_f2[inmo][i][(j*3)+k][l] << endl;
          }
          fs.close();

          fs2 << "a <- read.csv(file = \"" << "Figs/"+quecosa[i]+"/d"+bas.inttostring(dista[inmo])+"_"+losfen[j]+prop[k]+".txt\", head=FALSE, sep =\"\\t\")\n";
          wilcoxon(fs2, "Figs/"+quecosa[i]+"/stats", 1, 2, 'g', "d"+bas.inttostring(dista[inmo])+"_"+losfen[j]+prop[k]+"F1eq", "d"+bas.inttostring(dista[inmo])+"_"+losfen[j]+prop[k]+"F2df");

        }
      }
    fs2.close();
    orden = "Rscript Figs/"+quecosa[i]+"/paR.sh";
    system(orden.c_str());
     
  
  //native
  bas.open_ofstream(fs2, "Figs/native/paR.sh");
  fs2 << "library(ggplot2)\n";
  j = 0;
      for (inmo=0; inmo < 5; inmo++) {
        bas.open_ofstream(fs, "Figs/native/d"+bas.inttostring(dista[inmo])+"_"+losfen[j]+".txt");
        for (l=0; l < cuantos; l++) {
          fs << native_f1[inmo][1-j][l] << "\t" << native_f2[inmo][1-j][l] << endl;
        }
        fs.close();

        fs2 << "a <- read.csv(file = \"" << "Figs/native/d"+bas.inttostring(dista[inmo])+"_"+losfen[j]+".txt\", head=FALSE, sep =\"\\t\")\n";
        wilcoxon(fs2, "Figs/native/stats", 1, 2, 'g', "d"+bas.inttostring(dista[inmo])+"_"+losfen[j]+"F1eq", "d"+bas.inttostring(dista[inmo])+"_"+losfen[j]+"F2df");

      }
  fs2.close();
  orden = "Rscript Figs/native/paR.sh";
  system(orden.c_str());
     
  jacta.close_rng();
  return 0;
}


void impomat(ifstream& fli, string archi, double **mat, int colfil, int rofil) {
  int i, j;
  fli.open(archi.c_str());
  for (i = 0; i < rofil; i++)
    for (j = 0; j < colfil; j++)
      fli >> mat[j][i];
  fli.close();
  return;
  
}

void impovec(ifstream& fli, string archi, int *vec, int tam) {
  int i;
  fli.open(archi.c_str());
  for (i = 0; i < tam; i++)
    fli >> vec[i];
  fli.close();
  return;
}

void impovec(ifstream& fli, string archi, double *vec, int tam) {
  int i;
  fli.open(archi.c_str());
  for (i = 0; i < tam; i++)
    fli >> vec[i];
  fli.close();
  return;
}


void wilcoxon(ofstream& flusal, string archstat, int colx, int coly, char tgol, string descx, string descy) { //primera que que segunda
    string alti;
  if ((tgol == 't') || (tgol =='T'))
    alti = "two.sided";
  else {
    if ((tgol =='G') ||(tgol=='g'))
      alti = "greater";
    else {
      if ((tgol =='L') ||(tgol=='l'))
        alti = "less";
    }
  }
  flusal << "wtest <- wilcox.test(a$V" << colx << ",a$V" << coly << ", paired=TRUE, alternative=\"" << alti << "\")\n";
  flusal << "letr <- \"Wilcoxon: " << descx << " " << alti << " than " << descy << "\"\n";
  flusal << "capture.output(letr, file=\"" << archstat << ".txt\", append = TRUE)\n";
  flusal << "capture.output(wtest, file=\"" << archstat << ".txt\", append = TRUE)\n";
  flusal << "letr <- \"****************************************\"\n";
  flusal << "capture.output(letr, file=\"" << archstat << ".txt\", append = TRUE)\n";
  flusal << endl;
  return;
}
