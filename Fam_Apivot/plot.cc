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

void uu(string dirarr, string quefamF1oF2, double **midp, double *dis, double **pa, double **ma, double ***off, int **nat, int nudats, int nodes);
void uudu(string dirarr, string queparoff, int fenpb0a1, int fenob0a1, double **midp, double **pa, double **ma, double ***off, int **nat, int nudats, int nodes);
void uudutu(string dirarr, string quep, int fenpb0a1, int fenob0a1, double **queparent, double ***off, int **nat, int nudats, int nodes, double incre);
void uudutucu(string dirarr, string mmmms, int fenpb0a1, int fenob0a1, double **queparent, double **off, int nudats, int nodes, double incre);
void uudutucd(string dirarr, int fenpb0a1, int fenob0a1, double **queparent, int **nat, int nudats, int nodes, double incre);
void ud(double ***off1, double ***off2, double **wif1, double **wif2, int **nat1, int **nat2, int nudats, int nodes, double *dis1, double *dis2);
void uddd(string dirarr, double *dis1, double *dis2, int nudats);
void uddu(int qfen0B_1A, string dirarr, string descdat, double ***off1, double ***off2, double **wif1, double **wif2, int **nat1, int **nat2, int nudats, int nodes);
void uddutu(int qfen0B_1A, string mmmms, string dirarr, double **off1, double **off2, double **wife1, double **wife2, int nudats, int nodes);
void uddutd(int qfen0B_1A, string dirarr, double **wife1, double **wife2, int **nat1, int **nat2, int nudats, int nodes);
void impomat(ifstream& fli, string archi, double **mat, int colfil, int rofil);
void impovec(ifstream& fli, string archi, int *vec, int tam);
void impovec(ifstream& fli, string archi, double *vec, int tam);
void wilcoxon(ofstream& flusal, string archstat, int colx, int coly, char tgol, string descx, string descy);
void skahis_idl(ofstream& flusal, string archpdf, int colx, int divx, int coly, int divy, string labx, string laby, double alfa);
void cajaplo2(ofstream& flusal, string archpdf, int colx, int divx, int coly, int divy, string labx, string laby, double umbral);

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
  int i, j;
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
  for (i = 0; i < 27; i++) {
    orden = "cat "+losdir[i]+"1.txt "+losdir[i]+"2.txt > "+losdir[i]+"all.txt";
    system(orden.c_str());
    for (j=3; j <= nufiles; j++) {
      orden = "cat "+losdir[i]+bas.inttostring(j)+".txt >> "+losdir[i]+"all.txt";
      system(orden.c_str());
    }
  }
  double **male, **wife1, **wife2, **midp1, **midp2;
  bas.create_array(male, 6, cuantos);
  bas.create_array(wife1, 6, cuantos);
  bas.create_array(wife2, 6, cuantos);
  bas.create_array(midp1, 6, cuantos);
  bas.create_array(midp2, 6, cuantos);
  double *dist1, *dist2;
  dist1 = new double[cuantos];
  dist2 = new double[cuantos];
  double ***off_f1, ***off_f2;
  //orden es: mean, median, max, min, sd, 1q, 3q,
  bas.create_array(off_f1, 7, 6, cuantos);
  bas.create_array(off_f2, 7, 6, cuantos);
  int **native_f1, **native_f2;
  //orden es A,B,C
  bas.create_array(native_f1, 3, cuantos);
  bas.create_array(native_f2, 3, cuantos);
  
  impomat(fe, losdir[0]+"all.txt", male, 6, cuantos);
  impomat(fe, losdir[1]+"all.txt", wife1, 6, cuantos);
  impomat(fe, losdir[2]+"all.txt", wife2, 6, cuantos);
  impovec(fe, losdir[3]+"all.txt", dist1, cuantos);
  impovec(fe, losdir[4]+"all.txt", dist2, cuantos);
  for (i=0; i< 7; i++) {
    impomat(fe, losdir[5+i]+"all.txt", off_f1[i], 6, cuantos);
    impomat(fe, losdir[15+i]+"all.txt", off_f2[i], 6, cuantos);
  }
  for (i=0; i<3; i++) {
    impovec(fe, losdir[12+i]+"all.txt", native_f1[i], cuantos);
    impovec(fe, losdir[22+i]+"all.txt", native_f2[i], cuantos);
  }
  impomat(fe, losdir[25]+"all.txt", midp1, 6, cuantos);
  impomat(fe, losdir[26]+"all.txt", midp2, 6, cuantos);
  
  uu("Figs/", "F1/", midp1, dist1, male, wife1, off_f1, native_f1, cuantos, nodes);
  uu("Figs/", "F2/", midp2, dist2, male, wife2, off_f2, native_f2, cuantos, nodes);
  ud(off_f1, off_f2, wife1, wife2, native_f1, native_f2, cuantos, nodes, dist1, dist2);
     
     
  jacta.close_rng();
  return 0;
}

void uu(string dirarr, string quefamF1oF2, double **midp, double *dis, double **pa, double **ma, double ***off, int **nat, int nudats, int nodes) {//dirarr quefam/
  ofstream flusal;
  string orden;
  string eldir = dirarr+"BewteenGens"+quefamF1oF2;
  orden = "mkdir "+eldir;
  system(orden.c_str());
  uudu(eldir, "ParB_offB/", 0, 0, midp, pa, ma, off, nat, nudats, nodes);
  return;
}

void uudu(string dirarr, string queparoff, int fenpb0a1, int fenob0a1, double **midp, double **pa, double **ma, double ***off, int **nat, int nudats, int nodes) {
  ofstream flusal;
  string orden;
  string eldir = dirarr+queparoff; /////
  orden = "mkdir "+eldir;
  system(orden.c_str());
  
  uudutu(eldir, "midp/", fenpb0a1, fenob0a1, midp, off, nat, nudats, nodes, 0.5);
  return;
}

void uudutu(string dirarr, string quep, int fenpb0a1, int fenob0a1, double **queparent, double ***off, int **nat, int nudats, int nodes, double incre) {
  int i = 0;
  ofstream flusal;
  string orden;
  string eldir = dirarr+quep; /////
  orden = "mkdir "+eldir;
  system(orden.c_str());
  string dirmmmms[5] = {"mean/", "median/", "max/", "min/", "sd/"};
  i = 0;
  uudutucu(eldir, dirmmmms[i], fenpb0a1, fenob0a1, queparent, off[i], nudats, nodes, incre);
  uudutucd(eldir, fenpb0a1, fenob0a1, queparent, nat, nudats, nodes, incre);
  return;
}

void uudutucu(string dirarr, string mmmms, int fenpb0a1, int fenob0a1, double **queparent, double **off, int nudats, int nodes, double incre) {
  int i;
  ofstream flusal;
  string orden;
  string eldir = dirarr+mmmms; /////
  orden = "mkdir "+eldir;
  system(orden.c_str());
  string fenpar = "B";
  if (fenpb0a1==1)
    fenpar = "A";
  string fenoff = "B";
  if (fenob0a1==1)
    fenoff = "A";
  
  double losminpar[3];
  losminpar[1] = nodes-1;
  double *vecor;
  vecor = new double[nudats];
  string traits[3] = {"n1", "n2", "mu"};
  double losdiv[3] = {1, (nodes*(nodes-1))/2, 1};
  orden = eldir+"description.txt";
  flusal.open(orden.c_str());
  int j,k,lmin, lmax;
  for (i=0; i < 3; i++) {
    flusal << "**" << fenpar << " " << traits[i] << "in parent:\n";
    flusal << "Mean: " << bas.get_mean(queparent[(fenpb0a1*3)+i], nudats)/losdiv[i] << endl;
    flusal << "Std. Dev.: " << bas.get_sample_stddev(queparent[(fenpb0a1*3)+i], nudats)/losdiv[i] << endl;
    bas.sort(queparent[(fenpb0a1*3)+i], vecor, nudats);
    flusal << "Median: " << bas.get_midpoint(vecor, nudats)/losdiv[i] << endl;
    flusal << "Minimum: " << vecor[0]/losdiv[i] << endl;
    flusal << "Q1: " << bas.get_q1(vecor, nudats)/losdiv[i] << endl;
    flusal << "Q3: " << bas.get_q3(vecor, nudats)/losdiv[i] << endl;
    flusal << "Maximum: " << vecor[nudats-1]/losdiv[i] << endl << endl;
    if (i != 1)
      losminpar[i] = vecor[0];
    flusal << "**" << fenoff << " " << traits[i] << "in offspring:\n";
    flusal << "Mean: " << bas.get_mean(off[(fenob0a1*3)+i], nudats)/losdiv[i] << endl;
    flusal << "Std. Dev.: " << bas.get_sample_stddev(off[(fenob0a1*3)+i], nudats)/losdiv[i] << endl;
    bas.sort(off[(fenob0a1*3)+i], vecor, nudats);
    flusal << "Median: " << bas.get_midpoint(vecor, nudats)/losdiv[i] << endl;
    flusal << "Minimum: " << vecor[0]/losdiv[i] << endl;
    flusal << "Q1: " << bas.get_q1(vecor, nudats)/losdiv[i] << endl;
    flusal << "Q3: " << bas.get_q3(vecor, nudats)/losdiv[i] << endl;
    flusal << "Maximum: " << vecor[nudats-1]/losdiv[i] << endl << endl;
    flusal << endl;
  }
  flusal.close();
  delete [] vecor;
  
  if (mmmms == "mean/") {
    int contminpar[3];
    bas.fillv0(contminpar, 3);
    
    bas.open_ofstream(flusal, eldir+"data_minfact.txt");
    for (i=0; i < nudats; i++) {
      for (j=0; j < 3; j++)
      flusal << queparent[(fenpb0a1*3)+j][i] << "\t";
      for (j=0; j < 3; j++)
      flusal << off[(fenob0a1*3)+j][i] << "\t";
      for (j=0; j < 3; j++) {
        if (queparent[(fenpb0a1*3)+j][i] <= losminpar[j]) {
          flusal << 0 << "\t";
          contminpar[j]++;
        }
        else
          flusal << 1 << "\t";
      }
      flusal << endl;
    }
    flusal.close();
    
    bas.open_ofstream(flusal, eldir+"count_lowestclass.txt");
    flusal << "Number of data points in lowest class (n2 <= 2/nodes):\n";
    for (j=0; j < 3; j++) {
      flusal << "Parental " << fenpar << " " << traits[j] << ": " << contminpar[j] << endl;
    }
    flusal.close();
    ofstream flusald;
    double ****pahimm;
    pahimm = new double***[3];
    for (i=0; i < 3; i++) {
      pahimm[i] = new double**[3];
      for (j=0; j < 3; j++) {
        pahimm[i][j] = new double*[2];
        pahimm[i][j][0] = new double[contminpar[i]];
        pahimm[i][j][1] = new double[nudats-contminpar[i]];
      }
    }
    
    
    for (i=0; i < 3; i++) {
      bas.open_ofstream(flusal, eldir+"Min_par"+fenpar+"_"+traits[i]+"_off"+fenoff+"n1n2mu.txt");
      bas.open_ofstream(flusald, eldir+"Rest_par"+fenpar+"_"+traits[i]+"_off"+fenoff+"n1n2mu.txt");
      lmin = 0;
      lmax = 0;
      for (j=0; j < nudats; j++) {
        if (queparent[(fenpb0a1*3)+i][j] <= losminpar[i]) {
          for (k=0; k < 3; k++) {
            flusal << off[(fenob0a1*3)+k][j] << "\t";
            pahimm[i][k][0][lmin] = off[(fenob0a1*3)+k][j];
          }
          lmin++;
          flusal << endl;
        }
        else {
          for (k=0; k < 3; k++) {
            flusald << off[(fenob0a1*3)+k][j] << "\t";
            pahimm[i][k][1][lmax] = off[(fenob0a1*3)+k][j];
          }
          lmax++;
          flusald << endl;
        }
      }
      flusal.close();
      flusald.close();
    }
        
    orden = eldir+"description_bin.txt";
    flusal.open(orden.c_str());
    for (i=0; i < 3; i++) {
      if ((contminpar[i] > 0) && (contminpar[i] < nudats)) {
        vecor = new double[contminpar[i]];
        for (j=0; j< 3; j++) {
          flusal << "****Parents with min value (n2 <= 2/nodes) of " << fenpar << traits[i] << "have offspring with " << fenoff << traits[j] << " value of:\n";
          
          flusal << "Mean: " << bas.get_mean(pahimm[i][j][0], contminpar[i])/losdiv[i] << endl;
          flusal << "Std. Dev.: " << bas.get_sample_stddev(pahimm[i][j][0], contminpar[i])/losdiv[i] << endl;
          bas.sort(pahimm[i][j][0], vecor, contminpar[i]);
          flusal << "Median: " << bas.get_midpoint(vecor, contminpar[i])/losdiv[i] << endl;
          flusal << "Minimum: " << vecor[0]/losdiv[i] << endl;
          flusal << "Q1: " << bas.get_q1(vecor, contminpar[i])/losdiv[i] << endl;
          flusal << "Q3: " << bas.get_q3(vecor, contminpar[i])/losdiv[i] << endl;
          flusal << "Maximum: " << vecor[contminpar[i]-1]/losdiv[i] << endl << endl;
        }
        delete [] vecor;
        vecor = new double[nudats-contminpar[i]];
        for (j=0; j< 3; j++) {
          flusal << "****Parents with greater than min value (n2 > 2/nodes) of " << fenpar << traits[i] << "have offspring with " << fenoff << traits[j] << " value of:\n";
          
          flusal << "Mean: " << bas.get_mean(pahimm[i][j][1], nudats-contminpar[i])/losdiv[i] << endl;
          flusal << "Std. Dev.: " << bas.get_sample_stddev(pahimm[i][j][1], nudats-contminpar[i])/losdiv[i] << endl;
          bas.sort(pahimm[i][j][1], vecor, nudats-contminpar[i]);
          flusal << "Median: " << bas.get_midpoint(vecor, nudats-contminpar[i])/losdiv[i] << endl;
          flusal << "Minimum: " << vecor[0]/losdiv[i] << endl;
          flusal << "Q1: " << bas.get_q1(vecor, nudats-contminpar[i])/losdiv[i] << endl;
          flusal << "Q3: " << bas.get_q3(vecor, nudats-contminpar[i])/losdiv[i] << endl;
          flusal << "Maximum: " << vecor[(nudats-contminpar[i])-1]/losdiv[i] << endl << endl;
        }
        delete [] vecor;
      }
    }
    flusal.close();
    for (i=0; i < 3; i++) {
      for (j=0; j < 3; j++) {
        for (k=0; k<2; k++) {
          delete [] pahimm[i][j][k];
        }
        delete [] pahimm[i][j];
      }
      delete [] pahimm[i];
    }
    delete [] pahimm;
    
    
    string letpar[3] = {"Parental access to \'*italic("+fenpar+")*\' through plasticity", "Parental access to \'*italic("+fenpar+")*\' through plasticity", "Parental access to \'*italic("+fenpar+")*\' through mutation"};
    string letoff[3] = {"Offspring\\\'s access to \'*italic("+fenoff+")*\' through plasticity", "Offspring\\\'s access to \'*italic("+fenoff+")*\' through plasticity", "Offspring\\\'s access to \'*italic("+fenoff+")*\' through mutation"};
    
    bas.open_ofstream(flusal, eldir+"paR2.sh");
    flusal << "library(ggplot2)\n";
    flusal << "a <- read.csv(file = \"" << eldir << "data_minfact" << ".txt\", head=FALSE, sep =\"\\t\")\n";
    for (i=0; i < 3; i++) {
      if ((contminpar[i] > 0) && (contminpar[i] < nudats)) {
        flusal << "a$V" << 7+i << " <- factor(a$V" << 7+i << ")\n";
        for (j=0; j < 3; j++)
        cajaplo2(flusal, eldir+"bp_par"+fenpar+traits[i]+"_off"+fenoff+traits[j], 7+i, losdiv[i], j+4, losdiv[j], letpar[i], letoff[j], losminpar[i]);
      }
    }
    for (i=0; i < 3; i++) {
      if ((contminpar[i] > 0) && (contminpar[i] < nudats)) {
        flusal << "m <- read.csv(file = \"" << eldir+"Min_par"+fenpar+"_"+traits[i]+"_off"+fenoff+"n1n2mu" << ".txt\", head=FALSE, sep =\"\\t\")\n";
        flusal << "r <- read.csv(file = \"" << eldir+"Rest_par"+fenpar+"_"+traits[i]+"_off"+fenoff+"n1n2mu" << ".txt\", head=FALSE, sep =\"\\t\")\n";
        for (j = 0; j < 3; j++) {
          flusal << "wtest <- wilcox.test(m$V" << 1+j << ",r$V" << 1+j << ", paired=FALSE, alternative=\"less\")\n";
          flusal << "letr <- \"Mann-Whitney: Offspring's " << fenoff << traits[j]  << " with minimum parental " << fenpar << traits[i] << " less than with higher values of parental " << fenpar << traits[i] << "\"" << endl;
          flusal << "capture.output(letr, file=\"" << eldir << "stats2.txt\", append = TRUE)\n";
          flusal << "capture.output(wtest, file=\"" << eldir << "stats2.txt\", append = TRUE)\n";
          flusal << "letr <- \"****************************************\"\n";
          flusal << "capture.output(letr, file=\"" << eldir << "stats2.txt\", append = TRUE)\n";
          flusal << endl;
        }
      }
    }
    flusal.close();
    orden = "Rscript "+eldir+"paR2.sh";
    system(orden.c_str());
  }
  return;
}

void uudutucd(string dirarr, int fenpb0a1, int fenob0a1, double **queparent, int **nat, int nudats, int nodes, double incre) {//circir
  int i;
  ofstream flusal;
  string orden;
  string eldir = dirarr+"native/"; /////
  orden = "mkdir "+eldir;
  system(orden.c_str());
  
  double min = 1000, max = 0;
  orden = eldir+"data.txt";
  flusal.open(orden.c_str());
  for (i=0; i < nudats; i++) {//parent(n1, n2, mu) off(nat)
    flusal << queparent[(fenpb0a1*3)][i] << "\t" << queparent[(fenpb0a1*3)+1][i] << "\t" << queparent[(fenpb0a1*3)+2][i] << "\t" << nat[1-fenob0a1][i] << endl;
    if (queparent[(fenpb0a1*3)][i] < min)
      min = queparent[(fenpb0a1*3)][i];
    if (queparent[(fenpb0a1*3)][i] > max)
      max = queparent[(fenpb0a1*3)][i];
  }
  flusal.close();
  string fenpar = "B";
  if (fenpb0a1==1)
    fenpar = "A";
  string fenoff = "B";
  if (fenob0a1==1)
  fenoff = "A";
  
  double losminpar[3];
  losminpar[1] = nodes-1;
  double *vecor;
  vecor = new double[nudats];
  int *vecori;
  vecori = new int[nudats];
  string traits[3] = {"n1", "n2", "mu"};
  double losdiv[3] = {1, (nodes*(nodes-1))/2, 1};
  orden = eldir+"description.txt";
  flusal.open(orden.c_str());
  int j,lmin,lmax;
  for (i=0; i < 3; i++) {
    flusal << "**" << fenpar << " " << traits[i] << "in parent:\n";
    flusal << "Mean: " << bas.get_mean(queparent[(fenpb0a1*3)+i], nudats)/losdiv[i] << endl;
    flusal << "Std. Dev.: " << bas.get_sample_stddev(queparent[(fenpb0a1*3)+i], nudats)/losdiv[i] << endl;
    bas.sort(queparent[(fenpb0a1*3)+i], vecor, nudats);
    flusal << "Median: " << bas.get_midpoint(vecor, nudats)/losdiv[i] << endl;
    flusal << "Minimum: " << vecor[0]/losdiv[i] << endl;
    flusal << "Q1: " << bas.get_q1(vecor, nudats)/losdiv[i] << endl;
    flusal << "Q3: " << bas.get_q3(vecor, nudats)/losdiv[i] << endl;
    flusal << "Maximum: " << vecor[nudats-1]/losdiv[i] << endl << endl;
    if (i != 1)
      losminpar[i] = vecor[0];
  }
  delete [] vecor;
  flusal << "**" << fenoff << " native in offspring:\n";
  flusal << "Mean: " << bas.get_mean(nat[1-fenob0a1], nudats) << endl;
  flusal << "Std. Dev.: " << bas.get_sample_stddev(nat[1-fenob0a1], nudats) << endl;
  bas.sort(nat[1-fenob0a1], vecori, nudats);
  flusal << "Median: " << bas.get_midpoint(vecori, nudats) << endl;
  flusal << "Minimum: " << vecori[0] << endl;
  flusal << "Q1: " << bas.get_q1(vecori, nudats) << endl;
  flusal << "Q3: " << bas.get_q3(vecori, nudats) << endl;
  flusal << "Maximum: " << vecori[nudats-1] << endl << endl;
  flusal << endl;
  flusal.close();
  delete [] vecori;
  
  int contminpar[3];
  bas.fillv0(contminpar, 3);
  
  bas.open_ofstream(flusal, eldir+"data_minfact.txt");
  for (i=0; i < nudats; i++) {
    for (j=0; j < 3; j++) {
      flusal << queparent[(fenpb0a1*3)+j][i] << "\t";
    }
    flusal << nat[1-fenob0a1][i] << "\t";
    for (j=0; j < 3; j++) {
      if (queparent[(fenpb0a1*3)+j][i] <= losminpar[j]) {
        flusal << 0 << "\t";
        contminpar[j]++;
      }
      else
        flusal << 1 << "\t";
    }
    flusal << endl;
  }
  flusal.close();
  bas.open_ofstream(flusal, eldir+"count_lowestclass.txt");
  flusal << "Number of data points in lowest class (n2 <= 2/nodes):\n";
  for (j=0; j < 3; j++) {
    flusal << "Parental " << fenpar << " " << traits[j] << ": " << contminpar[j] << endl;
  }
  flusal.close();
  
  int ***pamm;
  pamm = new int**[3];
  for (i=0; i < 3; i++) {
    pamm[i] = new int*[2];
    pamm[i][0] = new int[contminpar[i]];
    pamm[i][1] = new int[nudats-contminpar[i]];
  }
  
  ofstream flusald;
  for (i=0; i < 3; i++) {
    bas.open_ofstream(flusal, eldir+"Min_par"+fenpar+"_"+traits[i]+"_off"+fenoff+"native.txt");
    bas.open_ofstream(flusald, eldir+"Rest_par"+fenpar+"_"+traits[i]+"_off"+fenoff+"native.txt");
    lmin = 0;
    lmax = 0;
    for (j=0; j < nudats; j++) {
      if (queparent[(fenpb0a1*3)+i][j] <= losminpar[i]) {
        flusal << nat[1-fenob0a1][j] << "\t";
        flusal << endl;
        pamm[i][0][lmin] = nat[1-fenob0a1][j];
        lmin++;
      }
      else {
        flusald << nat[1-fenob0a1][j] << "\t";
        flusald << endl;
        pamm[i][1][lmax] = nat[1-fenob0a1][j];
        lmax++;
      }
    }
    flusal.close();
    flusald.close();
  }
  orden = eldir+"description_bin.txt";
  flusal.open(orden.c_str());
  for (i=0; i < 3; i++) {
    if ((contminpar[i] > 0) && (contminpar[i] < nudats)) {
      vecori = new int[contminpar[i]];
      flusal << "****Parents with min value (n2 <= 2/nodes) of " << fenpar << traits[i] << "have offspring with native " << fenoff << " value of:\n";
      flusal << "Mean: " << bas.get_mean(pamm[i][0], contminpar[i])/losdiv[i] << endl;
      flusal << "Std. Dev.: " << bas.get_sample_stddev(pamm[i][0], contminpar[i])/losdiv[i] << endl;
      bas.sort(pamm[i][0], vecori, contminpar[i]);
      flusal << "Median: " << bas.get_midpoint(vecori, contminpar[i])/losdiv[i] << endl;
      flusal << "Minimum: " << vecori[0]/losdiv[i] << endl;
      flusal << "Q1: " << bas.get_q1(vecori, contminpar[i])/losdiv[i] << endl;
      flusal << "Q3: " << bas.get_q3(vecori, contminpar[i])/losdiv[i] << endl;
      flusal << "Maximum: " << vecori[contminpar[i]-1]/losdiv[i] << endl << endl;
      delete [] vecori;
      
      vecori = new int[nudats-contminpar[i]];
      flusal << "****Parents with greater than min value (n2 > 2/nodes) of " << fenpar << traits[i] << "have offspring with native " << fenoff << " value of:\n";
      flusal << "Mean: " << bas.get_mean(pamm[i][1], nudats-contminpar[i])/losdiv[i] << endl;
      flusal << "Std. Dev.: " << bas.get_sample_stddev(pamm[i][1], nudats-contminpar[i])/losdiv[i] << endl;
      bas.sort(pamm[i][1], vecori, nudats-contminpar[i]);
      flusal << "Median: " << bas.get_midpoint(vecori, nudats-contminpar[i])/losdiv[i] << endl;
      flusal << "Minimum: " << vecori[0]/losdiv[i] << endl;
      flusal << "Q1: " << bas.get_q1(vecori, nudats-contminpar[i])/losdiv[i] << endl;
      flusal << "Q3: " << bas.get_q3(vecori, nudats-contminpar[i])/losdiv[i] << endl;
      flusal << "Maximum: " << vecori[(nudats-contminpar[i])-1]/losdiv[i] << endl << endl;
      delete [] vecori;
    }
  }
  flusal.close();
  for (i=0; i<3; i++) {
    for (j=0; j<2; j++) {
      delete [] pamm[i][j];
    }
    delete [] pamm[i];
  }
  delete [] pamm;
  
  string letpar[3] = {"Parental access to \'*italic("+fenpar+")*\' through plasticity", "Parental access to \'*italic("+fenpar+")*\' through plasticity", "Parental access to \'*italic("+fenpar+")*\' through mutation"};
  
  bas.open_ofstream(flusal, eldir+"paR2.sh");
  flusal << "library(ggplot2)\n";
  flusal << "a <- read.csv(file = \"" << eldir << "data_minfact" << ".txt\", head=FALSE, sep =\"\\t\")\n";
  for (i=0; i < 3; i++) {
    if ((contminpar[i] > 0) && (contminpar[i] < nudats)) {
      flusal << "a$V" << 5+i << " <- factor(a$V" << 5+i << ")\n";
      cajaplo2(flusal, eldir+"bp_par"+fenpar+traits[i]+"_off"+fenoff+"native", 5+i, losdiv[i], 4, 1, letpar[i], "Offspring with native \'*italic("+fenpar+")*\'", losminpar[i]);
    }
  }
  for (i=0; i < 3; i++) {
    if ((contminpar[i] > 0) && (contminpar[i] < nudats)) {
      flusal << "m <- read.csv(file = \"" << eldir+"Min_par"+fenpar+"_"+traits[i]+"_off"+fenoff+"native" << ".txt\", head=FALSE, sep =\"\\t\")\n";
      flusal << "r <- read.csv(file = \"" << eldir+"Rest_par"+fenpar+"_"+traits[i]+"_off"+fenoff+"native" << ".txt\", head=FALSE, sep =\"\\t\")\n";
      flusal << "wtest <- wilcox.test(m$V" << 1 << ",r$V" << 1 << ", paired=FALSE, alternative=\"less\")\n";
      flusal << "letr <- \"Mann-Whitney: Offspring's " << fenoff << " native with minimum parental " << fenpar << traits[i] << " less than with higher values of parental " << fenpar << traits[i] << "\"" << endl;
      flusal << "capture.output(letr, file=\"" << eldir << "stats2.txt\", append = TRUE)\n";
      flusal << "capture.output(wtest, file=\"" << eldir << "stats2.txt\", append = TRUE)\n";
      flusal << "letr <- \"****************************************\"\n";
      flusal << "capture.output(letr, file=\"" << eldir << "stats2.txt\", append = TRUE)\n";
      flusal << endl;
    }
  }
  flusal.close();
  orden = "Rscript "+eldir+"paR2.sh";
  system(orden.c_str());
  return;
}

void ud(double ***off1, double ***off2, double **wif1, double **wif2, int **nat1, int **nat2, int nudats, int nodes, double *dis1, double *dis2) {
  system("mkdir Figs/BetweenFamilies");
  uddu(0, "Figs/BetweenFamilies/", "Alldata", off1, off2, wif1, wif2, nat1, nat2, nudats, nodes);
  uddd("Figs/BetweenFamilies/", dis1, dis2, nudats);

  int i,j,k;
  double ***off1men, ***off2men, **wif1men, **wif2men;
  int **nat1men, **nat2men, nudatsmen=0, vanmen=0;

  for (i=0; i < nudats; i++) {
    if (dis1[i] >= dis2[i])
      nudatsmen++;
  }

  off1men = new double**[5];
  off2men = new double**[5];
  for (i=0; i < 5; i++) {
    off1men[i] = new double*[6];
    off2men[i] = new double*[6];
    for (j=0; j < 6; j++) {
      off1men[i][j] = new double[nudatsmen];
      off2men[i][j] = new double[nudatsmen];
    }
  }
  wif1men = new double*[6];
  wif2men = new double*[6];
  nat1men = new int*[3];
  nat2men = new int*[3];
  for (i = 0; i < 6; i++) {
    wif1men[i] = new double[nudatsmen];
    wif2men[i] = new double[nudatsmen];
  }
  for (i = 0; i < 3; i++) {
    nat1men[i] = new int[nudatsmen];
    nat2men[i] = new int[nudatsmen];
  }


  for (i=0; i < nudats; i++) {
    if (dis1[i] >= dis2[i]) {
      for (j=0; j<6; j++) {
        wif2men[j][vanmen] = wif2[j][i];
        wif1men[j][vanmen] = wif1[j][i];
        for (k=0; k<5; k++) {
          off2men[k][j][vanmen] = off2[k][j][i];
          off1men[k][j][vanmen] = off1[k][j][i];
        }
      }
      for (j=0; j<3; j++) {
        nat2men[j][vanmen] = nat2[j][i];
        nat1men[j][vanmen] = nat1[j][i];
      }
      vanmen++;
    }
  }

  cout << nudatsmen << " " << vanmen << endl;

  uddu(0, "Figs/BetweenFamilies/", "Dist_men", off1men, off2men, wif1men, wif2men, nat1men, nat2men, nudatsmen, nodes);

  for (j=0; j < 6; j++) {
    delete [] wif2men[j];
    delete [] wif1men[j];
  }
  for (j=0; j < 3; j++) {
    delete [] nat2men[j];
    delete [] nat1men[j];
  }
  delete [] nat2men;
  delete [] nat1men;
  delete [] wif2men;
  delete [] wif1men;

  for (i = 0; i < 5; i++) {
    for (j=0; j < 6; j++) {
      delete [] off1men[i][j];
      delete [] off2men[i][j];
    }
    delete [] off1men[i];
    delete [] off2men[i];
  }
  delete [] off1men;
  delete [] off2men;

  return;
}

void uddd(string dirarr, double *dis1, double *dis2, int nudats) {
  Alea cion(1);
  Basics bax(cion);
  ofstream flusal;
  string orden;
  string eldir = dirarr+"dis2_dis1/";
  orden = "mkdir "+eldir;
  system(orden.c_str());
  orden = eldir+"data.txt";
  int i;
  flusal.open(orden.c_str());
  for (i=0; i < nudats; i++)
    flusal << dis2[i] << "\t" << dis1[i] << endl;
  flusal.close();
  
  
  orden = eldir+"paR.sh";
  flusal.open(orden.c_str());
  flusal << "library(ggplot2)\n";
  flusal << "a <- read.csv(file = \""<< eldir << "data" << ".txt\", head=FALSE, sep =\"\\t\")\n";
  wilcoxon(flusal, eldir+"stats", 2, 1, 'g', "dist_fam1_pair", "dist_fam2_pair");
  
  flusal << "pdf(file=\""<< eldir << "his_dis2_dis1.pdf\", width=7.5,height=5.32)\n";
  flusal << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
  flusal << "ggplot(a, aes(x=V2, fill=V1)) + geom_histogram(color=\"black\", alpha=0.5, position=\"identity\", binwidth=2) + scale_fill_grey(labels = c(expression(\'A\'%*%\'AB\'),expression(\'A\'%*%\'A\')))+  labs(x=expression(\'Genetic distance\'), y = expression(\'Network pairs\'))  + theme_classic() + theme(text = element_text(size=14), legend.position = c(.9, .9), legend.title = element_blank())\n";
  flusal << "dev.off()\n";
  
  flusal.close();
  orden = "Rscript "+eldir+"paR.sh";
  system(orden.c_str());
  
  double mientras;
  double *acom;
  acom = new double[nudats];
  
  orden = eldir+"desc_dis_f1_f2.txt";
  flusal.open(orden.c_str());
  flusal << "*Distance between AxAB parents:\n";
  mientras = bax.get_mean(dis2, nudats);
  flusal << "Mean: " << mientras << endl;
  flusal << "Std. Dev: " << bax.get_sample_stddev(dis2, nudats, mientras) << endl;
  bax.sort(dis2, acom, nudats);
  flusal << "Min: " << acom[0] << endl;
  flusal << "Max: " << acom[nudats-1] << endl;
  flusal << "Median: " << bax.get_midpoint(acom, nudats) << endl;
  flusal << "1st quartile: " << bax.get_q1(acom, nudats) << endl;
  flusal << "3rd quartile: " << bax.get_q3(acom, nudats) << endl << endl << endl;
  
  flusal << "*Distance between AxA parents:\n";
  mientras = bax.get_mean(dis1, nudats);
  flusal << "Mean: " << mientras << endl;
  flusal << "Std. Dev: " << bax.get_sample_stddev(dis1, nudats, mientras) << endl;
  bax.sort(dis1, acom, nudats);
  flusal << "Min: " << acom[0] << endl;
  flusal << "Max: " << acom[nudats-1] << endl;
  flusal << "Median: " << bax.get_midpoint(acom, nudats) << endl;
  flusal << "1st quartile: " << bax.get_q1(acom, nudats) << endl;
  flusal << "3rd quartile: " << bax.get_q3(acom, nudats) << endl << endl << endl;
  flusal.close();
  delete [] acom;
  
  
  return;
}

void uddu(int qfen0B_1A, string dirarr, string descdat, double ***off1, double ***off2, double **wif1, double **wif2, int **nat1, int **nat2, int nudats, int nodes) {
  int i = 0;
  string dirmmmms[5] = {"mean/", "median/", "max/", "min/", "sd/"};
  string elfen = "A";
  if (qfen0B_1A==0)
    elfen = "B";
  string eldir = dirarr+descdat+"Phen"+elfen+"/";
  string orden = "mkdir "+eldir;
  system(orden.c_str());
  ofstream flusal;
  bas.open_ofstream(flusal, eldir+"datapoints.txt");
  flusal << nudats << endl;
  flusal.close();
  
  i=0;
    uddutu(qfen0B_1A, dirmmmms[i], eldir, off1[i], off2[i], wif1, wif2, nudats, nodes);
  uddutd(qfen0B_1A, eldir, wif1, wif2, nat1, nat2, nudats, nodes);
  
  return;
}

void uddutu(int qfen0B_1A, string mmmms, string dirarr, double **off1, double **off2, double **wife1, double **wife2, int nudats, int nodes) { //deoff pasar rebcorrecta
  ofstream flusal;
  string orden;
  string eldir = dirarr+mmmms;
  orden = "mkdir "+eldir;
  system(orden.c_str());
  int i, j, indfen, indnat;
  indfen =3*qfen0B_1A;
  indnat = 1-qfen0B_1A;
  string elfen = "A";
  if (qfen0B_1A==0)
    elfen = "B";
  
  
  
  //descripcion
  double *vecor;
  vecor = new double[nudats];
  string traits[3] = {"n1", "n2", "mu"};
  double losdiv[3] = {1, (nodes*(nodes-1))/2, 1};
  orden = eldir+"description.txt";
  flusal.open(orden.c_str());
  for (i=0; i < 3; i++) {
    flusal << "**" << elfen << traits[i] << " in offspring of family 1 (equal parents)\n";
    flusal << "Mean: " << bas.get_mean(off1[i+indfen], nudats)/losdiv[i] << endl;
    flusal << "Std. Dev.: " << bas.get_sample_stddev(off1[i+indfen], nudats)/losdiv[i] << endl;
    bas.sort(off1[i+indfen], vecor, nudats);
    flusal << "Median: " << bas.get_midpoint(vecor, nudats)/losdiv[i] << endl;
    flusal << "Minimum: " << vecor[0]/losdiv[i] << endl;
    flusal << "Q1: " << bas.get_q1(vecor, nudats)/losdiv[i] << endl;
    flusal << "Q3: " << bas.get_q3(vecor, nudats)/losdiv[i] << endl;
    flusal << "Maximum: " << vecor[nudats-1]/losdiv[i] << endl << endl;
    
    flusal << "**" << elfen << traits[i] << " in offspring of family 2 (different parents)\n";
    flusal << "Mean: " << bas.get_mean(off2[i+indfen], nudats)/losdiv[i] << endl;
    flusal << "Std. Dev.: " << bas.get_sample_stddev(off2[i+indfen], nudats)/losdiv[i] << endl;
    bas.sort(off2[i+indfen], vecor, nudats);
    flusal << "Median: " << bas.get_midpoint(vecor, nudats)/losdiv[i] << endl;
    flusal << "Minimum: " << vecor[0]/losdiv[i] << endl;
    flusal << "Q1: " << bas.get_q1(vecor, nudats)/losdiv[i] << endl;
    flusal << "Q3: " << bas.get_q3(vecor, nudats)/losdiv[i] << endl;
    flusal << "Maximum: " << vecor[nudats-1]/losdiv[i] << endl << endl << endl;
  }
  flusal.close();
  delete [] vecor;
  //hastaaqui
  
  orden = eldir+"data.txt";
  flusal.open(orden.c_str());
  for (i = 0; i < nudats; i++) {//off1(pn1 pn2 mu) off2(pn1 pn2 mu) w1-w2(pn1 pn2 mu) off1-off2(pn1 pn2 mu)
    for (j = 0; j < 3; j++)
      flusal << off1[j+indfen][i] << "\t";
    for (j = 0; j < 3; j++)
      flusal << off2[j+indfen][i] << "\t";
    for (j = 0; j < 3; j++) {
      flusal << wife1[j+indfen][i]-wife2[j+indfen][i] << "\t";
    }
    for (j = 0; j < 3; j++) {
      flusal << off1[j+indfen][i]-off2[j+indfen][i] << "\t";
    }
    flusal << endl;
  }
  flusal.close();
  orden = eldir+"paR.sh";
  flusal.open(orden.c_str());
  flusal << "library(ggplot2)\n";
  flusal << "a <- read.csv(file = \""<< eldir << "data" << ".txt\", head=FALSE, sep =\"\\t\")\n";
  skahis_idl(flusal, eldir+elfen+"off2n1_off1n1", 4, 1, 1, 1, "A\'%*%\'AB offspring\\'s access to \'*italic("+elfen+")*\' through plasticity", "A\'%*%\'A offspring\\'s access to \'*italic("+elfen+")*\' through plasticity", 0.15);

  skahis_idl(flusal, eldir+elfen+"off2n2_off1n2", 5, ((nodes*(nodes-1))/2.0), 2, ((nodes*(nodes-1))/2.0), "A\'%*%\'AB offspring\\'s access to \'*italic("+elfen+")*\' through plasticity", "A\'%*%\'A offspring\\'s access to \'*italic("+elfen+")*\' through plasticity", 0.15);
  skahis_idl(flusal, eldir+elfen+"off2mu_off1mu", 6, 1, 3, 1, "A\'%*%\'AB offspring\\'s access to \'*italic("+elfen+")*\' through mutation", "A\'%*%\'A offspring\\'s access to \'*italic("+elfen+")*\' through mutation", 0.15);
  wilcoxon(flusal, eldir+"stats", 1, 4, 'l', elfen+"_off_f1_n1", elfen+"_off_f2_n1");

  wilcoxon(flusal, eldir+"stats", 2, 5, 'l', elfen+"_off_f1_n2", elfen+"_off_f2_n2");
  wilcoxon(flusal, eldir+"stats", 3, 6, 'l', elfen+"_off_f1_mu", elfen+"_off_f2_mu");
  
  flusal.close();
  orden = "Rscript "+eldir+"paR.sh";
  system(orden.c_str());
  
  return;
}

void uddutd(int qfen0B_1A, string dirarr, double **wife1, double **wife2, int **nat1, int **nat2, int nudats, int nodes) { //dirarr con /terminal
  ofstream flusal;
  string orden;
  string eldir = dirarr+"native/";
  orden = "mkdir "+eldir;
  system(orden.c_str());
  orden = eldir+"data.txt";
  int i, j, indfen, indnat;
  indfen =3*qfen0B_1A;
  indnat = 1-qfen0B_1A;
  string elfen = "A";
  if (qfen0B_1A==0)
    elfen = "B";

  double minwifrespn1 = 1000, maxwifrespn1 = 0;
  
  flusal.open(orden.c_str());
  for (i = 0; i < nudats; i++) {//w1(pn1 pn2 mu) w2(pn1 pn2 mu) w1-w2(pn1 pn2 mu) nat1 nat2 (nat1-nat2)
    for (j = 0; j < 3; j++)
      flusal << wife1[j+indfen][i] << "\t";
    for (j = 0; j < 3; j++)
      flusal << wife2[j+indfen][i] << "\t";
    for (j = 0; j < 3; j++) {
      flusal << wife1[j+indfen][i]-wife2[j+indfen][i] << "\t";
      if (j == 0) {
        if (wife1[j+indfen][i]-wife2[j+indfen][i] > maxwifrespn1)
          maxwifrespn1 = wife1[j+indfen][i]-wife2[j+indfen][i];
        if (wife1[j+indfen][i]-wife2[j+indfen][i] < minwifrespn1)
          minwifrespn1 = wife1[j+indfen][i]-wife2[j+indfen][i];
      }
    }
    flusal << nat1[indnat][i] << "\t" << nat2[indnat][i] << "\t" << nat1[indnat][i]-nat2[indnat][i] << endl;
  }
  flusal.close();
  
  //descripcion
  int *vecor;
  vecor = new int[nudats];
  orden = eldir+"description.txt";
  flusal.open(orden.c_str());
  flusal << "**" << elfen << " native in offspring of family 1 (equal parents)\n";
  flusal << "Mean: " << bas.get_mean(nat1[indnat], nudats) << endl;
  flusal << "Std. Dev.: " << bas.get_sample_stddev(nat1[indnat], nudats) << endl;
  bas.sort(nat1[indnat], vecor, nudats);
  flusal << "Median: " << bas.get_midpoint(vecor, nudats) << endl;
  flusal << "Minimum: " << vecor[0] << endl;
  flusal << "Q1: " << bas.get_q1(vecor, nudats) << endl;
  flusal << "Q3: " << bas.get_q3(vecor, nudats) << endl;
  flusal << "Maximum: " << vecor[nudats-1] << endl << endl;
  
  flusal << "**" << elfen << " native in offspring of family 2 (different parents)\n";
  flusal << "Mean: " << bas.get_mean(nat2[indnat], nudats) << endl;
  flusal << "Std. Dev.: " << bas.get_sample_stddev(nat2[indnat], nudats) << endl;
  bas.sort(nat2[indnat], vecor, nudats);
  flusal << "Median: " << bas.get_midpoint(vecor, nudats) << endl;
  flusal << "Minimum: " << vecor[0] << endl;
  flusal << "Q1: " << bas.get_q1(vecor, nudats) << endl;
  flusal << "Q3: " << bas.get_q3(vecor, nudats) << endl;
  flusal << "Maximum: " << vecor[nudats-1] << endl << endl;
  flusal.close();
  delete [] vecor;
  
  orden = eldir+"paR.sh";
  flusal.open(orden.c_str());
  flusal << "library(ggplot2)\n";
  flusal << "a <- read.csv(file = \""<< eldir << "data" << ".txt\", head=FALSE, sep =\"\\t\")\n";
  skahis_idl(flusal, eldir+elfen+"_natf2_natf1", 11, 1, 10, 1, "A\'%*%\'AB offspring with native \'*italic("+elfen+")*\'", "A\'%*%\'A offspring with native \'*italic("+elfen+")*\'", 0.15);
  wilcoxon(flusal, eldir+"stats", 10, 11, 'l', elfen+"_off_f1_native", elfen+"_off_f2_native");
  flusal.close();
  orden = "Rscript "+eldir+"paR.sh";
  system(orden.c_str());
  return;
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

void skahis_idl(ofstream& flusal, string archpdf, int colx, int divx, int coly, int divy, string labx, string laby, double alfa) { //alfa default 0.15
  flusal << "pdf(file=\"" << archpdf << ".pdf\", width=7.5,height=5.32)\n";
  flusal << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
  flusal << "p <- ggplot(a, aes(x=V" << colx;
  if (divx != 1)
    flusal << "/" << divx;
  flusal << ", y=V" << coly;
  if (divy != 1)
    flusal << "/" << divy;
  flusal << ")) + geom_point(alpha=" << alfa << ") + labs(x=expression(\'" << labx << "\'), y=expression(\'" << laby << "\')) + theme_classic() + geom_abline(intercept = 0, slope = 1) + theme(text=element_text(size=14))\n";
  flusal << "ggExtra::ggMarginal(p, type = \"histogram\", size=8)\n";
  flusal << "dev.off()\n";
  flusal << endl;
  return;
}

void cajaplo2(ofstream& flusal, string archpdf, int colx, int divx, int coly, int divy, string labx, string laby, double umbral) {
  flusal << "pdf(file=\"" << archpdf << ".pdf\", width=7.5,height=5.32)\n";
  flusal << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
  flusal << "ggplot(a, aes(x=V" << colx;
  flusal << ", y=V" << coly;
  if (divy != 1)
    flusal << "/" << divy;
  flusal << ")) + geom_point(size=1, shape=16, alpha=0.2, colour=\"black\", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour=\"grey20\")+ scale_x_discrete(labels=c(";
  if (divx > 1)
    flusal << "expression(\"\"<=" << umbral << "/" << divx << ")";
  else
    flusal << umbral;
  flusal << ",expression(\"\">";
  if (divx > 1)
    flusal << umbral << "/" << divx << ")";
  else
    flusal << umbral << ")";
  flusal << ")) + labs(x=expression(\'" << labx << "\'), y = expression(\'" << laby << "\')) + stat_summary(geom = \"point\", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))\n";
  flusal << "dev.off()\n";
  flusal << endl;
  return;
}
