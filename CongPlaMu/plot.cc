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

void cajaplo2(ofstream& flusal, string archpdf, int colx, int divx, int coly, int divy, string labx, string laby, double umbral);

int main(int argc, char **argv) {
  int cuantos = atoi(argv[1]);
  int nodes = atoi(argv[2]);
  ifstream fe;
  ofstream fs, fsd;
  int seed = 1;
  Alea jacta(seed);
  Basics bas(jacta);
  int i, j,k,l;

  string hazesto = "mkdir Figs";
  system(hazesto.c_str());
  
  double **datos;
  bas.create_array(datos, 6, cuantos);

  bas.open_ifstream(fe, "Results/CongPlaMu.txt");
  for (i=0; i<cuantos; i++)
    for (j=0; j < 6; j++)
      fe >> datos[j][i];
  fe.close();
  
  
  string losfen[2] = {"B", "A"};
  string traits[3] = {"n1", "n2", "mu"};
  
  double losmin[6];
  losmin[1] = nodes-1;
  losmin[4] = nodes-1;
  
  double *vecor;
  vecor = new double[cuantos];
  
  int *contmin;
  contmin = new int[6];
  bas.fillv0(contmin, 6);
  
  bas.open_ofstream(fs, "Figs/desc_data.txt");
  i=0;
    for (j=0; j < 3; j++) {
      bas.fillv0(vecor, cuantos);
      fs << "**Phenotype " << losfen[i] << "; trait: " << traits[j] << endl;
      fs << "Mean: " << bas.get_mean(datos[(i*3)+j], cuantos) << endl;
      fs << "Std. Dev.: " << bas.get_sample_stddev(datos[(i*3)+j], cuantos) << endl;
      bas.sort(datos[(i*3)+j], vecor, cuantos);
      fs << "Median: " << bas.get_midpoint(vecor, cuantos) << endl;
      fs << "Minimum: " << vecor[0] << endl;
      fs << "Q1: " << bas.get_q1(vecor, cuantos) << endl;
      fs << "Q3: " << bas.get_q3(vecor, cuantos) << endl;
      fs << "Maximum: " << vecor[cuantos-1] << endl << endl;
      if (j != 1)
        losmin[(i*3)+j] = vecor[0];
    }
  fs.close();
  
  bas.open_ofstream(fs, "Figs/data_minfact.txt");
  for (i=0; i < cuantos; i++) {
    for (j=0; j < 6; j++)
      fs << datos[j][i] << "\t";
    for (j=0; j < 6; j++) {
      if (datos[j][i] <= losmin[j]) {
        fs << 0 << "\t";
        contmin[j]++;
      }
      else
        fs << 1 << "\t";
    }
    fs << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Figs/count_lowestclass.txt");
  fs << "Number of data points in lowest class (n2 <= 2/nodes):\n";
  i=0;
    for (j=0; j < 3; j++)
      fs << losfen[i] << " " << traits[j] << ": " << contmin[(i*3)+j] << endl;
  fs.close();
  
  double ****minifentra, ****resifentra; //2,3,6
  int **contminft, **contresft;
  minifentra = new double***[2];
  resifentra = new double***[2];
  contminft = new int*[2];
  contresft = new int*[2];
  for (i=0; i < 2; i++) {
    minifentra[i] = new double**[3];
    resifentra[i] = new double**[3];
    contminft[i] = new int[3];
    contresft[i] = new int[3];
    for (j=0; j < 3; j++) {
      minifentra[i][j] = new double*[6];
      resifentra[i][j] = new double*[6];
      for (k = 0; k < 6; k++) {
        minifentra[i][j][k] = new double[contmin[(i*3)+j]];
        resifentra[i][j][k] = new double[cuantos-contmin[(i*3)+j]];
      }
      contminft[i][j] = 0;
      contresft[i][j] = 0;
    }
  }
  
  i=0;
    for (j=0; j < 3; j++) {
      bas.open_ofstream(fs, "Figs/Min_"+losfen[i]+"_"+traits[j]+".txt");
      bas.open_ofstream(fsd, "Figs/Rest_"+losfen[i]+"_"+traits[j]+".txt");
      for (k=0; k < cuantos; k++) {
        if (datos[(i*3)+j][k] <= losmin[(i*3)+j]) {
          for (l=0; l < 6; l++) {
            fs << datos[l][k] << "\t";
            minifentra[i][j][l][contminft[i][j]] = datos[l][k];
          }
          fs << endl;
          contminft[i][j]++;
        }
        else {
          for (l=0; l < 6; l++) {
            fsd << datos[l][k] << "\t";
            resifentra[i][j][l][contresft[i][j]] = datos[l][k];
          }
          fsd << endl;
          contresft[i][j]++;
        }
      }
      fsd.close();
      fs.close();
    }
  delete [] vecor;
  
  bas.open_ofstream(fs, "Figs/desc_data_2cat.txt");
  i=0;
    fs << "#################################\n";
    for (j=0; j < 3; j++) {
      if (contminft[i][j] > 0) {
        fs << "****Phenotype " << losfen[i] << "; minimum trait: " << traits[j] << endl;
        fs << "******************************\n";
        vecor = new double[contminft[i][j]];
        bas.fillv0(vecor, contminft[i][j]);
        for (k = 0; k < 3; k++) {
          fs << "++Trait: " << traits[k] << endl;
          fs << "Mean: " << bas.get_mean(minifentra[i][j][(i*3)+k], contminft[i][j]) << endl;
          fs << "Std Dev: " << bas.get_sample_stddev(minifentra[i][j][(i*3)+k], contminft[i][j]) << endl;
          bas.sort(minifentra[i][j][(i*3)+k], vecor,contminft[i][j]);
          fs << "Median: " << bas.get_midpoint(vecor, contminft[i][j]) << endl;
          fs << "Minimum: " << vecor[0] << endl;
          fs << "Q1: " << bas.get_q1(vecor, contminft[i][j]) << endl;
          fs << "Q3: " << bas.get_q3(vecor, contminft[i][j]) << endl;
          fs << "Maximum: " << vecor[contminft[i][j]-1] << endl << endl;
        }
        delete [] vecor;
        fs << endl << endl;
      }
      if (contresft[i][j] > 0) {
        fs << "****Phenotype " << losfen[i] << "; NOT minimum trait: " << traits[j] << endl;
        fs << "******************************\n";
        vecor = new double[contresft[i][j]];
        bas.fillv0(vecor, contresft[i][j]);
        for (k = 0; k < 3; k++) {
          fs << "++Trait: " << traits[k] << endl;
          fs << "Mean: " << bas.get_mean(resifentra[i][j][(i*3)+k], contresft[i][j]) << endl;
          fs << "Std Dev: " << bas.get_sample_stddev(resifentra[i][j][(i*3)+k], contresft[i][j]) << endl;
          bas.sort(resifentra[i][j][(i*3)+k], vecor,contresft[i][j]);
          fs << "Median: " << bas.get_midpoint(vecor, contresft[i][j]) << endl;
          fs << "Minimum: " << vecor[0] << endl;
          fs << "Q1: " << bas.get_q1(vecor, contresft[i][j]) << endl;
          fs << "Q3: " << bas.get_q3(vecor, contresft[i][j]) << endl;
          fs << "Maximum: " << vecor[contresft[i][j]-1] << endl << endl;
        }
        delete [] vecor;
      }
    }
  

  fs.close();
  
  
  double divx, divy;
  string letreros[6] = {"Access to \'*italic(B)*\' through plasticity", "Access to \'*italic(B)*\' through plasticity", "Access to \'*italic(B)*\' through mutations", "Access to \'*italic(A)*\' through plasticity", "Access to \'*italic(A)*\' through plasticity", "Access to \'*italic(A)*\' through mutations"};
  bas.open_ofstream(fs, "Figs/paR2.sh");
  fs << "library(ggplot2)\n";
  fs << "a <- read.csv(file = \"Figs/data_minfact" << ".txt\", head=FALSE, sep =\"\\t\")\n";
  i=0;
    for (j=0; j < 3; j++) {
      fs << "a$V" << 7+(i*3)+j << " <- factor(a$V" << 7+(i*3)+j << ")\n";
      if (j == 1)
        divx = (nodes*(nodes-1))/2.0;
      else
        divx = 1;
      for (k=0; k < 3; k++) {
        if (k != j) {
          if (k == 1)
            divy = (nodes*(nodes-1))/2.0;
          else
            divy = 1;
          cajaplo2(fs, "Figs/2cat"+losfen[i]+"_"+traits[j]+"_"+traits[k], 7+(i*3)+j, divx, k+1, divy, letreros[(i*3)+j], letreros[(i*3)+k], losmin[(i*3)+j]);
        }
      }
    }
  i=0;
    for (j=0; j < 3; j++) {
      fs << "m <- read.csv(file = \"Figs/Min_"+losfen[i]+"_"+traits[j] << ".txt\", head=FALSE, sep =\"\\t\")\n";
      fs << "r <- read.csv(file = \"Figs/Rest_"+losfen[i]+"_"+traits[j] << ".txt\", head=FALSE, sep =\"\\t\")\n";
      for (k=0; k < 3; k++) {
        fs << "wtest <- wilcox.test(m$V" << 1+(i*3)+k << ",r$V" << 1+(i*3)+k << ", paired=FALSE, alternative=\"less\")\n";
        fs << "letr <- \"Mann-Whitney: Phen" << losfen[i] << " " << traits[k] << " with minimum " << traits[j] << " less than with higher values of " << traits[j] << "\"" << endl;
        fs << "capture.output(letr, file=\"Figs/stats2.txt\", append = TRUE)\n";
        fs << "capture.output(wtest, file=\"Figs/stats2.txt\", append = TRUE)\n";
        fs << "letr <- \"****************************************\"\n";
        fs << "capture.output(letr, file=\"Figs/stats2.txt\", append = TRUE)\n";
        fs << endl;
      }

    }
  fs.close();
  
  string orden = "Rscript Figs/paR2.sh";
  system(orden.c_str());
  
  
  
  jacta.close_rng();
  return 0;
}

void cajaplo2(ofstream& flusal, string archpdf, int colx, int divx, int coly, int divy, string labx, string laby, double umbral) {
  flusal << "pdf(file=\"" << archpdf << ".pdf\", width=7.5,height=5.32)\n";
  flusal << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
  flusal << "ggplot(a, aes(x=V" << colx;
  flusal << ", y=V" << coly;
  if (divy != 1)
    flusal << "/" << divy;
  flusal << ")) + geom_point(size=1, shape=16, alpha=0.2, colour=\"black\", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour=\"grey20\")+ scale_x_discrete(labels=c(";
  if ((colx == 8) || (colx == 11))
    flusal << "expression(\"\"<=" << umbral << "/" << divx << ")";
  else
    flusal << umbral;
  flusal << ",expression(\"\">";
  if ((colx == 8) || (colx == 11))
    flusal << umbral << "/" << divx << ")";
  else
    flusal << umbral << ")";
  flusal << ")) + labs(x=expression(\'" << labx << "\'), y = expression(\'" << laby << "\')) + stat_summary(geom = \"point\", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))\n";
  flusal << "dev.off()\n";
  flusal << endl;
  return;
}
