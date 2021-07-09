#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include "alea.h"
#include "basics.h"

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

//se corre desde donde vas a crear los directorios

int main(int argc, char **argv)
{
  int nusim = atoi(argv[1]);
  
  Alea jacta(1);
  Basics bas(jacta);
  ofstream fs;
  ifstream fe;
  int i, j;
  double gar;
  double **meanaccpperAfromB, **meanaccpredAfromB;
  bas.create_array(meanaccpperAfromB, 3, nusim);
  bas.create_array(meanaccpredAfromB, 3, nusim);
  
  for (i = 0; i < nusim; i++) {
    bas.open_ifstream(fe, "Results/pn1/meanaccpperAfromB/"+bas.inttostring(i+1)+".txt");
    for (j=0; j < 3; j++) {
      fe >> gar;
      fe >> meanaccpperAfromB[j][i];
    }
    fe.close();
    bas.open_ifstream(fe, "Results/pn1/meanaccpredAfromB/"+bas.inttostring(i+1)+".txt");
    for (j=0; j < 3; j++) {
      fe >> gar;
      fe >> meanaccpredAfromB[j][i];
    }
    fe.close();
  }
  
  int van550per = 0, van580per = 0, van5080per = 0;
  int van550red = 0, van580red = 0, van5080red = 0;
  for (i = 0; i < nusim; i++) {
    if (meanaccpperAfromB[0][i] >= 0) {
      if (meanaccpperAfromB[1][i] >= 0) {
        van550per++;
      }
      if (meanaccpperAfromB[2][i] >= 0) {
        van580per++;
      }
    }
    if ((meanaccpperAfromB[1][i] >= 0) && (meanaccpperAfromB[2][i] >= 0)) {
      van5080per++;
    }
    
    if (meanaccpredAfromB[0][i] >= 0) {
      if (meanaccpredAfromB[1][i] >= 0) {
        van550red++;
      }
      if (meanaccpredAfromB[2][i] >= 0) {
        van580red++;
      }
    }
    if ((meanaccpredAfromB[1][i] >= 0) && (meanaccpredAfromB[2][i] >= 0)) {
      van5080red++;
    }
  }
  
  bas.open_ofstream(fs, "Graphs/NforaccfromB.txt");
  fs << "**Pert:\n";
  fs << "van550: " << van550per << endl;
  fs << "van580: " << van580per << endl;
  fs << "van5080: " << van5080per << endl << endl;
  
  fs << "**Networks:\n";
  fs << "van550: " << van550red << endl;
  fs << "van580: " << van580red << endl;
  fs << "van5080: " << van5080red << endl << endl;
  fs.close();
  
  double **per550;
  bas.create_array(per550, 2, van550per);
  j = 0;
  for (i=0; i < nusim; i++) {
    if ((meanaccpperAfromB[0][i] >= 0) && (meanaccpperAfromB[1][i] >= 0)) {
      per550[0][j] = meanaccpperAfromB[0][i];
      per550[1][j] = meanaccpperAfromB[1][i];
      j++;
    }
  }
  
  double **per580;
  bas.create_array(per580, 2, van580per);
  j = 0;
  for (i=0; i < nusim; i++) {
    if ((meanaccpperAfromB[0][i] >= 0) && (meanaccpperAfromB[2][i] >= 0)) {
      per580[0][j] = meanaccpperAfromB[0][i];
      per580[1][j] = meanaccpperAfromB[2][i];
      j++;
    }
  }
  
  double **per5080;
  bas.create_array(per5080, 2, van5080per);
  j = 0;
  for (i=0; i < nusim; i++) {
    if ((meanaccpperAfromB[1][i] >= 0) && (meanaccpperAfromB[2][i] >= 0)) {
      per5080[0][j] = meanaccpperAfromB[1][i];
      per5080[1][j] = meanaccpperAfromB[2][i];
      j++;
    }
  }
  //
  double **red550;
  bas.create_array(red550, 2, van550red);
  j = 0;
  for (i=0; i < nusim; i++) {
    if ((meanaccpredAfromB[0][i] >= 0) && (meanaccpredAfromB[1][i] >= 0)) {
      red550[0][j] = meanaccpredAfromB[0][i];
      red550[1][j] = meanaccpredAfromB[1][i];
      j++;
    }
  }
  
  double **red580;
  bas.create_array(red580, 2, van580red);
  j = 0;
  for (i=0; i < nusim; i++) {
    if ((meanaccpredAfromB[0][i] >= 0) && (meanaccpredAfromB[2][i] >= 0)) {
      red580[0][j] = meanaccpredAfromB[0][i];
      red580[1][j] = meanaccpredAfromB[2][i];
      j++;
    }
  }
  
  double **red5080;
  bas.create_array(red5080, 2, van5080red);
  j = 0;
  for (i=0; i < nusim; i++) {
    if ((meanaccpredAfromB[1][i] >= 0) && (meanaccpredAfromB[2][i] >= 0)) {
      red5080[0][j] = meanaccpredAfromB[1][i];
      red5080[1][j] = meanaccpredAfromB[2][i];
      j++;
    }
  }
  bas.open_ofstream(fs, "Graphs/stats2.txt");
  fs << "**Pert:\n";
  fs << "*5 vs 50:\n";
  fs << "Mean 5: " << bas.get_mean(per550[0], van550per) << endl;
  fs << "SD 5: " << bas.get_sample_stddev(per550[0], van550per) << endl;
  fs << "Mean 50: " << bas.get_mean(per550[1], van550per) << endl;
  fs << "SD 50: " << bas.get_sample_stddev(per550[1], van550per) << endl;
  fs << endl;
  fs << "*5 vs 80:\n";
  fs << "Mean 5: " << bas.get_mean(per580[0], van580per) << endl;
  fs << "SD 5: " << bas.get_sample_stddev(per580[0], van580per) << endl;
  fs << "Mean 80: " << bas.get_mean(per580[1], van580per) << endl;
  fs << "SD 80: " << bas.get_sample_stddev(per580[1], van580per) << endl;
  fs << endl;
  fs << "*50 vs 80:\n";
  fs << "Mean 50: " << bas.get_mean(per5080[0], van5080per) << endl;
  fs << "SD 50: " << bas.get_sample_stddev(per5080[0], van5080per) << endl;
  fs << "Mean 80: " << bas.get_mean(per5080[1], van5080per) << endl;
  fs << "SD 80: " << bas.get_sample_stddev(per5080[1], van5080per) << endl;
  fs << endl;
  fs << endl;
  
  fs << "**Red:\n";
  fs << "*5 vs 50:\n";
  fs << "Mean 5: " << bas.get_mean(red550[0], van550red) << endl;
  fs << "SD 5: " << bas.get_sample_stddev(red550[0], van550red) << endl;
  fs << "Mean 50: " << bas.get_mean(red550[1], van550red) << endl;
  fs << "SD 50: " << bas.get_sample_stddev(red550[1], van550red) << endl;
  fs << endl;
  fs << "*5 vs 80:\n";
  fs << "Mean 5: " << bas.get_mean(red580[0], van580red) << endl;
  fs << "SD 5: " << bas.get_sample_stddev(red580[0], van580red) << endl;
  fs << "Mean 80: " << bas.get_mean(red580[1], van580red) << endl;
  fs << "SD 80: " << bas.get_sample_stddev(red580[1], van580red) << endl;
  fs << endl;
  fs << "*50 vs 80:\n";
  fs << "Mean 50: " << bas.get_mean(red5080[0], van5080red) << endl;
  fs << "SD 50: " << bas.get_sample_stddev(red5080[0], van5080red) << endl;
  fs << "Mean 80: " << bas.get_mean(red5080[1], van5080red) << endl;
  fs << "SD 80: " << bas.get_sample_stddev(red5080[1], van5080red) << endl;
  fs << endl;
  fs.close();
  
  bas.open_ofstream(fs, "Graphs/accAperfromB_5vs50.txt");
  for (i=0; i < van550per; i++) {
    fs << per550[0][i] << "\t" << per550[1][i] << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Graphs/accAperfromB_5vs80.txt");
  for (i=0; i < van580per; i++) {
    fs << per580[0][i] << "\t" << per580[1][i] << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Graphs/accAperfromB_50vs80.txt");
  for (i=0; i < van5080per; i++) {
    fs << per5080[0][i] << "\t" << per5080[1][i] << endl;
  }
  fs.close();
  
  bas.open_ofstream(fs, "Graphs/accAredfromB_5vs50.txt");
  for (i=0; i < van550red; i++) {
    fs << red550[0][i] << "\t" << red550[1][i] << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Graphs/accAredfromB_5vs80.txt");
  for (i=0; i < van580red; i++) {
    fs << red580[0][i] << "\t" << red580[1][i] << endl;
  }
  fs.close();
  bas.open_ofstream(fs, "Graphs/accAredfromB_50vs80.txt");
  for (i=0; i < van5080red; i++) {
    fs << red5080[0][i] << "\t" << red5080[1][i] << endl;
  }
  fs.close();
  
  
  bas.open_ofstream(fs, "Graphs/paR2.sh");
  ///////////
  //iNTRO
  
  fs << "equis <- read.csv(file=\"Graphs/accAperfromB_5vs50.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "wtest <- wilcox.test(equis$V1, equis$V2, paired=TRUE, alternative=\"greater\")\n";
  fs << "a <- \"Wilcoxon test, 5pc vs 50pc, for access (pert) to GAP A through plasticity in networks with default B:\"\n";
  fs << "capture.output(a, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  fs << "capture.output(wtest, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  
  
  fs << "equis <- read.csv(file=\"Graphs/accAperfromB_5vs80.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "wtest <- wilcox.test(equis$V1, equis$V2, paired=TRUE, alternative=\"greater\")\n";
  fs << "a <- \"Wilcoxon test, 5pc vs 80pc, for access (pert) to GAP A through plasticity in networks with default B:\"\n";
  fs << "capture.output(a, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  fs << "capture.output(wtest, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  
  
  fs << "equis <- read.csv(file=\"Graphs/accAperfromB_50vs80.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "wtest <- wilcox.test(equis$V1, equis$V2, paired=TRUE, alternative=\"greater\")\n";
  fs << "a <- \"Wilcoxon test, 50pc vs 80pc, for access (pert) to GAP A through plasticity in networks with default B:\"\n";
  fs << "capture.output(a, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  fs << "capture.output(wtest, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  
  //
  fs << "equis <- read.csv(file=\"Graphs/accAredfromB_5vs50.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "wtest <- wilcox.test(equis$V1, equis$V2, paired=TRUE, alternative=\"greater\")\n";
  fs << "a <- \"Wilcoxon test, 5pc vs 50pc, for access (red) to GAP A through plasticity in networks with default B:\"\n";
  fs << "capture.output(a, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  fs << "capture.output(wtest, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  
  
  fs << "equis <- read.csv(file=\"Graphs/accAredfromB_5vs80.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "wtest <- wilcox.test(equis$V1, equis$V2, paired=TRUE, alternative=\"greater\")\n";
  fs << "a <- \"Wilcoxon test, 5pc vs 80pc, for access (red) to GAP A through plasticity in networks with default B:\"\n";
  fs << "capture.output(a, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  fs << "capture.output(wtest, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  
  
  fs << "equis <- read.csv(file=\"Graphs/accAredfromB_50vs80.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "wtest <- wilcox.test(equis$V1, equis$V2, paired=TRUE, alternative=\"greater\")\n";
  fs << "a <- \"Wilcoxon test, 50pc vs 80pc, for access (red) to GAP A through plasticity in networks with default B:\"\n";
  fs << "capture.output(a, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  fs << "capture.output(wtest, file=\"Graphs/stats2.txt\", append=TRUE)\n";
  
  
  ////
  
  fs.close();
  
  bas.run_command("Rscript Graphs/paR2.sh");
  
  
  
  
  
  
  return 0;
}

