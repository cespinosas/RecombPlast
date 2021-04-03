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
  int AganaBp = 0, AganaCp = 0, Aignoranciap = 0;
  int AganaB50 = 0, AganaC50 = 0, Aignorancia50 = 0;

  int ABganaBp = 0, ABganaCp = 0, ABignoranciap = 0;
  int ABganaB50 = 0, ABganaC50 = 0, ABignorancia50 = 0;

  int *Atp, *ABtp, *At50, *ABt50;
  
  int **AtpBC, **At50BC, **ABtpBC, **ABt50BC;
  int i;
  bas.create_array(AtpBC, 2, nusim);
  bas.create_array(At50BC, 2, nusim);
  bas.create_array(ABtpBC, 2, nusim);
  bas.create_array(ABt50BC, 2, nusim);
  for (i=0; i < nusim; i++) {
    bas.open_ifstream(fe, "Results/sinplas/tpB/"+bas.inttostring(i+1)+"tpB.txt");
    fe >> AtpBC[0][i];
    fe.close();
    bas.open_ifstream(fe, "Results/sinplas/t50B/"+bas.inttostring(i+1)+"t50B.txt");
    fe >> At50BC[0][i];
    fe.close();
    bas.open_ifstream(fe, "Results/sinplas/tpC/"+bas.inttostring(i+1)+"tpC.txt");
    fe >> AtpBC[1][i];
    fe.close();
    bas.open_ifstream(fe, "Results/sinplas/t50C/"+bas.inttostring(i+1)+"t50C.txt");
    fe >> At50BC[1][i];
    fe.close();
    
    bas.open_ifstream(fe, "Results/pn1/tpB/"+bas.inttostring(i+1)+"tpB.txt");
    fe >> ABtpBC[0][i];
    fe.close();
    bas.open_ifstream(fe, "Results/pn1/t50B/"+bas.inttostring(i+1)+"t50B.txt");
    fe >> ABt50BC[0][i];
    fe.close();
    bas.open_ifstream(fe, "Results/pn1/tpC/"+bas.inttostring(i+1)+"tpC.txt");
    fe >> ABtpBC[1][i];
    fe.close();
    bas.open_ifstream(fe, "Results/pn1/t50C/"+bas.inttostring(i+1)+"t50C.txt");
    fe >> ABt50BC[1][i];
    fe.close();
  }
//  int AganaBp = 0, AganaCp = 0, Aignoranciap = 0;
//  int AganaB50 = 0, AganaC50 = 0, Aignorancia50 = 0;
//
//  int ABganaBp = 0, ABganaCp = 0, ABignoranciap = 0;
//  int ABganaB50 = 0, ABganaC50 = 0, ABignorancia50 = 0;
  for (i = 0; i < nusim; i++) {
    if ((AtpBC[0][i] < 0) && (AtpBC[1][i] < 0))
      Aignoranciap++;
    else {
      if (((AtpBC[0][i] < 0) && (AtpBC[1][i] >= 0)) || ((AtpBC[0][i] >= 0) && (AtpBC[1][i] >= 0) && (AtpBC[1][i] <  AtpBC[0][i])))
        AganaCp++;
      else {
        if (((AtpBC[1][i] < 0) && (AtpBC[0][i] >= 0)) || ((AtpBC[1][i] >= 0) && (AtpBC[0][i] >= 0) && (AtpBC[0][i] <  AtpBC[1][i])))
          AganaBp++;
        else {
          if ((AtpBC[1][i] >= 0) && (AtpBC[0][i] >= 0) && (AtpBC[0][i] == AtpBC[1][i])) {
            cout << "ah jijo\n";
            exit(1);
          }
        }
      }
    }
    
    if ((At50BC[0][i] < 0) && (At50BC[1][i] < 0))
      Aignorancia50++;
    else {
      if (((At50BC[0][i] < 0) && (At50BC[1][i] >= 0)) || ((At50BC[0][i] >= 0) && (At50BC[1][i] >= 0) && (At50BC[1][i] <  At50BC[0][i])))
        AganaC50++;
      else {
        if (((At50BC[1][i] < 0) && (At50BC[0][i] >= 0)) || ((At50BC[1][i] >= 0) && (At50BC[0][i] >= 0) && (At50BC[0][i] <  At50BC[1][i])))
          AganaB50++;
        else {
          if ((At50BC[1][i] >= 0) && (At50BC[0][i] >= 0) && (At50BC[0][i] == At50BC[1][i])) {
            cout << "ah jijo\n";
            exit(1);
          }
        }
      }
    }
    
    
    if ((ABtpBC[0][i] < 0) && (ABtpBC[1][i] < 0))
      ABignoranciap++;
    else {
      if (((ABtpBC[0][i] < 0) && (ABtpBC[1][i] >= 0)) || ((ABtpBC[0][i] >= 0) && (ABtpBC[1][i] >= 0) && (ABtpBC[1][i] <  ABtpBC[0][i])))
        ABganaCp++;
      else {
        if (((ABtpBC[1][i] < 0) && (ABtpBC[0][i] >= 0)) || ((ABtpBC[1][i] >= 0) && (ABtpBC[0][i] >= 0) && (ABtpBC[0][i] <  ABtpBC[1][i])))
          ABganaBp++;
        else {
          if ((ABtpBC[1][i] >= 0) && (ABtpBC[0][i] >= 0) && (ABtpBC[0][i] == ABtpBC[1][i])) {
            cout << "ah jijo\n";
            exit(1);
          }
        }
      }
    }
    
    if ((ABt50BC[0][i] < 0) && (ABt50BC[1][i] < 0))
      ABignorancia50++;
    else {
      if (((ABt50BC[0][i] < 0) && (ABt50BC[1][i] >= 0)) || ((ABt50BC[0][i] >= 0) && (ABt50BC[1][i] >= 0) && (ABt50BC[1][i] <  ABt50BC[0][i])))
        ABganaC50++;
      else {
        if (((ABt50BC[1][i] < 0) && (ABt50BC[0][i] >= 0)) || ((ABt50BC[1][i] >= 0) && (ABt50BC[0][i] >= 0) && (ABt50BC[0][i] <  ABt50BC[1][i])))
          ABganaB50++;
        else {
          if ((ABt50BC[1][i] >= 0) && (ABt50BC[0][i] >= 0) && (ABt50BC[0][i] == ABt50BC[1][i])) {
            cout << "ah jijo\n";
            exit(1);
          }
        }
      }
    }
  }
  
//  int *Atp, *ABtp, *At50, *ABt50;
//
//  int **AtpBC, **At50BC, **ABtpBC, **ABt50BC;
  
  int j;
  j = 0;
  int ambasA = 0;
  Atp = new int[AganaBp+AganaCp];
  for (i = 0; i < nusim; i++) {
    if ((AtpBC[0][i] >= 0) && (AtpBC[1][i] >= 0)) {
      ambasA++;
      if (AtpBC[0][i] < AtpBC[1][i])
        Atp[j] = AtpBC[0][i];
      else
        Atp[j] = AtpBC[1][i];
      j++;
    }
    if ((AtpBC[0][i] >= 0) && (AtpBC[1][i] < 0)) {
      Atp[j] = AtpBC[0][i];
      j++;
    }
    if ((AtpBC[1][i] >= 0) && (AtpBC[0][i] < 0)) {
      Atp[j] = AtpBC[1][i];
      j++;
    }
  }
  
  int ambasAB = 0;
  j = 0;
  ABtp = new int[ABganaBp+ABganaCp];
  for (i = 0; i < nusim; i++) {
    if ((ABtpBC[0][i] >= 0) && (ABtpBC[1][i] >= 0)) {
      ambasAB++;
      if (ABtpBC[0][i] < ABtpBC[1][i])
        ABtp[j] = ABtpBC[0][i];
      else
        ABtp[j] = ABtpBC[1][i];
      j++;
    }
    if ((ABtpBC[0][i] >= 0) && (ABtpBC[1][i] < 0)) {
      ABtp[j] = ABtpBC[0][i];
      j++;
    }
    if ((ABtpBC[1][i] >= 0) && (ABtpBC[0][i] < 0)) {
      ABtp[j] = ABtpBC[1][i];
      j++;
    }
  }
  
  
  j = 0;
  At50 = new int[AganaB50+AganaC50];
  for (i = 0; i < nusim; i++) {
    if ((At50BC[0][i] >= 0) && (At50BC[1][i] >= 0)) {
      if (At50BC[0][i]  < At50BC[1][i])
        At50[j] = At50BC[0][i];// - AtpBC[0][i];
      else
        At50[j] = At50BC[1][i];// - AtpBC[1][i];
      j++;
    }
    if ((At50BC[0][i] >= 0) && (At50BC[1][i] < 0)) {
      At50[j] = At50BC[0][i];// - AtpBC[0][i];
      j++;
    }
    if ((At50BC[1][i] >= 0) && (At50BC[0][i] < 0)) {
      At50[j] = At50BC[1][i];// - AtpBC[1][i];
      j++;
    }
  }
  
  
  j = 0;
  ABt50 = new int[AganaB50+AganaC50];
  for (i = 0; i < nusim; i++) {
    if ((ABt50BC[0][i] >= 0) && (ABt50BC[1][i] >= 0)) {
      if (ABt50BC[0][i] < ABt50BC[1][i])
        ABt50[j] = ABt50BC[0][i];// - ABtpBC[0][i];
      else
        ABt50[j] = ABt50BC[1][i];// - ABtpBC[1][i];
      j++;
    }
    if ((ABt50BC[0][i] >= 0) && (ABt50BC[1][i] < 0)) {
      ABt50[j] = ABt50BC[0][i];// - ABtpBC[0][i];
      j++;
    }
    if ((ABt50BC[1][i] >= 0) && (ABt50BC[0][i] < 0)) {
      ABt50[j] = ABt50BC[1][i];// - ABtpBC[1][i];
      j++;
    }
  }
  
//  int AganaBp = 0, AganaCp = 0, Aignoranciap = 0;
//  int AganaB50 = 0, AganaC50 = 0, Aignorancia50 = 0;
//
//  int ABganaBp = 0, ABganaCp = 0, ABignoranciap = 0;
//  int ABganaB50 = 0, ABganaC50 = 0, ABignorancia50 = 0;
  
  bas.open_ofstream(fs, "Figs/stats.txt");
  fs << "** From A:\n";
  fs << "* Primer genotipo:\n";
  fs << "gana B: " << AganaBp << endl;
  fs << "gana C: " << AganaCp << endl;
  fs << "gana ignorancia: " << Aignoranciap << endl;
  fs << "Ambas: " << ambasA << endl;
  fs << "Promedio tiempo: " << bas.get_mean(Atp, AganaCp+AganaBp) << endl;
  fs << "Desv std: " << bas.get_sample_stddev(Atp, AganaCp+AganaBp) << endl << endl;
  
  
  fs << "* Final (total):\n";
  fs << "gana B: " << AganaB50 << endl;
  fs << "gana C: " << AganaC50 << endl;
  fs << "gana ignorancia: " << Aignorancia50 << endl;
  fs << "Promedio tiempo: " << bas.get_mean(At50, AganaC50+AganaB50) << endl;
  fs << "Desv std: " << bas.get_sample_stddev(At50, AganaC50+AganaB50) << endl << endl << endl;
  
  
  
  fs << "** From AB:\n";
  fs << "* Primer genotipo:\n";
  fs << "gana B: " << ABganaBp << endl;
  fs << "gana C: " << ABganaCp << endl;
  fs << "gana ignorancia: " << ABignoranciap << endl;
  fs << "Ambas: " << ambasAB << endl;
  fs << "Promedio tiempo: " << bas.get_mean(ABtp, ABganaCp+ABganaBp) << endl;
  fs << "Desv std: " << bas.get_sample_stddev(ABtp, ABganaCp+ABganaBp) << endl << endl;
  
  
  fs << "* Final (total):\n";
  fs << "gana B: " << ABganaB50 << endl;
  fs << "gana C: " << ABganaC50 << endl;
  fs << "gana ignorancia: " << ABignorancia50 << endl;
  fs << "Promedio tiempo: " << bas.get_mean(ABt50, ABganaC50+ABganaB50) << endl;
  fs << "Desv std: " << bas.get_sample_stddev(ABt50, ABganaC50+ABganaB50) << endl;
  
  fs.close();
  
  bas.open_ofstream(fs, "Figs/Atp.txt");
  for (i = 0; i < (AganaCp+AganaBp); i++)
    fs << Atp[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, "Figs/At50.txt");
  for (i = 0; i < (AganaC50+AganaB50); i++)
    fs << At50[i] << endl;
  fs.close();
  
  
  bas.open_ofstream(fs, "Figs/ABtp.txt");
  for (i = 0; i < (ABganaCp+ABganaBp); i++)
    fs << ABtp[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, "Figs/ABt50.txt");
  for (i = 0; i < (ABganaC50+ABganaB50); i++)
    fs << ABt50[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, "Figs/pacajaprimer.txt");
  for (i = 0; i < (AganaCp+AganaBp); i++)
    fs << Atp[i] << "\tA\n" ;
  for (i = 0; i < (ABganaCp+ABganaBp); i++)
    fs << ABtp[i] << "\tAB\n" ;
  fs.close();
  
  bas.open_ofstream(fs, "Figs/pacaja50.txt");
  for (i = 0; i < (AganaC50+AganaB50); i++)
    fs << At50[i] << "\tA\n" ;
  for (i = 0; i < (ABganaC50+ABganaB50); i++)
    fs << ABt50[i] << "\tAB\n" ;
  fs.close();
  
  
  
  
  
  
  bas.open_ofstream(fs, "Figs/paysfiga.tex");
  fs << "\\documentclass{standalone}\n";
  fs << "\\usepackage{tikz}\n";
  fs << "\\usetikzlibrary{arrows,decorations.pathmorphing,decorations.shapes,backgrounds,fit,snakes,shapes,positioning}\n";
  fs << "\\usepackage{times}\n";
  fs << "\\usepackage[T1]{fontenc}\n";
  fs << "\\renewcommand{\\familydefault}{\\sfdefault}\n";
  fs << "\\begin{document}\n";
  fs << "\\begin{tikzpicture}[rectangle, minimum width=0.33cm, minimum height=0.333cm, bla/.style={draw, fill=white},neg/.style={draw, fill=black}, gri/.style={draw, fill=black!50}, gena/.style={font=\\bf\\fontsize{8.04}{9.67}\\selectfont}, pac/.style={font=\\fontsize{9}{7.035} \\selectfont}]\n\n";
  
  fs << "\\node at (5,5) {Optimum first becoming a native phenotype};\n";
  fs << "\\filldraw[fill=black] (2, 2.5) -- +(90:1.5) arc (90:" << ((AganaBp/(1.0*nusim))*360)+90 << ":1.5) -- cycle;\n";
  fs << "\\filldraw[fill=white] (2, 2.5) -- +(" << ((AganaBp/(1.0*nusim))*360)+90 << ":1.5) arc (" << ((AganaBp/(1.0*nusim))*360)+90 << ":" << ((AganaBp/(1.0*nusim))*360)+90+((AganaCp/(1.0*nusim))*360) << ":1.5) -- cycle;\n";
  if (Aignoranciap > 0)
    fs << "\\filldraw[fill=black!50] (2, 2.5) -- +(" << ((AganaBp/(1.0*nusim))*360)+90+((AganaCp/(1.0*nusim))*360) << ":1.5) arc (" << ((AganaBp/(1.0*nusim))*360)+90+((AganaCp/(1.0*nusim))*360) << ":" << ((AganaBp/(1.0*nusim))*360)+90+((AganaCp/(1.0*nusim))*360)+((Aignoranciap/(1.0*nusim))*360) << ":1.5) -- cycle;\n";
  fs << "\\node[pac] at (2, 0.5) {Founder: $A$};\n";
  
  fs << "\\node[neg] at (4, 4) {};\n";
  fs << "\\node[pac] at (5.2, 4) {Phenotype $B$};\n";
  fs << "\\node[bla] at (4, 3.5) {};\n";
  fs << "\\node[pac] at (5.2, 3.5) {Phenotype $C$};\n";
  if (Aignoranciap > 0) {
    fs << "\\node[gri] at (4, 3) {};\n";
    fs << "\\node[pac] at (4.65, 3) {None};\n";
  }
  
  fs << "\\filldraw[fill=black] (8, 2.5) -- +(90:1.5) arc (90:" << ((ABganaBp/(1.0*nusim))*360)+90 << ":1.5) -- cycle;\n";
  fs << "\\filldraw[fill=white] (8, 2.5) -- +(" << ((ABganaBp/(1.0*nusim))*360)+90 << ":1.5) arc (" << ((ABganaBp/(1.0*nusim))*360)+90 << ":" << ((ABganaBp/(1.0*nusim))*360)+90+((ABganaCp/(1.0*nusim))*360) << ":1.5) -- cycle;\n";
  if (ABignoranciap > 0)
    fs << "\\filldraw[fill=black!50] (8, 2.5) -- +(" << ((ABganaBp/(1.0*nusim))*360)+90+((ABganaCp/(1.0*nusim))*360) << ":1.5) arc (" << ((ABganaBp/(1.0*nusim))*360)+90+((ABganaCp/(1.0*nusim))*360) << ":" << ((ABganaBp/(1.0*nusim))*360)+90+((ABganaCp/(1.0*nusim))*360)+((ABignoranciap/(1.0*nusim))*360) << ":1.5) -- cycle;\n";
  fs << "\\node[pac] at (8, 0.5) {Founder: $AB$};\n";
  fs << "\\end{tikzpicture}\n";
  fs << "\\end{document}\n";
  fs.close();
  system("pdflatex -output-directory Figs/ Figs/paysfiga.tex");
  
  
  bas.open_ofstream(fs, "Figs/paysega.tex");
  fs << "\\documentclass{standalone}\n";
  fs << "\\usepackage{tikz}\n";
  fs << "\\usetikzlibrary{arrows,decorations.pathmorphing,decorations.shapes,backgrounds,fit,snakes,shapes,positioning}\n";
  fs << "\\usepackage{times}\n";
  fs << "\\usepackage[T1]{fontenc}\n";
  fs << "\\renewcommand{\\familydefault}{\\sfdefault}\n";
  fs << "\\begin{document}\n";
  fs << "\\begin{tikzpicture}[rectangle, minimum width=0.33cm, minimum height=0.333cm, bla/.style={draw, fill=white},neg/.style={draw, fill=black}, gri/.style={draw, fill=black!50}, gena/.style={font=\\bf\\fontsize{8.04}{9.67}\\selectfont}, pac/.style={font=\\fontsize{9}{7.035} \\selectfont}]\n\n";
  
  fs << "\\node at (5,5) {Prevalent genetically assimilated optimum};\n";
  fs << "\\filldraw[fill=black] (2, 2.5) -- +(90:1.5) arc (90:" << ((AganaB50/(1.0*nusim))*360)+90 << ":1.5) -- cycle;\n";
  fs << "\\filldraw[fill=white] (2, 2.5) -- +(" << ((AganaB50/(1.0*nusim))*360)+90 << ":1.5) arc (" << ((AganaB50/(1.0*nusim))*360)+90 << ":" << ((AganaB50/(1.0*nusim))*360)+90+((AganaC50/(1.0*nusim))*360) << ":1.5) -- cycle;\n";
  if (Aignorancia50 > 0)
    fs << "\\filldraw[fill=black!50] (2, 2.5) -- +(" << ((AganaB50/(1.0*nusim))*360)+90+((AganaC50/(1.0*nusim))*360) << ":1.5) arc (" << ((AganaB50/(1.0*nusim))*360)+90+((AganaC50/(1.0*nusim))*360) << ":" << ((AganaB50/(1.0*nusim))*360)+90+((AganaC50/(1.0*nusim))*360)+((Aignorancia50/(1.0*nusim))*360) << ":1.5) -- cycle;\n";
  fs << "\\node[pac] at (2, 0.5) {Founder: $A$};\n";
  
  fs << "\\node[neg] at (4, 4) {};\n";
  fs << "\\node[pac] at (5.2, 4) {Phenotype $B$};\n";
  fs << "\\node[bla] at (4, 3.5) {};\n";
  fs << "\\node[pac] at (5.2, 3.5) {Phenotype $C$};\n";
  if (Aignorancia50 > 0) {
    fs << "\\node[gri] at (4, 3) {};\n";
    fs << "\\node[pac] at (4.65, 3) {None};\n";
  }
  
  fs << "\\filldraw[fill=black] (8, 2.5) -- +(90:1.5) arc (90:" << ((ABganaB50/(1.0*nusim))*360)+90 << ":1.5) -- cycle;\n";
  fs << "\\filldraw[fill=white] (8, 2.5) -- +(" << ((ABganaB50/(1.0*nusim))*360)+90 << ":1.5) arc (" << ((ABganaB50/(1.0*nusim))*360)+90 << ":" << ((ABganaB50/(1.0*nusim))*360)+90+((ABganaC50/(1.0*nusim))*360) << ":1.5) -- cycle;\n";
  if (ABignorancia50 > 0)
    fs << "\\filldraw[fill=black!50] (8, 2.5) -- +(" << ((ABganaB50/(1.0*nusim))*360)+90+((ABganaC50/(1.0*nusim))*360) << ":1.5) arc (" << ((ABganaB50/(1.0*nusim))*360)+90+((ABganaC50/(1.0*nusim))*360) << ":" << ((ABganaB50/(1.0*nusim))*360)+90+((ABganaC50/(1.0*nusim))*360)+((ABignorancia50/(1.0*nusim))*360) << ":1.5) -- cycle;\n";
  fs << "\\node[pac] at (8, 0.5) {Founder: $AB$};\n";
  fs << "\\end{tikzpicture}\n";
  fs << "\\end{document}\n";
  fs.close();
  system("pdflatex -output-directory Figs/ Figs/paysega.tex");
  
  
  int elpri;
  bas.open_ofstream(fs, "Figs/paR.sh");


  fs << "a <- \"\\n\\n\\n***A founder:\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs << "a <- \"**Probability of obtaining a greater deviation from expected (0.5) frequency of times that B appears first (two-tailed test), according to a binomial distribution:\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs << "a <- \"Number of succesful simulations: " << AganaBp+AganaCp << "\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
//  prob = 0;
  if (AganaBp > AganaCp)
    elpri = AganaBp;
  else
    elpri = AganaCp;
  fs << " bt <- binom.test(" << elpri << ", " << AganaBp+AganaCp << ", 0.5,alternative=\"two.sided\")\n";
  fs << "capture.output(bt, file=\"Figs/stats.txt\", append=TRUE)\n";


  
  fs << "a <- \"\\n\\n**Probability of obtaining a greater deviation from expected (0.5) frequency of times that B evolves (two-tailed test), according to a binomial distribution:\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs << "a <- \"Number of succesful simulations: " << AganaB50+AganaC50 << "\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";

//  prob = 0;
  if (AganaB50 > AganaC50)
    elpri = AganaB50;
  else
    elpri = AganaC50;
  fs << " bt <- binom.test(" << elpri << ", " << AganaB50+AganaC50 << ", 0.5,alternative=\"two.sided\")\n";
  fs << "capture.output(bt, file=\"Figs/stats.txt\", append=TRUE)\n";

 
  
  fs << "a <- \"\\n\\n\\n***AB founder:\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs << "a <- \"**Probability of obtaining a greater fraction of times that B appears first (one-tailed test), according to a binomial distribution:\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs << "a <- \"Number of succesful simulations: " << ABganaBp+ABganaCp << "\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  if (ABganaBp > ABganaCp)
    elpri = ABganaBp;
  else
    elpri = ABganaCp;
  fs << " bt <- binom.test(" << elpri << ", " << ABganaBp+ABganaCp << ", 0.5,alternative=\"two.sided\")\n";
  fs << "capture.output(bt, file=\"Figs/stats.txt\", append=TRUE)\n";
  
  fs << "a <- \"\\n\\n**Probability of obtaining a greater fraction of times that B evolves (one-tailed test), according to a binomial distribution:\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs << "a <- \"Number of succesful simulations: " << ABganaB50+ABganaC50 << "\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  if (ABganaB50 > ABganaC50)
    elpri = ABganaB50;
  else
    elpri = ABganaC50;
  fs << " bt <- binom.test(" << elpri << ", " << ABganaB50+ABganaC50 << ", 0.5,alternative=\"two.sided\")\n";
  fs << "capture.output(bt, file=\"Figs/stats.txt\", append=TRUE)\n";
  
  
  fs << "library(ggplot2)\n";
  fs << "a <- read.csv(file = \"Figs/pacajaprimer" << ".txt\", head=FALSE, sep =\"\\t\")\n";
  fs << "a$V2 <- factor(a$V2)\n";
  fs << "pdf(file=\"Figs/tfirstGA.pdf\", width=7.5,height=5.32)\n";
  fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
  fs << "ggplot(a, aes(x=V2, y=V1)) + geom_point(size=1, shape=16, alpha=0.2, colour=\"black\", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour=\"grey20\")+labs(x=expression(\'Founder\'), y = expression(\'First appearance of optimum as native phenotype\')) + stat_summary(geom = \"point\", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))\n";
  fs << "dev.off()\n";
  
  fs << "a <- read.csv(file = \"Figs/pacaja50" << ".txt\", head=FALSE, sep =\"\\t\")\n";
  fs << "a$V2 <- factor(a$V2)\n";
  fs << "pdf(file=\"Figs/tuntil50.pdf\", width=7.5,height=5.32)\n";
  fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
  fs << "ggplot(a, aes(x=V2, y=V1)) + geom_point(size=1, shape=16, alpha=0.2, colour=\"black\", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour=\"grey20\")+labs(x=expression(\'Founder\'), y = expression(\'Time until prevalence of assimilated optimum\')) + stat_summary(geom = \"point\", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))\n";
  fs << "dev.off()\n";
  
  
  fs << "equis <- read.csv(file=\"Figs/Atp.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "ye <- read.csv(file=\"Figs/ABtp.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "wtest <- wilcox.test(equis$V1, ye$V1, paired=FALSE, alternative=\"greater\")\n";
  fs << "a <- \"Mann-Whitney for time of appearance (one-tailed; time when founder A > time when founder AB):\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs << "capture.output(wtest, file=\"Figs/stats.txt\", append=TRUE)\n";

  
  fs << "equis <- read.csv(file=\"Figs/At50.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "ye <- read.csv(file=\"Figs/ABt50.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
  fs << "wtest <- wilcox.test(equis$V1, ye$V1, paired=FALSE, alternative=\"greater\")\n";
  fs << "a <- \"Mann-Whitney for time until dominance (one-tailed; time when founder A > time when founder AB):\"\n";
  fs << "capture.output(a, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs << "capture.output(wtest, file=\"Figs/stats.txt\", append=TRUE)\n";
  fs.close();
  
  bas.run_command("Rscript Figs/paR.sh");
  
  
  return 0;
}

