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
    int maxgen = atoi(argv[2]);
    
    Alea jacta(1);
    Basics bas(jacta);
    ofstream fs;
    ifstream fe;
    system("mkdir Graphs");
    bas.open_ofstream(fs, "Graphs/maxgen.txt");
    fs << maxgen << endl;
    fs.close();
    int AganaBp = 0, AganaCp = 0, Aignoranciap = 0;
    int AganaB5 = 0, AganaC5 = 0, Aignorancia5 = 0;
    int AganaB50 = 0, AganaC50 = 0, Aignorancia50 = 0;
    
    int ABganaBp = 0, ABganaCp = 0, ABignoranciap = 0;
    int ABganaB5 = 0, ABganaC5 = 0, ABignorancia5 = 0;
    int ABganaB50 = 0, ABganaC50 = 0, ABignorancia50 = 0;
    int ABganaB80 = 0, ABganaC80 = 0, ABignorancia80 = 0;
    
    int *Atp, *ABtp, *At50, *ABt50;
    int *At5, *ABt5, *ABt80; //, *At80
    
    int **AtpBC, **At50BC, **ABtpBC, **ABt50BC;
    int **At5BC, **At80BC, **ABt5BC, **ABt80BC;
    int i, j;
    bas.create_array(AtpBC, 2, nusim);
    bas.create_array(At50BC, 2, nusim);
    bas.create_array(ABtpBC, 2, nusim);
    bas.create_array(ABt50BC, 2, nusim);
    
    bas.create_array(At5BC, 2, nusim);
    bas.create_array(At80BC, 2, nusim);
    bas.create_array(ABt5BC, 2, nusim);
    bas.create_array(ABt80BC, 2, nusim);
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
        
        bas.open_ifstream(fe, "Results/sinplas/t5B/"+bas.inttostring(i+1)+"t5B.txt");
        fe >> At5BC[0][i];
        fe.close();
        bas.open_ifstream(fe, "Results/sinplas/t80B/"+bas.inttostring(i+1)+"t80B.txt");
        fe >> At80BC[0][i];
        fe.close();
        bas.open_ifstream(fe, "Results/sinplas/t5C/"+bas.inttostring(i+1)+"t5C.txt");
        fe >> At5BC[1][i];
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
        
        bas.open_ifstream(fe, "Results/pn1/t5B/"+bas.inttostring(i+1)+"t5B.txt");
        fe >> ABt5BC[0][i];
        fe.close();
        bas.open_ifstream(fe, "Results/pn1/t80B/"+bas.inttostring(i+1)+"t80B.txt");
        fe >> ABt80BC[0][i];
        fe.close();
        bas.open_ifstream(fe, "Results/pn1/t5C/"+bas.inttostring(i+1)+"t5C.txt");
        fe >> ABt5BC[1][i];
        fe.close();
        bas.open_ifstream(fe, "Results/pn1/t80C/"+bas.inttostring(i+1)+"t80C.txt");
        fe >> ABt80BC[1][i];
        fe.close();
    }
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
        //
        if ((At5BC[0][i] < 0) && (At5BC[1][i] < 0))
            Aignorancia5++;
        else {
            if (((At5BC[0][i] < 0) && (At5BC[1][i] >= 0)) || ((At5BC[0][i] >= 0) && (At5BC[1][i] >= 0) && (At5BC[1][i] <  At5BC[0][i])))
                AganaC5++;
            else {
                if (((At5BC[1][i] < 0) && (At5BC[0][i] >= 0)) || ((At5BC[1][i] >= 0) && (At5BC[0][i] >= 0) && (At5BC[0][i] <  At5BC[1][i])))
                    AganaB5++;
                else {
                    if ((At5BC[1][i] >= 0) && (At5BC[0][i] >= 0) && (At5BC[0][i] == At5BC[1][i])) {
                        cout << "ah jijo\n";
                        exit(1);
                    }
                }
            }
        }
        //
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
        //
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
        //
        if ((ABt5BC[0][i] < 0) && (ABt5BC[1][i] < 0))
            ABignorancia5++;
        else {
            if (((ABt5BC[0][i] < 0) && (ABt5BC[1][i] >= 0)) || ((ABt5BC[0][i] >= 0) && (ABt5BC[1][i] >= 0) && (ABt5BC[1][i] <  ABt5BC[0][i])))
                ABganaC5++;
            else {
                if (((ABt5BC[1][i] < 0) && (ABt5BC[0][i] >= 0)) || ((ABt5BC[1][i] >= 0) && (ABt5BC[0][i] >= 0) && (ABt5BC[0][i] <  ABt5BC[1][i])))
                    ABganaB5++;
                else {
                    if ((ABt5BC[1][i] >= 0) && (ABt5BC[0][i] >= 0) && (ABt5BC[0][i] == ABt5BC[1][i])) {
                        cout << "ah jijo\n";
                        exit(1);
                    }
                }
            }
        }//
        
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
        //
        if ((ABt80BC[0][i] < 0) && (ABt80BC[1][i] < 0))
            ABignorancia80++;
        else {
            if (((ABt80BC[0][i] < 0) && (ABt80BC[1][i] >= 0)) || ((ABt80BC[0][i] >= 0) && (ABt80BC[1][i] >= 0) && (ABt80BC[1][i] <  ABt80BC[0][i])))
                ABganaC80++;
            else {
                if (((ABt80BC[1][i] < 0) && (ABt80BC[0][i] >= 0)) || ((ABt80BC[1][i] >= 0) && (ABt80BC[0][i] >= 0) && (ABt80BC[0][i] <  ABt80BC[1][i])))
                    ABganaB80++;
                else {
                    if ((ABt80BC[1][i] >= 0) && (ABt80BC[0][i] >= 0) && (ABt80BC[0][i] == ABt80BC[1][i])) {
                        cout << "ah jijo\n";
                        exit(1);
                    }
                }
            }
        }
        
        
    }
    int ambasA = 0;
    int ambasAB = 0;
    
    for (i=0; i < nusim; i++) {
        
        if ((AtpBC[0][i] >= 0) && (AtpBC[1][i] >= 0))
            ambasA++;
        if ((ABtpBC[0][i] >= 0) && (ABtpBC[1][i] >= 0))
            ambasAB++;
        for (j=0; j < 2; j++) {
            if (AtpBC[j][i] < 0)
                AtpBC[j][i] = maxgen;
            if (At5BC[j][i] < 0)
                At5BC[j][i] = maxgen;
            if (At50BC[j][i] < 0)
                At50BC[j][i] = maxgen;
            
            if (ABtpBC[j][i] < 0)
                ABtpBC[j][i] = maxgen;
            if (ABt5BC[j][i] < 0)
                ABt5BC[j][i] = maxgen;
            if (ABt50BC[j][i] < 0)
                ABt50BC[j][i] = maxgen;
            
            if (ABt80BC[j][i] < 0)
                ABt80BC[j][i] = maxgen;
        }
    }
    
    Atp = new int[nusim];
    At5 = new int[nusim];
    At50 = new int[nusim];
    
    ABtp = new int[nusim];
    ABt5 = new int[nusim];
    ABt50 = new int[nusim];
    ABt80 = new int[nusim];
    
    for (i=0; i < nusim; i++) {
        if (AtpBC[0][i] < AtpBC[1][i])
            Atp[i] = AtpBC[0][i];
        else
            Atp[i] = AtpBC[1][i];
        if (At5BC[0][i] < At5BC[1][i])
            At5[i] = At5BC[0][i];
        else
            At5[i] = At5BC[1][i];
        if (At50BC[0][i] < At50BC[1][i])
            At50[i] = At50BC[0][i];
        else
            At50[i] = At50BC[1][i];
        
        if (ABtpBC[0][i] < ABtpBC[1][i])
            ABtp[i] = ABtpBC[0][i];
        else
            ABtp[i] = ABtpBC[1][i];
        if (ABt5BC[0][i] < ABt5BC[1][i])
            ABt5[i] = ABt5BC[0][i];
        else
            ABt5[i] = ABt5BC[1][i];
        if (ABt50BC[0][i] < ABt50BC[1][i])
            ABt50[i] = ABt50BC[0][i];
        else
            ABt50[i] = ABt50BC[1][i];

        if (ABt80BC[0][i] < ABt80BC[1][i])
            ABt80[i] = ABt80BC[0][i];
        else
            ABt80[i] = ABt80BC[1][i];
    }
    
    int *vo;
    vo = new int[nusim];
    
    
    bas.open_ofstream(fs, "Graphs/stats.txt");
    fs << "** From A:\n";
    fs << "* Primer genotipo:\n";
    fs << "gana B: " << AganaBp << endl;
    fs << "gana C: " << AganaCp << endl;
    fs << "gana ignorancia: " << Aignoranciap << endl;
    fs << "Ambas: " << ambasA << endl;
    fs << "Promedio tiempo: " << bas.get_mean(Atp, nusim) << endl;
    fs << "Desv std: " << bas.get_sample_stddev(Atp, nusim) << endl;
    bas.sort(Atp, vo, nusim);
    fs << "Mediana: " << bas.get_midpoint(vo, nusim) << endl;
    fs << "Q1: " <<  bas.get_q1(vo, nusim) << endl;
    fs << "Q3: " <<  bas.get_q3(vo, nusim) << endl << endl;
    
    
    fs << "* 5pc:\n";
    fs << "gana B: " << AganaB5 << endl;
    fs << "gana C: " << AganaC5 << endl;
    fs << "gana ignorancia: " << Aignorancia5 << endl;
    fs << "Promedio tiempo: " << bas.get_mean(At5, nusim) << endl;
    fs << "Desv std: " << bas.get_sample_stddev(At5, nusim) << endl;
    bas.sort(At5, vo, nusim);
    fs << "Mediana: " << bas.get_midpoint(vo, nusim) << endl;
    fs << "Q1: " <<  bas.get_q1(vo, nusim) << endl;
    fs << "Q3: " <<  bas.get_q3(vo, nusim) << endl << endl;
    
    fs << "* 50pc:\n";
    fs << "gana B: " << AganaB50 << endl;
    fs << "gana C: " << AganaC50 << endl;
    fs << "gana ignorancia: " << Aignorancia50 << endl;
    fs << "Promedio tiempo: " << bas.get_mean(At50, nusim) << endl;
    fs << "Desv std: " << bas.get_sample_stddev(At50, nusim) << endl;
    bas.sort(At50, vo, nusim);
    fs << "Mediana: " << bas.get_midpoint(vo, nusim) << endl;
    fs << "Q1: " <<  bas.get_q1(vo, nusim) << endl;
    fs << "Q3: " <<  bas.get_q3(vo, nusim) << endl << endl << endl;
    
    
    fs << "** From AB:\n";
    fs << "* Primer genotipo:\n";
    fs << "gana B: " << ABganaBp << endl;
    fs << "gana C: " << ABganaCp << endl;
    fs << "gana ignorancia: " << ABignoranciap << endl;
    fs << "Ambas: " << ambasAB << endl;
    fs << "Promedio tiempo: " << bas.get_mean(ABtp, nusim) << endl;
    fs << "Desv std: " << bas.get_sample_stddev(ABtp, nusim) << endl;
    bas.sort(ABtp, vo, nusim);
    fs << "Mediana: " << bas.get_midpoint(vo, nusim) << endl;
    fs << "Q1: " <<  bas.get_q1(vo, nusim) << endl;
    fs << "Q3: " <<  bas.get_q3(vo, nusim) << endl << endl;
    
    
    fs << "* 5pc:\n";
    fs << "gana B: " << ABganaB5 << endl;
    fs << "gana C: " << ABganaC5 << endl;
    fs << "gana ignorancia: " << ABignorancia5 << endl;
    fs << "Promedio tiempo: " << bas.get_mean(ABt5, nusim) << endl;
    fs << "Desv std: " << bas.get_sample_stddev(ABt5, nusim) << endl;
    bas.sort(ABt5, vo, nusim);
    fs << "Mediana: " << bas.get_midpoint(vo, nusim) << endl;
    fs << "Q1: " <<  bas.get_q1(vo, nusim) << endl;
    fs << "Q3: " <<  bas.get_q3(vo, nusim) << endl << endl;
    
    
    fs << "* 50pc:\n";
    fs << "gana B: " << ABganaB50 << endl;
    fs << "gana C: " << ABganaC50 << endl;
    fs << "gana ignorancia: " << ABignorancia50 << endl;
    fs << "Promedio tiempo: " << bas.get_mean(ABt50, nusim) << endl;
    fs << "Desv std: " << bas.get_sample_stddev(ABt50, nusim) << endl;
    bas.sort(ABt50, vo, nusim);
    fs << "Mediana: " << bas.get_midpoint(vo, nusim) << endl;
    fs << "Q1: " <<  bas.get_q1(vo, nusim) << endl;
    fs << "Q3: " <<  bas.get_q3(vo, nusim) << endl << endl;
    
    
    fs << "* 80pc:\n";
    fs << "gana B: " << ABganaB80 << endl;
    fs << "gana C: " << ABganaC80 << endl;
    fs << "gana ignorancia: " << ABignorancia80 << endl;
    fs << "Promedio tiempo: " << bas.get_mean(ABt80, nusim) << endl;
    fs << "Desv std: " << bas.get_sample_stddev(ABt80, nusim) << endl;
    bas.sort(ABt80, vo, nusim);
    fs << "Mediana: " << bas.get_midpoint(vo, nusim) << endl;
    fs << "Q1: " <<  bas.get_q1(vo, nusim) << endl;
    fs << "Q3: " <<  bas.get_q3(vo, nusim) << endl << endl << endl;
    
    
    fs.close();
    
    ///ELECCIONES
    bas.open_ofstream(fs, "Graphs/AganaCp.txt");
    fs << AganaCp << endl ;
    fs.close();
    bas.open_ofstream(fs, "Graphs/AganaBp.txt");
    fs << AganaBp << endl ;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/AganaC5.txt");
    fs << AganaC5 << endl ;
    fs.close();
    bas.open_ofstream(fs, "Graphs/AganaB5.txt");
    fs << AganaB5 << endl ;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/AganaC50.txt");
    fs << AganaC50 << endl ;
    fs.close();
    bas.open_ofstream(fs, "Graphs/AganaB50.txt");
    fs << AganaB50 << endl ;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/ABganaCp.txt");
    fs << ABganaCp << endl ;
    fs.close();
    bas.open_ofstream(fs, "Graphs/ABganaBp.txt");
    fs << ABganaBp << endl ;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/ABganaC5.txt");
    fs << ABganaC5 << endl ;
    fs.close();
    bas.open_ofstream(fs, "Graphs/ABganaB5.txt");
    fs << ABganaB5 << endl ;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/ABganaC50.txt");
    fs << ABganaC50 << endl ;
    fs.close();
    bas.open_ofstream(fs, "Graphs/ABganaB50.txt");
    fs << ABganaB50 << endl ;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/ABganaC80.txt");
    fs << ABganaC80 << endl ;
    fs.close();
    bas.open_ofstream(fs, "Graphs/ABganaB80.txt");
    fs << ABganaB80 << endl ;
    fs.close();
    
    //
    bas.open_ofstream(fs, "Graphs/Atp.txt");
    for (i = 0; i < nusim; i++)
    fs << Atp[i] << endl;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/At5.txt");
    for (i = 0; i < nusim; i++)
    fs << At5[i] << endl;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/At50.txt");
    for (i = 0; i < nusim; i++)
    fs << At50[i] << endl;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/ABtp.txt");
    for (i = 0; i < nusim; i++)
    fs << ABtp[i] << endl;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/ABt5.txt");
    for (i = 0; i < nusim; i++)
    fs << ABt5[i] << endl;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/ABt50.txt");
    for (i = 0; i < nusim; i++)
    fs << ABt50[i] << endl;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/ABt80.txt");
    for (i = 0; i < nusim; i++)
    fs << ABt80[i] << endl;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/pacajaprimer.txt");
    for (i = 0; i < nusim; i++)
    fs << Atp[i] << "\tA\n" ;
    for (i = 0; i < nusim; i++)
    fs << ABtp[i] << "\tAB\n" ;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/pacaja5.txt");
    for (i = 0; i < nusim; i++)
    fs << At5[i] << "\tA\n" ;
    for (i = 0; i < nusim; i++)
    fs << ABt5[i] << "\tAB\n" ;
    fs.close();
    
    bas.open_ofstream(fs, "Graphs/pacaja50.txt");
    for (i = 0; i < nusim; i++)
    fs << At50[i] << "\tA\n" ;
    for (i = 0; i < nusim; i++)
    fs << ABt50[i] << "\tAB\n" ;
    fs.close();
    
    
    
    
    bas.open_ofstream(fs, "Graphs/paysfiga.tex");
    fs << "\\documentclass{standalone}\n";
    fs << "\\usepackage{tikz}\n";
    fs << "\\usetikzlibrary{arrows,decorations.pathmorphing,decorations.shapes,backgrounds,fit,snakes,shapes,positioning}\n";
    fs << "\\usepackage{times}\n";
    fs << "\\usepackage[T1]{fontenc}\n";
    fs << "\\renewcommand{\\familydefault}{\\sfdefault}\n";
    fs << "\\begin{document}\n";
    fs << "\\begin{tikzpicture}[rectangle, minimum width=0.33cm, minimum height=0.333cm, bla/.style={draw, fill=white},neg/.style={draw, fill=black}, gri/.style={draw, fill=black!50}, gena/.style={font=\\bf\\fontsize{8.04}{9.67}\\selectfont}, pac/.style={font=\\fontsize{9}{7.035} \\selectfont}]\n\n";
    
    fs << "\\node at (5,5) {First genetically assimilated optimum};\n";// Optimum first becoming a native phenotype};\n";
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
    system("pdflatex -output-directory Graphs/ Graphs/paysfiga.tex");
    
    ///50
    bas.open_ofstream(fs, "Graphs/pays50.tex");
    fs << "\\documentclass{standalone}\n";
    fs << "\\usepackage{tikz}\n";
    fs << "\\usetikzlibrary{arrows,decorations.pathmorphing,decorations.shapes,backgrounds,fit,snakes,shapes,positioning}\n";
    fs << "\\usepackage{times}\n";
    fs << "\\usepackage[T1]{fontenc}\n";
    fs << "\\renewcommand{\\familydefault}{\\sfdefault}\n";
    fs << "\\begin{document}\n";
    fs << "\\begin{tikzpicture}[rectangle, minimum width=0.33cm, minimum height=0.333cm, bla/.style={draw, fill=white},neg/.style={draw, fill=black}, gri/.style={draw, fill=black!50}, gena/.style={font=\\bf\\fontsize{8.04}{9.67}\\selectfont}, pac/.style={font=\\fontsize{9}{7.035} \\selectfont}]\n\n";
    
    fs << "\\node at (5,5) {Prevalent genetically assimilated optimum (50\\%)};\n";
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
    system("pdflatex -output-directory Graphs/ Graphs/pays50.tex");
    //
    
    int elpri;
    bas.open_ofstream(fs, "Graphs/paR.sh");
    
    //iNTRO
    fs << "a <- \"\\n\\n\\n***A founder:\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    //ap
    fs << "a <- \"**Probability of obtaining a greater deviation from expected (0.5) frequency of times that B appears first (two-tailed test), according to a binomial distribution:\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "a <- \"Number of succesful simulations: " << AganaBp+AganaCp << "\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    if (AganaBp > AganaCp)
        elpri = AganaBp;
    else
        elpri = AganaCp;
    fs << " bt <- binom.test(" << elpri << ", " << AganaBp+AganaCp << ", 0.5,alternative=\"two.sided\")\n";
    fs << "capture.output(bt, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
    //A5
    fs << "a <- \"\\n\\n**Probability of obtaining a greater deviation from expected (0.5) frequency of times that B evolves (5pc) (two-tailed test), according to a binomial distribution:\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "a <- \"Number of succesful simulations: " << AganaB5+AganaC5 << "\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    if (AganaB5 > AganaC5)
        elpri = AganaB5;
    else
        elpri = AganaC5;
    fs << " bt <- binom.test(" << elpri << ", " << AganaB5+AganaC5 << ", 0.5,alternative=\"two.sided\")\n";
    fs << "capture.output(bt, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
    
    //A50
    fs << "a <- \"\\n\\n**Probability of obtaining a greater deviation from expected (0.5) frequency of times that B evolves (50pc) (two-tailed test), according to a binomial distribution:\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "a <- \"Number of succesful simulations: " << AganaB50+AganaC50 << "\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    if (AganaB50 > AganaC50)
        elpri = AganaB50;
    else
        elpri = AganaC50;
    fs << " bt <- binom.test(" << elpri << ", " << AganaB50+AganaC50 << ", 0.5,alternative=\"two.sided\")\n";
    fs << "capture.output(bt, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
  
    
    //INTRO
    fs << "a <- \"\\n\\n\\n***AB founder:\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    //ABP
    fs << "a <- \"**Probability of obtaining a greater fraction of times that B appears first (two-tailed test), according to a binomial distribution:\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "a <- \"Number of succesful simulations: " << ABganaBp+ABganaCp << "\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    if (ABganaBp > ABganaCp)
        elpri = ABganaBp;
    else
        elpri = ABganaCp;
    fs << " bt <- binom.test(" << elpri << ", " << ABganaBp+ABganaCp << ", 0.5,alternative=\"two.sided\")\n";
    fs << "capture.output(bt, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
    //AB5
    fs << "a <- \"\\n\\n**Probability of obtaining a greater fraction of times that B evolves (two-tailed test), according to a binomial distribution:\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "a <- \"Number of succesful simulations: " << ABganaB5+ABganaC5 << "\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    if (ABganaB5 > ABganaC5)
        elpri = ABganaB5;
    else
        elpri = ABganaC5;
    fs << " bt <- binom.test(" << elpri << ", " << ABganaB5+ABganaC5 << ", 0.5,alternative=\"two.sided\")\n";
    fs << "capture.output(bt, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
    //AB50
    fs << "a <- \"\\n\\n**Probability of obtaining a greater fraction of times that B evolves (two-tailed test), according to a binomial distribution:\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "a <- \"Number of succesful simulations: " << ABganaB50+ABganaC50 << "\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    if (ABganaB50 > ABganaC50)
        elpri = ABganaB50;
    else
        elpri = ABganaC50;
    fs << " bt <- binom.test(" << elpri << ", " << ABganaB50+ABganaC50 << ", 0.5,alternative=\"two.sided\")\n";
    fs << "capture.output(bt, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
    
    
    fs << "library(ggplot2)\n";
    fs << "a <- read.csv(file = \"Graphs/pacajaprimer" << ".txt\", head=FALSE, sep =\"\\t\")\n";
    fs << "a$V2 <- factor(a$V2)\n";
    fs << "pdf(file=\"Graphs/tfirstGA.pdf\", width=7.5,height=5.32)\n";
    fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
    fs << "ggplot(a, aes(x=V2, y=V1)) + geom_point(size=1, shape=16, alpha=0.2, colour=\"black\", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour=\"grey20\")+labs(x=expression(\'Founder\'), y = expression(\'First genetically assimilated optimum\')) + stat_summary(geom = \"point\", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))\n";
    fs << "dev.off()\n";
    
    fs << "a <- read.csv(file = \"Graphs/pacaja5" << ".txt\", head=FALSE, sep =\"\\t\")\n";
    fs << "a$V2 <- factor(a$V2)\n";
    fs << "pdf(file=\"Graphs/tuntil5.pdf\", width=7.5,height=5.32)\n";
    fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
    fs << "ggplot(a, aes(x=V2, y=V1)) + geom_point(size=1, shape=16, alpha=0.2, colour=\"black\", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour=\"grey20\")+labs(x=expression(\'Founder\'), y = expression(\'Time until prevalence (0.05) of assimilated optimum\')) + stat_summary(geom = \"point\", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))\n";
    fs << "dev.off()\n";
    
    fs << "a <- read.csv(file = \"Graphs/pacaja50" << ".txt\", head=FALSE, sep =\"\\t\")\n";
    fs << "a$V2 <- factor(a$V2)\n";
    fs << "pdf(file=\"Graphs/tuntil50.pdf\", width=7.5,height=5.32)\n";
    fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
    fs << "ggplot(a, aes(x=V2, y=V1)) + geom_point(size=1, shape=16, alpha=0.2, colour=\"black\", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour=\"grey20\")+labs(x=expression(\'Founder\'), y = expression(\'Time until prevalence (0.50) of assimilated optimum\')) + stat_summary(geom = \"point\", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))\n";
    fs << "dev.off()\n";
    
    
    fs << "equis <- read.csv(file=\"Graphs/Atp.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
    fs << "ye <- read.csv(file=\"Graphs/ABtp.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
    fs << "wtest <- wilcox.test(equis$V1, ye$V1, paired=FALSE, alternative=\"greater\")\n";
    fs << "a <- \"Mann-Whitney for time of appearance (one-tailed; time when founder A > time when founder AB):\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "capture.output(wtest, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
    
    fs << "equis <- read.csv(file=\"Graphs/At5.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
    fs << "ye <- read.csv(file=\"Graphs/ABt5.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
    fs << "wtest <- wilcox.test(equis$V1, ye$V1, paired=FALSE, alternative=\"greater\")\n";
    fs << "a <- \"Mann-Whitney for time until prevalence (0.05) (one-tailed; time when founder A > time when founder AB):\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "capture.output(wtest, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
    fs << "equis <- read.csv(file=\"Graphs/At50.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
    fs << "ye <- read.csv(file=\"Graphs/ABt50.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
    fs << "wtest <- wilcox.test(equis$V1, ye$V1, paired=FALSE, alternative=\"greater\")\n";
    fs << "a <- \"Mann-Whitney for time until prevalence (0.50) (one-tailed; time when founder A > time when founder AB):\"\n";
    fs << "capture.output(a, file=\"Graphs/stats.txt\", append=TRUE)\n";
    fs << "capture.output(wtest, file=\"Graphs/stats.txt\", append=TRUE)\n";
    
    
    fs.close();
    
    bas.run_command("Rscript Graphs/paR.sh");
    
    return 0;
}

