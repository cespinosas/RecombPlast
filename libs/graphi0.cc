#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "graphi0.h"
#include <string>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ifstream;
using std::istream;
using std::string;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;

GraphI0::GraphI0()
{
  ionary = false;
}

GraphI0::GraphI0(Alea& jacta)
{
  start_rng(jacta);
  ionary = false;
}

void GraphI0::start_rng(Alea& jacta)
{
  est = jacta;
  basic.start_rng(est);
}

GraphI0::GraphI0(int n, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, dirnw);
}

void GraphI0::make_nw(int n, bool dirnw)
{
  int i;
  size = n;
  nw = new int*[size];
  for (i = 0; i < size; i++) 
		nw[i] = new int[size];
  basic.fillmat0(nw, size, size);
  ionary = false;
  directed = dirnw;
  set_default_exclusive_vars();
  edo = new int[size];
  ima = new int[size];
}

void GraphI0::copy(GraphI0 &templ)
{
  int n = templ.number_of_nodes();
	bool dirnw = templ.is_directed();
  make_nw(n, dirnw);
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      nw[i][j] = templ.weight(j, i);
}

GraphI0::GraphI0(int n, double c, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, c, dirnw);
}

void GraphI0::make_nw(int n, double c, bool dirnw)
{
  if (c >= 1) {
		cout << "[Error]: You are trying to construct an overconnected network. GraphI0::make_nw.\n";
		exit(1);
	}
  int i, j,newval;
  double nure;
  bool pos;
  make_nw(n, dirnw);
	if (directed) {
		for (i = 0; i < size; i++)
			for(j = 0; j < size; j++) {
				nure = est.randreal();
				if (nure <= c) {
					pos = est.toss();
					if (pos)
            newval = 1;
          else
            newval = -1;
          change_interaction(j,i,newval);
				}
			}
	}
	else {
		for (i=0; i<size; i++)
			for (j=0; j<=i; j++) {
				nure = est.randreal();
				if (nure <= c) {
					pos = est.toss();
					if (pos)
            newval = 1;
          else
            newval = -1;
          change_interaction_undir(j,i,newval);
				}
			}
	}
}

GraphI0::GraphI0(int n, int e, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, e, dirnw);
}

void GraphI0::make_nw(int n, int e, bool dirnw)
{
  if (e > (n*n)) {
		cout << "[Error]: You are trying to construct an overconnected network. GraphI0::make_nw.\n";
		exit(1);
	}
  int i, j, k,newval;
  make_nw(n, dirnw);
  bool pos;
  for (i = 0; i < e; i++) {
		do {
			j = est.randint(0, size);
			k = est.randint(0, size);
		} while (nw[j][k] != 0);
		pos = est.toss();
		if (pos)
			newval = 1;
		else 
			newval = -1;
    if (directed)
      change_interaction(k,j,newval);
    else
      change_interaction_undir(k,j,newval);
	}
}


int GraphI0::rnd_nw_with_2fp1tn1(int n, int e, int *fp1, int *fp2, int *ic1) {
  int conta = 0;
  bool otra = false;
  int i,j,k;
  int *unvec;
  unvec = new int[n];
  for (i = 0; i < n; i++)
    unvec[i] = ic1[i];
  
  int **quiecua, *cuandecual;
  /*    ci  a   b   nomb
   0:   -   -   -   T-
   1:   -   -   +   B+
   2:   -   +   -   A+
   3:   -   +   +   C-
   4:   +   -   -   C+
   5:   +   -   +   A-
   6:   +   +   -   B-
   7:   +   +   +   T+
   */
  basic.create_array(quiecua, 8, n);
  basic.fillmatm1(quiecua, 8, n);
  cuandecual = new int[8];
  basic.fillv0(cuandecual, 8);
  for (i=0; i < n; i++) {
    j = (((ic1[i]+1)/2)*4) + (((fp1[i]+1)/2)*2) + ((fp2[i]+1)/2);
    quiecua[j][cuandecual[j]] = i;
    cuandecual[j]++;
  }
  int swi;
  int vanin=0;
  make_nw(n,true);
  for (i = 0; i < cuandecual[7]; i++) {
    for (j=0; j < cuandecual[7]; j++)
      force_interaction(quiecua[7][i], quiecua[7][j], 1);
    for (j=0; j < cuandecual[4]; j++)
      force_interaction(quiecua[7][i], quiecua[4][j], -1);
    for (j=0; j < cuandecual[3]; j++)
      force_interaction(quiecua[7][i], quiecua[3][j], 1);
  }
  vanin += cuandecual[7]*(cuandecual[7]+cuandecual[4]+cuandecual[3]);
  if (cuandecual[1] > 0) {
    swi =quiecua[1][0];
    for (i = 0; i < cuandecual[1]; i++)
      force_interaction(swi, quiecua[1][i], 1);
    for (i = 0; i < cuandecual[5]; i++)
      force_interaction(swi, quiecua[5][i], 1);
    for (i = 0; i < cuandecual[2]; i++)
      force_interaction(swi, quiecua[2][i], -1);
    for (i = 0; i < cuandecual[6]; i++)
      force_interaction(swi, quiecua[6][i], -1);
    vanin+= (cuandecual[1]+cuandecual[5]+cuandecual[2]+cuandecual[6]);
  }
  if (cuandecual[2] > 0) {
    swi = quiecua[2][0];
    for (i = 1; i < cuandecual[2]; i++)
      force_interaction(swi, quiecua[2][i], 1);
    for (i = 0; i < cuandecual[5]; i++)
      force_interaction(swi, quiecua[5][i], -1);
    vanin += ((cuandecual[2]-1) + cuandecual[5]);
    if (cuandecual[7] > 0) {
      force_interaction(quiecua[7][0], swi, 1);
      vanin++;
    }
    if (cuandecual[5] > 0) {
      force_interaction(quiecua[5][0], swi, -1);
      vanin++;
    }
    if (cuandecual[6] > 0) {
      force_interaction(quiecua[6][0], swi, 1);
      vanin++;
    }
  }
  for (i = vanin; i < e ; i++) { //7,0,-1 ; 0,X,X
    conta = 0;
    if (est.toss()) {
      do {
        j = est.randint(0, cuandecual[7]);
        k = est.randint(0, cuandecual[0]);
        conta++;
      } while((weight(quiecua[7][j], quiecua[0][k]) != 0) && (conta < 20));
      if (conta < 20) {
        force_interaction(quiecua[7][j], quiecua[0][k], (-1));
      }
      
    } else {
      do {
        j = est.randint(0, cuandecual[0]);
        k = est.randint(0, n);
        conta++;
      } while((weight(quiecua[0][j], k) != 0) && (conta < 20));
      if (conta < 20) {
        if (est.toss())
          force_interaction(quiecua[0][j], k, (-1));
        else
          force_interaction(quiecua[0][j], k, 1);
      }
    }
    if (conta >= 20)
      i--;
  }
  set_as_state(ic1);
  find_an_attractor();
  conta = 1;
  if (atsize != 1)
    otra = true;
  else {
    if (!basic.eqvec(attractor[0], size, fp1, size))
      otra = true;
  }
  if (!otra) {
    cout << "fenotipo A desde ci, bien\n";
    otra = true;
    for (i = 0; i < n; i++) {
      clear_attractor();
      unvec[i]*=(-1);
      set_as_state(unvec);
      find_an_attractor();
      if (atsize == 1) {
        if (basic.eqvec(attractor[0], size, fp2, size)) {
          otra = false;
          cout << "ahorasi\n";
        }
      }
      unvec[i]*=(-1);
      if (!otra)
        break;
    }
  }
  if (otra) {
    conta = -1;
    cout << "[Error]: The program did not produce the expected network.\n";
    exit(1);
  }
  return conta;
}

void GraphI0::rwalk_1traj(int n, int e, int *fp, int *ic, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  rnd_nw_with_2fp1tn1(n, e, fp, fp, ic);
  clear_attractor();
  int h,i,j,k,l,sig, vval, numutr,mr, enarch;
  ofstream fs;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphI0::rwalk_1fp.\n";
    exit(1);
  }
  enarch = cuan/numarch;
  bool sube, estasi;
  for (h = 0; h < numarch; h++) {
    basic.open_ofstream(fs, nomarch+basic.inttostring(h+archdesde)+".txt");
    for (i=0; i<enarch; i++) {
      do {
        mr = n*n*numut;
        numutr = est.randint(int(mr*8/10), int(mr*12/10));
        for (j=0; j<numutr; j++) {
          if (number_of_edges() == e+mame)
            sube = false;
          else {
            if (number_of_edges() == e-mame)
              sube = true;
            else {
              sube = false;
              if (est.toss())
                sube = true;
            }
          }
          if (!sube) {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)==0);
            vval = weight(k,l);
            force_interaction(k,l,0);
            set_as_state(ic);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp, size)))
              force_interaction(k,l,vval);
          } else {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)!=0);
            vval = weight(k,l);
            if (est.toss())
              sig = 1;
            else
              sig = -1;
            force_interaction(k,l,sig);
            set_as_state(ic);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp, size)))
              force_interaction(k,l,vval);
          }
        }
        estasi = true;
        if (estri)
          if (number_of_edges() != e)
            estasi = false;
      } while(!estasi);
      export_nw(fs);
      fs << "*\n";
    }
    fs.close();
  }
}

void GraphI0::rwalk_2fp1t(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  rnd_nw_with_2fp1tn1(n, e, fp1, fp2, ic1);
  clear_attractor();
  int h,i,j,k,l,sig, vval, numutr,mr, enarch;
  ofstream fs;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphI0::rwalk_2fp1t.\n";
    exit(1);
  }
  enarch = cuan/numarch;
  bool sube, estasi;
  for (h = 0; h < numarch; h++) {
    basic.open_ofstream(fs, nomarch+basic.inttostring(h+archdesde)+".txt");
    for (i=0; i<enarch; i++) {
      do {
        mr = n*n*numut;
        numutr = est.randint(int(mr*8/10), int(mr*12/10));
        for (j=0; j<numutr; j++) {
          if (number_of_edges() == e+mame)
            sube = false;
          else {
            if (number_of_edges() == e-mame)
              sube = true;
            else {
              sube = false;
              if (est.toss())
                sube = true;
            }
          }
          if (!sube) {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)==0);
            vval = weight(k,l);
            force_interaction(k,l,0);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              set_as_state(fp2);
              find_an_attractor();
              if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp2, size)))
                force_interaction(k,l,vval);
              else
                clear_attractor();
            }
          } else {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)!=0);
            vval = weight(k,l);
            if (est.toss())
              sig = 1;
            else
              sig = -1;
            force_interaction(k,l,sig);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              set_as_state(fp2);
              find_an_attractor();
              if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp2, size)))
                force_interaction(k,l,vval);
              else
                clear_attractor();
            }
          }
        }
        estasi = true;
        if (estri)
          if (number_of_edges() != e)
            estasi = false;
      } while(!estasi);
      export_nw(fs);
      fs << "*\n";
    }
    fs.close();
  }
}

void GraphI0::rwalk_2fp1tn1(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  rnd_nw_with_2fp1tn1(n, e, fp1, fp2, ic1);
  clear_attractor();
  int h,i,j,k,l,sig, vval, numutr,mr, enarch,ii;
  int *ic2;
  ic2 = new int[n];
  ofstream fs;
  bool ucon;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphI0::rwalk_2fp1tn1.\n";
    exit(1);
  }
  enarch = cuan/numarch;
  bool sube, estasi;
  for (h = 0; h < numarch; h++) {
    basic.open_ofstream(fs, nomarch+basic.inttostring(h+archdesde)+".txt");
    for (i=0; i<enarch; i++) {
      do {
        mr = n*n*numut;
        numutr = est.randint(int(mr*8/10), int(mr*12/10));
        for (j=0; j<numutr; j++) {
          if (number_of_edges() == e+mame)
            sube = false;
          else {
            if (number_of_edges() == e-mame)
              sube = true;
            else {
              sube = false;
              if (est.toss())
                sube = true;
            }
          }
          if (!sube) {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)==0);
            vval = weight(k,l);
            force_interaction(k,l,0);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              for (ii=0; ii<size; ii++)
                ic2[ii] = ic1[ii];
              ucon = false;
              for (ii=0; ii<size; ii++) {
                ic2[ii]*=(-1);
                set_as_state(ic2);
                find_an_attractor();
                ic2[ii]*=(-1);
                if ((attractor_size() == 1) && (basic.eqvec(attractor[0], size, fp2, size))) {
                  clear_attractor();
                  ucon = true;
                  break;
                }
                clear_attractor();
              }
              if (!ucon)
                force_interaction(k,l,vval);
            }
          } else {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)!=0);
            vval = weight(k,l);
            if (est.toss())
              sig = 1;
            else
              sig = -1;
            force_interaction(k,l,sig);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              for (ii=0; ii<size; ii++)
                ic2[ii] = ic1[ii];
              ucon = false;
              for (ii=0; ii<size; ii++) {
                ic2[ii]*=(-1);
                set_as_state(ic2);
                find_an_attractor();
                ic2[ii]*=(-1);
                if ((attractor_size() == 1) && (basic.eqvec(attractor[0], size, fp2, size))) {
                  clear_attractor();
                  ucon = true;
                  break;
                }
                clear_attractor();
              }
              if (!ucon)
                force_interaction(k,l,vval);
            }
          }
        }
        estasi = true;
        if (estri)
          if (number_of_edges() != e)
            estasi = false;
      } while(!estasi);
      export_nw(fs);
      fs << "*\n";
    }
    fs.close();
  }
  delete [] ic2;
}

void GraphI0::rwalk_2fp1tn2(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  rnd_nw_with_2fp1tn1(n, e, fp1, fp2, ic1);
  clear_attractor();
  int h,i,j,k,l,sig, vval, numutr,mr, enarch,ii,jj;
  int *ic2;
  ic2 = new int[n];
  ofstream fs;
  bool ucon;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphI0::rwalk_2fp1tn2.\n";
    exit(1);
  }
  enarch = cuan/numarch;
  bool sube, estasi;
  for (h = 0; h < numarch; h++) {
    basic.open_ofstream(fs, nomarch+basic.inttostring(h+archdesde)+".txt");
    for (i=0; i<enarch; i++) {
      do {
        mr = n*n*numut;
        numutr = est.randint(int(mr*8/10), int(mr*12/10));
        for (j=0; j<numutr; j++) {
          if (number_of_edges() == e+mame)
            sube = false;
          else {
            if (number_of_edges() == e-mame)
              sube = true;
            else {
              sube = false;
              if (est.toss())
                sube = true;
            }
          }
          if (!sube) {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)==0);
            vval = weight(k,l);
            force_interaction(k,l,0);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              for (ii=0; ii<size; ii++)
                ic2[ii] = ic1[ii];
              ucon = false;
              for (ii=0; ii<size; ii++) {
                ic2[ii]*=(-1);
                for (jj = ii; jj < size; jj++) {
                  if (jj != ii)
                    ic2[jj] *= (-1);
                  set_as_state(ic2);
                  find_an_attractor();
                  if (jj != ii)
                    ic2[jj] *= (-1);
                  if ((attractor_size() == 1) && (basic.eqvec(attractor[0], size, fp2, size))) {
                    clear_attractor();
                    ucon = true;
                    break;
                  }
                  clear_attractor();
                }
                ic2[ii]*=(-1);
                if (ucon)
                  break;
              }
              if (!ucon)
                force_interaction(k,l,vval);
            }
          } else {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)!=0);
            vval = weight(k,l);
            if (est.toss())
              sig = 1;
            else
              sig = -1;
            force_interaction(k,l,sig);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              for (ii=0; ii<size; ii++)
                ic2[ii] = ic1[ii];
              ucon = false;
              for (ii=0; ii<size; ii++) {
                ic2[ii]*=(-1);
                for (jj = ii; jj < size; jj++) {
                  if (jj != ii)
                    ic2[jj] *= (-1);
                  set_as_state(ic2);
                  find_an_attractor();
                  if (jj != ii)
                    ic2[jj] *= (-1);
                  if ((attractor_size() == 1) && (basic.eqvec(attractor[0], size, fp2, size))) {
                    clear_attractor();
                    ucon = true;
                    break;
                  }
                  clear_attractor();
                }
                ic2[ii]*=(-1);
                if (ucon)
                  break;
              }
              if (!ucon)
                force_interaction(k,l,vval);
            }
          }
        }
        estasi = true;
        if (estri)
          if (number_of_edges() != e)
            estasi = false;
      } while(!estasi);
      export_nw(fs);
      fs << "*\n";
    }
    fs.close();
  }
  delete [] ic2;
}

void GraphI0::clear()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] nw[i];
  delete [] nw;
  delete [] ima;
  delete [] edo;
	
  if (ionary)
    clear_dict();
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yadeg)
    clear_deg();
  if (yadegdist)
		clear_degdist();
  if (yatra)
    clear_attractor();
	if (yae)
		clear_edge_count();
  if (yatras)
    clear_attractors();
  if (yamadya)
    clear_adj();
}

void GraphI0::clear_edge_count()
{
	numofe = 0;
	yae = false;
}

void GraphI0::clear_loops()
{
  delete [] loops;
  totalloops = 0;
  int i;
  for (i = 0; i < size; i++)
    delete [] guada[i];
  delete [] guada;
  floops = false;
}

void GraphI0::clear_comp()
{
  delete [] components;
  yacomp = false;
}

void GraphI0::clear_bc()
{
  delete [] bc;
  bcya = false;
}

void GraphI0::clear_scc()
{
  delete [] scc;
  cfc = false;
  numbscc = 0;
}

void GraphI0::clear_adj() {
  int i;
  for (i = 0; i < size; i++)
    delete [] matadya[i];
  delete [] matadya;
  yamadya = false;
}

void GraphI0::clear_attractor()
{
  int i;
  for (i = 0; i < atsize; i++)
    delete [] attractor[i];
  delete [] attractor;
  atsize = 0;
  palen = 0;
  yatra = false;
}

void GraphI0::clear_attractors()
{
  int i,j;
  for (i=0; i < numatrs; i++) {
    for (j=0; j< atsizes[i]; j++)
      delete [] attractors[i][j];
    delete [] attractors[i];
  }
  delete [] attractors;
  delete [] atsizes;
  numatrs = 0;
  yatras = false;
}

void GraphI0::clear_dima()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] dima[i];
  delete [] dima;
  dist = false;
}

void GraphI0::clear_ds()
{
  int i;
  for (i = 0; i < size; i++)
    downstream[i].clear();
  delete [] downstream;
  infl = false;
}

void GraphI0::clear_dict()
{
  delete [] dict;
  ionary = false;
}

void GraphI0::clear_deg()
{
  delete [] degree;
  delete [] outdegree;
  delete [] indegree;
  yadeg = false;
	if (yadegdist)
		clear_degdist();
}

void GraphI0::clear_degdist()
{
	delete [] degdist;
	delete [] odegdist;
	delete [] idegdist;
	yadegdist = false;
}

//access
int GraphI0::number_of_edges()
{
  if (!yae) {
		numofe = 0;
		int i, j;
		if (directed) {
			for (i = 0; i < size; i++)
				for (j = 0; j < size; j++)
					if (nw[i][j] != 0)
						numofe++;
		}
		else {
			for (i = 0; i < size; i++)
				for (j = 0; j <= i; j++)
					if (nw[i][j] != 0)
						numofe++;
		}
		yae = true;
	}
  return numofe;
}

int GraphI0::number_of_nodes()
{
  return size;
}

int GraphI0::weight(int source, int target)
{
  return nw[target][source];
}

bool GraphI0::is_directed()
{
  return directed;
}

bool GraphI0::attractor_exists()
{
  return yatra;
}

bool GraphI0::attractors_exist()
{
  return yatras;
}

//modification
void GraphI0::force_interaction(int source, int target, int value)
{
	if (directed) {
		if (yae) {
			if ((value == 0) && (nw[target][source] != 0))
				numofe--;
			else 
				if ((value != 0) && (nw[target][source] == 0))
					numofe++;
		}
		if (yadeg) {
			if ((value == 0) && (nw[target][source] != 0)) {
				outdegree[source]--;
				indegree[target]--;
				degree[source]--;
				degree[target]--;
			}
			else {
				if ((value != 0) && (nw[target][source] == 0)) {
					outdegree[source]++;
					indegree[target]++;
					degree[source]++;
					degree[target]++;
				}
			}
		}
		nw[target][source] = value;
	}
	else {
		cout << "[Error]: GraphI0::force_interaction does not work for undirected graphs.\n";
		exit(1);
  }
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yamadya)
    clear_adj();
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  if (yadegdist)
    clear_degdist();
}

void GraphI0::change_interaction(int source, int target, int value)
{
  if (directed) {
    if (value == nw[target][source]) {
      cout << "[Error]: GraphI0::change_interaction changes to stay the same.\n";
      exit(1);
    }
    if (yae) {
      if ((value == 0) && (nw[target][source] != 0))
        numofe--;
      else
        if ((value != 0) && (nw[target][source] == 0))
          numofe++;
    }
    if (yadeg) {
      if ((value == 0) && (nw[target][source] != 0)) {
        outdegree[source]--;
        indegree[target]--;
        degree[source]--;
        degree[target]--;
      }
      else {
        if ((value != 0) && (nw[target][source] == 0)) {
          outdegree[source]++;
          indegree[target]++;
          degree[source]++;
          degree[target]++;
        }
      }
    }
    nw[target][source] = value;
  }
  else {
    cout << "[Error]: GraphI0::change_interaction does not work for undirected graphs.\n";
    exit(1);
  }
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yamadya)
    clear_adj();
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  if (yadegdist)
    clear_degdist();
}

void GraphI0::change_interaction_undir(int source, int target, int value)
{
  if (directed) {
    cout << "[Error]: GraphI0::change_interaction_undir does not work for directed graphs.\n";
    exit(1);
  }
  if (value == nw[target][source]) {
    cout << "[Error]: GraphI0::change_interaction changes to stay the same.";
    exit(1);
  }
  if (yae) {
    if ((value==0) && (nw[target][source] != 0))
      numofe--;
    else
      if ((value != 0) && (nw[target][source] == 0))
        numofe++;
  }
  if (yadeg) {
    if ((value==0) && (nw[target][source] != 0)) {
      outdegree[source]--;
      outdegree[target]--;
      indegree[target]--;
      indegree[source]--;
      degree[source]--;
      degree[target]--;
    }
    else
      if ((value!=0) && (nw[target][source] == 0)) {
        outdegree[source]++;
        outdegree[target]++;
        indegree[target]++;
        indegree[source]++;
        degree[source]++;
        degree[target]++;
      }
  }
  nw[target][source] = value;
  nw[source][target] = value;
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yamadya)
    clear_adj();
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  if (yadegdist)
    clear_degdist();
}
//exclusive
//import
void GraphI0::get_dir_nw_from_file(int nn, istream& en, int e)
{
  make_nw(nn, true);
  int i, j, k, l;
  for (l = 0; l < e; l++) {
    en >> j;
    en >> i;
    en >> k;
    if (((k==1) || (k==(-1))) && (i < size) && (j < size) && (i >= 0) && (j >= 0))
      nw[i][j] = k;
    else {
      cout << "[Error]: This is not an int network. GraphI0::get_dir_nw_from_file.\n";
      exit(1);
    }
  }
  set_default_exclusive_vars();
}

void GraphI0::get_dir_nw_from_file(int nn, string arch)
{
	ifstream sal;
	basic.open_ifstream(sal, arch);
	get_dir_nw_from_file(nn, sal);	
	sal.close();
}

void GraphI0::get_dir_nw_from_file(int nn, istream& en)
{
  make_nw(nn, true);
  int i, j, k;
  while (en >> j) {
		en >> i;
		en >> k;
		if (((k==1) || (k==(-1))) && (i < size) && (j < size) && (i >= 0) && (j >= 0))
			nw[i][j] = k;
		else {
			cout << "[Error]: This is not an int network. GraphI0::get_dir_nw_from_file.\n";
			exit(1);
		}
	}
  set_default_exclusive_vars();
}

//Export
void GraphI0::export_nw(string arch) {
  ofstream fs;
  basic.open_ofstream(fs, arch);
  export_nw(fs);
  fs.close();
}

void GraphI0::export_nw(ostream& fs) {
  int i, j;
  if (directed) {
    for (i = 0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i) != 0)
          fs << j << "\t" << i << "\t" << weight(j, i) << endl;
  }
  else {
    for (i = 0; i < size; i++)
      for (j=0; j <= i; j++)
        if (weight(j,i) != 0)
          fs << j << "\t" << i << "\t" << weight(j, i) << endl;
  }
}

//status
//analyses
int GraphI0::distance_from_nw(GraphI0& otra) {
  int res = 0;
  int i, j;
  for (i=0; i < size; i++)
    for (j = 0; j < size; j++)
      res += abs(weight(i,j) - otra.weight(i,j));
  return res;
}

bool GraphI0::equal_adjmat(GraphI0& otra) {
  bool res = true;
  int i, j;
  if ((otra.number_of_nodes() != size) || (otra.number_of_edges() != number_of_edges()) || (otra.is_directed() != is_directed()))
    res = false;
  else {
    for (i = 0; i < size; i++) {
      if (!res)
        break;
      for (j = 0; j < size; j++) {
        if (((weight(i,j) == 0) && (otra.weight(i, j) != 0)) || ((weight(i,j) != 0) && (otra.weight(i, j) == 0))) {
          res = false;
          break;
        }
      }
    }
  }
  return res;
}


bool GraphI0::equal_nw(GraphI0 &templ)
{
  int i, j;
  bool res = true;
  if ((templ.number_of_nodes() != size) || (templ.is_directed() != is_directed()))
    res = false;
  else {
    for (i = 0; i < size; i++) {
			for (j=0; j<size; j++)
				if (templ.weight(j, i) != weight(j, i)) {
					res = false;
					break;
				}
			if (!res)
				break;
		}
  }
  return res;
}

void GraphI0::synchrony()
{
  int i, j, sum;
  for (i = 0; i < size; i++) {
		sum = 0;
		for (j = 0; j < size; j++) {
			if ((weight(j,i) != 0) && (edo[j] > 0))
        sum = sum + (edo[j]*weight(j,i));
		}
		if (sum > 0)
			ima[i] = 1;
		if (sum < 0)
			ima[i] = -1;
		if (sum == 0)
			ima[i] = edo[i]; //it was 0 before 2017
	}
  for (i = 0; i < size; i++)
    edo[i] = ima[i];
}

void GraphI0::set_as_state(int *vec)
{
  int i;
  for (i=0; i< size; i++)
    edo[i] = vec[i];
}

void GraphI0::find_an_attractor()
{
  if (yatra) {
    cout << "[Error]: Attractor already found when GraphI0::find_an_attractor was called.\n";
    exit(1);
  }
  int i, j, k, l;
  int maxsta = int(pow(2.0, size));
  int **trajectory;
  trajectory = new int*[maxsta];
  for (i = 0; i < maxsta; i++) 
		trajectory[i] = new int[size];
  k = 0;
  do {
    for (i = 0; i < size; i++)
      trajectory[k][i] = edo[i];
    synchrony();
    k++;
    if (basic.eqvec(edo, size, trajectory[k-1], size)) {
			l = k-1;
			break;
		}
    l = basic.vecinmat(trajectory, k, size, edo, size);
  } while (l==(-1));
  attractor = new int*[k-l];
  for (i= 0; i<(k-l); i++) 
		attractor[i] = new int[size];
  for (i = 0; i < (k-l); i++)
		for (j = 0; j < size; j++)
			attractor[i][j] = trajectory[i+l][j];
  atsize = k-l;
  palen = l;
  for (i = 0; i < maxsta; i++)
    delete [] trajectory[i];
  delete [] trajectory;
  yatra = true;
}

void GraphI0::find_attractors(int **ic, int nic)
{
  if (yatras) {
    cout << "[Error]: Attractors already found when GraphI0::find_attractors was called.\n";
    exit(1);
  }
  numatrs = nic;
  attractors = new int**[numatrs];
  atsizes = new int[numatrs];
  int i,j,k;
  for (i=0; i<nic; i++) {
    set_as_state(ic[i]);
    find_an_attractor();
    atsizes[i] = atsize;
    attractors[i] = new int*[atsizes[i]];
    for (j=0; j<atsizes[i]; j++)
      attractors[i][j] = new int[size];
    for (j=0; j<atsizes[i]; j++)
      for (k=0; k< size; k++)
        attractors[i][j][k] = attractor_element(j,k);
    clear_attractor();
  }
  yatras=true;
}

int GraphI0::attractor_size()
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphI0::attractor_size was called.\n";
    exit(1);
  }
  return atsize;
}

int GraphI0::attractor_element(int row, int node)
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphI0::attractor_element was called.\n";
    exit(1);
  }
  if ((row >=atsize)||(row < 0)) {
    cout << "[Error]: Attractor has no " << row << "th row. GraphI0::attractor_element.\n";
    exit(1);
  }
  if ((node >= size) || (node < 0)) {
    cout << "[Error]: Node " << node << " does not exist. GraphI0::attractor_element.\n";
  }
  return attractor[row][node];
}

void GraphI0::print_attractor(ostream& sal)
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphI0::print_attractor was called.\n";
    exit(1);
  }
  int i,j;
  for (i=0; i < atsize; i++) {
    for (j=0; j<size; j++)
      sal << attractor_element(i,j) << " ";
    sal << endl;
  }
}

int GraphI0::attractor_size(int wh)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphI0::attractor_size was called.\n";
    exit(1);
  }
  return atsizes[wh];
}

int GraphI0::attractor_element(int wh, int row, int node)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphI0::attractor_element was called.\n";
    exit(1);
  }
  if ((row >=atsizes[wh])||(row < 0)) {
    cout << "[Error]: Attractor has no " << row << "th row. GraphI0::attractor_element.\n";
    exit(1);
  }
  if ((node >= size) || (node < 0)) {
    cout << "[Error]: Node " << node << " does not exist. GraphI0::attractor_element.\n";
  }
  return attractors[wh][row][node];
}

void GraphI0::print_attractor(int wh, ostream& sal)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphI0::print_attractor was called.\n";
    exit(1);
  }
  int i,j;
  for (i=0; i < attractor_size(wh); i++) {
    for (j=0; j<size; j++)
      sal << attractor_element(wh,i,j) << " ";
    sal << endl;
  }
}

int GraphI0::number_of_attractors()
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphI0::number_of_attractors was called.\n";
    exit(1);
  }
  return numatrs;
}

void GraphI0::penetrance_n1_2p(int *cior, int *pA, int *pB, int& countA, int& countB) {
  int i;
  countA = 0;
  countB = 0;
  int *vt;
  vt = new int[size];
  for (i=0; i < size; i++)
    vt[i] = cior[i];
  for (i=0; i < size; i++) {
    vt[i]*=(-1);
    set_as_state(vt);
    find_an_attractor();
    if (attractor_size() == 1) {
      if (basic.eqvec(attractor[0], size, pA, size))
        countA++;
      else if (basic.eqvec(attractor[0], size, pB, size))
        countB++;
    }
    clear_attractor();
    vt[i]*=(-1);
  }
  return;
}

void GraphI0::penetrance_n2_ex_2p(int *cior, int *pA, int *pB, int& countA, int& countB) {
  countA = 0;
  countB = 0;
  int i,j;
  int *vt;
  vt = new int[size];
  for (i=0; i < size; i++)
    vt[i] = cior[i];
  for (i=0; i < (size-1); i++) {
    vt[i]*=(-1);
    for (j = i+1; j < size; j++) {
      vt[j] *=(-1);
      set_as_state(vt);
      find_an_attractor();
      if (attractor_size() == 1) {
        if (basic.eqvec(attractor[0], size, pA, size))
          countA++;
        else if (basic.eqvec(attractor[0], size, pB, size))
          countB++;
      }
      clear_attractor();
      vt[j] *=(-1);
    }
    vt[i]*=(-1);
  }
  return;
}

void GraphI0::mutational_access_2p(int *ci, int *pA, int *pB, int& countA, int& countB) {
  countA = 0;
  countB = 0;
  int i,j;
  int vval;
  for (i=0; i < size; i++) {
    for (j=0; j < size; j++) {
      vval = weight(i,j);
      if (vval!=0) {
        force_interaction(i,j,0);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1) {
          if (basic.eqvec(attractor[0], size, pA, size))
            countA++;
          else if (basic.eqvec(attractor[0], size, pB, size))
            countB++;
        }
        clear_attractor();
      } else {
        force_interaction(i,j,1);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1) {
          if (basic.eqvec(attractor[0], size, pA, size))
            countA++;
          else if (basic.eqvec(attractor[0], size, pB, size))
            countB++;
        }
        clear_attractor();
        force_interaction(i,j,-1);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1) {
          if (basic.eqvec(attractor[0], size, pA, size))
            countA++;
          else if (basic.eqvec(attractor[0], size, pB, size))
            countB++;
        }
        clear_attractor();
      }
      force_interaction(i,j,vval);
    }
  }
  return;
}

//for evolution
void GraphI0::build_offspring(GraphI0 &mother, GraphI0 &father, bool *frommom) {
  if (mother.number_of_nodes() != father.number_of_nodes()) {
    cout << "[Error]: Parents with different genome size. GraphI::build_offspring.\n";
    exit(1);
  }
  if ((!mother.is_directed()) || (!father.is_directed())) {
    cout << "[Error]: One parent is not a directed network. GraphI::build_offspring.\n";
    exit(1);
  }
  int i, j, n;
  n = father.number_of_nodes();
  make_nw(n, true);
  for (i=0; i < size; i++) {
    if (frommom[i]) {
      for (j= 0; j < size; j++)
        force_interaction(j,i,mother.weight(j,i));
    } else {
      for (j= 0; j < size; j++)
        force_interaction(j,i,father.weight(j,i));
    }
  }
  return;
}

void GraphI0::mate(GraphI0 &mother, GraphI0 &father)
{
  if (mother.number_of_nodes() != father.number_of_nodes()) {
		cout << "[Error]: Parents with different genome size. GraphI0::mate.\n";
		exit(1);
	}
	if (mother.is_directed() != father.is_directed()) {
		cout << "[Error]: One parent is a directed network while the other is an undirected network. GraphI0::mate.\n";
		exit(1);
	}
  int i, j;
  bool paoma;
  int n = father.number_of_nodes();
  make_nw(n, father.is_directed());
	if (directed) {
		for (i = 0; i < size; i++) {
			paoma = est.toss();
			if (paoma) {
				for (j = 0; j < size; j++)
          force_interaction(j,i,father.weight(j,i));
			}
			else {
				for (j = 0; j < size; j++)
          force_interaction(j,i,mother.weight(j,i));
			}
		}
	}
	else {
		for (i=0; i<size; i++)
			for (j=0; j<=i; j++) {
				paoma = est.toss();
				if (paoma)
          force_interaction(j,i,father.weight(j,i));
				else
          force_interaction(j,i,mother.weight(j,i));
			}
	}
}
	

void GraphI0::mutate(double con) {
  int i = est.randint(0, size);
  mutate_gene(i, con);
}

void GraphI0::mutate_gene(int gene, double con) {
  if (!directed) {
    cout << "[Error]: GraphI0::mutate_gene is not defined for undirected networks.\n";
    exit(1);
  }
  bool sube = false;
  int j;
  j = est.randint(0, size);
  int oldval = weight(j,gene);
  if (est.randreal() < con)
    sube = true;
  if (sube) {
    if (oldval == 0) {
      if (est.toss())
        force_interaction(j,gene,1);
      else
        force_interaction(j,gene,-1);
    }
  }
  else {
    if (oldval != 0)
      force_interaction(j, gene, 0);
  }
}

void GraphI0::consider_mutation(double muratepg, double con) {
  int i;
  for (i = 0; i < size; i++)
    if (est.randreal() < muratepg)
      mutate_gene(i, con);
}

//private

void GraphI0::set_default_exclusive_vars()
{
  infl = false;
  dist = false;
  yacomp = false;
  cfc = false;
  bcya = false;
  floops = false;
  yae = false;
  yadeg = false;
  yatra = false;
	yadegdist=false;
  yatras = false;
  yamadya = false;

}
