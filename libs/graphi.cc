#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "graphi.h"
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

GraphI::GraphI()
{
  ionary = false;
}

GraphI::GraphI(Alea& jacta) //creates instance of class GraphI and assigns rng to jacta
{
  start_rng(jacta);
  ionary = false;
}

void GraphI::start_rng(Alea& jacta) //assigns rng to jacta
{
  est = jacta;
  basic.start_rng(est);
}

GraphI::GraphI(int n, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, dirnw);
}

void GraphI::make_nw(int n, bool dirnw) //creates empty directed/undirected network
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

void GraphI::make_robust_nw(int n, int e, int *cia, int *pa, int *cib, int *pb) {
  int nurcom=0,nurdif=0,comd=0,difd=0; //counter for genes in each class
  int i,j,l;
  for (i=0; i < n; i++) {
    if (((cia[i]*cib[i])>0) && ((pa[i]*pb[i])>0) && ((cia[i]*pa[i])>0))
      nurcom++;
    if (((cia[i]*cib[i])<0) && ((pa[i]*pb[i])<0) && ((cia[i]*pa[i])>0))
      nurdif++;
    if (pa[i]==pb[i])
      comd++;
    else
      difd++;
  }
  if ((nurcom==0) || ((nurdif==0)&&(difd>0))) {
    cout << "[Error]: Impossible network for GraphI::make_robust_nw.\n";
    exit(1);
  }
  int *vrcom, *vrdif, mc=0, md=0, *igudown, *difdown, cid=0, cdd = 0;
  vrcom = new int[nurcom];
  vrdif = new int[nurdif];
  igudown = new int[comd];
  difdown = new int[difd];
  for (i=0; i < n; i++) {
    if (((cia[i]*cib[i])>0) && ((pa[i]*pb[i])>0) && ((cia[i]*pa[i])>0)) {
      vrcom[mc] = i;
      mc++;
    }
    if (((cia[i]*cib[i])<0) && ((pa[i]*pb[i])<0) && ((cia[i]*pa[i])>0)) {
      vrdif[md] = i;
      md++;
    }
    if ((pa[i]*pb[i])>0) {
      igudown[cid] = i;
      cid++;
    }
    if ((pa[i]*pb[i])<0) {
      difdown[cdd] = i;
      cdd++;
    }
  }
  int maxd = nurdif*difd;
  int maxc = nurcom*comd;
  int k;//,minmax, maxmax;
  bool usobc;
  if (maxd > maxc)
    usobc = false;
  else
    usobc = true;
  
  if ((comd != n) && (difd != n)) {
    if (nurdif < nurcom)
      k = nurdif;
    else
      k = nurcom;
  } else {
    if (nurdif==0)
      k = nurcom;
    else
      k = nurdif;    
  }
  while (e < (k*n)) {
    k--;
  }
  if (k < 1) {
    cout << "[Error]: Impossible network for GraphI::make_robust_nw.\n";
    exit(1);
  }
  int val, maxinbue, vanin, conta=0;
  bool yast = false;
  do {
    make_nw(n, true);
    for (l = 0; l < k; l++) {
      for (i=0; i < n; i++) {
        if ((pa[i]*pb[i])>0) {
          do {
            j = est.randint(0,nurcom);
          } while (weight(vrcom[j],i) != 0);
          if ((pa[vrcom[j]]*pa[i])>0)
            val = 1;
          else
            val = -1;
          force_interaction(vrcom[j],i,val);
        } else {
          do {
            j = est.randint(0,nurdif);
          } while (weight(vrdif[j],i) != 0);
          if ((pa[vrdif[j]]*pa[i])>0)
            val = 1;
          else
            val = -1;
          force_interaction(vrdif[j],i,val);
        }
      }
    }
    vanin = k*n;
    if (e < (maxd+maxc))
      maxinbue = e;
    else
      maxinbue = maxc+maxd;
    if (usobc) {
      for (i = vanin; i < maxinbue; i++) {
        do {
          j = est.randint(0,nurcom);
          l = est.randint(0,comd);
        } while (weight(vrcom[j],igudown[l]) != 0) ;
        if ((pa[vrcom[j]]*pa[igudown[l]]) > 0)
          val = 1;
        else
          val = -1;
        force_interaction(vrcom[j],igudown[l],val);
      }
    } else {
      for (i = vanin; i < maxinbue; i++) {
        do {
          j = est.randint(0,nurdif);
          l = est.randint(0,difd);
        } while (weight(vrdif[j],difdown[l]) != 0) ;
        if ((pa[vrdif[j]]*pa[difdown[l]]) > 0)
          val = 1;
        else
          val = -1;
        force_interaction(vrdif[j],difdown[l],val);
      }
    }
    vanin = maxinbue;
    for (i = vanin; i < e; i++) {
      do {
        j = est.randint(0, n);
        l = est.randint(0, n);
      } while (weight(j,l) != 0);
      if (est.toss())
        val = 1;
      else
        val = -1;
      force_interaction(j,l,val);
    }
    yast = true;
    set_as_state(cia);
    find_an_attractor();
    if (attractor_size() != 1)
      yast = false;
    else {
      if (!basic.eqvec(pa, n, attractor[0], n))
        yast = false;
      else {
        clear_attractor();
        set_as_state(cib);
        find_an_attractor();
        if (attractor_size() != 1)
          yast = false;
        else
          if (!basic.eqvec(pb, n, attractor[0], n))
            yast = false;
        clear_attractor();
      }
    }
    if (!yast)
      clear();
    cout << ++conta << endl;
    if (conta > 200) {
      cout << "[Error]: GraphI::make_robust_nw surpassed limit number of tries.\n";
      exit(1);
    }
  }while(!yast);
  delete [] vrcom;
  delete [] vrdif;
  delete [] igudown;
  delete [] difdown;
}

void GraphI::copy(GraphI &templ)
{
  int n = templ.number_of_nodes();
	bool dirnw = templ.is_directed();
  make_nw(n, dirnw);
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      nw[i][j] = templ.weight(j, i);
}

GraphI::GraphI(int n, double c, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, c, dirnw);
}

void GraphI::make_nw(int n, double c, bool dirnw)
{
  if (c >= 1) {
		cout << "[Error]: You are trying to construct an overconnected network. GraphI::make_nw.\n";
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

GraphI::GraphI(int n, int e, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, e, dirnw);
}

void GraphI::make_nw(int n, int e, bool dirnw)
{
  if (e > (n*n)) {
		cout << "[Error]: You are trying to construct an overconnected network. GraphI::make_nw.\n";
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

void GraphI::rwalk_1traj(int n, int e, int *fp, int *ic, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  make_robust_nw(n, e, ic, fp, ic, fp); //nuevo
  int h,i,j,k,l,sig, vval, numutr,mr, enarch;
  ofstream fs;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphI::rwalk_1traj.\n";
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

void GraphI::rwalk_2fp1t(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  make_robust_nw(n, e, ic1, fp1, fp2, fp2); //nuevo
  int h,i,j,k,l,sig, vval, numutr,mr, enarch;
  ofstream fs;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphI::rwalk_2fp1t.\n";
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

void GraphI::rwalk_2fp1tn1(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  
  int h,i,j,k,l,sig, vval, numutr,mr, enarch,ii;
  int *ic2;
  ic2 = new int[n];
  for (ii=0; ii<n; ii++) //nuevo
    ic2[ii] = ic1[ii];
  for (i=0; i < n; i++) {
    if ((ic1[i]==fp1[i]) && (fp1[i] != fp2[i])) {
      ic2[i] *= (-1);
      break;
    }
  }
  if (i==n) {
    cout << "[Error]: Impossible initial network for GraphI::rwalk_2fp1tn1.\n";
    exit(1);
  }
  make_robust_nw(n,e, ic1, fp1, ic2, fp2);//hastaaqui
  
  
  ofstream fs;
  bool ucon;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphI::rwalk_2fp1tn1.\n";
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
                ic2[ii] = ic1[ii]; //corrected
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
                ic2[ii] = ic1[ii]; //corrected
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

void GraphI::rwalk_2fp1tn2(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {

  int h,i,j,k,l,sig, vval, numutr,mr, enarch,ii,jj;
  int *ic2;
  ic2 = new int[n];
  int uy = 0;//nuevo
  for (ii=0; ii<n; ii++)
    ic2[ii] = ic1[ii];
  for (i=0; i < n; i++) {
    if ((ic1[i]==fp1[i]) && (fp1[i] != fp2[i])) {
      ic2[i] *= (-1);
      if (uy == 0)
        uy++;
      if (uy > 0)
        break;
    }
  }
  if (i==n) {
    cout << "[Error]: Impossible initial network for GraphI::rwalk:2fp1tn2.\n";
    exit(1);
  }
  make_robust_nw(n,e, ic1, fp1, ic2, fp2);//hastaaqui
  
  
  
  ofstream fs;
  bool ucon;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphI::rwalk_2fp1tn2.\n";
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
                ic2[ii] = ic1[ii]; //corrected
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

void GraphI::clear()
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

void GraphI::clear_edge_count()
{
	numofe = 0;
	yae = false;
}

void GraphI::clear_loops()
{
  delete [] loops;
  totalloops = 0;
  int i;
  for (i = 0; i < size; i++)
    delete [] guada[i];
  delete [] guada;
  floops = false;
}

void GraphI::clear_comp()
{
  delete [] components;
  yacomp = false;
}

void GraphI::clear_bc()
{
  delete [] bc;
  bcya = false;
}

void GraphI::clear_scc()
{
  delete [] scc;
  cfc = false;
  numbscc = 0;
}

void GraphI::clear_adj() {
  int i;
  for (i = 0; i < size; i++)
    delete [] matadya[i];
  delete [] matadya;
  yamadya = false;
}

void GraphI::clear_attractor()
{
  int i;
  for (i = 0; i < atsize; i++)
    delete [] attractor[i];
  delete [] attractor;
  atsize = 0;
  palen = 0;
  yatra = false;
}

void GraphI::clear_attractors()
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

void GraphI::clear_dima()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] dima[i];
  delete [] dima;
  dist = false;
}

void GraphI::clear_ds()
{
  int i;
  for (i = 0; i < size; i++)
    downstream[i].clear();
  delete [] downstream;
  infl = false;
}

void GraphI::clear_dict()
{
  delete [] dict;
  ionary = false;
}

void GraphI::clear_deg()
{
  delete [] degree;
  delete [] outdegree;
  delete [] indegree;
  yadeg = false;
	if (yadegdist)
		clear_degdist();
}

void GraphI::clear_degdist()
{
	delete [] degdist;
	delete [] odegdist;
	delete [] idegdist;
	yadegdist = false;
}

//access
int GraphI::number_of_edges()
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

int GraphI::number_of_nodes() 
{
  return size;
}

int GraphI::weight(int source, int target)
{
  return nw[target][source];
}

bool GraphI::is_directed()
{
  return directed;
}

bool GraphI::attractor_exists()
{
  return yatra;
}

bool GraphI::attractors_exist()
{
  return yatras;
}

//modification
void GraphI::force_interaction(int source, int target, int value)
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
		cout << "[Error]: GraphI::force_interaction does not work for undirected graphs.\n";
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

void GraphI::change_interaction(int source, int target, int value)
{
  if (directed) {
    if (value == nw[target][source]) {
      cout << "[Error]: GraphI::change_interaction changes to stay the same.\n";
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
    cout << "[Error]: GraphI::change_interaction does not work for undirected graphs.\n";
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

void GraphI::change_interaction_undir(int source, int target, int value)
{
  if (directed) {
    cout << "[Error]: GraphI::change_interaction_undir does not work for directed graphs.\n";
    exit(1);
  }
  if (value == nw[target][source]) {
    cout << "[Error]: GraphI::change_interaction changes to stay the same.";
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
void GraphI::get_dir_nw_from_file(int nn, istream& en, int e) //get directed network from file
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
      cout << "[Error]: This is not an int network. GraphI::get_dir_nw_from_file.\n";
      exit(1);
    }
  }
  set_default_exclusive_vars();
}

void GraphI::get_dir_nw_from_file(int nn, string arch)
{
	ifstream sal;
	basic.open_ifstream(sal, arch);
	get_dir_nw_from_file(nn, sal);	
	sal.close();
}

void GraphI::get_dir_nw_from_file(int nn, istream& en)
{
  make_nw(nn, true);
  int i, j, k;
  while (en >> j) {
		en >> i;
		en >> k;
		if (((k==1) || (k==(-1))) && (i < size) && (j < size) && (i >= 0) && (j >= 0))
			nw[i][j] = k;
		else {
			cout << "[Error]: This is not an int network. GraphI::get_dir_nw_from_file.\n";
			exit(1);
		}
	}
  set_default_exclusive_vars();
}

//Export
void GraphI::export_nw(string arch) {
  ofstream fs;
  basic.open_ofstream(fs, arch);
  export_nw(fs);
  fs.close();
}

void GraphI::export_nw(ostream& fs) {
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
int GraphI::distance_from_nw(GraphI& otra) {
  int res = 0;
  int i, j;
  for (i=0; i < size; i++)
    for (j = 0; j < size; j++)
      res += abs(weight(i,j) - otra.weight(i,j));
  return res;
}


bool GraphI::equal_adjmat(GraphI& otra) {
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


bool GraphI::equal_nw(GraphI &templ)
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

void GraphI::synchrony() //one step in network dynamics
{
  int i, j, sum;
  for (i = 0; i < size; i++) {
		sum = 0;
		for (j = 0; j < size; j++) {
			if (weight(j,i) != 0)
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

void GraphI::set_as_state(int *vec)
{
  int i;
  for (i=0; i< size; i++)
    edo[i] = vec[i];
}

void GraphI::find_an_attractor()
{
  if (yatra) {
    cout << "[Error]: Attractor already found when GraphI::find_an_attractor was called.\n";
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

void GraphI::find_attractors(int **ic, int nic)
{
  if (yatras) {
    cout << "[Error]: Attractors already found when GraphI::find_attractors was called.\n";
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

int GraphI::attractor_size()
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphI::attractor_size was called.\n";
    exit(1);
  }
  return atsize;
}

int GraphI::attractor_element(int row, int node)
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphI::attractor_element was called.\n";
    exit(1);
  }
  if ((row >=atsize)||(row < 0)) {
    cout << "[Error]: Attractor has no " << row << "th row. GraphI::attractor_element.\n";
    exit(1);
  }
  if ((node >= size) || (node < 0)) {
    cout << "[Error]: Node " << node << " does not exist. GraphI::attractor_element.\n";
  }
  return attractor[row][node];
}

void GraphI::print_attractor(ostream& sal)
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphI::print_attractor was called.\n";
    exit(1);
  }
  int i,j;
  for (i=0; i < atsize; i++) {
    for (j=0; j<size; j++)
      sal << attractor_element(i,j) << " ";
    sal << endl;
  }
}

int GraphI::attractor_size(int wh)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphI::attractor_size was called.\n";
    exit(1);
  }
  return atsizes[wh];
}

int GraphI::attractor_element(int wh, int row, int node)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphI::attractor_element was called.\n";
    exit(1);
  }
  if ((row >=atsizes[wh])||(row < 0)) {
    cout << "[Error]: Attractor has no " << row << "th row. GraphI::attractor_element.\n";
    exit(1);
  }
  if ((node >= size) || (node < 0)) {
    cout << "[Error]: Node " << node << " does not exist. GraphI::attractor_element.\n";
  }
  return attractors[wh][row][node];
}

void GraphI::print_attractor(int wh, ostream& sal)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphI::print_attractor was called.\n";
    exit(1);
  }
  int i,j;
  for (i=0; i < attractor_size(wh); i++) {
    for (j=0; j<size; j++)
      sal << attractor_element(wh,i,j) << " ";
    sal << endl;
  }
}

int GraphI::number_of_attractors()
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphI::number_of_attractors was called.\n";
    exit(1);
  }
  return numatrs;
}

void GraphI::penetrance_n1_2p(int *cior, int *pA, int *pB, int& countA, int& countB) {
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

void GraphI::penetrance_n2_ex_2p(int *cior, int *pA, int *pB, int& countA, int& countB) {
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

void GraphI::mutational_access_2p(int *ci, int *pA, int *pB, int& countA, int& countB) {
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
void GraphI::build_offspring(GraphI &mother, GraphI &father, bool *frommom) {
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

void GraphI::mate(GraphI &mother, GraphI &father)
{
  if (mother.number_of_nodes() != father.number_of_nodes()) {
		cout << "[Error]: Parents with different genome size. GraphI::mate.\n";
		exit(1);
	}
	if (mother.is_directed() != father.is_directed()) {
		cout << "[Error]: One parent is a directed network while the other is an undirected network. GraphI::mate.\n";
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

void GraphI::mutate(double con) {
  int i = est.randint(0, size);
  mutate_gene(i, con);
}

void GraphI::mutate_gene(int gene, double con) {
  if (!directed) {
    cout << "[Error]: GraphI::mutate_gene is not defined for undirected networks.\n";
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

void GraphI::consider_mutation(double muratepg, double con) {
  int i;
  for (i = 0; i < size; i++)
    if (est.randreal() < muratepg)
      mutate_gene(i, con);
}


//private

void GraphI::set_default_exclusive_vars() //sets default variables
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
