#ifndef GRAPHI_H
#define GRAPHI_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"

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

class GraphI
{
public:
  GraphI();
  GraphI(Alea& jacta);
  void start_rng(Alea& jacta);
  GraphI(int n, Alea& jacta, bool dirnw);
  void make_nw(int n, bool dirnw);
  void make_robust_nw(int n, int e, int *cia, int *pa, int *cib, int *pb);
  void copy(GraphI &templ);
  GraphI(int n, double c, Alea& jacta, bool dirnw);
  void make_nw(int n, double c, bool dirnw);
  GraphI(int n, int e, Alea& jacta, bool dirnw);
  void make_nw(int n, int e, bool dirnw);
  
  
  void rwalk_1traj(int n, int e, int *fp, int *ic, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  void rwalk_2fp1t(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  void rwalk_2fp1tn1(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  void rwalk_2fp1tn2(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  
  void clear();
	void clear_edge_count();
  void clear_loops();
  void clear_comp();
  void clear_adj();
  void clear_scc();
  void clear_bc();
  void clear_dima();
  void clear_attractor();
  void clear_attractors();
  void clear_ds();
  void clear_dict();
  void clear_deg();
	void clear_degdist();
  
  //access
  int number_of_edges();
  int number_of_nodes();
  int weight(int source, int target);
  bool is_directed();
	bool is_name(string elna);
  bool attractor_exists();
  bool attractors_exist();
  
  //modification
  void force_interaction(int source, int target, int value);
  void change_interaction(int source, int target, int value);
  void change_interaction_undir(int source, int target, int value);
  
	//exclusive
  //import
  void get_dir_nw_from_file(int nn, istream& en, int e); //nuevisimo
  void get_dir_nw_from_file(int nn, string arch);
  void get_dir_nw_from_file(int nn, istream& en);
	
  //export
  void export_nw(string arch);
  void export_nw(ostream& fs);

  //status
	
  //analyses
  int distance_from_nw(GraphI& otra);
  bool equal_adjmat(GraphI& otra);
  bool equal_nw(GraphI &templ);
  
  void synchrony();
  void set_as_state(int *vec);
  void find_an_attractor();
  void find_attractors(int **ic, int nic);
  int attractor_size();
  int attractor_element(int row, int node);
  void print_attractor(ostream& sal);
  int attractor_size(int wh);
  int attractor_element(int wh, int row, int node);
  void print_attractor(int wh, ostream& sal);
  int number_of_attractors();
  void penetrance_n1_2p(int *cior, int *pA, int *pB, int& countA, int& countB);
  void penetrance_n2_ex_2p(int *cior, int *pA, int *pB, int& countA, int& countB);
  void mutational_access_2p(int *ci, int *pA, int *pB, int& countA, int& countB);

  //for evolution
  void build_offspring(GraphI &mother, GraphI &father, bool *frommom);
  void mate(GraphI &mother, GraphI &father);
  void mutate(double con); //201707-\infty; It is already decided that mutation takes place. a single mutation somewhere. (mutation may be without any effect).
  void mutate_gene(int gene, double con);
  void consider_mutation(double muratepg, double con); //probably multiple mutations per genotype
	

  
private:
  int size;
  int **nw;
  int **dima;
  int **matadya;
  bool yamadya;
  Alea est;
  string *dict;
  bool ionary;
  Basics basic;
  bool directed;
  //exclusive
  set<int> *downstream;
  bool infl;
  bool dist;
  bool yacomp;
  int numbcomp;
  int *components;
  int *scc;
  int numbscc;
  bool cfc;
  double *bc;
  bool bcya;
  int *loops;
  int **guada;
  int totalloops;
  bool floops;
  int numofe;
  bool yae;
  int *outdegree;
  int *indegree;
  int *degree;
  bool yadeg;
	int *degdist; // falta poner en graphi y graphc
	int *odegdist;
	int *idegdist;
	bool yadegdist; // falta poner en graphi y graphc... recordar que falta tambien en inicializadas
  int *edo;
  int *ima;
  int **attractor;
  int atsize;
  int palen;
  bool yatra;
  int ***attractors;
  int *atsizes;
  int numatrs;
  bool yatras;
  
	void set_default_exclusive_vars();
};


#endif

