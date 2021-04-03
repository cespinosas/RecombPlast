#ifndef BASICS_H
#define BASICS_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <set>
#include <list>
#include <string>
#include "alea.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <sys/stat.h>
#include <dirent.h>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;

class Basics
{
public:
	Basics();
	Basics(Alea& jacta);
	void start_rng(Alea& jacta);
	
	//vectors
	bool eqvec(int vec1[], int s1, int vec2[], int s2);
	bool eqvec(double vec1[], int s1, double vec2[], int s2);
	bool eqvec(bool vec1[], int s1, bool vec2[], int s2);
  bool eqvec(char vec1[], int s1, char vec2[], int s2);
  bool eqvec(string vec1[], int s1, string vec2[], int s2);
  
	int difsinvecs(int vec1[], int s1, int vec2[], int s2);
	int difsinvecs(double vec1[], int s1, double vec2[], int s2);
	int difsinvecs(char vec1[], int s1, char vec2[], int s2);
	int difsinvecs(bool vec1[], int s1, bool vec2[], int s2);
	
	void fillv0(int vec[], int s);
	void fillv0(double vec[], int s);
	void fillv0(bool vec[], int s);
  void fillv0(string vec[], int s);
	
	void fillv1(int vec[], int s);
	void fillv1(double vec[], int s);
	void fillv1(bool vec[], int s);
  
	void fillvm1(int vec[], int s);
	void fillvm1(double vec[], int s);
	
	bool last_vector(bool *ve, int siz); //checks if vector is full of true values
	bool last_vector(int *ve, int siz, int max);
	
	void next_in_vector(bool *ve, int siz);
	void next_in_vector(int *ve, int siz, int min, int max, int step);
	
	//vector formatting
	void printv(ostream& sal, int vec[], int s);
	void printv(ostream& sal, double vec[], int s);
	void printv(ostream& sal, bool vec[], int s);
	void printv(ostream& sal, string vec[], int s);
	
	//vectors and matrices
	int vecinmat(int** mat, int rows, int cols, int vec[], int s);
	int vecinmat(double** mat, int rows, int cols, double vec[], int s);
	int vecinmat(bool** mat, int rows, int cols, bool vec[], int s);
	int vecinmat(string** mat, int rows, int cols, string vec[], int s);
	
	//matrices
	void fillmat0(int** mat, int rows, int cols);
	void fillmat0(double** mat, int rows, int cols);
	void fillmat0(bool** mat, int rows, int cols);
	
	void fillmatm1(int** mat, int rows, int cols);
	void fillmatm1(double** mat, int rows, int cols);
	bool is_symmetric(int** mat, int rows, int cols);
	bool is_symmetric(double** mat, int rows, int cols);
	bool is_symmetric(bool** mat, int rows, int cols);
	bool is_symmetric(string** mat, int rows, int cols);
	
	bool last_matrix(bool **ve, int rows, int cols); //checks whether matrix is full of true values
	bool last_matrix(int **ve, int rows, int cols, int max);
	
	void next_in_matrix(bool **ve, int rows, int cols);
	void next_in_matrix(int **ve, int rows, int cols, int min, int max, int step);
  
  bool eqmatrix(int **mat1, int r1, int c1, int **mat2, int r2, int c2);
  bool eqmatrix(double **mat1, int r1, int c1, double **mat2, int r2, int c2);
  bool eqmatrix(char **mat1, int r1, int c1, char **mat2, int r2, int c2);
  bool eqmatrix(bool **mat1, int r1, int c1, bool **mat2, int r2, int c2);
  bool eqmatrix(string **mat1, int r1, int c1, string **mat2, int r2, int c2);
  
	//matrix formatting
  void printm(ostream& sal, int** mat, int rows, int cols);
  void printm(ostream& sal, double** mat, int rows, int cols);
  void printm(ostream& sal, bool** mat, int rows, int cols);
  void printm(ostream& sal, string** mat, int rows, int cols);
  
  void printm_latex(ostream& sal, int** mat, int rows, int cols);
  void printm_latex(ostream& sal, double** mat, int rows, int cols);
	
	//others
  void run_command(string cuerda);
  
  void create_array(int** &arr, int rows, int cols);
  void create_array(bool** &arr, int rows, int cols);
  void create_array(double** &arr, int rows, int cols);
  void create_array(char** &arr, int rows, int cols);
  void create_array(string** &arr, int rows, int cols);
  
  void create_array(int*** &arr, int slices, int rows, int cols);
  void create_array(bool*** &arr, int slices, int rows, int cols);
  void create_array(double*** &arr, int slices, int rows, int cols);
  void create_array(char*** &arr, int slices, int rows, int cols);
  void create_array(string*** &arr, int slices, int rows, int cols);
	
	
	void open_ifstream(ifstream& fe, string nomb);
	void open_ofstream(ofstream& fs, string nomb);
	void open_ofstream_to_append(ofstream& fs, string nomb);
	
  bool dir_exists(string dir);
	
	//strings
	string inttostring(int num);
	
	//statistics
	double get_mean(int* vec, int tam);
	double get_mean(double* vec, int tam);
  double get_mean(int* vec, int tam, set<int> &cua);
  double get_mean(double* vec, int tam, set<int> &cua);
  double get_mean_by_col(int **mat, int rows, int cols, int col);
  double get_mean_by_col(double **mat, int rows, int cols, int col);
	double get_sample_variance(int* vec, int tam);
	double get_sample_variance(double* vec, int tam);
  double get_sample_variance(int* vec, int tam, set<int> &cua);
  double get_sample_variance(double* vec, int tam, set<int> &cua);
  double get_sample_variance(int* vec, int tam, double ave);
  double get_sample_variance(double* vec, int tam, double ave);
  double get_sample_variance(int* vec, int tam, double ave, set<int> &cua);
  double get_sample_variance(double* vec, int tam, double ave, set<int> &cua);
  double get_pop_variance(int* vec, int tam);
  double get_pop_variance(double* vec, int tam);
  double get_pop_variance(int* vec, int tam, set<int> &cua);
  double get_pop_variance(double* vec, int tam, set<int> &cua);
  double get_pop_variance(int* vec, int tam, double ave);
  double get_pop_variance(double* vec, int tam, double ave);
  double get_pop_variance(int* vec, int tam, double ave, set<int> &cua);
  double get_pop_variance(double* vec, int tam, double ave, set<int> &cua);
  double get_sample_stddev(int* vec, int tam);
  double get_sample_stddev(double* vec, int tam);
	double get_sample_stddev(int* vec, int tam, set<int> &cua);
	double get_sample_stddev(double* vec, int tam, set<int> &cua);
  double get_sample_stddev(int* vec, int tam, double ave);
  double get_sample_stddev(double* vec, int tam, double ave);
  double get_sample_stddev(int* vec, int tam, double ave, set<int> &cua);
  double get_sample_stddev(double* vec, int tam, double ave, set<int> &cua);
	double get_pop_stddev(int* vec, int tam);
	double get_pop_stddev(double* vec, int tam);
  double get_pop_stddev(int* vec, int tam, set<int> &cua);
  double get_pop_stddev(double* vec, int tam, set<int> &cua);
  double get_pop_stddev(int* vec, int tam, double ave);
  double get_pop_stddev(double* vec, int tam, double ave);
  double get_pop_stddev(int* vec, int tam, double ave, set<int> &cua);
  double get_pop_stddev(double* vec, int tam, double ave, set<int> &cua);
	
	double find_min(double* vec, int tam);
	int find_min(int* vec, int tam);
  
	double find_max(double* vec, int tam);
  int find_max(int* vec, int tam);
	void sort(double* vec, double* nvec, int tam);
	void sort(int* vec, int* nvec, int tam);
  
  double get_midpoint(int* vo, int tam); //mediana en vector ordenado
  double get_midpoint(double* vo, int tam); //mediana en vector ordenado
  double get_q1(int* vo, int tam);
  double get_q1(double* vo, int tam);
  double get_q3(int* vo, int tam);
  double get_q3(double* vo, int tam);
  
  
	//para gsl y distribucion normal: si z > 37.515, entonces p < 2.62271 \times 10^{-308}
	Alea est;
	
private:
};

#endif



