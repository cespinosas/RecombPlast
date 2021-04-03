#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <set>
#include <list>
#include <string>
#include "alea.h"
#include "basics.h"
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

Basics::Basics()
{
}

Basics::Basics(Alea& jacta) //creates instance of class Basics and assigns rng to jacta
{
	start_rng(jacta);
}

void Basics::start_rng(Alea& jacta) //assigns rng to jacta
{
  est = jacta;
}

//vectors
bool Basics::eqvec(int vec1[], int s1, int vec2[], int s2) //are two vectors equal
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i]) {
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(double vec1[], int s1, double vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i]) {
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(bool vec1[], int s1, bool vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])	{
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(char vec1[], int s1, char vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])	{
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(string vec1[], int s1, string vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])	{
				res = false;
				break;
			}
  }
  return res;
}

int Basics::difsinvecs(int vec1[], int s1, int vec2[], int s2) //number of differences between two vectors
{
  int i;
  int dif = 0;
  if (s1 != s2)
  {
    cout << "[Error]: You can not count differences between vectors with different sizes using Basics::difsinvecs.\n";
    exit(1);
  }
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])
				dif++;
  }
  return dif;
}

int Basics::difsinvecs(double vec1[], int s1, double vec2[], int s2) //number of differences between two vectors
{
  int i;
  int dif = 0;
  if (s1 != s2)
  {
    cout << "[Error]: You can not count differences between vectors with different sizes using Basics::difsinvecs.\n";
    exit(1);
  }
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])
				dif++;
  }
  return dif;
}

int Basics::difsinvecs(char vec1[], int s1, char vec2[], int s2)
{
  int i;
  int dif = 0;
  if (s1 != s2)
  {
		cout << "[Error]: You can not count differences between vectors with different sizes using Basics::difsinvecs.\n";
		exit(1);
  }
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])
				dif++;
  }
  return dif;
}

int Basics::difsinvecs(bool vec1[], int s1, bool vec2[], int s2)
{
  int i;
  int dif = 0;
  if (s1 != s2)
  {
    cout << "[Error]: You can not count differences between vectors with different sizes using Basics::difsinvecs.\n";
    exit(1);
  }
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])
				dif++;
  }
  return dif;
}

void Basics::fillv0(int vec[], int s) //fills 1d array with 0s
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 0;
}

void Basics::fillv0(double vec[], int s) //fills 1d array with 0s
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 0;
}

void Basics::fillv0(bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = false;
}

void Basics::fillv0(string vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = "";  
}

void Basics::fillv1(int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 1;
}

void Basics::fillv1(double vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 1;
}

void Basics::fillv1(bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = true;
}

void Basics::fillvm1(int vec[], int s) //fills vector with -1
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = -1;
}

void Basics::fillvm1(double vec[], int s) //fills vector with -1
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = -1;
}

bool Basics::last_vector(bool *ve, int siz)
{
	bool res = true;
	int i;
	for (i=0; i<siz; i++)
		if (!ve[i]){
			res = false;
			break;
		}
	return res;
}

bool Basics::last_vector(int *ve, int siz, int max)
{
	bool res = true;
	int i;
	for (i=0; i<siz; i++)
		if (ve[i] != max){
			res = false;
			break;
		}
	return res;
}

void Basics::next_in_vector(bool *ve, int siz)
{
	int i;
	if (last_vector(ve, siz)) {
		cout << "[Error]: Last vector already reached in Basics::next_in_vector.\n";
		exit(1);
	}
	for (i=0; i<siz; i++)
		if (!ve[i]) {
			ve[i] = true;
			break;
		}
	fillv0(ve, i);
}

void Basics::next_in_vector(int *ve, int siz, int min, int max, int step)
{
	int i,j;
	if (last_vector(ve, siz, max)) {
		cout << "[Error]: Last vector already reached in Basics::next_in_vector.\n";
		exit(1);
  }
  for (i=0; i<siz; i++)
		if (ve[i] < max) {
			ve[i] = ve[i]+step;
			break;
		}
  for (j=0; j<i; j++)
    ve[j]=min;
}

//vector formatting
void Basics::printv(ostream& sal, int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    sal << vec[i] << "\t";
}

void Basics::printv(ostream& sal, double vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    sal << vec[i] << "\t";
}

void Basics::printv(ostream& sal, bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    sal << vec[i] << "\t";
}

void Basics::printv(ostream& sal, string vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    sal << vec[i] << "\t";
}

//vectors and matrices
int Basics::vecinmat(int** mat, int rows, int cols, int vec[], int s) //vector appears as matrix row
{
  int i;
  int res = -1;
  if (cols != s) {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols))	{
				res = i;
				break;
			}
  }
  return res;
}


int Basics::vecinmat(double** mat, int rows, int cols, double vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s) {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols)){
				res = i;
				break;
			}
  }
  return res;
}

int Basics::vecinmat(bool** mat, int rows, int cols, bool vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s)   {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols))	{
				res = i;
				break;
			}
  }
  return res;
}

int Basics::vecinmat(string** mat, int rows, int cols, string vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s)   {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols))	{
				res = i;
				break;
			}
  }
  return res;
}

//matrices
void Basics::fillmat0(int** mat, int rows, int cols) //fills 2d array with 0s
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmat0(double** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmat0(bool** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmatm1(int** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillvm1(mat[i], cols);
}

void Basics::fillmatm1(double** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillvm1(mat[i], cols);
}

bool Basics::is_symmetric(int** mat, int rows, int cols)
{
  int i, j;
  bool res = true;
  if (rows != cols)
    res = false;
  else {
    for (i = 0; i < rows; i++)
      for (j = 0; j < i; j++)
				if (mat[i][j] != mat[j][i]) {
					res = false;
					break;
				}
  }
  return res;
}

bool Basics::is_symmetric(double** mat, int rows, int cols)
{
  int i, j;
  bool res = true;
  if (rows != cols)
    res = false;
  else {
    for (i = 0; i < rows; i++)
      for (j = 0; j < i; j++)
				if (mat[i][j] != mat[j][i])	  {
					res = false;
					break;
				}
  }
  return res;
}

bool Basics::is_symmetric(bool** mat, int rows, int cols)
{
  int i, j;
  bool res = true;
  if (rows != cols)
    res = false;
  else {
    for (i = 0; i < rows; i++)
      for (j = 0; j < i; j++)
				if (mat[i][j] != mat[j][i])	  {
					res = false;
					break;
				}
  }
  return res;
}

bool Basics::is_symmetric(string** mat, int rows, int cols)
{
  int i, j;
  bool res = true;
  if (rows != cols)
    res = false;
  else {
    for (i = 0; i < rows; i++)
      for (j = 0; j < i; j++)
				if (mat[i][j] != mat[j][i])	  {
					res = false;
					break;
				}
  }
  return res;
}

bool Basics::last_matrix(bool **ve, int rows, int cols)
{
	bool res = true;
	int i;
	for (i=0; i< rows; i++)
		if (!last_vector(ve[i], cols)) {
			res = false;
			break;
		}
	return res;
}

bool Basics::last_matrix(int **ve, int rows, int cols, int max)
{
	bool res = true;
	int i;
	for (i=0; i< rows; i++)
		if (!last_vector(ve[i], cols, max)) {
			res = false;
			break;
		}
	return res;
}

void Basics::next_in_matrix(bool **ve, int rows, int cols)
{
	int i,j;
	if (last_matrix(ve, rows, cols)) {
		cout << "[Error]: Last matrix already reached in Basics::next_in_matrix.\n";
		exit(1);
	}
	for (i=0; i<rows; i++)
		if (!last_vector(ve[i], cols)) {
			next_in_vector(ve[i], cols);
			break;
		}
	for (j=0; j<i; j++)
		fillv0(ve[j], cols);
}


void Basics::next_in_matrix(int **ve, int rows, int cols, int min, int max, int step)
{
	int i,j,k;
	if (last_matrix(ve, rows, cols, max)) {
		cout << "[Error]: Last matrix already reached in Basics::next_in_matrix.\n";
		exit(1);
	}
	for (i=0; i<rows; i++)
		if (!last_vector(ve[i], cols, max)) {
			next_in_vector(ve[i], cols, min, max, step);
			break;
		}
	for (j=0; j<i; j++)
		for (k=0; k < cols; k++)
			ve[j][k] = min;
}

bool Basics::eqmatrix(int **mat1, int r1, int c1, int **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::eqmatrix(double **mat1, int r1, int c1, double **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::eqmatrix(char **mat1, int r1, int c1, char **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::eqmatrix(bool **mat1, int r1, int c1, bool **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::eqmatrix(string **mat1, int r1, int c1, string **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

//matrix formatting
void Basics::printm(ostream& sal, int** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)    {
		printv(sal, mat[i], cols);
		sal << endl;
	}
}

void Basics::printm(ostream& sal, double** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)    {
		printv(sal, mat[i], cols);
		sal << endl;
	}
}

void Basics::printm(ostream& sal, bool** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)    {
		printv(sal, mat[i], cols);
		sal << endl;
	}
}

void Basics::printm(ostream& sal, string** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)    {
		printv(sal, mat[i], cols);
		sal << endl;
	}
}

void Basics::printm_latex(ostream& sal, int** mat, int rows, int cols)
{
  int i,j;
  sal << "\\left[ \\begin{array}{";
  for (i= 0; i < cols; i++)
    sal << "c";
  sal << "}\n";
  for (i=0; i < rows; i++) {
    sal << mat[i][0];
    for (j=1; j < cols; j++)
      sal << "&" << mat[i][j];
    if (i < (rows-1))
      sal << "\\\\ \n";
  }
  sal << "\n \\end{array} \\right]\n";
}

void Basics::printm_latex(ostream& sal, double** mat, int rows, int cols)
{
  int i,j;
  sal << "\\left[ \\begin{array}{";
  for (i= 0; i < cols; i++)
    sal << "c";
  sal << "}\n";
  for (i=0; i < rows; i++) {
    sal << mat[i][0];
    for (j=1; j < cols; j++)
      sal << "&" << mat[i][j];
    if (i < (rows-1))
      sal << "\\\\ \n";
  }
  sal << "\n \\end{array} \\right]\n";
}

//others
void Basics::run_command(string cuerda){
  system(cuerda.c_str());
  return;
}

void Basics::create_array(int** &arr, int rows, int cols) {
  int i;
  arr = new int*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new int[cols];
  return;
}

void Basics::create_array(bool** &arr, int rows, int cols) {
  int i;
  arr = new bool*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new bool[cols];
  return;
}

void Basics::create_array(double** &arr, int rows, int cols) {
  int i;
  arr = new double*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new double[cols];
  return;
}

void Basics::create_array(char** &arr, int rows, int cols) {
  int i;
  arr = new char*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new char[cols];
  return;
}

void Basics::create_array(string** &arr, int rows, int cols) {
  int i;
  arr = new string*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new string[cols];
  return;
}

void Basics::create_array(int*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new int**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new int*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new int[cols];
  }
  return;
}

void Basics::create_array(bool*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new bool**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new bool*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new bool[cols];
  }
  return;
}

void Basics::create_array(double*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new double**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new double*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new double[cols];
  }
  return;
}

void Basics::create_array(char*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new char**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new char*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new char[cols];
  }
  return;
}

void Basics::create_array(string*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new string**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new string*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new string[cols];
  }
  return;
}

void Basics::open_ifstream(ifstream& fe, string nomb) //opens input file stream
{
	fe.open(nomb.c_str());
	if (!fe.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ifstream.\n";
		exit(1);
	}
}

void Basics::open_ofstream(ofstream& fs, string nomb) //opens output file stream
{
	fs.open(nomb.c_str());
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ofstream.\n";
		exit(1);
	}
}

void Basics::open_ofstream_to_append(ofstream& fs, string nomb)
{
	fs.open(nomb.c_str(), ios::app);
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ofstream.\n";
		exit(1);
	}
}

bool Basics::dir_exists(string dir)
{
  struct stat buf;
  int esmu = stat(dir.c_str(), & buf);
  if ((esmu >= 0) && (S_ISDIR(buf.st_mode)))
    return true;
  else
    return false;
}

//strings

string Basics::inttostring(int num)
{
	char buff[50];
	int j;
	j = sprintf(buff, "%d", num);
	if (j <0) {
		cout << "[Error]: Transformation of string to int was not possible using Basics::inttostring.\n";
		exit(1);
	}
	string res(buff);
	return res;
}

//statistics
double Basics::get_mean(int* vec, int tam)
{
  double res = 0;
  int i;
  for (i = 0; i < tam; i++)
    res = res + vec[i];
  res = res/(tam*1.0);
  return res;
}

double Basics::get_mean(double* vec, int tam)
{
  double res = 0;
  int i;
  for (i = 0; i < tam; i++)
    res = res + vec[i];
  res = res/(tam*1.0);
  return res;
}

double Basics::get_mean(int* vec, int tam, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += double(vec[*it]);
  res /= double(cua.size());
  return res;
}

double Basics::get_mean(double* vec, int tam, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += vec[*it];
  res /= double(cua.size());
  return res;
}

double Basics::get_mean_by_col(int **mat, int rows, int cols, int col)
{
  if (col >= cols) {
    cout << "[Error]: Column does not exist in matrix. Basics::get_mean_by_col.\n";
    exit(1);
  }
  double res = 0;
  int i;
  for (i= 0; i< rows; i++)
    res = res + double(mat[i][col]);
  res = res/double(rows);
  return res;
}

double Basics::get_mean_by_col(double **mat, int rows, int cols, int col)
{
  if (col >= cols) {
    cout << "[Error]: Column does not exist in matrix. Basics::get_mean_by_col.\n";
    exit(1);
  }
  double res = 0;
  int i;
  for (i= 0; i< rows; i++)
    res = res + mat[i][col];
  res = res/double(rows);
  return res;
}

double Basics::get_sample_variance(int* vec, int tam)
{
  double mean = get_mean(vec, tam);
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((mean - vec[i])*(mean - vec[i]));
  res = sum/(tam-1.0);
  return res;
}

double Basics::get_sample_variance(double* vec, int tam)
{
  double mean = get_mean(vec, tam);
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((mean - vec[i])*(mean - vec[i]));
  res = sum/(tam-1.0);
  return res;
}

double Basics::get_sample_variance(int* vec, int tam, set<int> &cua) {
  double ave = get_mean(vec, tam, cua);
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()-1.0);
  return res;
}

double Basics::get_sample_variance(double* vec, int tam, set<int> &cua) {
  double ave = get_mean(vec, tam, cua);
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()-1.0);
  return res;
}

double Basics::get_sample_variance(int* vec, int tam, double ave)
{
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((ave - vec[i])*(ave - vec[i]));
  res = sum/(tam-1.0);
  return res;
}

double Basics::get_sample_variance(double* vec, int tam, double ave)
{
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((ave - vec[i])*(ave - vec[i]));
  res = sum/(tam-1.0);
  return res;
}

double Basics::get_sample_variance(int* vec, int tam, double ave, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()-1.0);
  return res;
}

double Basics::get_sample_variance(double* vec, int tam, double ave, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()-1.0);
  return res;
}

double Basics::get_pop_variance(int* vec, int tam)
{
  double mean = get_mean(vec, tam);
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((mean - vec[i])*(mean - vec[i]));
  res = sum/(tam*1.0);
  return res;
}

double Basics::get_pop_variance(double* vec, int tam)
{
  double mean = get_mean(vec, tam);
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((mean - vec[i])*(mean - vec[i]));
  res = sum/(tam*1.0);
  return res;
}

double Basics::get_pop_variance(int* vec, int tam, set<int> &cua) {
  double ave = get_mean(vec, tam, cua);
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()*1.0);
  return res;
}

double Basics::get_pop_variance(double* vec, int tam, set<int> &cua) {
  double ave = get_mean(vec, tam, cua);
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()*1.0);
  return res;
}

double Basics::get_pop_variance(int* vec, int tam, double ave)
{
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((ave - vec[i])*(ave - vec[i]));
  res = sum/(tam*1.0);
  return res;
}

double Basics::get_pop_variance(double* vec, int tam, double ave)
{
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((ave - vec[i])*(ave - vec[i]));
  res = sum/(tam*1.0);
  return res;
}

double Basics::get_pop_variance(int* vec, int tam, double ave, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()*1.0);
  return res;
}

double Basics::get_pop_variance(double* vec, int tam, double ave, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()*1.0);
  return res;
}

double Basics::get_sample_stddev(int* vec, int tam)
{
  double var = get_sample_variance(vec, tam);
	double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(double* vec, int tam)
{
  double var = get_sample_variance(vec, tam);
	double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(int* vec, int tam, set<int> &cua)
{
  double var = get_sample_variance(vec, tam, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(double* vec, int tam, set<int> &cua)
{
  double var = get_sample_variance(vec, tam, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(int* vec, int tam, double ave)
{
  double var = get_sample_variance(vec, tam, ave);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(double* vec, int tam, double ave)
{
  double var = get_sample_variance(vec, tam, ave);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(int* vec, int tam, double ave, set<int> &cua)
{
  double var = get_sample_variance(vec, tam, ave, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(double* vec, int tam, double ave, set<int> &cua)
{
  double var = get_sample_variance(vec, tam, ave, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(int* vec, int tam)
{
  double var = get_pop_variance(vec, tam);
	double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(double* vec, int tam)
{
  double var = get_pop_variance(vec, tam);
	double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(int* vec, int tam, set<int> &cua)
{
  double var = get_pop_variance(vec, tam, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(double* vec, int tam, set<int> &cua)
{
  double var = get_pop_variance(vec, tam, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(int* vec, int tam, double ave)
{
  double var = get_pop_variance(vec, tam, ave);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(double* vec, int tam, double ave)
{
  double var = get_pop_variance(vec, tam, ave);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(int* vec, int tam, double ave, set<int> &cua)
{
  double var = get_pop_variance(vec, tam, ave, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(double* vec, int tam, double ave, set<int> &cua) {
  double var = get_pop_variance(vec, tam, ave, cua);
  double res = sqrt(var);
  return res;
}

double Basics::find_min(double* vec, int tam)
{
  double res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] < res)
      res = vec[i];
  return res;
}

int Basics::find_min(int* vec, int tam)
{
  int res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] < res)
      res = vec[i];
  return res;
}

double Basics::find_max(double* vec, int tam)
{
  double res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] > res)
      res = vec[i];
  return res;
}

int Basics::find_max(int* vec, int tam)
{
  int res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] > res)
      res = vec[i];
  return res;
}

void Basics::sort(int* vec, int* nvec, int tam)
{
	int min, max;
	int i, j, quedan, k;
	int* wov;
	wov = new int[tam];
	int* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
	max = find_max(wov, tam);
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
				i++;
				break;
			}
		}
		wov2 = new int[quedan-1];
		for (k = 0; k < j; k++)
			wov2[k] = wov[k];
		for (k=j; k < (quedan-1); k++)
			wov2[k] = wov[k+1];
		delete [] wov;
		quedan--;
		wov = new int[quedan];
		for (k=0; k < quedan; k++)
			wov[k] = wov2[k];
		delete [] wov2;
	}
	if (min != max) {
		cout << "[Error]: sort does not converge in Basics::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}



void Basics::sort(double* vec, double* nvec, int tam)
{
	double min, max;
	int i, j, quedan, k;
	double* wov;
	wov = new double[tam];
	double* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
	max = find_max(wov, tam);
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
				i++;
				break;
			}
		}
		wov2 = new double[quedan-1];
		for (k = 0; k < j; k++)
			wov2[k] = wov[k];
		for (k=j; k < (quedan-1); k++)
			wov2[k] = wov[k+1];
		delete [] wov;
		quedan--;
		wov = new double[quedan];
		for (k=0; k < quedan; k++)
			wov[k] = wov2[k];
		delete [] wov2;
	}
	if (min != max) {
		cout << "[Error]: sort does not converge in Basics::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}

double Basics::get_midpoint(int* vo, int tam)//vector debe estar ordenado
{
  double res;
  if ((tam %2)==0)
    res = (vo[tam/2] + vo[(tam/2)-1])/2.0;
  else
    res = double(vo[(tam-1)/2]);
  return res;
}

double Basics::get_midpoint(double* vo, int tam)//vector debe estar ordenado
{
  double res;
  if ((tam %2)==0)
    res = (vo[tam/2] + vo[(tam/2)-1])/2.0;
  else
    res = vo[(tam-1)/2];
  return res;
}

double Basics::get_q1(int* vo, int tam)
{
  double res;
  if ((tam%2)==0)
    res = get_midpoint(vo, (tam/2));
  else
    res = get_midpoint(vo, ((tam-1)/2));
  return res;
}

double Basics::get_q1(double* vo, int tam)
{
  double res;
  if ((tam%2)==0)
    res = get_midpoint(vo, (tam/2));//no incluye mediana
  else
    res = get_midpoint(vo, ((tam-1)/2));
  return res;
}

double Basics::get_q3(int* vo, int tam)
{
  double res;
  if ((tam%2)==0) {
    if ((tam%4)==0)
      res = (vo[((tam*3)/4)] + vo[((tam*3)/4)-1])/2.0;
    else
      res = vo[((3*tam)-2)/4];
  }
  else {
    if ((((tam-1)/2)%2)==0)
      res = (vo[(tam-1)*(3/4)] + vo[((tam-1)*(3/4))+1])/2.0;
    else
      res = vo[((3*tam)-1)/4];
  }
  return res;
}

double Basics::get_q3(double* vo, int tam)
{
  double res;
  if ((tam%2)==0) {
    if ((tam%4)==0)
      res = (vo[((tam*3)/4)] + vo[((tam*3)/4)-1])/2.0;
    else
      res = vo[((3*tam)-2)/4];
  }
  else {
    if ((((tam-1)/2)%2)==0)
      res = (vo[(tam-1)*(3/4)] + vo[((tam-1)*(3/4))+1])/2.0;
    else
      res = vo[((3*tam)-1)/4];
  }
  return res;
}
