//This program finds the initial conditions by solving the ODE's in Goldbeter et al 2007.

#include <boost/array.hpp>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "libs/horse.h"

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;
using std::cout;
using std::cin;
using std::endl;
using std::setprecision;
using std::fixed;
using std::setw;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::string;

//Constants
const double Vs1 = 4.97;
const double Mf = 1.5;
const double V0 = 0.365;
const double Vsc = 7.1;
const double Kd1 = 1;
const double Kd2 = 0.28;
const double Kd3 = 1;
const double Kd4 = 1;
const double Kd5 = 0;
const double Ks2 = 1;
const double Ks3 = 1;
const double Ka = 0.2;
const double Ki = 0.25;
const double n = 2;
const double m = 2;

const double A = 4.97;//Vs1
const double B = 0.365;//V0
const double C = 7.1;//Vsc
const double D = 0.04;//Ka^n = 0.2²
const double E = 0.28;//Kd2
const double F = 1.5;//ks3*Mf= 1*1.5
const double G = 0.0625;//Ki^m= 0.25²

typedef boost::array< double , 4 > state_type;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = A - ( x[2] * x[0] );
    dxdt[1] = B + C * ( pow(x[3],2) / ( D + pow(x[3],2) ) ) - x[1];
    dxdt[2] = x[1] - E * x[2];
    dxdt[3] = F * ( G / ( G + pow(x[0],2) ) ) - x[3];
}

void write_lorenz( const state_type &x , const double t )
{
    Horse hor;
    ofstream fo;
    hor.open_ofstream_new(fo, ""+ root +"/getci4_tray.txt");
    fo << t << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
    fo.close();
}

int main(int argc, char **argv){
    
    system("date");
    time_t rawtime;
    time ( &rawtime );
    
    string folder = argv[1];
    int seed = atoi(argv[2]);
    string root = argv[3];
    
    Horse hor(seed);
    ofstream fo;
    ifstream fi;
    
    double condi1, condi2, condi3, condi4;
    double alfa1=0.1, alfa2=0.1, alfa3=0.1, alfa4=0.2;
    //From figure 3 : "starting from various initial conditions for F, in the range 0–100 nM, and setting the other initial conditions as RA=MC=C= 0.1 nM."
    double t, x0, x1, x2, x3;
    string basura;
    double ratiox0, ratiox3, rho;
    
    do {
        //Select random initial conditions near alfa
        condi1 = alfa1*hor.randreal() + (alfa1*0.5);
        condi2 = alfa2*hor.randreal() + (alfa2*0.5);
        condi3 = alfa3*hor.randreal() + (alfa3*0.5);
        condi4 = alfa4*hor.randreal() + (alfa4*0.5);
        
        state_type x = { condi1 , condi2 , condi3, condi4 }; // initial conditions
        integrate( lorenz , x , 0.0 , 20.0 , 0.1 , write_lorenz );
        //Integrate parameters: system, initial conditions, initial t, final t, interval, observer.
        
        hor.open_ifstream(fi, ""+ root +"/getci4_tray.txt");
        fi >> t;
        fi >> x0;
        fi >> x1;
        fi >> x2;
        fi >> x3;
        fi.close();
        
        ratiox0 = x0 / (x0 + 1);
        ratiox3 = x3 / (x3 + 1);
        rho = ratiox3 / ratiox0;

    }while(rho > 0.5); // x3 over x0. Because of "values of rho larger or smaller than unity will reflect predominance of FGF over RA signaling, or of RA over FGF", we select a rho below 0.5 to ensure RA predominance over FGF.
    cout << "rho: " << rho << endl;
    //Check in Octave
    hor.open_ofstream_new(fo, ""+ root +"/octave_script.m");
    fo << "addpath(genpath('../'));\n";
    fo << "[t,y] = ode45(@(t,y) ra_fgf(t, y, "<< Vs1 << ", " << Mf << ", " << V0 << ", " << Vsc << ", " << kd1 << ", " << kd2 << ", " << kd3 << ", " << kd4 << ", " << kd5 << ", " << ks2 << ", " << ks3 << ", " << Ka << ", " << Ki << ", " << n << ", " << m <<"), [0 20],["<< condi1 <<"; "<< condi2 <<"; "<< condi3 <<"; "<< condi4 <<"]);\n";
    fo << "[m, n] = size(y);";
    fo << "ratiox0 = y(m,1) / (y(m,1) + 1);\n";
    fo << "ratiox3 = y(m,4) / (y(m,4) + 1);\n";
    fo << "rho = ratiox3 / ratiox0;\n";
    fo << "fileID = fopen('octave_res.txt','w');\n";
    fo << "fprintf(fileID,'%f %f %f %f %f', rho, y(m,1), y(m,2), y(m,3), y(m,4));\n";
    fo << "fclose(fileID);\n";
    fo << "exit;\n";
    fo.close();
    
    system("octave "+ root +"/octave_script.m");
    
    hor.open_ifstream(fi, ""+ root +"/octave_res.txt");
    fi >> rho;
    fi.close();
    cout << "rho_octave: " << rho << endl;
    if(rho > 0.5)
        cout << "Error in seed" << endl;
    else{
        //Save the initial conditions
        hor.open_ofstream(fo, ""+ root +"/S"+ folder +"/focal_ic.txt");
        fo << condi1 << endl;
        fo << condi2 << endl;
        fo << condi3 << endl;
        fo << condi4 << endl;
        fo.close();
    }
    
    hor.close_rng();
    
     cout << "------------Ejecución terminada!!-----------\n";
    
    cout << "Inicio: " << ctime (&rawtime) << endl;
    system("date");
    
    return 0;
}
