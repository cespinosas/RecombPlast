//This program finds the initial conditions by solving the ODE's developed by Gardner et al 2000.

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include "libs/horse.h"

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

int main(int argc, char **argv)
{
    int seed = atoi(argv[1]);
    Horse hor(seed);
    
    ofstream fs;
    ifstream fe;
    double y1_0, y2_0, y1f, y2f;
    double alfa1=50, alfa2=50;
    int i = 0;
    double disA, disB;
    do {
        cout << i << endl;
        y1_0 = alfa1*hor.randreal() + (alfa1*0.5);
        y2_0 = alfa2*hor.randreal() + (alfa2*0.5);
        hor.open_ofstream_new(fs, "paci.m");
        fs << "addpath(genpath(\'../\'));\n";
        fs << "[t,y] = ode45(@(t,y) toggle(t, y, " << alfa1 << ", " << alfa2 << ", 3, 3), [0 20],[" << y1_0 << "; " << y2_0 << "]);\n";//genotipo focal [50 50]
        fs << "[m, n] = size(y);\n";
        fs << "fileID = fopen(\'mientras.txt\',\'w\');\n";
        fs << "fprintf(fileID,\'%f %f\', y(m,1), y(m,2));\n";
        fs << "fclose(fileID);\n";
        fs << "exit;\n";
        fs.close();
        system("octave paci.m");
        
        hor.open_ifstream(fe, "mientras.txt");
        fe >> y1f;
        fe >> y2f;
        fe.close();
        i++;
        disA = sqrt(((y1f-alfa1)*(y1f-alfa1))+(y2f*y2f));
        disB = sqrt(((y2f-alfa2)*(y2f-alfa2))+(y1f*y1f));
        
    } while(disA >= disB);
    
    string cad;
    cad = "mkdir S"+hor.inttostring(seed);
    system(cad.c_str());
    hor.open_ofstream(fs, "S"+hor.inttostring(seed)+"/focal_ic.txt");
    fs << y1_0 << endl;
    fs << y2_0 << endl;
    fs.close();
    
    hor.close_rng();
    
    return 0;
}

