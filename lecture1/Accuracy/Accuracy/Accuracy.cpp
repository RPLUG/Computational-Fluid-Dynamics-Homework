#include <iostream>
#include "Calcu.h"
#include <utility>
#include <math.h>
#include "para.h"
int main(){
    string path = GetCurrentPath();
    string soutionoutpath = path + "\\SolutionResult.txt";
    string parapath = path + "\\para.ini";
    string wavepath= path + "\\WaveResult.txt";
    /*Already find the best solution as a6=-9.59 when deltax = 0.001*/
    int Mesh_Number;
    int iterationtimes;
    double deltaT;
    vector<double> A6;
    bool build = false;
    para_read(parapath, &Mesh_Number, &deltaT, &iterationtimes, A6);
    double deltaX = (2 * 4*atan(1)) / Mesh_Number;
    cout << "The input para are as follows" << endl;
    cout << "Mesh_Number: " << Mesh_Number << endl;
    cout << "Iterationtimes: " << iterationtimes << endl;
    cout << "DeltaT: " << deltaT << endl;
    cout << "DeltaX: " << deltaX << endl;
    cout << "Finding the Optimized solution" << endl;
    pair<vector<KData>, double> st = Draw(A6,deltaX, soutionoutpath,build);
    if (build == false) {
        cout << "find falult" << endl;
        abort();
    }

    result_out_put(soutionoutpath,st,deltaX);
   Solution X(st, deltaX, deltaT, iterationtimes, Mesh_Number);
    X.Getpara();
    X.ExtractSolINIT();
    X.loop();
    X.out_put(wavepath);
    return 0;
}