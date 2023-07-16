#pragma once
#include <vector>
#include <math.h>
#include <string>
#include <conio.h>
#include <fstream>
#include <direct.h>
#include <iomanip>
#include <iostream>
using namespace std;
struct KData
{   
    double alpha;
    double Kr;
    double Ki;
    KData(){
        alpha=0;
        Kr=0;
        Ki=0;
    }
};
pair<vector<KData>, double> Draw(const vector<double>& A6, double& deltax, const string& outpath, bool& build);
double ExperationForKr(const double& alpha, const double& deltax, double& a6);
double ExperationForKi(const double& alpha, const double& deltax, double& a6);
string GetCurrentPath();


void result_out_put(const string& outpath, pair<vector<KData>, double> st,const double& deltax);

class Solution {
private:
    pair<vector<KData>, double> st;
    bool build;
    double deltaX;
    double deltaT;
    int iterationtimes;
    int Mesh_Number;
    vector<double> middle, storage, para,pre;
public:
    Solution(pair<vector<KData>, double> para, double deltaX,double deltaT, int iterationtimes,int Mesh_Number);

    void Getpara();

    double LU(int& j);

    void ExtractSolINIT();

    void loop();

    void RungeKutta();

    void out_put(string& wavepath);
};
