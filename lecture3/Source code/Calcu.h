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
//求最优解
pair<vector<KData>, double> Draw(const vector<double>& A6, double& deltax, const string& outpath, bool& build);
//得到ki和kr
double ExperationForKr(const double& alpha, const double& deltax, double& a6);
double ExperationForKi(const double& alpha, const double& deltax, double& a6);
string GetCurrentPath();

//输出最优解结果
void result_out_put(const string& outpath, pair<vector<KData>, double> st,const double& deltax);

class Solution {
private:
    pair<vector<KData>, double> st;
    bool build;
    double deltaX;
    double deltaT;
    int iterationtimes;//迭代次数
    int Mesh_Number;//网格数
    //RungeKutta用的临时存储
    vector<double> middle, storage, para,pre;
public:
    Solution(pair<vector<KData>, double> para, double deltaX,double deltaT, int iterationtimes,int Mesh_Number);
    //获得差分格式的各个参数
    void Getpara();
    //顾名思义
    double LU(int& j);
    //精确解
    void ExtractSolINIT();
    //主体循环
    void loop();
    //顾名思义
    void RungeKutta();
    //结果输出
    void out_put(string& wavepath);
};
