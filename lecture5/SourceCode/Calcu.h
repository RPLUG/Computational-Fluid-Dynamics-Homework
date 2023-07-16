#pragma once
#include <vector>
#include <math.h>
#include <string>
#include <conio.h>
#include <fstream>
#include <direct.h>
#include <iomanip>
#include <iostream>
#include "BasicData.h"

using namespace std;

string GetCurrentPath();

class Solution {
private:
	static constexpr double Gamma = 1.4;

	inline void UpData(vector<UFData>& Source);

	inline void ExactSol(int& Tcase);

	void GetUAVData();

	void GetModLambda();

	void Get_S_Inv_S();

	vector<vector<double>> _SASU(int i);

	void JacobinA_U();

	void WenoC(const int& i, const char ch);

	void Method();

	void GetHalfPoint();

	void RungeKutta();

	void out_put(const string& rote);

	vector<vector<double>> Matrix_Multiplication(const vector<vector<double>>& MatrixA, const vector<vector<double>>& MatrixB);

	int Mesh_Number, PIODT, PIEM;

	string outpath;

	double deltaX;

	int Ite_Tim;

	double Ite_Step, Beg_pos, End_pos;

	bool Exc;
	vector<Data> Storage;//this storage is for point in 1/2;

	vector<UFData> StorageUF, StorageV; //this storage is for point in 1.2.3.4;

	vector<vector<vector<double>>> Inv_S, _S, ModLambda;//in this means ¡Ä

	vector<temp> Iteres;//just for rk iteration

	int testcase;//which case do you want to choose 
public:
	Solution(const int& Number_Mesh, const double& Beg_pos, const double& End_pos, const int& IterationTimes, const double& step, const string& rote, const bool& Ex, const int& Tcase,const string& path);

	void Calcu();

};
