#pragma once
#include <vector>
#include <stdio.h>
#include <math.h>
#include "DataS.h"
#include <fstream>

using namespace std;


class MeshGenerate
{
public:
	MeshGenerate(para para);
	void Mesh_Init();
	~MeshGenerate();
	void SloveAB(double& A, double& B, const double& hb, const double& he,const int& n);
	void Print_res();
	void LUSGS();
	double LUSGSCalcu(int& i, int& j,const char& Direction);	
	void JKB();
	void JKBCalcu(int& i, int& j);
private:
	struct Coordinate {
		double x;
		double y;
		bool boundary;
		double prex;
		double prey;
		Coordinate() {
			x = 0;
			y = 0;
			prex = 0;
			prey = 0;
			boundary = 0;
		}
	};
	string Rote;
	const int Row,Col;
	const int BEnum;
	const double Lside = 0;
	const double Rside = 0;
	const double FlyBe = 0;
	const double center = 0;
	const double H1b,H1e;
	const double H2b,H2e;
	const double H3b,H3e;
	const double HenStart;
	const double MaxIteration;
	const int Tole;
	double A1 = 0,B1 = 0;
	double B2 = 0,A2 = 0;
	double B3 = 0,A3 = 0;
	vector<vector<Coordinate>> Org_Coor;
};


string GetCurrentPath();