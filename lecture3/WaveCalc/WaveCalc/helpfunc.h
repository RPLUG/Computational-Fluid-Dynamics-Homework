#pragma once
#include <math.h>
#include <iostream>
#include <vector>
#include "para.h"
#include <string>
#include <direct.h>
#include <windows.h>
#include <fstream>
#include <iomanip>
using namespace std;

string GetCurrentPath();

void SetCurrentPath(const string& path);

class NewTon {
	//this function only used in fomation like x^n
private:
	//index 1 ->left,index 2->right
	double Init_Time;
	double p_, p1, p2;//initial value
	double u1, u2, ui;//initial value 
	double rou1, rou2, roui;
	double roustar1, roustar2;//rou in region (3)->roustar1, (4)->roustar2
	double leftu, rightu;//the speed in region (3) & (4),vakue are same
	double Z1, Z2;//speed of shock wave 
	double Z1head, Z1tail;//if left side is expension wave
	double Z2head, Z2tail;//if right side is expension wave
	double c1star, c2star;//sound speed in region (3)& (4)
	static constexpr double gamma = 1.4;//for air gamma is 1.4
	double elip;//precision for NewTon
	bool rightex = true, leftex = true;//¼¤²¨true£¬ÅòÕÍfalse
	double c1, c2, ci;//c2->right side,  front of the wave 
	static constexpr double divnumb = 2000;//iteration times
	bool succ;//Weather successfully slove p*
	double Pos_Left_Shok_End, Pos_Right_Shok_End;
	static constexpr int Exp_Wave_Data_Density = 100;//data density
	string outroute;
	struct ExpWavePro
	{
		double Positon;
		double C_in_Exp;
		double P_in_Exp;
		double U_in_Exp;
		double rou_in_Exp;
	};
	vector<ExpWavePro> storage;//to store the data in expwave

private:

	pair<double, double> GetEqua(const double& p_, const double& pi, const double& roui, const double& ui, const double& ci);

	pair<double, double> GetDEqua(const double& p_, const double& pi, const double& roui, const double& ui, const double& ci);

	double GetA(const double& p_, const double& pi, const double& roui, const double& ci);

	void ExpWavePro(const double& Zhead, const double& Ztail);

public:
	NewTon(const double& time, const double& p1, const double& p2, const double& rou1, const double& rou2, const double& u1, const double& u2, const double& elip,const string& rote);

	double NewTonCalcu();

	pair<double, double> GetExpre(const double& p_, const double& pi, const double& roui, const double& ui, const double& ci, bool& status);

	void GetProperity(const double & p);

	void Result_out_put();

	~NewTon();
}; 
