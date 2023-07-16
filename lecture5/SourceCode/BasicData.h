#pragma once
#include <vector>
using namespace std;

struct BasicData
{
	double RouL;
	double RouR;
	double HL;
	double HR;
	double EL;
	double ER;
	double PR;
	double PL;
	double UL;
	double UR;

	BasicData() {
		RouL = 0;
		RouR = 0;
		HL = 0;
		HR = 0;
		EL = 0;
		ER = 0;
		PR = 0;
		PL = 0;
		UL = 0;
		UR = 0;
	}
};


struct UAVData
{
	double U_;
	double c_;
	double Rou_;
	double H_;
	double P_;
	UAVData() {
		U_ = 0;
		c_ = 0;
		Rou_ = 0;
		H_ = 0;
		P_ = 0;
	}
};

struct HarfFData//用于记录半点处的值
{
	double F1;
	double F2;
	double F3;
	HarfFData() {
		F1 = 0;
		F2 = 0;
		F3 = 0;
	}
};

struct HarF_RES//用于记录半点处的值
{
	double F1;
	double F2;
	double F3;
	HarF_RES() {
		F1 = 0;
		F2 = 0;
		F3 = 0;
	}
};

struct WENO {
	double IS0;
	double IS1;
	double IS2;
	double Alpha0;
	double Alpha1;
	double Alpha2;
	double Omega0;
	double Omega1;
	double Omega2;
	WENO() {
		IS0 = 0;
		IS1 = 0;
		IS2 = 0;
		Alpha0 = 0;
		Alpha1 = 0;
		Alpha2 = 0;
		Omega0 = 0;
		Omega1 = 0;
		Omega2 = 0;
	}
};

struct Data //半点值，包含当前点的左右情况，以及当前点的平均值
{
	HarF_RES FRES;
	//WENO WeL, WeR;
	HarfFData HarFL,HarFR;
	BasicData BaData;
	UAVData U;
};


struct UFData//网格点上的值
{
	double Rou;
	double U;
	double P;
	double C;
	double E;
	double RouU;
	double RouE;
	double F1;
	double F2;
	double F3;
	UFData() {
		
		Rou = 0;
		U = 0;
		P = 0;
		C = 0;
		E = 0;
		RouU = 0;
		RouE = 0;
		F1 = 0;
		F2 = 0;
		F3 = 0;
	}
};

struct temp
{
	double F1;
	double F2;
	double F3;
	temp() {
		F1 = 0;
		F2 = 0;
		F3 = 0;
	}
};