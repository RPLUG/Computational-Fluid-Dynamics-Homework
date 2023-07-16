#pragma once
#include <string.h>
#include <string>
using namespace std;
struct NData {
	//ÿ��������ٶȺ�omega
	double u;
	double v;
	double omega;
	//�����ñ���������Ħ�
	double Delt_P_OME_X;//��omega in U directionǰ��һ�׵���
	double Delt_P_OME_Y;//��omega in V direction
	double Delt_M_OME_X;//��omega in U direction���һ�׵���
	double Delt_M_OME_Y;//��omega in V direction
	//������������
	double Delt_OME_X_SEC;//��omega in U direction���׵�
	double Delt_OME_Y_SEC;//��omega in V direction

	double Phy;
	NData() {
		u = 0;
		v = 0;
		omega = 0;
		Delt_P_OME_X = 0;
		Delt_P_OME_Y = 0;
		Delt_M_OME_X = 0;
		Delt_M_OME_Y = 0;
		Delt_OME_X_SEC = 0;
		Delt_OME_Y_SEC = 0;
		Phy = 0;
	}
};
struct Data {
	double x;
	double y;
	double P;
	double Lu;
	NData Pre, Cur;
	Data() {
		x = 0;
		y = 0;
		P = 0;
		Lu = 0;;
	}
};
struct Para
{
	int Xnumber;
	int Ynumber;
	double DeltX;
	double DeltY;
	double Ite_Step;
	double End_Time;
	string out_route;
	bool Print;
	bool Cal_P_Flag;
	int Inter_ite_time;
	double Re;
	Para() {
		Xnumber = 0;
		Ynumber = 0;
		Ite_Step = 0;
		End_Time = 0;
		Inter_ite_time = 0;
		Re = 0;
		out_route = "";
		Print = false;
		Cal_P_Flag = false;
		DeltY = 0;
	}
};