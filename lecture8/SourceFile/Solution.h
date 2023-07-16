#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <direct.h>
#include "BasicData.h"
#include <fstream>
using namespace std;
class Solution_SCF
{
public:
	Solution_SCF(Para para);
	void Print_Res(string Out_route);
	void Calcu();
	void RK();
	void Cal_Omega();
	void SolveP();
	void StorageInit();
	void Cal_Relative_Number(); 
	void SolveUV();
	void check();
	~Solution_SCF();

private:
	const int Col = 0;
	const int Row = 0;
	const double Ite_Step = 0;
	const double End_Time = 0;
	const double DeltX;
	const double DeltY;
	const double DeltXX;
	const double DeltYY;
	const string Out_route;
	const bool PrintFlag;
	const bool Cal_P_Flag;
	const double Re;
	const double UU;
	const int Inter_ite_time;
	vector<vector<Data>> Storage;
};

Solution_SCF::Solution_SCF(Para para) 
	:Col(para.Xnumber + 5), Row(para.Ynumber + 5), Ite_Step(para.Ite_Step), End_Time(para.End_Time),
	Out_route(para.out_route), PrintFlag(para.Print), DeltX(para.DeltX), DeltY(para.DeltY),
	Cal_P_Flag(para.Cal_P_Flag), Re(para.Re),Inter_ite_time(para.Inter_ite_time),
	DeltXX(para.DeltX*para.DeltX),DeltYY(para.DeltY*para.DeltY),UU(400)
{
	cout << "Now you Are running a demo, if you want to designate more parameters, plz try to change source code" << endl;
	double time = 0;
	Storage.resize(Row, vector<Data>(Col));
	StorageInit();
	Calcu();
	int count = 0;
	while (time < End_Time)
	{
		time += Ite_Step;
		RK();
		if (count%100 == 0) {
			cout << "Now Iteration number is " << count <<" time step are " << time<< endl;
		}
		count++;
	}
	
	if (Cal_P_Flag) {
		SolveP();
	}
	if (PrintFlag) {
		Print_Res(Out_route);
	}
	else {
		cout << "Waring, Out put neglected" << endl;
	}
}

//初始化网格数据
void Solution_SCF::StorageInit() {
	//网格坐标给定
	for (int i = 0; i < Row; i++) {
		for (int j = 0; j < Col; j++) {
			Storage[i][j].x = j * DeltX;
			Storage[i][j].y = i * DeltY;
		}
	}

	//下表面
	for (int i = 0; i < 2; i++) {
		for (int j = 2; j < Col-3; j++) {
			Storage[i][j].Cur.u = Storage[i][j].Pre.u = 0;
			Storage[i][j].Cur.v = Storage[i][j].Pre.v = 0;
			Storage[i][j].Cur.Phy = 0;
		}
	}
	//for (int i = 0; i < 2; i++) {
	//	for (int j = 2; j < Col - 2; j++) {
	//		Cal_Omega(i, j);
	//	}
	//}
	
	//左表面
	for (int j = 0; j < 2; j++) {
		for (int i = 2; i < Row - 3; i++) {
			Storage[i][j].Cur.u = Storage[i][j].Pre.u = 0;
			Storage[i][j].Cur.v = Storage[i][j].Pre.v = 0;
			Storage[i][j].Cur.Phy = 0; 
		}
	}
	//for (int j = 0; j < 2; j++) {
	//	for (int i = 2; i < Row - 2; i++) {
	//		Cal_Omega(i, j);
	//	}
	//}
	
	//右表面
	for (int j = Col - 1, i = 2; i < Row - 3; i++) {
		Storage[i][j].Cur.u = Storage[i][j].Pre.u = 0;
		Storage[i][j].Cur.v = Storage[i][j].Pre.v = 0;
		Storage[i][j].Cur.Phy = 0;
	}
	//右侧剩下的两个计算
	for (int j = Col - 3; j < Col-1; j++) {
		for (int i = 2; i < Row - 3; i++) {
			Storage[i][j].Cur.u = Storage[i][j].Pre.u = 0;
			Storage[i][j].Cur.v = Storage[i][j].Pre.v = 0;
			Storage[i][j].Cur.Phy = 0;
		}
	}	
	//for (int j = Col - 3; j < Col - 1; j++) {
	//	for (int i = 2; i < Row - 3; i++) {
	//		Cal_Omega(i, j);
	//	}
	//}
	
	//上表面最上侧的不计算omega
	for (int j = 2, i = Row - 1; j < Col-3; j++) {
		Storage[i][j].Cur.u = Storage[i][j].Pre.u = UU;
		Storage[i][j].Cur.v = Storage[i][j].Pre.v = 0;
		Storage[i][j].Cur.Phy = 0;
	}
	//上表面
	for (int i = Row - 3; i < Row-1; i++) {
		for (int j = 2; j < Col - 3; j++) {
			Storage[i][j].Cur.u = Storage[i][j].Pre.u = UU;
			Storage[i][j].Cur.v = Storage[i][j].Pre.v = 0;
			Storage[i][j].Cur.Phy = 0;
		}
	}
	//for (int i = Row - 3; i < Row - 1; i++) {
	//	for (int j = 2; j < Col - 3; j++) {
	//		Cal_Omega(i, j);
	//	}
	//}
}

//求解i.j位置的omega
void Solution_SCF::Cal_Omega() {
	for (int i = Row - 3; i >= 2; i--) {
		for (int j = Col - 3; j >= 2; j--) {
			Storage[i][j].Cur.omega = (Storage[i+1][j].Cur.u - Storage[i][j].Cur.u) / ( DeltY)
				- (Storage[i][j+1].Cur.v - Storage[i][j].Cur.v) / ( DeltX);
			//Storage[i][j].Cur.omega = (Storage[i - 2][j].Cur.u - 6 * Storage[i - 1][j].Cur.u + 3 * Storage[i][j].Cur.u + 2 * Storage[i + 1][j].Cur.u) / (6 * DeltY)
			//	- (Storage[i][j - 2].Cur.u - 6 * Storage[i][j - 1].Cur.u + 3 * Storage[i][j].Cur.u + 2 * Storage[i][j + 1].Cur.u) / (6 * DeltX);
		}
	}
}

//求解一阶导，二阶导
void Solution_SCF::Cal_Relative_Number() {
	for (int i = 2; i < Row - 3; i++) {
		for (int j = 2; j < Col - 3; j++) {
			//对流项的离散采用三阶迎风，黏性项离散采用中心差分
			//后差
			Storage[i][j].Cur.Delt_M_OME_X = (Storage[i][j - 2].Cur.omega - 6 * Storage[i][j - 1].Cur.omega + 3 * Storage[i][j].Cur.omega + 2 * Storage[i][j + 1].Cur.omega) / (6 * DeltX);
			Storage[i][j].Cur.Delt_M_OME_Y = (Storage[i - 2][j].Cur.omega - 6 * Storage[i - 1][j].Cur.omega + 3 * Storage[i][j].Cur.omega + 2 * Storage[i + 1][j].Cur.omega) / (6 * DeltY);
			//前插
			Storage[i][j].Cur.Delt_P_OME_X = (-Storage[i][j + 2].Cur.omega + 6 * Storage[i][j + 1].Cur.omega - 3 * Storage[i][j].Cur.omega - 2 * Storage[i][j - 1].Cur.omega) / (6 * DeltX);
			Storage[i][j].Cur.Delt_P_OME_Y = (-Storage[i + 2][j].Cur.omega + 6 * Storage[i + 1][j].Cur.omega - 3 * Storage[i][j].Cur.omega - 2 * Storage[i - 1][j].Cur.omega) / (6 * DeltY);
	
			//黏性用的二阶采用中心差分
			Storage[i][j].Cur.Delt_OME_X_SEC = (Storage[i][j + 1].Cur.omega + Storage[i][j - 1].Cur.omega - 2 * Storage[i][j].Cur.omega) / (DeltXX);
			Storage[i][j].Cur.Delt_OME_Y_SEC = (Storage[i + 1][j].Cur.omega + Storage[i - 1][j].Cur.omega - 2 * Storage[i][j].Cur.omega) / (DeltYY);
			}
	}
	
}

//单轮计算
void Solution_SCF::Calcu() {
	//check();
	SolveUV();
	Cal_Omega();
	Cal_Relative_Number();
	//求解
	for (int i = 2; i < Row - 3; i++) {
		for (int j = 2; j < Col - 3; j++) {
			double up = (Storage[i][j].Cur.u + abs(Storage[i][j].Cur.u)) / 2;
			double ud = (Storage[i][j].Cur.u - abs(Storage[i][j].Cur.u)) / 2;
			double vp = (Storage[i][j].Cur.v + abs(Storage[i][j].Cur.v)) / 2;
			double vd = (Storage[i][j].Cur.v - abs(Storage[i][j].Cur.v)) / 2;
			Storage[i][j].Lu = -((up * Storage[i][j].Cur.Delt_M_OME_X + ud * Storage[i][j].Cur.Delt_P_OME_X)
								   + (vp * Storage[i][j].Cur.Delt_M_OME_Y + vd * Storage[i][j].Cur.Delt_P_OME_Y))
				+ (1.0 / Re) * (Storage[i][j].Cur.Delt_OME_X_SEC + Storage[i][j].Cur.Delt_OME_Y_SEC);
		}
	}
}

//RK推进
void Solution_SCF::RK() {
	for (int i = 2; i < Row - 3; i++) {
		for (int j = 2; j < Col - 3; j++) {
			Storage[i][j].Cur.omega = Storage[i][j].Cur.omega + Ite_Step * Storage[i][j].Lu;
		}
	}
	int count = 0;
	//Jocabi迭代，迭代次数由参数给定
	while (count++ < Inter_ite_time)
	{
		for (int i = 2; i < Row - 3; i++) {
			for (int j = 2; j < Col - 3; j++) {
				Storage[i][j].Cur.Phy = -(
					Storage[i][j].Cur.omega * DeltXX * DeltYY
					- DeltXX * Storage[i + 1][j].Pre.Phy
					- DeltXX * Storage[i - 1][j].Pre.Phy
					- DeltYY * Storage[i][j + 1].Pre.Phy
					- DeltYY * Storage[i][j - 1].Pre.Phy)
					/ ((2*DeltXX + 2*DeltYY));
			}
		}
		for (int i = 2; i < Row-3; i++) {
			for (int j = 2; j < Col-3; j++) {
				Storage[i][j].Pre.Phy = Storage[i][j].Cur.Phy;
			}
		}
	}
	Calcu();
}

//求解完成后计算速度
void Solution_SCF::SolveUV() {
	for (int i = 2; i < Row - 3; i++) {
		for (int j = 2; j < Col - 3; j++) {
			Storage[i][j].Cur.u = (Storage[i + 1][j].Cur.Phy - Storage[i - 1][j].Cur.Phy) / (2 * DeltY);
			Storage[i][j].Cur.v = -(Storage[i][j + 1].Cur.Phy - Storage[i][j - 1].Cur.Phy) / (2 * DeltX);

		}
	}
};

//析构函数
Solution_SCF::~Solution_SCF()
{
	cout << "Calculation End" << endl;
}

//待完成，计算求解压力
void Solution_SCF::SolveP() {
	for (int i = 0; i < Row; i++) {
		for (int j = 0; j < Col; j++) {
			//Through using Solution Method like LU-SGS.ADI,Gauss-seild et.al to solve this equation
			//Equation in the powerprint Chapter 11, Page 12
			Storage[i][j].P = 1;//remain to be done
		}
	}
}

//打印结果
void Solution_SCF::Print_Res(string Out_route) {
	ofstream File(Out_route+to_string(Re)+"XCoordinate" + ".txt");
	for (int j = 2,i=2; j < Col - 3; j++) {
		File << Storage[i][j].x<<"\n";
	}
	File.close();

	ofstream File1(Out_route + "Re" + to_string(Re) + "YCoordinate" + ".txt");
	for (int j = 2, i = 2; i < Row - 3; i++) {
		File1 << Storage[i][j].y << "\n";
	}
	File1.close();

	ofstream File2(Out_route + "Re" + to_string(Re) + "U" + ".txt");
	for (int i = 2; i < Row - 3; i++) {
		for (int j = 2; j < Col - 3; j++) {
			File2 << Storage[i][j].Cur.u << " ";
		}
		File2 << "\n";
	}
	File2.close();

	ofstream File3(Out_route + "Re" + to_string(Re) + "V" + ".txt");
	for (int i = 2; i < Row - 3; i++) {
		for (int j = 2; j < Col - 3; j++) {
			File3 << Storage[i][j].Cur.v << " ";
		}
		File3 <<"\n";
	}
	File3.close();
	cout << "Out Put Success" << endl;
	cout << "Attention, the out put file is in matrix foam, Using Matlab to deal the data" << endl;
	cout << " You can use matlab code like blow ↓" << endl;
	cout << "clear" << endl;
	cout << "clc" << endl;
	cout << "X = readmatrix('ResultRe400.000000X.txt')" << endl;
	cout << "Y = readmatrix('ResultRe400.000000Y.txt')" << endl;
	cout << "U = readmatrix('ResultRe400.000000U.txt')" << endl;
	cout << "V = readmatrix('ResultRe400.000000V.txt')" << endl;
	cout << "streamslice(X, Y, U, V)" << endl;
	return;
}

void Solution_SCF::check() {
	for (int i = 0; i < Row; i++) {
		for (int j = 0; j < Col; j++) {
			if (isnan(Storage[i][j].Cur.Delt_M_OME_X) || isnan(Storage[i][j].Cur.Delt_M_OME_Y) ||
				isnan(Storage[i][j].Cur.Delt_OME_X_SEC) || isnan(Storage[i][j].Cur.Delt_OME_Y_SEC) ||
				isnan(Storage[i][j].Cur.Delt_P_OME_X) || isnan(Storage[i][j].Cur.Delt_P_OME_Y) ||
				isnan(Storage[i][j].Cur.omega) || isnan(Storage[i][j].Cur.Phy) ||
				isnan(Storage[i][j].Cur.u) || isnan(Storage[i][j].Cur.v)) {
				cout << endl;
			}
		}
	}
}
string GetCurrentPath() {
	char curr_path[1024];
	char* nul;
	nul = _getcwd(curr_path, 1024);
	string res = curr_path;
	return res;
}