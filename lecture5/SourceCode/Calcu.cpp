#include "Calcu.h"

inline void Solution::UpData(vector<UFData>& Source) {
	//这里已知的是Rou，Rou*U和E，该函数用来更新剩余的点值
	for (int i = 0; i < Source.size(); i++) {
		Source[i].U = Source[i].RouU / Source[i].Rou;
		Source[i].RouE = Source[i].Rou * Source[i].E;
		Source[i].P = (Gamma - 1) * ((Source[i].E) - 0.5 * (Source[i].Rou * Source[i].U * Source[i].U));
		Source[i].C = sqrt(Gamma * Source[i].P / Source[i].Rou);
		Source[i].F1 = Source[i].RouU;
		Source[i].F2 = Source[i].Rou * Source[i].U * Source[i].U + Source[i].P;
		Source[i].F3 = (Source[i].E + Source[i].P) * Source[i].U;
	}
}

inline void Solution::ExactSol(int& Tcase) {
	if (Tcase == 1) {
		for (int i = 0; i < PIEM * 0.1; i++) {
			StorageUF[i].U = 2.629;
			StorageUF[i].P = 10.333;
			StorageUF[i].Rou = 3.857;
			StorageUF[i].C = sqrt(Gamma * StorageUF[i].P / StorageUF[i].Rou);
			StorageUF[i].E = StorageUF[i].P / (Gamma - 1) + 0.5 * StorageUF[i].Rou * StorageUF[i].U * StorageUF[i].U;
			StorageUF[i].RouU = StorageUF[i].Rou * StorageUF[i].U;
			StorageUF[i].RouE = StorageUF[i].E * StorageUF[i].Rou;
			StorageUF[i].F1 = StorageUF[i].RouU;
			StorageUF[i].F2 = StorageUF[i].Rou * StorageUF[i].U * StorageUF[i].U + StorageUF[i].P;
			StorageUF[i].F3 = (StorageUF[i].E + StorageUF[i].P) * StorageUF[i].U;
		}
		for (int i = PIEM * 0.1; i < PIEM; i++) {
			StorageUF[i].U = 0;
			StorageUF[i].P = 1;
			StorageUF[i].Rou = 1 + 0.2 * sin(5 * deltaX * i);
			StorageUF[i].C = sqrt(Gamma * StorageUF[i].P / StorageUF[i].Rou);
			StorageUF[i].E = StorageUF[i].P / (Gamma - 1) + 0.5 * StorageUF[i].Rou * StorageUF[i].U * StorageUF[i].U;
			StorageUF[i].RouU = StorageUF[i].Rou * StorageUF[i].U;
			StorageUF[i].RouE = StorageUF[i].Rou * StorageUF[i].E;
			StorageUF[i].F1 = StorageUF[i].RouU;
			StorageUF[i].F2 = StorageUF[i].Rou * StorageUF[i].U * StorageUF[i].U + StorageUF[i].P;
			StorageUF[i].F3 = (StorageUF[i].E + StorageUF[i].P) * StorageUF[i].U;
		}
	}
	else if (Tcase == 2) {
		for (int i = 0; i < PIEM / 2; i++) {
			StorageUF[i].U = 0;
			StorageUF[i].P = 1;
			StorageUF[i].Rou = 1;
			StorageUF[i].C = sqrt(Gamma * StorageUF[i].P / StorageUF[i].Rou);
			StorageUF[i].E = StorageUF[i].P / (Gamma - 1) + 0.5 * StorageUF[i].Rou * StorageUF[i].U * StorageUF[i].U;
			StorageUF[i].RouU = StorageUF[i].Rou * StorageUF[i].U;
			StorageUF[i].RouE = StorageUF[i].E * StorageUF[i].Rou;
			StorageUF[i].F1 = StorageUF[i].RouU;
			StorageUF[i].F2 = StorageUF[i].Rou * StorageUF[i].U * StorageUF[i].U + StorageUF[i].P;
			StorageUF[i].F3 = (StorageUF[i].E + StorageUF[i].P) * StorageUF[i].U;
		}
		for (int i = PIEM / 2; i < PIEM; i++) {
			StorageUF[i].U = 0;
			StorageUF[i].P = 0.1;
			StorageUF[i].Rou = 0.125;
			StorageUF[i].C = sqrt(Gamma * StorageUF[i].P / StorageUF[i].Rou);
			StorageUF[i].E = StorageUF[i].P / (Gamma - 1) + 0.5 * StorageUF[i].Rou * StorageUF[i].U * StorageUF[i].U;
			StorageUF[i].RouU = StorageUF[i].Rou * StorageUF[i].U;
			StorageUF[i].RouE = StorageUF[i].Rou * StorageUF[i].E;
			StorageUF[i].F1 = StorageUF[i].RouU;
			StorageUF[i].F2 = StorageUF[i].Rou * StorageUF[i].U * StorageUF[i].U + StorageUF[i].P;
			StorageUF[i].F3 = (StorageUF[i].E + StorageUF[i].P) * StorageUF[i].U;
		}
	}

};

void Solution::GetUAVData() {
	for (int i = 0; i < PIODT; i++) {
		double sqrouL = sqrt(Storage[i].BaData.RouL);
		double sqrouR = sqrt(Storage[i].BaData.RouR);
		Storage[i].U.Rou_ = pow(((sqrouL + sqrouR) / 2), 2);
		Storage[i].U.U_ = (sqrouL * Storage[i].BaData.UL + sqrouR * Storage[i].BaData.UR) / (sqrouL + sqrouR);
		Storage[i].U.H_ = (sqrouL * Storage[i].BaData.HL + sqrouR * Storage[i].BaData.HR) / (sqrouL + sqrouR);
		double uu = pow(Storage[i].U.U_, 2);
		Storage[i].U.c_ = sqrt((Gamma - 1) * (Storage[i].U.H_ - 0.5 * uu));
		Storage[i].U.P_ = ((Gamma - 1) / (Gamma)) * (Storage[i].U.Rou_ * Storage[i].U.H_ - 0.5 * Storage[i].U.Rou_ * uu);
	}
	return;
};

void Solution::GetModLambda() {
	vector<vector<double>> res(3, vector<double>(3, 0));
	double tolerence = 0.0001;
	for (int i = 0; i < PIODT; i++) {
		res[0][0] = fabs(Storage[i].U.U_) > tolerence ? fabs(Storage[i].U.U_) : ((Storage[i].U.U_) * (Storage[i].U.U_) + tolerence * tolerence) / (2 * tolerence);
		res[1][1] = fabs(Storage[i].U.U_ - Storage[i].U.c_) > tolerence ? fabs(Storage[i].U.U_ - Storage[i].U.c_) : (((Storage[i].U.U_ - Storage[i].U.c_) * (Storage[i].U.U_ - Storage[i].U.c_) + tolerence * tolerence) / (2 * tolerence));
		res[2][2] = fabs(Storage[i].U.U_ + Storage[i].U.c_) > tolerence ? fabs(Storage[i].U.U_ + Storage[i].U.c_) : (((Storage[i].U.U_ + Storage[i].U.c_) * (Storage[i].U.U_ + Storage[i].U.c_) + tolerence * tolerence) / (2 * tolerence));
		ModLambda[i] = res;
	}
	return;
};

void Solution::Get_S_Inv_S() {
	vector<vector<double>> res2(3, vector<double>(3, 0));
	vector<vector<double>> res(3, vector<double>(3, 0));
	for (int i = 0; i < PIODT; i++) {
		double cc = pow(Storage[i].U.c_, 2);
		double uu = pow(Storage[i].U.U_, 2);
		double u = Storage[i].U.U_;
		double c = Storage[i].U.c_;
		double uc = Storage[i].U.U_ * Storage[i].U.c_;
		res[0][0] = (uu / 2) - (cc / (Gamma - 1));
		res[0][1] = -u;
		res[0][2] = 1.0;
		res[1][0] = -u - ((Gamma - 1) / c) * (uu / 2);
		res[1][1] = 1 + ((Gamma - 1) / c) * u;
		res[1][2] = -(Gamma - 1) / c;
		res[2][0] = -u + ((Gamma - 1) / c) * (uu / 2);
		res[2][1] = 1 - ((Gamma - 1) / c) * u;
		res[2][2] = (Gamma - 1) / c;
		_S[i] = res;
		res2[0][0] = -(Gamma - 1) / cc;
		res2[0][1] = -1.0 / (2 * c);
		res2[0][2] = 1.0 / (2 * c);
		res2[1][0] = (-(Gamma - 1.0) * u) / cc;
		res2[1][1] = -(u - c) / (2 * c);
		res2[1][2] = (c + u) / (2 * c);
		res2[2][0] = (-(Gamma - 1) / cc) * (uu / 2);
		res2[2][1] = (-1.0 / (2 * c)) * ((uu / 2) + (cc / (Gamma - 1)) - (uc));
		res2[2][2] = (1.0 / (2 * c)) * ((uu / 2) + (cc / (Gamma - 1)) + (uc));
		Inv_S[i] = res2;
	}
};

vector<vector<double>> Solution::_SASU(int i) {
	vector<vector<double>> res;
	res = Matrix_Multiplication(Inv_S[i], ModLambda[i]);
	res = Matrix_Multiplication(res, _S[i]);
	double _rou = Storage[i].BaData.RouR - Storage[i].BaData.RouL;
	double _rouU = Storage[i].BaData.RouR * Storage[i].BaData.UR - Storage[i].BaData.RouL * Storage[i].BaData.UL;
	double _E = Storage[i].BaData.ER - Storage[i].BaData.EL;
	vector<vector<double>> temp = { {_rou},{_rouU},{_E} };
	res = Matrix_Multiplication(res, temp);
	return res;
}

void Solution::WenoC(const int& i, const char ch) {
	/*if (i < 0 || i > PIEM - 2) {
		cout << "wrong input i in func WenoC" << endl;
	}*/
	bool flag0 = false;
	bool flagM1 = false;
	bool flagM2 = false;
	bool flagM3 = false;
	bool flagP1 = false;
	bool flagP2 = false;
	bool flagP3 = false; {}
	double vp1, vp2, vp3, v0, vm1, vm2, IS0, IS1, IS2, Alpha0, Alpha1, Alpha2, Omega0, Omega1, Omega2, f0, f1, f2;
	if (ch == 'R') {
		flag0 = (i<0 || i>PIEM - 1) ? false : true;
		flagM1 = (i - 1 < 0) ? false : true;
		flagP1 = (i > PIEM - 2) ? false : true;
		flagP2 = (i > PIEM - 3) ? false : true;
		flagP3 = (i > PIEM - 4) ? false : true;

		/*Calcu Rou*/
		if (flag0) {
			v0 = StorageUF[i].Rou;
		}
		if (flagM1) {
			vm1 = StorageUF[i - 1].Rou;
		}
		if (flagP1) {
			vp1 = StorageUF[i + 1].Rou;
		}
		if (flagP2) {
			vp2 = StorageUF[i + 2].Rou;
		}
		if (flagP3) {
			vp3 = StorageUF[i + 3].Rou;
		}
		if (flagP1 && flag0 && flagM1) {
			IS0 = (13.0 / 12) * pow((vp1 - 2 * v0 + vm1), 2) + 0.25 * pow((3 * vp1 - 4 * v0 + vm1), 2);
			Alpha0 = 0.3 / pow((10E-6 + IS0), 2);
			f0 = (1.0 / 3) * vp1 + (5.0 / 6) * v0 - (1.0 / 6) * vm1;
		}
		else {
			IS0 = 0;
			Alpha0 = 0;
			f0 = 0;
		}
		if (flagP2 && flagP1 && flag0) {
			IS1 = (13.0 / 12) * pow((vp2 - 2 * vp1 + v0), 2) + 0.25 * pow((vp2 - v0), 2);
			Alpha1 = 0.6 / pow((10E-6 + IS1), 2);
			f1 = (-1.0 / 6) * vp2 + (5.0 / 6) * vp1 + (1.0 / 3) * v0;
		}
		else {
			IS1 = 0;
			Alpha1 = 0;
			f1 = 0;
		}
		if (flagP3 && flagP2 && flagP1) {
			IS2 = (13.0 / 12) * pow((vp3 - 2 * vp2 + vp1), 2) + 0.25 * pow((vp3 - 4 * vp2 + 3 * vp1), 2);
			Alpha2 = 0.1 / pow((10E-6 + IS2), 2);
			f2 = (1.0 / 3) * vp3 - (7.0 / 6) * vp2 + (11.0 / 6) * vp1;
		}
		else {
			IS2 = 0;
			Alpha2 = 0;
			f2 = 0;
		}
		Omega0 = Alpha0 / (Alpha0 + Alpha1 + Alpha2);
		Omega1 = Alpha1 / (Alpha0 + Alpha1 + Alpha2);
		Omega2 = Alpha2 / (Alpha0 + Alpha1 + Alpha2);
		Storage[i].BaData.RouR = Omega0 * f0 + Omega1 * f1 + Omega2 * f2;

		/*Calcu U*/
		if (flag0) {
			v0 = StorageUF[i].U;
		}
		if (flagM1) {
			vm1 = StorageUF[i - 1].U;
		}
		if (flagP1) {
			vp1 = StorageUF[i + 1].U;
		}
		if (flagP2) {
			vp2 = StorageUF[i + 2].U;
		}
		if (flagP3) {
			vp3 = StorageUF[i + 3].U;
		}
		if (flagM1 && flag0 && flagP1) {
			IS0 = (13.0 / 12) * pow((vp1 - 2 * v0 + vm1), 2) + 0.25 * pow((3 * vp1 - 4 * v0 + vm1), 2);
			Alpha0 = 0.3 / pow((10E-6 + IS0), 2);
			f0 = (1.0 / 3) * vp1 + (5.0 / 6) * v0 - (1.0 / 6) * vm1;
		}
		else {
			IS0 = 0;
			Alpha0 = 0;
			f0 = 0;
		}
		if (flagP2 && flagP1 && flag0) {
			IS1 = (13.0 / 12) * pow((vp2 - 2 * vp1 + v0), 2) + 0.25 * pow((vp2 - v0), 2);
			Alpha1 = 0.6 / pow((10E-6 + IS1), 2);
			f1 = (-1.0 / 6) * vp2 + (5.0 / 6) * vp1 + (1.0 / 3) * v0;
		}
		else {
			IS1 = 0;
			Alpha1 = 0;
			f1 = 0;
		}
		if (flagP3 && flagP2 && flagP1) {
			IS2 = (13.0 / 12) * pow((vp3 - 2 * vp2 + vp1), 2) + 0.25 * pow((vp3 - 4 * vp2 + 3 * vp1), 2);
			Alpha2 = 0.1 / pow((10E-6 + IS2), 2);
			f2 = (1.0 / 3) * vp3 - (7.0 / 6) * vp2 + (11.0 / 6) * vp1;
		}
		else {
			IS2 = 0;
			Alpha2 = 0;
			f2 = 0;
		}
		Omega0 = Alpha0 / (Alpha0 + Alpha1 + Alpha2);
		Omega1 = Alpha1 / (Alpha0 + Alpha1 + Alpha2);
		Omega2 = Alpha2 / (Alpha0 + Alpha1 + Alpha2);

		Storage[i].BaData.UR = Omega0 * f0 + Omega1 * f1 + Omega2 * f2;

		/* Calcu (E+p)*u */
		if (flag0) {
			v0 = (StorageUF[i].E);
		}
		if (flagM1) {
			vm1 = (StorageUF[i - 1].E);
		}
		if (flagP1) {
			vp1 = (StorageUF[i + 1].E);
		}
		if (flagP2) {
			vp2 = (StorageUF[i + 2].E);
		}
		if (flagP3) {
			vp3 = (StorageUF[i + 3].E);
		}
		if (flagP1 && flag0 && flagM1) {
			IS0 = (13.0 / 12) * pow((vp1 - 2 * v0 + vm1), 2) + 0.25 * pow((3 * vp1 - 4 * v0 + vm1), 2);
			Alpha0 = 0.3 / pow((10E-6 + IS0), 2);
			f0 = (1.0 / 3) * vp1 + (5.0 / 6) * v0 - (1.0 / 6) * vm1;
		}
		else {
			IS0 = 0;
			Alpha0 = 0;
			f0 = 0;
		}
		if (flagP2 && flagP1 && flag0) {
			IS1 = (13.0 / 12) * pow((vp2 - 2 * vp1 + v0), 2) + 0.25 * pow((vp2 - v0), 2);
			Alpha1 = 0.6 / pow((10E-6 + IS1), 2);
			f1 = (-1.0 / 6) * vp2 + (5.0 / 6) * vp1 + (1.0 / 3) * v0;
		}
		else {
			IS1 = 0;
			Alpha1 = 0;
			f1 = 0;
		}
		if (flagP3 && flagP2 && flagP1) {
			IS2 = (13.0 / 12) * pow((vp3 - 2 * vp2 + vp1), 2) + 0.25 * pow((vp3 - 4 * vp2 + 3 * vp1), 2);
			Alpha2 = 0.1 / pow((10E-6 + IS2), 2);
			f2 = (1.0 / 3) * vp3 - (7.0 / 6) * vp2 + (11.0 / 6) * vp1;
		}
		else {
			IS2 = 0;
			Alpha2 = 0;
			f2 = 0;
		}
		Omega0 = Alpha0 / (Alpha0 + Alpha1 + Alpha2);
		Omega1 = Alpha1 / (Alpha0 + Alpha1 + Alpha2);
		Omega2 = Alpha2 / (Alpha0 + Alpha1 + Alpha2);

		Storage[i].BaData.ER = Omega0 * f0 + Omega1 * f1 + Omega2 * f2;
		/*updata UR,RouR....etl*/
		Storage[i].BaData.PR = (Gamma - 1) * (Storage[i].BaData.ER - 0.5 * Storage[i].BaData.RouR * Storage[i].BaData.UR * Storage[i].BaData.UR);
		Storage[i].BaData.HR = (Storage[i].BaData.ER + Storage[i].BaData.PR) / Storage[i].BaData.RouR;
		Storage[i].HarFL.F1 = Storage[i].BaData.UR * Storage[i].BaData.RouR;
		Storage[i].HarFL.F2 = Storage[i].BaData.UR * Storage[i].BaData.UR * Storage[i].BaData.RouR + Storage[i].BaData.PR;
		Storage[i].HarFL.F3 = (Storage[i].BaData.ER + Storage[i].BaData.PR) * Storage[i].BaData.UR;

	}

	else if (ch == 'L') {
		flag0 = (i<0 || i>PIEM - 1) ? false : true;
		flagM1 = (i < 1) ? false : true;
		flagM2 = (i < 2) ? false : true;
		flagP1 = (i > PIEM - 2) ? false : true;
		flagP2 = (i > PIEM - 3) ? false : true;
		flagP3 = (i > PIEM - 4) ? false : true;
		/*Calcu Rou*U*/
		if (flag0) {
			v0 = StorageUF[i].Rou;
		}
		if (flagM1) {
			vm1 = StorageUF[i - 1].Rou;
		}
		if (flagP1) {
			vp1 = StorageUF[i + 1].Rou;
		}
		if (flagP2) {
			vp2 = StorageUF[i + 2].Rou;
		}
		if (flagM2) {
			vm2 = StorageUF[i - 2].Rou;
		}
		if (flagP3) {
			vp3 = StorageUF[i + 3].Rou;
		}
		if (flagP1 && flag0 && flagP2) {
			IS0 = (13.0 / 12) * pow((v0 - 2 * vp1 + vp2), 2) + 0.25 * pow((3 * v0 - 4 * vp1 + vp2), 2);
			Alpha0 = 0.3 / pow((10E-6 + IS0), 2);
			f0 = (1.0 / 3) * v0 + (5.0 / 6) * vp1 - (1.0 / 6) * vp2;
		}
		else {
			IS0 = 0;
			Alpha0 = 0;
			f0 = 0;
		}
		if (flagM1 && flag0 && flagP1) {
			IS1 = (13.0 / 12) * pow((vm1 - 2 * v0 + vp1), 2) + 0.25 * pow((vm1 - vp1), 2);
			Alpha1 = 0.6 / pow((10E-6 + IS1), 2);
			f1 = (-1.0 / 6) * vm1 + (5.0 / 6) * v0 + (1.0 / 3) * vp1;
		}
		else {
			IS1 = 0;
			Alpha1 = 0;
			f1 = 0;
		}
		if (flag0 && flagM2 && flagM1) {
			IS2 = (13.0 / 12) * pow((vm2 - 2 * vm1 + v0), 2) + 0.25 * pow((vm2 - 4 * vm1 + 3 * v0), 2);
			Alpha2 = 0.1 / pow((10E-6 + IS2), 2);
			f2 = (1.0 / 3) * vm2 - (7.0 / 6) * vm1 + (11.0 / 6) * v0;
		}
		else {
			IS2 = 0;
			Alpha2 = 0;
			f2 = 0;
		}
		Omega0 = Alpha0 / (Alpha0 + Alpha1 + Alpha2);
		Omega1 = Alpha1 / (Alpha0 + Alpha1 + Alpha2);
		Omega2 = Alpha2 / (Alpha0 + Alpha1 + Alpha2);
		Storage[i].BaData.RouL = Omega0 * f0 + Omega1 * f1 + Omega2 * f2;

		/*Calcu Rou*U*U+p*/
		if (flag0) {
			v0 = StorageUF[i].U;
		}
		if (flagM1) {
			vm1 = StorageUF[i - 1].U;
		}
		if (flagP1) {
			vp1 = StorageUF[i + 1].U;
		}
		if (flagP2) {
			vp2 = StorageUF[i + 2].U;
		}
		if (flagP3) {
			vp3 = StorageUF[i + 3].U;
		}
		if (flag0 && flagP1 && flagP2) {
			IS0 = (13.0 / 12) * pow((v0 - 2 * vp1 + vp2), 2) + 0.25 * pow((3 * v0 - 4 * vp1 + vp2), 2);
			Alpha0 = 0.3 / pow((10E-6 + IS0), 2);
			f0 = (1.0 / 3) * v0 + (5.0 / 6) * vp1 - (1.0 / 6) * vp2;
		}
		else {
			IS0 = 0;
			Alpha0 = 0;
			f0 = 0;
		}
		if (flag0 && flagP1 && flagM1) {
			IS1 = (13.0 / 12) * pow((vm1 - 2 * v0 + vp1), 2) + 0.25 * pow((vm1 - vp1), 2);
			Alpha1 = 0.6 / pow((10E-6 + IS1), 2);
			f1 = (-1.0 / 6) * vm1 + (5.0 / 6) * v0 + (1.0 / 3) * vp1;
		}
		else {
			IS1 = 0;
			Alpha1 = 0;
			f1 = 0;
		}
		if (flagM2 && flagM1 && flag0) {
			IS2 = (13.0 / 12) * pow((vm2 - 2 * vm1 + v0), 2) + 0.25 * pow((vm2 - 4 * vm1 + 3 * v0), 2);
			Alpha2 = 0.1 / pow((10E-6 + IS2), 2);
			f2 = (1.0 / 3) * vm2 - (7.0 / 6) * vm1 + (11.0 / 6) * v0;
		}
		else {
			IS2 = 0;
			Alpha2 = 0;
			f2 = 0;
		}
		Omega0 = Alpha0 / (Alpha0 + Alpha1 + Alpha2);
		Omega1 = Alpha1 / (Alpha0 + Alpha1 + Alpha2);
		Omega2 = Alpha2 / (Alpha0 + Alpha1 + Alpha2);

		Storage[i].BaData.UL = Omega0 * f0 + Omega1 * f1 + Omega2 * f2;

		/*Calcu Rou*U*/
		if (flag0) {
			v0 = (StorageUF[i].E);
		}
		if (flagM1) {
			vm1 = (StorageUF[i - 1].E);
		}
		if (flagP1) {
			vp1 = (StorageUF[i + 1].E);
		}
		if (flagP2) {
			vp2 = (StorageUF[i + 2].E);
		}
		if (flagP3) {
			vp3 = (StorageUF[i + 3].E);
		}
		if (flag0 && flagP1 && flagP2) {
			IS0 = (13.0 / 12) * pow((v0 - 2 * vp1 + vp2), 2) + 0.25 * pow((3 * v0 - 4 * vp1 + vp2), 2);
			Alpha0 = 0.3 / pow((10E-6 + IS0), 2);
			f0 = (1.0 / 3) * v0 + (5.0 / 6) * vp1 - (1.0 / 6) * vp2;
		}
		else {
			IS0 = 0;
			Alpha0 = 0;
			f0 = 0;
		}
		if (flag0 && flagM1 && flagP1) {
			IS1 = (13.0 / 12) * pow((vm1 - 2 * v0 + vp1), 2) + 0.25 * pow((vm1 - vp1), 2);
			Alpha1 = 0.6 / pow((10E-6 + IS1), 2);
			f1 = (-1.0 / 6) * vm1 + (5.0 / 6) * v0 + (1.0 / 3) * vp1;
		}
		else {
			IS1 = 0;
			Alpha1 = 0;
			f1 = 0;
		}
		if (flagM1 && flagM2 && flag0) {
			IS2 = (13.0 / 12) * pow((vm2 - 2 * vm1 + v0), 2) + 0.25 * pow((vm2 - 4 * vm1 + 3 * v0), 2);
			Alpha2 = 0.1 / pow((10E-6 + IS2), 2);
			f2 = (1.0 / 3) * vm2 - (7.0 / 6) * vm1 + (11.0 / 6) * v0;
		}
		else {
			IS2 = 0;
			Alpha2 = 0;
			f2 = 0;
		}
		Omega0 = Alpha0 / (Alpha0 + Alpha1 + Alpha2);
		Omega1 = Alpha1 / (Alpha0 + Alpha1 + Alpha2);
		Omega2 = Alpha2 / (Alpha0 + Alpha1 + Alpha2);

		Storage[i].BaData.EL = Omega0 * f0 + Omega1 * f1 + Omega2 * f2;

		Storage[i].BaData.PL = (Gamma - 1) * (Storage[i].BaData.EL - 0.5 * Storage[i].BaData.RouL * Storage[i].BaData.UL * Storage[i].BaData.UL);
		Storage[i].BaData.HL = (Storage[i].BaData.EL + Storage[i].BaData.PL) / Storage[i].BaData.RouL;
		Storage[i].HarFL.F1 = Storage[i].BaData.UL * Storage[i].BaData.RouL;
		Storage[i].HarFL.F2 = Storage[i].BaData.UL * Storage[i].BaData.UL * Storage[i].BaData.RouL + Storage[i].BaData.PL;
		Storage[i].HarFL.F3 = (Storage[i].BaData.EL + Storage[i].BaData.PL) * Storage[i].BaData.UL;
	}
}

void Solution::Method()//使用点值更新半点值 
{
	for (int i = 0; i < PIODT; i++) {
		WenoC(i, 'L');
		for (int i = 0; i < PIODT; i++) {
			WenoC(i, 'R');
		}
	}
}

void Solution::JacobinA_U() {
	Get_S_Inv_S();
	GetModLambda();
}

void Solution::GetHalfPoint()//用来得到半点处结果的值并更新ite
{
	Method();
	GetUAVData();
	JacobinA_U();
	for (int i = 0; i < PIODT; i++) {//更新半点结果
		vector<vector<double>> temp = _SASU(i);
		Storage[i].FRES.F1 = 0.5 * (Storage[i].HarFR.F1 + Storage[i].HarFL.F1) - 0.5 * temp[0][0];
		Storage[i].FRES.F2 = 0.5 * (Storage[i].HarFR.F2 + Storage[i].HarFL.F2) - 0.5 * temp[1][0];
		Storage[i].FRES.F3 = 0.5 * (Storage[i].HarFR.F3 + Storage[i].HarFL.F3) - 0.5 * temp[2][0];
	}

	for (int i = 0; i < Mesh_Number - 2; i++) {
		Iteres[i].F1 = -(Storage[i + 1].FRES.F1 - Storage[i].FRES.F1) / deltaX;
		Iteres[i].F2 = -(Storage[i + 1].FRES.F2 - Storage[i].FRES.F2) / deltaX;
		Iteres[i].F3 = -(Storage[i + 1].FRES.F3 - Storage[i].FRES.F3) / deltaX;
	}
}

void Solution::RungeKutta() {
	vector<UFData> Utemp;//storageU是当前的
	vector<UFData> UtempN = StorageUF;
	vector<UFData> rec = StorageUF;//用来记录当前的点值

	for (int i = 1; i < Mesh_Number - 1; i++) {
		UtempN[i].Rou = rec[i].Rou + Ite_Step * Iteres[i - 1].F1;
		UtempN[i].RouU = rec[i].RouU + Ite_Step * Iteres[i - 1].F2;
		UtempN[i].E = rec[i].E + Ite_Step * Iteres[i - 1].F3;
	}
	UpData(UtempN);//用来更新点值
	Utemp = UtempN;
	StorageUF = UtempN;
	GetHalfPoint();
	/*updata*/
	for (int i = 1; i < Mesh_Number - 1; i++) {
		UtempN[i].Rou = (3.0 / 4) * rec[i].Rou + (0.25) * Utemp[i].Rou + 0.25 * Ite_Step * Iteres[i - 1].F1;
		UtempN[i].RouU = (3.0 / 4) * rec[i].RouU + (0.25) * Utemp[i].RouU + 0.25 * Ite_Step * Iteres[i - 1].F2;
		UtempN[i].E = (3.0 / 4) * rec[i].E + (0.25) * Utemp[i].E + 0.25 * Ite_Step * Iteres[i - 1].F3;
	}
	UpData(UtempN);//用来更新点值
	Utemp = UtempN;
	StorageUF = UtempN;
	GetHalfPoint();//用来得到半点处结果的值并更新ite

	for (int i = 1; i < Mesh_Number - 1; i++) {
		UtempN[i].Rou = (1.0 / 3) * rec[i].Rou + (2.0 / 3) * Utemp[i].Rou + (2.0 / 3) * Ite_Step * Iteres[i - 1].F1;
		UtempN[i].RouU = (1.0 / 3) * rec[i].RouU + (2.0 / 3) * Utemp[i].RouU + (2.0 / 3) * Ite_Step * Iteres[i - 1].F2;
		UtempN[i].E = (1.0 / 3) * rec[i].E + (2.0 / 3) * Utemp[i].E + (2.0 / 3) * Ite_Step * Iteres[i - 1].F3;
	}
	UpData(UtempN);
	/*updata*/
	StorageUF = UtempN;
	return;
};

vector<vector<double>> Solution::Matrix_Multiplication(const vector<vector<double>>& MatrixA, const vector<vector<double>>& MatrixB) {
	int rowA = static_cast<int> (MatrixA.size());
	int colA = static_cast<int> (MatrixA[0].size());
	int rowB = static_cast<int> (MatrixB.size());
	int colB = static_cast<int> (MatrixB[0].size());
	if (colA != rowB) {
		cout << "Those two Matrix are not match" << endl;
		abort();
	}
	vector<vector<double>> res(rowA, vector<double>(colB, 0));
	for (int i = 0; i < rowA; i++) {
		for (int j = 0; j < colB; j++) {
			for (int k = 0; k < colA; k++) {
				res[i][j] += MatrixA[i][k] * MatrixB[k][j];
			}
		}
	}
	return res;
};

void Solution::Calcu() {
	int count = 0;
	ExactSol(testcase);//initiate 
	while (count++ < Ite_Tim) {
		GetHalfPoint();
		RungeKutta();
		cout << "Iteration step now are " << count << " remain " << Ite_Tim - count << endl;
	}
	out_put(outpath);
	return;
}

void Solution::out_put(const string& rote) {
	ofstream file(rote);
	string line;
	file << "The result as follows" << endl;
	file << "position     rou     u     p" << endl;
	int count = 0;
	for (auto& elem : this->StorageUF) {
		file << setprecision(3) << Beg_pos + this->deltaX * count << "  " << elem.Rou << "  " << elem.U << "  " << elem.P << endl;
		count++;
	}
	file.close();
	return;
};

Solution::Solution(const int& Number_Mesh, const double& Beg_pos, const double& End_pos, const int& IterationTimes, const double& step, const string& rote, const bool& Ex, const int& Tcase,const string& out_rote) {
	/*Basic para*/
	this->deltaX = (End_pos - Beg_pos) / Number_Mesh;
	this->Mesh_Number = Number_Mesh;
	this->Beg_pos = Beg_pos;
	this->End_pos = End_pos;
	this->Ite_Step = step;
	this->Ite_Tim = IterationTimes;
	this->Exc = Ex;
	this->outpath = out_rote;
	this->testcase = Tcase;
	/*Control the point number in storage*/
	/*Because of the Calculation requirement, the point in j+1/2 need point in j-1,j,j+1,j+2*/
	/*Thus, the point in j,have to be larger than in j+1/2*/
	PIODT = Mesh_Number;
	PIEM = Mesh_Number + 1;

	/*Init Ex matrix and storage structure*/
	Storage.resize(PIODT); Inv_S.resize(PIODT); _S.resize(PIODT); ModLambda.resize(PIODT);
	StorageUF.resize(PIEM);
	/*Mesh_Number+1 means the froentest and backest point are fake*/
	Iteres.resize(Mesh_Number - 1);
	/*Init_Val*/
};

string GetCurrentPath() {
	char curr_path[1024];
	char* nul;
	nul = _getcwd(curr_path, 1024);
	string res = curr_path;
	return res;
}

