#include "GeneF.h"
MeshGenerate::MeshGenerate(para para) :Rote(para.rote), Col(para.col), Row(para.row), BEnum(para.BEnum), Rside(para.Rside),
Lside(para.Lside), FlyBe(para.FlyBe), H1b(para.H1b), H1e(para.H1e), H2b(para.H2b), H2e(para.H2e), Tole(para.Tole),
center(para.center), HenStart(para.HenStart), H3b(para.H3b), H3e(para.H3e), MaxIteration(para.MaxIteration)
{
	Org_Coor.resize(Row, vector<Coordinate>(Col-1));
	Mesh_Init();
	//LUSGS();
	JKB();
}

void MeshGenerate::SloveAB(double& A, double& B,const double& hb,const double& he, const int& n) {
	//二分查根
	double ans = 100;
	double tp1 = 0;//给个大致范围就行
	double tp2 = 100;
	while (abs(0 - ans) > 1e-10) {
		double middle = tp1 + (tp2 - tp1) / 2;
		ans = he * pow(e, (middle / (n - 1))) - hb * pow(e, middle) - he + hb;
		if (ans >0 ) {
			tp1 = middle;
		}
		else {
			tp2 = middle;
		}
	}
	B = tp2;
	A = he / (pow(e, B) - 1);
	return;
}

MeshGenerate::~MeshGenerate()
{
}

void MeshGenerate::Mesh_Init() {
	//Row shoule be even
	//int mid = ((Col - 2 - BEnum) + (BEnum - 1 + 1) - 1) / 2;
	//计算放缩因子
	SloveAB(A1, B1, H1b, H1e, Row);
	cout << "竖直方向比例因子A1 = " << A1 << "B1 = " << B1 << endl;
	SloveAB(A2, B2, H2b, H2e, (Col) / 2);
	cout << "水平方向比例因子A2 = " << A2 << "B2 = " << B2 << endl;
	SloveAB(A3, B3, H3b, H3e, (Col) / 2);
	cout << "水平方向比例因子A3 = " << A3 << "B3 = " << B3 << endl;
	double s = 0;
	double dis = 0;
	double dist = 0;
	//初始化翼型上沿
	for (int j = (Col) / 2 - 1, i = 0, k = 1; j <= Col - 2; j++, k++) {
		s = double(k-1) / (double (Col)/2-1);//上侧分割25层；
		dis = A2 * (pow(e, B2 * s) - 1);
		if(j<=Col-1-BEnum){
			Org_Coor[i][j].x = FlyBe+dis;
			Org_Coor[i][j].y = center+0.1781 * sqrt(dis) - 0.0756 * dis - 0.2122 * pow(dis, 2) + 0.1705 * pow(dis, 3) - 0.0609 * pow(dis, 4);
			Org_Coor[i][j].prex = Org_Coor[i][j].x;
			Org_Coor[i][j].prey = Org_Coor[i][j].y;
			if (Org_Coor[i][j].y < 0) {
				Org_Coor[i][j].y = center;
				Org_Coor[i][j].prey = Org_Coor[i][j].y;
			}
		}
		else if (j >= Col - BEnum) {
			Org_Coor[i][j].x = FlyBe + dis;
			Org_Coor[i][j].y = center;
			Org_Coor[i][j].prex = Org_Coor[i][j].x;
			Org_Coor[i][j].prey = Org_Coor[i][j].y;
		}
		Org_Coor[i][j].boundary = 1;
	}
	
	//初始化翼型下沿
	for (int j = Col/2-2, i = 0, k = 2; j >=0 ; k++, j--) {
		s = double(k - 1) / (double(Col) / 2-1);
		dis = A2 * (pow(e, B2 * s) - 1);
		if(j>=BEnum-1){
			Org_Coor[i][j].x = FlyBe + dis;
			Org_Coor[i][j].y = center + -(0.1781 * sqrt(dis) - 0.0756 * dis - 0.2122 * pow(dis, 2) + 0.1705 * pow(dis, 3) - 0.0609 * pow(dis, 4));
			Org_Coor[i][j].prex = Org_Coor[i][j].x;
			Org_Coor[i][j].prey = Org_Coor[i][j].y;
			if (Org_Coor[i][j].y > 0) {
				Org_Coor[i][j].y = center;
				Org_Coor[i][j].prey = Org_Coor[i][j].y;
			}
		}
		else if (j <BEnum-1) {
			Org_Coor[i][j].x = FlyBe + dis;
			Org_Coor[i][j].y = center;
			Org_Coor[i][j].prex = Org_Coor[i][j].x;
			Org_Coor[i][j].prey = Org_Coor[i][j].y;
		}
		Org_Coor[i][j].boundary = 1;
	}
	//右侧沿
	
	for (int i = 1, j = Col - 2, k = 2; i < Row-1; k++, i++) {
		Org_Coor[i][j].x = Rside;
		s = double(k - 1) / (Row - 2);
		dist = A1 * (pow(e, B1* s) - 1);
		Org_Coor[i][j].y = center + dist;
		Org_Coor[i][j].boundary = 1;
		Org_Coor[i][j].prex = Org_Coor[i][j].x;
		Org_Coor[i][j].prey = Org_Coor[i][j].y;
	}
	
	//左侧沿
	dist = 0;
	for (int i = 1, j = 0, k = 2; i < Row-1; k++, i++) {
		Org_Coor[i][j].x = Rside;
		s = double(k - 1) / (Row - 2);
		dist = A1 * (pow(e, B1 * s) - 1);
		Org_Coor[i][j].y = center + -dist;
		Org_Coor[i][j].boundary = 1;
		Org_Coor[i][j].prex = Org_Coor[i][j].x;
		Org_Coor[i][j].prey = Org_Coor[i][j].y;
	}
	
	//外测沿上半部分
	for (int i = Row - 1, j = (Col) / 2 - 1, k = 1; j <= Col - 2; j++, k++) {
		s = double(k - 1) / (double(Col) / 2 - 1);
		dist = A3 * (pow(e, B3 * s) - 1);
		Org_Coor[i][j].x = Lside+dist;
		Org_Coor[i][j].y = center + sqrt(1 - (Org_Coor[i][j].x - HenStart) * (Org_Coor[i][j].x - HenStart));
		Org_Coor[i][j].prex = Org_Coor[i][j].x;
		Org_Coor[i][j].prey = Org_Coor[i][j].y;
		if (Org_Coor[i][j].x>=HenStart) {
			Org_Coor[i][j].y = 1;
			Org_Coor[i][j].prey = Org_Coor[i][j].y;
		}
		Org_Coor[i][j].boundary = 1;
	}
	//外测沿下半部分
	for (int i = Row - 1, j = Col / 2 - 2, k = 2; j >= 0; j--,k++) {
		s = double(k - 1) / (double(Col) / 2 - 1);
		dist = A3 * (pow(e, B3 * s) - 1);
		Org_Coor[i][j].x = Lside + dist;
		Org_Coor[i][j].y = center - sqrt(1 - (Org_Coor[i][j].x - HenStart) * (Org_Coor[i][j].x - HenStart));
		Org_Coor[i][j].prex = Org_Coor[i][j].x;
		Org_Coor[i][j].prey = Org_Coor[i][j].y;
		if (Org_Coor[i][j].x >=HenStart) {
			Org_Coor[i][j].y = -1;
			Org_Coor[i][j].prey = Org_Coor[i][j].y;
		}
		Org_Coor[i][j].boundary = 1;
	}
	Print_res();
	return; 
}

void MeshGenerate::LUSGS() {
	if (Org_Coor.size() == 0) {
		cout << "Matrix Size Wrong" << endl;
		abort();
	}
	int nr = Org_Coor.size();//row;
	int nc = Org_Coor[0].size();
	double Tolerence = 100;
	int time = 0;
	Tolerence = 0;
	//上扫过程
	int row = 0;
	int col = 0;
	int dir = 1;
	while (row < nr && col < nc) {
		int nerow = row + (dir == 1 ? -1 : 1);
		int necol = col + (dir == 1 ? 1 : -1);
		if (nerow < 0 || necol < 0 || nerow == nr || necol == nc) {
			if (dir == 1) {
				row += (col == nc - 1 ? 1 : 0);
				col += (col == nc - 1 ? 0 : 1);
			}
			else if (dir == 0) {
				col += (row == nr - 1 ? 1 : 0);
				row += (row == nr - 1 ? 0 : 1);
			}
			dir = 1 - dir;
		}
		else {
			row = nerow;
			col = necol;

		}
		//再次判断是否到达边界，因为这个while循环在外层控制，有可能内层已经达到边界但不能被while检测
		if (row == nr || col == nc) {
			break;
		}
		if (Org_Coor[row][col].boundary) {
			continue;
		}
		else {
			LUSGSCalcu(row, col, 'U');
			//cout << row << "  " << col << endl;;
		}
	}
	//下扫过程
	row = nr - 1;
	col = nc - 1;
	dir = 1;
	while (row >= 0 && col >= 0) {
		int nerow = row + (dir == 1 ? -1 : 1);
		int necol = col + (dir == 1 ? 1 : -1);
		if (nerow < 0 || necol < 0 || nerow == nr || necol == nc) {
			if (dir == 1) {
				col += (row == 0 ? -1 : 0);
				row += (row == 0 ? 0 : -1);
			}
			else if (dir == 0) {
				row += (col == 0 ? -1 : 0);
				col += (col == 0 ? 0 : -1);
			}
			dir = 1 - dir;
		}
		else {
			row = nerow;
			col = necol;
		}
		if (row == -1 || col == -1) {
			break;
		}
		if (Org_Coor[row][col].boundary) {
			continue;
		}
		else {
			//cout << row << "  " << col << endl;;
			Tolerence += LUSGSCalcu(row, col, 'D');
		}
	}

	Print_res();
}

void MeshGenerate::JKB() {
	int row = Org_Coor.size();
	int col = Org_Coor[0].size();
	for (int p = 0; p < MaxIteration; p++) {
		for (int i = 1; i < row-1; i++) {
			for (int j = 1; j < col-1; j++) {
				JKBCalcu(i, j);
			}
		}	
		//Synchronize
		for (auto& Block : Org_Coor) {
			for (auto& elem : Block) {
				elem.prex = elem.x;
				elem.prey = elem.y;
			}
		}
		if (p % 50 == 0) {
			cout << "Iteration Step Now At " << p << endl;
		}
	}
	Print_res();
}

double MeshGenerate::LUSGSCalcu(int& i,int& j, const char& direction) {
	if (i == 47 && j == 191) {
		cout << i << endl;
	}
		auto& Block = Org_Coor[i][j];
		double alpha = pow((Org_Coor[i][j + 1].x - Org_Coor[i][j - 1].x) / 2, 2) + pow((Org_Coor[i][j + 1].y - Org_Coor[i][j - 1].y) / 2, 2);
		double beta =  -(((Org_Coor[i + 1][j].x - Org_Coor[i - 1][j].x) / 2) * ((Org_Coor[i][j + 1].x - Org_Coor[i][j - 1].x) / 2) + 
					((Org_Coor[i + 1][j].y - Org_Coor[i - 1][j].y) / 2) * ((Org_Coor[i][j + 1].y - Org_Coor[i][j - 1].y) / 2));

		double gamma = pow((Org_Coor[i + 1][j].x - Org_Coor[i - 1][j].x) / 2, 2) + pow((Org_Coor[i + 1][j].y - Org_Coor[i - 1][j].y) / 2, 2);
		double Dij = -2 * (alpha + gamma);// - 8 * beta;
		double Dijp1 = gamma;//+ 2 * beta;
		double Dijm1 = Dijp1;
		double Dip1j = alpha; //+ 2 * beta;
		double Dim1j = Dip1j;
		double tolx = 0;
		double toly = 0;
		if (direction == 'U') {
			//计算x坐标
			Block.x = -((Dim1j * Org_Coor[i - 1][j].x) + (Dijm1 * Org_Coor[i][j - 1].x)) / (Dij);
			//计算y坐标
			Block.y = -((Dim1j * Org_Coor[i - 1][j].y) + (Dijm1 * Org_Coor[i][j - 1].y)) / (Dij);
		}else if (direction == 'D') {
			tolx = Block.x;
			toly = Block.y;
			Block.x = (Dij * Block.x - Dip1j * Org_Coor[i + 1][j].x - Dijp1 * Org_Coor[i][j + 1].x) / Dij;
			Block.y = (Dij * Block.y - Dip1j * Org_Coor[i + 1][j].y - Dijp1 * Org_Coor[i][j + 1].y) / Dij;
			tolx = Block.x - tolx;
			toly = Block.y - toly;
		}		
		if(isnan(Block.x)|| isnan(Block.y)){
		cout << beta << endl;
		}
		return tolx * tolx + toly * toly;
}

void MeshGenerate::JKBCalcu(int& i, int& j) {
	double alpha = pow((Org_Coor[i][j + 1].x - Org_Coor[i][j - 1].x) / 2, 2) + pow((Org_Coor[i][j + 1].y - Org_Coor[i][j - 1].y) / 2, 2);
	double beta = 1;//-(((Org_Coor[i + 1][j].x - Org_Coor[i - 1][j].x) / 2) * ((Org_Coor[i][j + 1].x - Org_Coor[i][j - 1].x) / 2) +
		//((Org_Coor[i + 1][j].y - Org_Coor[i - 1][j].y) / 2) * ((Org_Coor[i][j + 1].y - Org_Coor[i][j - 1].y) / 2));

	double gamma = pow((Org_Coor[i + 1][j].x - Org_Coor[i - 1][j].x) / 2, 2) + pow((Org_Coor[i + 1][j].y - Org_Coor[i - 1][j].y) / 2, 2);
	double Dij = -2 * (alpha + gamma) - 8 * beta;
	double Dijp1 = gamma + 2 * beta;
	double Dijm1 = Dijp1;
	double Dip1j = alpha + 2 * beta;
	double Dim1j = Dip1j;
	Org_Coor[i][j].x = -((Dim1j * Org_Coor[i - 1][j].prex) + (Dip1j * Org_Coor[i + 1][j].prex) + (Dijp1 * Org_Coor[i][j + 1].prex) + (Dijm1 * Org_Coor[i][j - 1].prex)) / (Dij);
	Org_Coor[i][j].y = -((Dim1j * Org_Coor[i - 1][j].prey) + (Dip1j * Org_Coor[i + 1][j].prey) + (Dijp1 * Org_Coor[i][j + 1].prey) + (Dijm1 * Org_Coor[i][j - 1].prey)) / (Dij);
}

void MeshGenerate::Print_res() {
	ofstream file(Rote + "test.txt");
	for (int i = 0; i < Row; i++) {
		for (int j = 0; j < Col - 1; j++) {
			file << Org_Coor[i][j].x << "  " << Org_Coor[i][j].y << endl;
		}
	}
	file.close();
}

string GetCurrentPath() {
	char curr_path[1024];
	char* nul;
	nul = _getcwd(curr_path, 1024);
	string res = curr_path;
	return res;
}