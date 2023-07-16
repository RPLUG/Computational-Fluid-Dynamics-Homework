#include "Calcu.h"
pair<vector<KData>, double> Draw(const vector<double>& A6, double& deltax, const string& outpath, bool& build) {
    KData* store = (KData*)malloc(sizeof(KData) * 1);
    if (store == nullptr) {
        cout << "malloc store crash" << endl;
        abort();
    }
    vector<KData> res;
    
    double a6res = 0;
    int times = 0;
    double para1 = 0x3f3f3f3f;
    double para2 = 0x3f3f3f3f;
    double para3 = 1;
    bool flag = true;
    for (double a6 = A6[0]; a6 < A6[1]; a6 += A6[2]) {
        //cout << "try times " << times << endl;
        flag = false;
        vector<KData> temp;
        double terminal = 4.0 / deltax;
        for (double k = 0; k < terminal; k += 0.2) {
            store->alpha = k * deltax;
            store->Ki = ExperationForKi(store->alpha, deltax, a6);
            store->Kr = ExperationForKr(store->alpha, deltax, a6);
            double judg = fabs(1 - (store->Ki / store->alpha));
            if (store->Kr < 0.02 && judg < 0.02) {
                if (store->alpha>para3){
                    flag = true;
                    para2 = judg;
                    build = true;
                }
            }
            temp.emplace_back(*store);      
        }
        if (flag) {
            res = temp;
            a6res = a6;
        }
        times++;
    }
    return {res,a6res};
};

double ExperationForKr(const double& alpha,const double& deltax, double& a6){
    double re = (-1.0 / 12 - (a6 * deltax)) * cos(-3.0 * alpha) + (5 * a6 * deltax + (1.0 / 2)) * cos(-2.0 * alpha) + (-3.0 / 2 - (a6 * deltax) * 10) * cos(-alpha) + (5.0 / 6 + 10.0 * (a6 * deltax)) + (-5.0 * a6 * deltax + 1.0 / 4) * cos(alpha) + (a6 * deltax) * cos(2 * alpha);
    //double re = -a6 * deltax / 10 * cos(3 * alpha) - (3.0 / 2) * a6 * deltax * cos(alpha) + 3.0 / 5 * a6 * deltax * cos(2 * alpha) + a6 * deltax;
    return re;
};

double ExperationForKi(const double& alpha,const double& deltax, double& a6){
    double re = (-1.0 / 12 + (a6 * deltax)) * sin(-3.0 * alpha) + (5 * a6 * deltax + (1.0 / 2)) * sin(-2.0 * alpha) + (-3.0 / 2 - (a6 * deltax) * 10) * sin(-alpha) + (-5.0 * a6 * deltax + 1.0 / 4) * sin(alpha) + (a6 * deltax) * sin(2 * alpha);
   // double re = a6 * deltax / 10 * sin(3 * alpha) + ((1.0 / 2) * a6 * deltax + 4.0 / 3) * sin(alpha) - (2.0 / 5 * a6 * deltax + 1.0 / 6) * sin(2 * alpha);
    return re;
};

void result_out_put(const string& outpath, pair<vector<KData>, double> st,const double& deltax){
    ofstream file(outpath);
    //cout << outpath;
    file<<"deltax is " << deltax<<endl;
    file << "When a6 = " << st.second << endl;
    file << "alpha        " << "Kr        " << "Ki        " << endl;
    int count=1;
    for(auto & Felem :st.first){
        file<<setprecision(4)<<Felem.alpha<<"    "<<Felem.Kr <<"    "<<Felem.Ki<<endl;
    }
    file.close();
    
    return;
}

Solution::Solution(pair<vector<KData>, double> para, double deltaX, double deltaT, int iterationtimes, int Mesh_Number) 
{
    this->deltaT = deltaT;
    this->deltaX = deltaX;
    this->build = false;
    this->iterationtimes = iterationtimes;
    this->st = para;
    this->Mesh_Number = Mesh_Number;
};

void Solution::loop() {
    int itm = this->iterationtimes;
    int n = this->Mesh_Number;
    int count = 0;
    while (count <= itm) {
        this->middle = this->storage;
        this->pre = this->storage;
        RungeKutta();
        count++;
    }
    return;
}

double Solution::LU(int& j) {
    int n = this->Mesh_Number;
    int j1 = j - 3 < 0 ? (n - 1 + (j + 1) - 3) : j - 3;
    int j2 = j - 2 < 0 ? (n - 1 + (j + 1) - 2) : j - 2;
    int j3 = j - 1 < 0 ? (n - 1 + (j + 1) - 1) : j - 1;
    int j4 = j + 1 > n - 1 ? ((j) - (n - 1)) : j + 1;
    int j5 = j + 2 > n - 1 ? ((j+1) - (n - 1)) : j + 2;
    return -(this->para[0] * this->pre[j1] + this->para[1] * this->pre[j2]
        + this->para[2] * this->pre[j3] + this->para[3] * this->pre[j]
        + this->para[4] * this->pre[j4] + this->para[5] * this->pre[j5]);
}

void Solution::RungeKutta() {
    //TVD Runge-Kutta
    int n = this->Mesh_Number;
    for (int i = 0; i < n; i++) {
        double u1 = this->storage[i] + this->deltaT * LU(i);
        this->middle[i] = u1;
    }
    this->pre = this->middle;
    for (int i = 0; i < n; i++) {
        double u2 = (3.0 / 4) * this->storage[i] + (1.0 / 4) * (this->pre[i] + deltaT * LU(i));
        this->middle[i] = u2;
    }
    this->pre = this->middle;
    for (int i = 0; i < n; i++) {
        double res = (1.0 / 3) * this->storage[i] + (2.0 / 3) * (this->pre[i] + deltaT * LU(i));
        this->middle[i] = res;
    }
    this->storage = this->middle;
    return ;
}

void Solution::ExtractSolINIT() {
    double deltax = this->deltaX;
    for (int i = 0; i < this->Mesh_Number; i++) {
        this->storage.emplace_back(sin(i * deltax));
    }
    return ;
}

void Solution::Getpara() {
    double k = this->st.second;
    double deltax = this->deltaX;
    this->para.emplace_back((( - 1.0 / (12 * deltax)) - k));;
    this->para.emplace_back(((1.0 / (2 * deltax)) + 5 * k));
    this->para.emplace_back(((-3.0 / (2 * deltax)) - (10 * k)));
    this->para.emplace_back(((5.0/(6*deltax))+(10*k)));
    this->para.emplace_back(((-5*k)+(1.0/(4*deltax))));
    this->para.emplace_back(k);
    for (auto & elem : para) {
        cout << elem << endl;
    }
}



void Solution::out_put(string& wavepath) {
    ofstream file(wavepath);
    file << "The following are the wave data" << endl;
    file << "Position    Value" << endl;
    int count = 0;
    for (auto& elem : this->storage) {
        file << setprecision(3)<<count*this->deltaX <<"    " << elem << endl;
        count++;
    }
    file.close();
    return;
}


string GetCurrentPath() {
	char curr_path[1024];
	_getcwd(curr_path, 1024);
	string res = curr_path;
	return res;
}

