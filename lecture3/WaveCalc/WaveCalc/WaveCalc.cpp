#include <iostream>
#include "helpfunc.h"
#include <vector>
#include <string>
#include "para.h"
#include <conio.h>
#include <fstream>
using namespace std;

int main(int argc, char** argv)
{
	string path = GetCurrentPath();
	string inpath = path + "\\para.ini";
	cout << "Launching dir: " << inpath << endl;
	string out = path + "\\result.txt";
	double p1, rou1, u1;
	double p2, rou2, u2;
	double elip;
	double time;
	para_read(inpath, &p1, &p2, &rou1, &rou2, &u1, &u2, &elip, &time);
	cout << "Initial conditon : " << endl;
	cout << "Pleft= " << p1 << ", " << "Rouleft= " << rou1 << ", " << "Uleft= " << u1 << endl;
	cout << "PRight= " << p2 << ", " << "Rouright= " << rou2 << ", " << "Uright= " << u2 << endl;
	cout << "Calculation begin" << endl;
	NewTon* test = new NewTon(time, p1, p2, rou1, rou2, u1, u2, elip, out);
	double newp_ = test->NewTonCalcu();
	test->GetProperity(newp_);
	test->Result_out_put();
	test->~NewTon();
	cout<<"The result are put out in current dir, name 'result.txt'"<<endl;
	cout << "Press anything to quit" << endl;
	_getch();
	return 0;
}


