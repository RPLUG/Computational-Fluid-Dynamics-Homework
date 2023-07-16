#include "para.h"
#include <vector>
#include <string>
void para_read(string& pararote, double* p1, double* p2, double* rou1, double* rou2, double* u1, double* u2, double* elip, double* time) {
	ifstream file(pararote);
	vector<double> ss;
	string rote;
	while (getline(file, rote)) {
		int pos = rote.find("=") + 1;
		while (rote[pos] == ' ') {
			pos++;
		}
		string data = rote.substr(pos);
		ss.emplace_back(atof(data.c_str()));
	}
	*p1 = ss[0];
	*p2 = ss[1];
	*rou1 = ss[2];
	*rou2 = ss[3];
	*u1 = ss[4];
	*u2 = ss[5];
	*elip = ss[6];
	*time = ss[7];
	file.close();
	return;
}