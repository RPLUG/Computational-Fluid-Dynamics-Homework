#include "para.h"
#include <vector>
#include <string>
void para_read(string& pararote, int* Mesh_Number, double* DeltaT, int* Iterationtimes, vector<double>& A6) {
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
	*Mesh_Number = ss[0];
	*DeltaT = ss[1];
	*Iterationtimes = ss[2];
	A6.push_back(ss[3]);
	A6.push_back(ss[4]);
	A6.push_back(ss[5]);
	file.close();
	return;
}