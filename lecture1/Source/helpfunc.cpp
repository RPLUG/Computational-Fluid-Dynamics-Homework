#include "helpfunc.h"
#include <vector>
using namespace std;
NewTon::NewTon(const double &time, const double &p1, const double &p2, const double &rou1, const double &rou2, const double &u1, const double &u2, const double &elip, const string &outrote)
{
	this->p1 = p1;
	this->p2 = p2;
	this->rou1 = rou1;
	this->rou2 = rou2;
	this->elip = elip;
	this->u1 = u1;
	this->u2 = u2;
	this->c1 = sqrt(this->gamma * p1 / rou1);
	this->c2 = sqrt(this->gamma * p2 / rou2);
	this->Init_Time = time;
	this->outroute = outrote;
}

double NewTon::NewTonCalcu()
{
	int count = 0;
	double p_ = (p1 * rou2 * c2 + p2 * rou1 * c1 + (u2 - u1) * (rou1 * c1 * rou2 * c2)) / (rou1 * c1 + rou2 * c2); // init data
	pair<double, double> expres1 = GetExpre(p_, p1, rou1, u1, c1, leftex);
	pair<double, double> expres2 = GetExpre(p_, p2, rou2, u2, c2, rightex);
	double newp_ = p_ - (expres1.first + expres2.first - u1 + u2) / (expres1.second + expres2.second);
	while ((newp_ - p_) > this->elip && count <= this->divnumb)
	{
		p_ = newp_;
		expres1 = GetExpre(p_, p1, rou1, u1, c1, leftex);
		expres2 = GetExpre(p_, p2, rou2, u2, c2, rightex);
		newp_ = p_ - (expres1.first + expres2.first - u1 + u2) / (expres2.second + expres2.second);
		++count;
		cout << "iterator times " << count << "  the value of p* in regin(3) is " << newp_ << endl;
	}

	if (count == divnumb && newp_ - p_ > elip)
	{
		cout << "iterator false" << endl;
		succ = 0;
	}
	else
	{
		cout << "Total iterator times: " << count << endl;
		succ = 1;
	}
	return newp_;
}

pair<double, double> NewTon::GetExpre(const double &p_, const double &pi, const double &roui, const double &ui, const double &ci, bool &status)
{
	double fx;
	double fxx;
	if (p_ >= pi)
	{
		pair<double, double> res = GetEqua(p_, pi, roui, ui, ci);
		fx = res.first;
		fxx = res.second;
		status = true;
	}
	else if (p_ < pi)
	{
		pair<double, double> res = GetDEqua(p_, pi, roui, ui, ci);
		fx = res.first;
		fxx = res.second;
		status = false;
	}
	return make_pair(fx, fxx);
}

double NewTon::GetA(const double &p_, const double &pi, const double &roui, const double &ci)
{
	double Ai = (roui * ci * sqrt(((p_ / pi) * (gamma + 1) / (2 * gamma)) + ((gamma - 1) / (2 * gamma))));
	return Ai;
}

pair<double, double> NewTon::GetEqua(const double &p_, const double &pi, const double &roui, const double &ui, const double &ci)
{
	// p_ >= pi
	double fx = (p_ - pi) / GetA(p_, pi, roui, ci);
	double fxx = ((gamma + 1) * (p_ / pi) + (3 * gamma - 1)) / ((4 * gamma * roui * ci) * pow(((((gamma + 1) / (2 * gamma)) * (p_ / pi)) + ((gamma - 1) / (2 * gamma))), 1.5));
	return make_pair(fx, fxx);
}

pair<double, double> NewTon::GetDEqua(const double &p_, const double &pi, const double &roui, const double &ui, const double &ci)
{
	// p_ < pi
	double fx = (2 * ci / (gamma - 1)) * (pow((p_ / pi), (gamma - 1) / (2 * gamma)) - 1);
	double fxx = ((ci) / (gamma * pi)) * pow(p_ / pi, (-gamma - 1) / (2 * gamma));
	return make_pair(fx, fxx);
}

void NewTon::GetProperity(const double &p)
{
	// region in (3) or (4)
	// storage contaion u*,rou1,rou2,
	p_ = p;
	double A1, A2;

	// in region (3) and (4),leftu=rightu, in there just to distinguish
	if (leftex && !rightex)
	{
		leftu = u1 - GetEqua(p, p1, rou1, u1, c1).first;
		rightu = u2 + GetDEqua(p, p2, rou2, u2, c2).first;

		/*Calculate the left properity*/
		A1 = GetA(p, p1, rou1, c1);
		Z1 = u1 - A1 / rou1;
		roustar1 = (rou1 * A1) / (A1 - rou1 * (u1 - leftu));
		Pos_Left_Shok_End = Z1 * Init_Time;

		/*Calculate the right properity*/
		c2star = c2 + (gamma - 1) * (u2 - rightu) / 2;
		roustar2 = gamma * p / pow(c2star, 2);
		Z2head = u2 + c2;
		Z2tail = rightu + c2star;

		/*Calculate the properity inside the expension wave*/
		ExpWavePro(Z2head, Z2tail);
	}
	else if (!leftex && rightex)
	{
		leftu = u1 - GetDEqua(p, p1, rou1, u1, c1).first;
		rightu = u2 + GetEqua(p, p2, rou2, u2, c2).first;

		/*Calculate the right properity*/
		A2 = GetA(p, p2, rou2, c2);
		roustar2 = (rou2 * A2) / (A2 + rou2 * (u2 - rightu));
		Z2 = u2 + A2 / rou2;
		Pos_Right_Shok_End = Z2 * Init_Time;
		/*Calculate the left properity*/
		c1star = c1 + (gamma - 1) * (u1 - leftu) / 2;
		roustar1 = gamma * p / pow(c1star, 2);
		Z1head = u1 - c1;
		Z1tail = leftu - c1star;
		/*Calculate the properity inside the expension wave*/
		ExpWavePro(Z1head, Z1tail);
	}
	else if (leftex && rightex)
	{
		leftu = u1 - GetEqua(p, p1, rou1, u1, c1).first;
		rightu = u2 + GetEqua(p, p2, rou2, u2, c2).first;
		A1 = GetA(p, p1, rou1, c1);
		A2 = GetA(p, p2, rou2, c2);
		roustar1 = (rou1 * A1) / (A1 - rou1 * (u1 - leftu));
		roustar2 = (rou2 * A2) / (A2 + rou2 * (u2 - rightu));
		Z1 = u1 - A1 / rou1;
		Z2 = u2 + A2 / rou2;
		roustar1 = (rou1 * A1) / (A1 - rou1 * (u1 - leftu));
		roustar2 = (rou2 * A2) / (A2 + rou2 * (u2 - rightu));
		Pos_Left_Shok_End = Z1 * Init_Time;
		Pos_Right_Shok_End = Z2 * Init_Time;
	}

	return;
}

void NewTon::ExpWavePro(const double &Zhead, const double &Ztail)
{
	double step = abs((abs(Zhead) * Init_Time - abs(Ztail) * Init_Time)) / Exp_Wave_Data_Density;
	storage.resize(Exp_Wave_Data_Density + 1);
	int count = 0;
	if (!rightex)
	{
		for (double position = Ztail * Init_Time; position >= Zhead * Init_Time; position += step)
		{
			storage[count].Positon = position;
			storage[count].C_in_Exp = ((gamma - 1) / (gamma + 1)) * ((position / Init_Time) - u2) + c2 * 2 / (gamma + 1);
			storage[count].P_in_Exp = p2 * pow(storage[count].C_in_Exp / c2, 2 * gamma / (gamma - 1));
			storage[count].rou_in_Exp = gamma * storage[count].P_in_Exp / pow(storage[count].C_in_Exp, 2);
			storage[count].U_in_Exp = (position / Init_Time) - storage[count].C_in_Exp;
			++count;
		}
	}

	if (!leftex)
	{
		for (double position = Ztail * Init_Time; position >= Zhead * Init_Time; position -= step)
		{
			storage[count].Positon = position;
			storage[count].C_in_Exp = ((gamma - 1) / (gamma + 1)) * (u1 - position / Init_Time) + c1 * 2 / (gamma + 1);
			storage[count].P_in_Exp = p1 * pow(storage[count].C_in_Exp / c1, 2 * gamma / (gamma - 1));
			storage[count].rou_in_Exp = gamma * storage[count].P_in_Exp / pow(storage[count].C_in_Exp, 2);
			storage[count].U_in_Exp = storage[count].C_in_Exp + position / Init_Time;
			++count;
		}
	}
}

void NewTon::Result_out_put()
{
	ofstream file(outroute);

	if (!file)
	{
		cout << "This file can't be opean" << endl;
		abort();
	}
	file.precision(5);
	if (!leftex && rightex)
	{
		file << "The left side is expension wave, and the right side is shock wave" << endl;
		file << "The following are their properity" << endl;
		file << "For shock wave, in region(3), we have Position and speed" << endl;
		file << Z2 * Init_Time << "    " << Z2 << endl;
		file << "For the gap we have rou*L, rou*R, u*, p* " << endl;
		file << roustar1 << "    " << roustar2 << "    " << rightu << "    " << p_ << endl;
		file << "For expension wave, we have wave's head " << Z1head * Init_Time << ", And the wave's tail " << Z1tail * Init_Time << endl;
		file << "Position             speed             rou             p" << endl;
		for (int i = 0; i < Exp_Wave_Data_Density; i++)
		{
			file << setiosflags(ios::left) << storage[i].Positon << "    " << storage[i].U_in_Exp << "    " << storage[i].rou_in_Exp << "   " << storage[i].P_in_Exp << endl;
		}
	}
	else if (leftex && !rightex)
	{
		file << "The right side is expension wave, and the left side is shok wave" << endl;
		file << "The following are their properity" << endl;
		file << "For shock wave, in region(3), we have Position and speed" << endl;
		file << Z1 * Init_Time << "    " << Z1 << endl;
		file << "For the gap we have rou*L, rou*R, u*, p* " << endl;
		file << roustar1 << "    " << roustar2 << "    " << rightu << "    " << p_ << endl;
		file << "For expension wave, we have wave's head " << Z2head * Init_Time << ", And the wave's tail " << Z2tail * Init_Time << endl;
		file << "Position             speed             rou             p" << endl;
		for (int i = 0; i < Exp_Wave_Data_Density; i++)
		{
			file << storage[i].Positon << "    " << storage[i].U_in_Exp << "    " << storage[i].rou_in_Exp << "   " << storage[i].P_in_Exp << endl;
		}
	}
	else
	{
		file << "The right and left side are both shok wave" << endl;
		file << "The following are their properity" << endl;
		file << "For left shock wave, in region(3), we have Position and speed" << endl;
		file << Z1 * Init_Time << "    " << Z1 << endl;
		file << "For right shock wave, in region(3), we have Position and speed" << endl;
		file << Z2 * Init_Time << "    " << Z2 << endl;
		file << "For the gap we have rou*L, rou*R, u*, p* " << endl;
		file << roustar1 << "    " << roustar2 << "    " << rightu << "    " << p_ << endl;
	}
}

NewTon::~NewTon()
{
	if (succ)
	{
		cout << "Calculation success" << endl;
	}
	else
	{
		cout << "Calculation failed" << endl;
	}
};
