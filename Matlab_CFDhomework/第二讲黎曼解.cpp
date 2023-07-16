# include<iostream>
# include<cmath>
# include<string>
# include<fstream>
# define ga 1.4
typedef double (*pf)(double, double, double, double, double);
double f(double p, double pi, double rhoi)
{
	double temp, ci;
	ci = sqrt(ga * pi / rhoi);
	if (p > pi)
		temp = (p - pi) / ((rhoi * ci) * sqrt((ga + 1) / (2 * ga) * (p / pi) + (ga - 1) / (2 * ga)));
	else
		temp = (2 * ci) / (ga - 1) * (pow((p / pi), (ga - 1) / (2 * ga)) - 1);
	return temp;
}
double F(double p, double p1, double rho1, double p2, double rho2)
{
	double temp;
	temp = f(p, p1, rho1) + f(p, p2, rho2);
	return temp;
}
double dF(pf F, double p, double p1,double rho1,double p2, double rho2)//����
{
	double temp,delta;
	delta = 0.00000001;
	temp = (F(p + delta,p1,rho1,p2,rho2) - F(p - delta, p1, rho1,p2,rho2)) / (2 * delta);
	return temp;
}
int main()
{
	std::ofstream OutFile("Test.txt");//���ù��캯������txt�ı������Ҵ򿪸��ı�����¼����
	double p1, p2, rho1, rho2, u1, u2,pk, p,x,t,u,rho;
	double u_star,u_head,u_tail,c1,c_starl,p_star,p_medium,rho_starl,rho_starr,u_medium,c_medium,rho_medium,x_head,x_tail,Z2,x_interval;
	t = 0.14;
	pk=2, p_star = 0.5;//��ѭ��һ������������������ֵ0.5
	u1 = 0, rho1 = 1, p1 = 1, u2 = 0, rho2 = 0.125, p2 = 0.1;
	c1 = sqrt(ga * p1 / rho1);
	while (abs(pk - p_star) > 0.000001)//���ѭ�����p*��ֵ
		{
		pk = p_star;
		p_star = pk - F(pk, p1, rho1, p2, rho2) / dF(F, pk, p1, rho1, p2, rho2);
		}
	rho_starl = rho1/pow(p1 / p_star , 1 / ga);//�ܶȼ��������ܶ�
	rho_starr = rho2 * ((ga - 1) * p2 + (ga + 1) * p_star) /( (ga + 1) * p2 + (ga - 1) * p_star);//�ܶȼ�����Ҳ��ܶ�
	c_starl = sqrt(ga * p_star / rho_starl);//�ܶȼ������ನ��
	u_star = u1 - f(p_star, p1, rho1);//���Ͳ��뼤��������
	u_head = u1 - c1;//�󲨲�ͷ�ٶ�
	u_tail = u_star - c_starl;//�󲨲�β�ٶ�
	x_head = u_head * t;//���Ͳ���ͷλ��
	x_tail = u_tail * t;//���Ͳ���βλ��
	Z2 = (rho1 * u2 - rho_starr * u_star) / (rho2 - rho_starr);//�Ҽ��������ٶ�
	x_interval = u_star * t;//�ܶȼ����λ��
	for (x=-1;x<=1;x=x+0.01)//�ٶȷֲ�
	{
		c_medium = (ga - 1) / (ga + 1) * (u1 - x / t) + 2 / (ga + 1) * c1;//���Ͳ��ڲ����ٷֲ�
		u_medium = c_medium + x / t;//���Ͳ��ڲ����ٶȷֲ�
		if (x < x_head)
			u = u1;
		else if (x_head <= x && x < x_tail)
			u = u_medium;
		else if (x < Z2 * t && x >= x_tail)
			u = u_star;
		else
			u = u2;
		std::cout<< "x=" << x << "     u=" << u << std::endl;
		OutFile <<  x <<" "<< u << std::endl;
		
	}
	for (x = -1; x <= 1; x = x + 0.01)//�ܶȷֲ�
	{
		c_medium = (ga - 1) / (ga + 1) * (u1 - x / t) + 2 / (ga + 1) * c1;
		p_medium = pow(p1*(c_medium/c1),2*ga/(ga-1) );//���Ͳ��ڲ�ѹǿ�ֲ�
		rho_medium = ga * p_medium / pow(c_medium, 2);//���Ͳ��ڲ��ܶȷֲ�
		if (x < x_head)
			rho = rho1;
		else if (x_head <= x && x < x_tail)
			rho = rho_medium;
		else if (x < x_interval * t && x >= x_tail)
			rho = rho_starl;
		else if (x <= x_interval * t && x <Z2*t )
			rho = rho_starr;
		else
			rho = rho2;
		std::cout << "x=" << x << "     ��=" << u << std::endl;
		OutFile << x << " " << rho << std::endl;
		
	}
	for (x = -1; x <= 1; x = x + 0.01)//ѹǿ�ֲ�
	{
		c_medium = (ga - 1) / (ga + 1) * (u1 - x / t) + 2 / (ga + 1) * c1;
		p_medium = pow(p1 * (c_medium / c1), 2 * ga / (ga - 1));//���Ͳ��ڲ�ѹǿ�ֲ�
		if (x < x_head)
			p = p1;
		else if (x_head <= x && x < x_tail)
			p = p_medium;
		else if (x < Z2 * t && x >= x_tail)
			p = p_star;
		else
			p = p2;
		std::cout << "x=" << x << "     p=" << u << std::endl;
		OutFile  << x << " " << p << std::endl;
	}
	return 0;
}