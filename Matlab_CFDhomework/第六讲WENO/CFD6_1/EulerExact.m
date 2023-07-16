%��[-0.5��0.5]��������⣬���ͳһ��������任��ƽ��
ga=1.4;
t = 0.14;
u1 = 0; rho1 = 1; p1 = 1; u2 = 0; rho2 = 0.125; p2 = 0.1;
c1 = sqrt(ga * p1 / rho1);
F=@(p_star)(f1(p_star, p1, rho1) + f1(p_star, p2, rho2));
   p_star=fzero(F,0.3);
	rho_starl = rho1/(p1 / p_star)^  (1 / ga);%�ܶȼ��������ܶ�
	rho_starr = rho2 * ((ga - 1) * p2 + (ga + 1) * p_star) /( (ga + 1) * p2 + (ga - 1) * p_star);%�ܶȼ�����Ҳ��ܶ�
	c_starl = sqrt(ga * p_star / rho_starl);%�ܶȼ������ನ��
	u_star = u1 - f1(p_star, p1, rho1);%���Ͳ��뼤��������
	u_head = u1 - c1;%�󲨲�ͷ�ٶ�
	u_tail = u_star - c_starl;%�󲨲�β�ٶ�
	x_head = u_head * t;%���Ͳ���ͷλ��
	x_tail = u_tail * t;%���Ͳ���βλ��
	Z2 = (rho1 * u2 - rho_starr * u_star) / (rho2 - rho_starr);%�Ҽ��������ٶ�
	x_interval = u_star * t;%�ܶȼ����λ��
    for i=1:101
        x(i)=-0.51+0.01*i;%�ٶȷֲ�
		c_medium = (ga - 1) / (ga + 1) * (u1 - x(i) / t) + 2 / (ga + 1) * c1;%���Ͳ��ڲ����ٷֲ�
		u_medium = c_medium + x(i) / t;%���Ͳ��ڲ����ٶȷֲ�
		if (x(i) < x_head)
			u(i) = u1; 
        elseif (x_head <= x(i) && x(i) < x_tail)
			u(i) = u_medium;   
        elseif (x(i) < Z2 * t && x(i) >= x_tail)
			u(i) = u_star;   
		else
			u(i) = u2;  
        end
    end
    i=0;
	  for i=1:101
        x(i)=-0.51+2/200*i;%�ܶȷֲ�
		c_medium = (ga - 1) / (ga + 1) * (u1 - x(i) / t) + 2 / (ga + 1) * c1;%���Ͳ��ڲ����ٷֲ�
		u_medium = c_medium + x(i) / t;%���Ͳ��ڲ����ٶȷֲ�
		p_medium = p1*(c_medium/c1)^(2*ga/(ga-1)) ;%���Ͳ��ڲ�ѹǿ�ֲ�
		rho_medium = ga * p_medium / c_medium^2;%���Ͳ��ڲ��ܶȷֲ�
		if (x(i) < x_head)
			rho(i) = rho1;
        elseif (x_head <= x(i) && x(i) < x_tail)
			rho(i) = rho_medium;
        elseif (x(i) < x_interval  && x(i) >= x_tail)
			rho(i) = rho_starl;
        elseif (x(i) >= x_interval  && x(i) <Z2*t )
			rho(i) = rho_starr;
		else
			rho(i) = rho2;
        end
    end
    i=0;	
    	  for i=1:101
        x(i)=-0.51+2/200*i;%�ٶȷֲ�
		c_medium = (ga - 1) / (ga + 1) * (u1 - x(i) / t) + 2 / (ga + 1) * c1;%���Ͳ��ڲ����ٷֲ�
		u_medium = c_medium + x(i) / t;%���Ͳ��ڲ����ٶȷֲ�
		p_medium = p1*(c_medium/c1)^(2*ga/(ga-1)) ;%���Ͳ��ڲ�ѹǿ�ֲ�
		rho_medium = ga * p_medium / c_medium^2;%���Ͳ��ڲ��ܶȷֲ�
	if (x(i) < x_head)
			p(i) = p1;
    elseif (x_head <= x(i) && x(i) < x_tail)
			p(i) = p_medium;
    elseif (x(i) < Z2 * t && x(i) >= x_tail)
			p(i) = p_star;
	else
			p(i) = p2;
    end
          end