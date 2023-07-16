%在[-0.5，0.5]区间上求解，最后统一进行坐标变换的平移
ga=1.4;
t = 0.14;
u1 = 0; rho1 = 1; p1 = 1; u2 = 0; rho2 = 0.125; p2 = 0.1;
c1 = sqrt(ga * p1 / rho1);
F=@(p_star)(f1(p_star, p1, rho1) + f1(p_star, p2, rho2));
   p_star=fzero(F,0.3);
	rho_starl = rho1/(p1 / p_star)^  (1 / ga);%密度间断面左侧密度
	rho_starr = rho2 * ((ga - 1) * p2 + (ga + 1) * p_star) /( (ga + 1) * p2 + (ga - 1) * p_star);%密度间断面右侧密度
	c_starl = sqrt(ga * p_star / rho_starl);%密度间断面左侧波速
	u_star = u1 - f1(p_star, p1, rho1);%膨胀波与激波间流速
	u_head = u1 - c1;%左波波头速度
	u_tail = u_star - c_starl;%左波波尾速度
	x_head = u_head * t;%膨胀波波头位置
	x_tail = u_tail * t;%膨胀波波尾位置
	Z2 = (rho1 * u2 - rho_starr * u_star) / (rho2 - rho_starr);%右激波传播速度
	x_interval = u_star * t;%密度间断面位置
    for i=1:101
        x(i)=-0.51+0.01*i;%速度分布
		c_medium = (ga - 1) / (ga + 1) * (u1 - x(i) / t) + 2 / (ga + 1) * c1;%膨胀波内部波速分布
		u_medium = c_medium + x(i) / t;%膨胀波内部的速度分布
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
        x(i)=-0.51+2/200*i;%密度分布
		c_medium = (ga - 1) / (ga + 1) * (u1 - x(i) / t) + 2 / (ga + 1) * c1;%膨胀波内部波速分布
		u_medium = c_medium + x(i) / t;%膨胀波内部的速度分布
		p_medium = p1*(c_medium/c1)^(2*ga/(ga-1)) ;%膨胀波内部压强分布
		rho_medium = ga * p_medium / c_medium^2;%膨胀波内部密度分布
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
        x(i)=-0.51+2/200*i;%速度分布
		c_medium = (ga - 1) / (ga + 1) * (u1 - x(i) / t) + 2 / (ga + 1) * c1;%膨胀波内部波速分布
		u_medium = c_medium + x(i) / t;%膨胀波内部的速度分布
		p_medium = p1*(c_medium/c1)^(2*ga/(ga-1)) ;%膨胀波内部压强分布
		rho_medium = ga * p_medium / c_medium^2;%膨胀波内部密度分布
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