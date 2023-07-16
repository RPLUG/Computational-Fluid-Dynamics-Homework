%作图
[rho_numerical,u_numerical,p_numerical] = conservation_to_physics(U);
EulerExact;
 plot(x,rho,'m','DisplayName','精确解密度');hold on;
 plot(x,u,'k','DisplayName','精确解速度');hold on;
 plot(x,p,'r','DisplayName','精确解压强');hold on; 
 plot(x,rho_numerical,'m+','DisplayName','数值解密度');hold on;
 plot(x,u_numerical,'k+','DisplayName','数值解速度');hold on;
 plot(x,p_numerical,'r+','DisplayName','数值解压强');hold on;
 legend('show');
 title('数值解与精确解对比图');xlabel('x');ylabel('各物理量值');

