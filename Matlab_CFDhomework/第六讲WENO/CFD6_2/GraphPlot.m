%作图
%求解一维Shu-Osher激波管问题
%计算参数设置，按原始文献参考值
global N;%节点数
N = 2001;
T = 1.8;
x = linspace(0,10,N);
dx = x(2) - x(1);
dt = 0.1*dx;
Tnumber = T/dt;
%物理参数
global gamma;
gamma = 1.4;
%设置初始条件
[rho0,u0,p0] = IC(x);
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2;  % 总能量密度
a0 = sqrt(gamma*p0./rho0);            % 声速
%开始计算
U_laste = physics_to_conservation(rho0,u0,p0);
for i = 1:Tnumber
 %时间离散使用三阶RK3格式
    U_laste = BC(U_laste,3);
    U1e = U_laste(:,4:end-3) - dt*FluxDiffRoe(U_laste,dx);
    U1e = BC(U1e,3);
    U2e = 0.75*U_laste(:,4:end-3) + 0.25*U1e(:,4:end-3) - 0.25*dt*FluxDiffRoe(U1e,dx);
    U2e = BC(U2e,3);
    Ue = 1/3*U_laste(:,4:end-3) + 2/3*U2e(:,4:end-3) - 2/3*dt*FluxDiffRoe(U2e,dx);
    %维护循环并作图
    U_laste = Ue;
end
[rho,u,p] = conservation_to_physics(Ue);
  plot(x,rho,'m','DisplayName','精确解密度');hold on;
  plot(x,u,'k','DisplayName','精确解速度');hold on;
  plot(x,p,'g','DisplayName','精确解压强');hold on; 
 legend('show');
 title('数值解与精确解对比图');xlabel('x');ylabel('各物理量值');

