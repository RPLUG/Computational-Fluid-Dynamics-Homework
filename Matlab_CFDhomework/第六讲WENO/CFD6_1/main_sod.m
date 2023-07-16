%求解一维sod激波管问题
clc; clear;
format long e;
%计算参数设置
global N;
N = 101;
T = 0.14;
x = linspace(0,1,N);
dx = x(2) - x(1);
dt = dx*0.1;
Tnumber = T/dt;
%物理参数
global gamma;
gamma = 1.4;
%设置初始条件
[rho0,u0,p0] = IC(x);
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2;  % 总能量密度
a0 = sqrt(gamma*p0./rho0);            % 声速
%开始计算
U_last = physics_to_conservation(rho0,u0,p0);
for i = 1:Tnumber
 %时间离散使用三阶RK3格式
    U_last = BC(U_last,3);
    U1 = U_last(:,4:end-3) - dt*FluxDiffRoe(U_last,dx);
    U1 = BC(U1,3);
    U2 = 0.75*U_last(:,4:end-3) + 0.25*U1(:,4:end-3) - 0.25*dt*FluxDiffRoe(U1,dx);
    U2 = BC(U2,3);
    U = 1/3*U_last(:,4:end-3) + 2/3*U2(:,4:end-3) - 2/3*dt*FluxDiffRoe(U2,dx);
    %维护循环并作图
    U_last = U;
end
    GraphPlot;

