%初始条件设置
function [rho,u,p] = IC(x)
         x1 = x<=1;
        x2 = x>1;
        rho = 3.857*x1 + (1+0.2*sin(5.0.*x)).*x2;
        u   = 2.629*x1+0*x2;
        p   = 10.333*x1 + 1*x2;     
end