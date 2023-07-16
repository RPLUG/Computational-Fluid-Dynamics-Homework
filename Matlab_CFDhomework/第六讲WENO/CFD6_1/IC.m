%初始条件设置
function [rho,u,p] = IC(x)
         x1 = x<=0.5;
        x2 = x>0.5;
        rho = x1 + 0.125*x2;
        u   = 0;
        p   = x1 + 0.1*x2;     
end