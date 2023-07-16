function [rho,u,v,p] = ctop2(U)

global gamma;

rho = U(1);
u = U(2)./U(1);
v = U(3)./U(1);
p = (U(4)./rho-0.5*(u.^2+v.^2))*(gamma+1)*rho;

end
