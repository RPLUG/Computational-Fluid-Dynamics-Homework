   
function U = physics_to_conservation(rho,vel,pre)
%将物理量转换成守恒量
global gamma;

%给出守恒量
U(1,:) = rho;
U(2,:) = rho.*vel;
U(3,:) = pre/(gamma-1) + 0.5*rho.*vel.^2;

end