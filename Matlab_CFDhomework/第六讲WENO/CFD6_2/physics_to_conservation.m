   
function U = physics_to_conservation(rho,vel,pre)
%��������ת�����غ���
global gamma;

%�����غ���
U(1,:) = rho;
U(2,:) = rho.*vel;
U(3,:) = pre/(gamma-1) + 0.5*rho.*vel.^2;

end