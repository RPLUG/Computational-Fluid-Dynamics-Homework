function [rho,vel,pre] = conservation_to_physics(U)
%�غ���ת��Ϊԭʼ����
global gamma;
rho = U(1,:);
vel = U(2,:)./U(1,:);
pre = (U(3,:)-0.5*U(2,:).^2./U(1,:))*(gamma-1);
end
