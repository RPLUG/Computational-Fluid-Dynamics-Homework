%�õ�һ����ά����洢�غ��� �����������ܶ�����������ܶȱ�ʾ������������ʦppt�Լ�toro����һ��
function U = ptoc2(rho,u,v,p)
global gamma;

%�����غ���
U(:,:,1) = rho;
U(:,:,2) = rho.*u;
U(:,:,3) = rho.*v;
U(:,:,4) =p/(gamma-1) + 0.5*rho.*(v.^2+u.^2);

end