%��ͼ
%���һάShu-Osher����������
%����������ã���ԭʼ���ײο�ֵ
global N;%�ڵ���
N = 2001;
T = 1.8;
x = linspace(0,10,N);
dx = x(2) - x(1);
dt = 0.1*dx;
Tnumber = T/dt;
%�������
global gamma;
gamma = 1.4;
%���ó�ʼ����
[rho0,u0,p0] = IC(x);
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2;  % �������ܶ�
a0 = sqrt(gamma*p0./rho0);            % ����
%��ʼ����
U_laste = physics_to_conservation(rho0,u0,p0);
for i = 1:Tnumber
 %ʱ����ɢʹ������RK3��ʽ
    U_laste = BC(U_laste,3);
    U1e = U_laste(:,4:end-3) - dt*FluxDiffRoe(U_laste,dx);
    U1e = BC(U1e,3);
    U2e = 0.75*U_laste(:,4:end-3) + 0.25*U1e(:,4:end-3) - 0.25*dt*FluxDiffRoe(U1e,dx);
    U2e = BC(U2e,3);
    Ue = 1/3*U_laste(:,4:end-3) + 2/3*U2e(:,4:end-3) - 2/3*dt*FluxDiffRoe(U2e,dx);
    %ά��ѭ������ͼ
    U_laste = Ue;
end
[rho,u,p] = conservation_to_physics(Ue);
  plot(x,rho,'m','DisplayName','��ȷ���ܶ�');hold on;
  plot(x,u,'k','DisplayName','��ȷ���ٶ�');hold on;
  plot(x,p,'g','DisplayName','��ȷ��ѹǿ');hold on; 
 legend('show');
 title('��ֵ���뾫ȷ��Ա�ͼ');xlabel('x');ylabel('��������ֵ');

