clc;clear;Mesh;
Re=400;nu=U*L/Re;dt=0.002;%��ŵ����ճ���Լ�ʱ�䲽��
NewIC;Tol=1;
it=0;
while(Tol>0.0002)
%% ���ٶȸ����ڲ�����
% ������
uomega_x=u(2:end-1,2:end-1).*(omega(2:end-1,3:end)-omega(2:end-1,1:end-2))/(2*dx);
vomega_y=v(2:end-1,2:end-1).*(omega(1:end-2,2:end-1)-omega(3:end,2:end-1))/(2*dy);
%ճ����
omega2_x2=(omega(2:end-1,3:end)+omega(2:end-1,1:end-2)-2*omega(2:end-1,2:end-1))/(dx^2);
omega2_y2=(omega(3:end,2:end-1)+omega(1:end-2,2:end-1)-2*omega(2:end-1,2:end-1))/(dy^2);
%������ʱ�䵼��
omega_t(2:end-1,2:end-1)=nu*(omega2_x2+omega2_y2)-uomega_x-vomega_y;
%��������
omega2(2:end-1,2:end-1)=omega_t(2:end-1,2:end-1)*dt+omega(2:end-1,2:end-1);
%% ����������������
for iterpsi=1:2%ÿ�ε������Ĵ����������������
psitempn=psi(1:end-2,2:end-1);
psitemps=psi(3:end,2:end-1);
psitempe=psi(2:end-1,3:end);
psitempw=psi(2:end-1,1:end-2);
psi(2:end-1,2:end-1)=1/(2*dx^2+2*dy^2)*(dx^2*(psitemps+psitempn)+dy^2*(psitempe+psitempw)-dx^2*dy^2*omega2(2:end-1,2:end-1));
end
%% �����������±߽�����
omega2(:,1)=2*psi(:,2)/dx^2;%w
omega2(:,end)=2*psi(:,end-1)/dx^2;%e
omega2(end,:)=2*psi(end-1,:)/dy^2;%s
omega2(1,2:end-1)=2*(psi(2,2:end-1)+u(1,2:end-1)*dy)/dy^2;%n
%% �������������ٶ�
u(2:end-1,2:end-1)=(psi(1:end-2,2:end-1)-psi(3:end,2:end-1))/2/dy;
v(2:end-1,2:end-1)=(psi(2:end-1,1:end-2)-psi(2:end-1,3:end))/2/dx;
%% �����о�
Tol=max(max((abs(omega2-omega))));
%% ��������
omega=omega2;
it=it+1;
end
%% ��ͼ
contour(real(XY),imag(XY),psi,[diag(psi(1:5:end,1:5:end));diag(psi(1:-5:end,end:5:1))],'-m');
 title('Re = 400')