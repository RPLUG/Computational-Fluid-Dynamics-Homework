clc
clear
Nmesh=20;%������
a4=2.4;%�ɱ����ֵ
dx=2*pi/Nmesh;dt=0.01;t1=20;t2=50;
x=[0:dx:2*pi-dx];
u_pre=sin(x);u5_pre=sin(x);u6=sin(x);%�������ʱ��u��sinxΪ��ֵ
u=[1:1:Nmesh];u_5=[1:1:Nmesh];u_x=[1:1:Nmesh];u_t=[1:1:Nmesh];%�������ʱ��u_x,u_t
K1=[1:1:Nmesh];K2=[1:1:Nmesh];%�������ʱRK���Ĳ���U1��U2
u_true1=sin(x-t1);%t=20��ȷ��
u_true2=sin(x-t2);%%t=50��ȷ��
t=0;%��ʼʱ��
while(t<t2)
    for i=1:Nmesh%���ѭ�������ʼ���ڵĸ���ֵ�͸���ռ�ʱ�䵼��,i����ռ�λ�õı�ţ��߽�����ʹ����ѭ������
    [u_x(i),u_t(i)]=fu(Nmesh,i,a4,dx,u_pre);%����һ��Un����ռ�ʱ�䵼��
    end
    for i=1:Nmesh%ʹ��RK����ʱ������ƽ���������u(i)����Ϊ��һʱ�䲽��ֵ
    K1(i)=u_pre(i)+u_t(i)*dt;%�ƽ���U1�������K��U1
    end
    for i=1:Nmesh
    [u_x(i),u_t(i)]=fu(Nmesh,i,a4,dx,K1);%����L(U1)
    end
    for i=1:Nmesh
    K2(i)=3/4*u_pre(i)+1/4*(K1(i)+u_t(i)*dt);%��U1��L(U1)���µõ�U2
    end
    for i=1:Nmesh
    [u_x(i),u_t(i)]=fu(Nmesh,i,a4,dx,K2);%����L(U2)
    end
    for i=1:Nmesh
    u(i)=1/3*u_pre(i)+2/3*(K2(i)+dt*u_t(i));%��U2��L(U2)���µõ�Un+1
    end
    u_pre=u;
    t=t+dt;
end
t=0;
while(t<t2)%�Ա���
    for i=1:Nmesh%���ѭ�������ʼ���ڵĸ���ֵ�͸���ռ�ʱ�䵼��,i����ռ�λ�õı�ţ��߽�����ʹ����ѭ������
    [u_x(i),u_t(i)]=fu2(Nmesh,i,dx,u5_pre);%����һ��Un����ռ�ʱ�䵼��
    end
    for i=1:Nmesh%ʹ��RK����ʱ������ƽ���������u(i)����Ϊ��һʱ�䲽��ֵ
    K1(i)=u5_pre(i)+u_t(i)*dt;%�ƽ���U1�������K��U1
    end
    for i=1:Nmesh
    [u_x(i),u_t(i)]=fu2(Nmesh,i,dx,K1);%����L(U1)
    end
    for i=1:Nmesh
    K2(i)=3/4*u5_pre(i)+1/4*(K1(i)+u_t(i)*dt);%��U1��L(U1)���µõ�U2
    end
    for i=1:Nmesh
    [u_x(i),u_t(i)]=fu2(Nmesh,i,dx,K2);%����L(U2)
    end
    for i=1:Nmesh
    u5(i)=1/3*u5_pre(i)+2/3*(K2(i)+dt*u_t(i));%��U2��L(U2)���µõ�Un+1
    end
    u5_pre=u5;
    t=t+dt;
end
t=0;

 plot(x,u,'m','DisplayName','�Ľ׾���u');hold on;
plot(x,u5,'--k','DisplayName','��׾���u');hold on;
  plot(x,u_true2,'--r','DisplayName','�ο���ȷ��');
  title('��ֵ���뾫ȷ��Ա�ͼ');xlabel('x');ylabel('��ֵ');
legend('show');
sumu=0;sumu5=0;
for n=1:Nmesh
    error=(u(n)-u_true2(n))^2*dx;
    sumu=sumu+error;
end
sum1=sqrt(sumu);
for n=1:Nmesh
    error=(u5(n)-u_true2(n))^2*dx;
    sumu5=sumu5+error;
end
sum2=sqrt(sumu5);
disp(sum1);
disp(sum2);
