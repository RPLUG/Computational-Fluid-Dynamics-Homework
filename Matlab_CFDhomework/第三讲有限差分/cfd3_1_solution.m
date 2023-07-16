clc
clear
Nmesh=20;%网格数
a4=2.4;%可变参数值
dx=2*pi/Nmesh;dt=0.01;t1=20;t2=50;
x=[0:dx:2*pi-dx];
u_pre=sin(x);u5_pre=sin(x);u6=sin(x);%储存迭代时的u，sinx为初值
u=[1:1:Nmesh];u_5=[1:1:Nmesh];u_x=[1:1:Nmesh];u_t=[1:1:Nmesh];%储存迭代时的u_x,u_t
K1=[1:1:Nmesh];K2=[1:1:Nmesh];%储存迭代时RK法的参数U1，U2
u_true1=sin(x-t1);%t=20精确解
u_true2=sin(x-t2);%%t=50精确解
t=0;%初始时刻
while(t<t2)
    for i=1:Nmesh%这个循环计算初始层内的各点值和各点空间时间导数,i代表空间位置的编号，边界条件使用了循环计算
    [u_x(i),u_t(i)]=fu(Nmesh,i,a4,dx,u_pre);%由上一步Un计算空间时间导数
    end
    for i=1:Nmesh%使用RK法对时间进行推进，最终让u(i)更新为下一时间步的值
    K1(i)=u_pre(i)+u_t(i)*dt;%推进到U1，这里的K是U1
    end
    for i=1:Nmesh
    [u_x(i),u_t(i)]=fu(Nmesh,i,a4,dx,K1);%计算L(U1)
    end
    for i=1:Nmesh
    K2(i)=3/4*u_pre(i)+1/4*(K1(i)+u_t(i)*dt);%由U1和L(U1)更新得到U2
    end
    for i=1:Nmesh
    [u_x(i),u_t(i)]=fu(Nmesh,i,a4,dx,K2);%计算L(U2)
    end
    for i=1:Nmesh
    u(i)=1/3*u_pre(i)+2/3*(K2(i)+dt*u_t(i));%由U2和L(U2)更新得到Un+1
    end
    u_pre=u;
    t=t+dt;
end
t=0;
while(t<t2)%对比用
    for i=1:Nmesh%这个循环计算初始层内的各点值和各点空间时间导数,i代表空间位置的编号，边界条件使用了循环计算
    [u_x(i),u_t(i)]=fu2(Nmesh,i,dx,u5_pre);%由上一步Un计算空间时间导数
    end
    for i=1:Nmesh%使用RK法对时间进行推进，最终让u(i)更新为下一时间步的值
    K1(i)=u5_pre(i)+u_t(i)*dt;%推进到U1，这里的K是U1
    end
    for i=1:Nmesh
    [u_x(i),u_t(i)]=fu2(Nmesh,i,dx,K1);%计算L(U1)
    end
    for i=1:Nmesh
    K2(i)=3/4*u5_pre(i)+1/4*(K1(i)+u_t(i)*dt);%由U1和L(U1)更新得到U2
    end
    for i=1:Nmesh
    [u_x(i),u_t(i)]=fu2(Nmesh,i,dx,K2);%计算L(U2)
    end
    for i=1:Nmesh
    u5(i)=1/3*u5_pre(i)+2/3*(K2(i)+dt*u_t(i));%由U2和L(U2)更新得到Un+1
    end
    u5_pre=u5;
    t=t+dt;
end
t=0;

 plot(x,u,'m','DisplayName','四阶精度u');hold on;
plot(x,u5,'--k','DisplayName','五阶精度u');hold on;
  plot(x,u_true2,'--r','DisplayName','参考精确解');
  title('数值解与精确解对比图');xlabel('x');ylabel('解值');
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
