clc;clear;format compact;close all;fclose all; 
m = 20; %y方向每单位长度网格数
n = 20; %x方向每单位长度网格数
%%
%绘制计算空间中的下方边界网格
x_plane=(exp(linspace(0,1,n))-1)/(exp(1)-1);%翼型按指数分布，钝头加密，横向布点
x_res=1+2.*(exp(linspace(1,3,2*n))-exp(1))/(exp(3)-exp(1));
x_res(1)=[];
x_r=[x_plane x_res];
x_l=fliplr(x_r);
x_l(end)=[];
x1=[x_l x_r];
y_r= [1i*(0.1781*sqrt(x_plane)-...
  0.0756*x_plane-0.2122*x_plane.^2+0.1705*x_plane.^3-...
     0.0609*x_plane.^4+0.0001) 0i*x_res];
 y_l=-fliplr(y_r);
y_l(end)=[];
y1=[y_l y_r];
xy1=x1+y1;
%%
%绘制计算空间中的上方边界网格
x_d=linspace(3,0,2*n);
x_u=fliplr(x_d);
x_ru=linspace(-2,0,n);
x_ru(end)=[];
x_rd=fliplr(x_ru);
x_rd(end)=[];
x_2=[x_d x_rd x_ru x_u];
y_d=-2i*linspace(1,1,2*n);
y_u=-fliplr(y_d);
y_ru=1i*sqrt((4-linspace(-2,0,n).^2));
y_ru(end)=[];
y_rd=-fliplr(y_ru);
y_rd(end)=[];
y_2=[y_d y_rd y_ru y_u];
xy2=x_2+y_2;
%%
%侧面边界布点，靠近中轴位置同样以指数函数加密
yR=-fliplr((1i*2*(exp(linspace(0,2,2*m))-1)/(exp(2)-1)))';
yR(end)=[];yR(1)=[];
yL=fliplr((1i*2*(exp(linspace(0,2,2*m))-1)/(exp(2)-1)))';
yL(end)=[];yL(1)=[];
xR=linspace(3,3,2*m-2)';xL=linspace(3,3,2*m-2)';
xyR=xR+yR;xyL=xL+yL;
%%
%创建各点迭代初值，注意只控制了C型网格的边界，内部点初值旨在保证各不相等以使迭代式分母不为0
xy=linspace(-1,2,6*n-3);
t=linspace(0,2,2*m)';
xy = repmat(xy, 2*m, 1) + t/(2)*(xy2 - xy1);
xy(1,:)=xy2;
xy(end,:)=xy1;
xy(2:end-1,1)=xyL;
xy(2:end-1,end)=xyR;
%%
% 用elliptic_mesh函数求解椭圆型方程，得到网格坐标
 xySmooth = elliptic_mesh(xy); 
% 绘图比较
figure 
subplot(1, 2, 1), hold on, axis equal 
plot(xy, 'k'), plot(xy.', 'k') 
title('初始网格') 
 subplot(1, 2, 2), hold on, axis equal 
 plot(xySmooth, 'm'), plot(xySmooth.', 'm') 
 title('结果网格')