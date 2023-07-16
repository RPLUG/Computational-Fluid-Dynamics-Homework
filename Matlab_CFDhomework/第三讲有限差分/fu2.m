function [ux,ut] = fu2(Nmesh,j,dx,u)%五阶精确格式
%计算空间点的函数
if(j>3&&j<Nmesh-1)
    ux=(-2*u(j-3)+15*u(j-2)-60*u(j-1)+20*u(j)+30*u(j+1)-3*u(j+2))/(60*dx);%空间导数差分
    ut=-ux;%时间导数差分
elseif(j==Nmesh-1)
    ux=(-2*u(j-3)+15*u(j-2)-60*u(j-1)+20*u(j)+30*u(j+1)-3*u(1))/(60*dx);%空间循环求导数差分
    ut=-ux;%时间导数差分
elseif(j==Nmesh)
    ux=(-2*u(j-3)+15*u(j-2)-60*u(j-1)+20*u(j)+30*u(1)-3*u(2))/(60*dx);%空间循环求导数差分
    ut=-ux;%时间导数差分
elseif(j==1)
       ux= (-2*u(Nmesh-2)+15*u(Nmesh-1)-60*u(Nmesh)+20*u(j)+30*u(j+1)-3*u(j+2))/(60*dx);%空间导数差分
        ut=-ux;%时间导数差分
elseif(j==2)
        ux=(-2*u(Nmesh-1)+15*u(Nmesh)-60*u(j-1)+20*u(j)+30*u(j+1)-3*u(j+2))/(60*dx);%空间导数差分
        ut=-ux;%时间导数差分
elseif(j==3)
        ux=(-2*u(Nmesh)+15*u(j-2)-60*u(j-1)+20*u(j)+30*u(j+1)-3*u(j+2))/(60*dx);%空间导数差分
        ut=-ux;%时间导数差分
end