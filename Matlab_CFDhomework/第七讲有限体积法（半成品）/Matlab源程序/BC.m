%% 设置虚网格原始物理量值，代表边界条件，虚网格的值由边界条件和插值方式决定
%拓扑关系上，虚网格矩阵与内部点矩阵对应，左w右e上n下s，行列序号与矩阵默认方式相同
%n，钝楔表面,固体壁面，壁温为来流温度的4.2倍
rho_dummyN(2,:)=0.5*(rho(1,:)+rho(2,:));rho_dummyN(1,:)=5*rho_dummyN(2,:)-4*rho(1,:);
p_dummyN(2,:)=0.5*(p(1,:)+p(2,:));p_dummyN(1,:)=5*p_dummyN(2,:)-4*p(1,:);
u_dummyN(2,:)=0.5*(u(2,:)-5*u(1,:));u_dummyN(1,:)=2*u(1,:)+5*u_dummyN(2,:);
v_dummyN(2,:)=0.5*(v(2,:)-5*v(1,:));v_dummyN(1,:)=2*v(1,:)+5*v_dummyN(2,:);
T_dummyN=8.4-T(1,:);
% 
% rho_dummyN(2,:)=rho(1,:);rho_dummyN(1,:)=rho(1,:);
% p_dummyN(2,:)=p(1,:);p_dummyN(1,:)=p(1,:);
% u_dummyN(2,:)=-u(1,:);u_dummyN(1,:)=-u(1,:);
% v_dummyN(2,:)=-v(1,:);v_dummyN(1,:)=-v(1,:);
% T_dummyN=4.2;
%s，超音速入口
rho_dummyS(1,:)=0.5*(rho(end-1,:)+rho(end,:));rho_dummyS(2,:)=5*rho_dummyS(1,:)-4*rho(end,:);
p_dummyS(1,:)=0.5*(p(end-1,:)+p(end,:));p_dummyS(2,:)=5*p_dummyS(1,:)-4*p(end,:);
u_dummyS(1,:)=0.5*(u(end-1,:)+6-5*u(end,:));u_dummyS(2,:)=-6+2*u(end,:)+5*u_dummyS(1,:);%入口u=1
v_dummyS(1,:)=0.5*(v(end-1,:)-5*v(end,:));v_dummyS(2,:)=2*u(end,:)+5*v_dummyS(1,:);%入口v=0
T_dummyS=2-T(end,:);
%  rho_dummyS=ones(2,N-1);
%  u_dummyS=ones(2,N-1);v_dummyS=ones(2,N-1);p_dummyS=ones(2,N-1);
%  T_dummyS=ones(1,N-1);
%e,w，超音速出口
rho_dummyW(:,2)=0.5*(rho(:,1)+rho(:,2));rho_dummyW(:,1)=5*rho_dummyW(:,2)-4*rho(:,1);
p_dummyW(:,2)=0.5*(p(:,1)+p(:,2));p_dummyW(:,1)=5*p_dummyW(:,2)-4*p(:,1);
u_dummyW(:,2)=0.5*(u(:,1)+u(:,2));u_dummyW(:,1)=5*u_dummyW(:,2)-4*u(:,1);
v_dummyW(:,2)=0.5*(v(:,1)+v(:,2));v_dummyW(:,1)=5*v_dummyW(:,2)-4*v(:,1);
T_dummyW=0.5*(T(:,1)+T(:,2));
rho_dummyE(:,1)=0.5*(rho(:,end)+rho(:,end-1));rho_dummyE(:,2)=5*rho_dummyE(:,1)-4*rho(:,end);
p_dummyE(:,1)=0.5*(p(:,end)+p(:,end-1));p_dummyE(:,2)=5*p_dummyE(:,1)-4*p(:,end);
u_dummyE(:,1)=0.5*(u(:,end)+u(:,end-1));u_dummyE(:,2)=5*u_dummyE(:,1)-4*u(:,end);
v_dummyE(:,1)=0.5*(v(:,end)+v(:,end-1));v_dummyE(:,2)=5*v_dummyE(:,1)-4*v(:,end);
T_dummyE=0.5*(T(:,end)+T(:,end-1));
