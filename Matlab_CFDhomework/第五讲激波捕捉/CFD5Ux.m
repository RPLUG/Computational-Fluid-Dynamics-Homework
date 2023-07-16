function [U_x,U_t] = CFD5Ux(U_former)

ga=1.4;dx=0.01;
rho_former=U_former(1,:);
u_former=U_former(2,:)./U_former(1,:);
p_former=(ga-1)*(U_former(3,:)-U_former(2,:).^2./(2*U_former(1,:)));
H_former=(U_former(3,:)+p_former)./rho_former;

for k=1:98%采用GVC格式，计算u_R和u_L，0.14时波传播不到边界，格点的值简化为初始值不变
    [rho_L(k),rho_R(k)]=GVC(k+1,rho_former);
    [u_L(k),u_R(k)]=GVC(k+1,u_former);
    [p_L(k),p_R(k)]=GVC(k+1,p_former);
    [H_L(k),H_R(k)]=GVC(k+1,H_former);
     UL(:,k)=[rho_L(k);rho_L(k)*u_L(k);rho_L(k)*H_L(k)-p_L(k)];
    UR(:,k)=[rho_R(k);rho_R(k)*u_R(k);rho_R(k)*H_R(k)-p_R(k)];
    fU_L(:,k)=[rho_L(k)* u_L(k);rho_L(k)* u_L(k)^2+p_L(k);u_L(k)*H_L(k)*rho_L(k)];
    fU_R(:,k)=[rho_R(k)* u_R(k);rho_R(k)* u_R(k)^2+p_R(k);u_R(k)*H_R(k)*rho_R(k)];
end
for k=1:98%计算roe平均值
    rho_ba(k)=((rho_L(k)^0.5+rho_R(k)^0.5)/2)^2;
    u_ba(k)=(rho_L(k)^0.5*u_L(k)+rho_R(k)^0.5*u_R(k))/(rho_L(k)^0.5+rho_R(k)^0.5);
    H_ba(k)=(rho_L(k)^0.5*H_L(k)+rho_R(k)^0.5*H_R(k))/(rho_L(k)^0.5+rho_R(k)^0.5);
    p_ba(k)=(ga-1)/ga*(rho_ba(k)*H_ba(k)-1/2*rho_ba(k)*u_ba(k)^2);
    c_ba(k)=sqrt((ga-1)*(H_ba(k)-0.5*u_ba(k)^2));
end
for k=1:98
S_inv(:,:,k)=[(1-ga)/c_ba(k)^2,-1/(2*c_ba(k)),1/(2*c_ba(k));
    (1-ga)*u_ba(k)/c_ba(k)^2,(c_ba(k)-u_ba(k))/(2*c_ba(k)),(c_ba(k)+u_ba(k))/(2*c_ba(k));
    (1-ga)*(u_ba(k)^2)/(2*c_ba(k)^2),(c_ba(k)^2/(1-ga)+c_ba(k)*u_ba(k)-u_ba(k)^2/2)/(2*c_ba(k)),(c_ba(k)^2/(ga-1)+c_ba(k)*u_ba(k)+u_ba(k)^2/2)/(2*c_ba(k))];
v_Lamda=[abs(u_ba(k)),abs(u_ba(k)-c_ba(k)),abs(u_ba(k)+c_ba(k))];
Lamda(:,:,k)=diag(v_Lamda);
S(:,:,k)=[u_ba(k)^2/2-c_ba(k)^2/(ga-1),-u_ba(k),1;-u_ba(k)-u_ba(k)^2*(ga-1)/(2*c_ba(k)),u_ba(k)*(ga-1)/c_ba(k)+1,(1-ga)/c_ba(k);u_ba(k)^2*(ga-1)/(2*c_ba(k))-u_ba(k),1-u_ba(k)*(ga-1)/c_ba(k),(ga-1)/c_ba(k)];
ABS_A(:,:,k)=S_inv(:,:,k)*Lamda(:,:,k)*S(:,:,k);
flux(:,k)=0.5*( fU_L(:,k)+fU_R(:,k))-0.5*ABS_A(:,:,k)*(UR(:,k)-UL(:,k));
end
for i=1:101%与101个格点对应的U的空间导数初值
U_x(:,i)=[0;0;0];
end
for k=1:97%计算空间导数值，同样认为边界格点值不变，故仅更新“中间”
U_x(:,k+2)=(flux(:,k+1)-flux(:,k))/dx;
end
U_t=-U_x;
end

