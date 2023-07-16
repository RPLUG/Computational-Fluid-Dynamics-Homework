function Flux = RoeFlux(UtE,UtW,M,N,nxE,nxW,nyE,nyW)
global gamma;
%% 计算各个roe平均值
Flux=UtE;%确定输出矩阵大小
utE=UtE(:,:,2).*nxE+UtE(:,:,3).*nyE;vtE=-UtE(:,:,2).*nxE+UtE(:,:,3).*nyE;
utW=UtW(:,:,2).*nxE+UtW(:,:,3).*nyE;vtW=-UtW(:,:,2).*nxE+UtW(:,:,3).*nyE;
rhotW=UtW(:,:,1);etW=UtW(:,:,4)./UtW(:,:,1);
ptW=(gamma - 1) *rhotW.* (etW - 0.5 * (utW .*utW + vtW.*vtW));
 HtW = etW + ptW./ rhotW;
rhotE=UtE(:,:,1);etE=UtE(:,:,4)./UtE(:,:,1);
ptE=(gamma - 1) *rhotE.* (etE - 0.5 * (utE .*utE + vtE.*vtE));
HtE = etE + ptE ./ rhotE;
rho_roe = ((sqrt(rhotW)+sqrt(rhotE))/2).^2;
u_roe = (sqrt(rhotW).*utW + sqrt(rhotE).*utE)./(sqrt(rhotW)+sqrt(rhotE));
v_roe = (sqrt(rhotW).*vtW + sqrt(rhotE).*vtE)./(sqrt(rhotW)+sqrt(rhotE));
H_roe = (sqrt(rhotW).*HtW + sqrt(rhotE).*HtE)./(sqrt(rhotW)+sqrt(rhotE));
p_roe = (gamma-1)/gamma*(rho_roe.*H_roe-0.5*rho_roe.*u_roe.^2);
c_roe = sqrt((gamma-1)*(H_roe-0.5*u_roe.^2));
vel_roe=sqrt(v_roe.^2+u_roe.^2);
%% 计算特征值，记特征值矩阵lamda,相似变换矩阵S，S_inv
lamda1=abs(u_roe-c_roe);
lamda2=abs(u_roe);lamda3=abs(u_roe);
lamda4=abs(u_roe+c_roe);
%EW，79*240，NS，80*239
for m=1:M
    for n=1:N
lamda=[lamda1(m,n),0,0,0;0,lamda2(m,n),0,0;0,0,lamda3(m,n),0;0,0,0,lamda4(m,n)];
S(1,1) = 1;
S(1,2) = 1;
S(1,3) = 0;
S(1,4) = 1;
S(2,1) = u_roe(m,n) - c_roe(m,n) ;
S(2,2) = u_roe(m,n);
S(2,3) = 0;
S(2,4) = u_roe(m,n) + c_roe(m,n);
S(3,1) = v_roe(m,n);
S(3,2) = v_roe(m,n);
S(3,3) = 1;
S(3,4) = v_roe(m,n);
S(4,1) = H_roe(m,n)-c_roe(m,n).*u_roe(m,n);
S(4,2) = 0.5*vel_roe(m,n).^2;
S(4,3) = v_roe(m,n);
S(4,4)=H_roe(m,n)+u_roe(m,n).*c_roe(m,n);

S_inv(1,1) = (- vel_roe(m,n).^2.*u_roe(m,n) - c_roe(m,n).*vel_roe(m,n).^2 + 2*c_roe(m,n).*u_roe(m,n).^2 + 2*H_roe(m,n).*u_roe(m,n) + 2*c_roe(m,n).*v_roe(m,n).^2)./(2*(- c_roe(m,n).*vel_roe(m,n).^2 + 2*H_roe(m,n).*c_roe(m,n)));
S_inv(1,2) = -(- vel_roe(m,n).^2 + 2*H_roe(m,n) + 2*c_roe(m,n).*u_roe(m,n))/(2*(- c_roe(m,n).*vel_roe(m,n).^2 + 2*H_roe(m,n).*c_roe(m,n)));
S_inv(1,3) = -v_roe(m,n)./(- vel_roe(m,n).^2 + 2*H_roe(m,n));
S_inv(1,4) = 1./(- vel_roe(m,n).^2 + 2*H_roe(m,n));
S_inv(2,1)= -(2*(u_roe(m,n).^2 + v_roe(m,n).^2 - H_roe(m,n)))./(- vel_roe(m,n).^2 + 2*H_roe(m,n));
S_inv(2,2) =(2*u_roe(m,n))./(- vel_roe(m,n).^2 + 2*H_roe(m,n));
S_inv(2,3) = (2*v_roe(m,n))./(- vel_roe(m,n).^2 + 2*H_roe(m,n));
S_inv(2,4) = -2./(- vel_roe(m,n).^2 + 2*H_roe(m,n));
S_inv(3,1) = -v_roe(m,n);
S_inv(3,2) = 0;
S_inv(3,3) = 1;
S_inv(3,4) = 0;
S_inv(4,1)=(vel_roe(m,n).^2.*u_roe(m,n) - c_roe(m,n).*vel_roe(m,n).^2 + 2*c_roe(m,n).*u_roe(m,n).^2 - 2*H_roe(m,n).*u_roe(m,n) + 2*c_roe(m,n).*v_roe(m,n).^2)./(2*(- c_roe(m,n).*vel_roe(m,n).^2 + 2*H_roe(m,n).*c_roe(m,n)));
S_inv(4,2)=-(vel_roe(m,n).^2 - 2*H_roe(m,n) + 2*c_roe(m,n).*u_roe(m,n))/(2*(- c_roe(m,n).*vel_roe(m,n).^2 + 2*H_roe(m,n).*c_roe(m,n)));
S_inv(4,3)=-v_roe(m,n)/(- vel_roe(m,n).^2 + 2*H_roe(m,n));
S_inv(4,4)=1./(- vel_roe(m,n).^2 + 2*H_roe(m,n));
A=S_inv*lamda*S;
UL=[UtW(m,n,1);UtW(m,n,2);UtW(m,n,3);UtW(m,n,4);];
UR=[UtE(m,n,1);UtE(m,n,2);UtE(m,n,3);UtE(m,n,4);];
[rhol,ul,vl,pl]=ctop2(UL);[rhor,ur,vr,pr]=ctop2(UR);
el=pl/(gamma+1)+0.5*(ul^2+vl^2);er=pr/(gamma+1)+0.5*(ur^2+vr^2);
FUL=[rhol*ul;rhol*ul*ul+pl;rhol*vl*ul;ul*(rhol*el+pl)];
FUR=[rhor*ur;rhor*ur*ur+pr;rhor*vr*ur;ur*(rhor*er+pr)];
tempU=0.5*(FUL+FUR)-0.5*A*(UR-UL);
Flux(m,n,1)=tempU(1);Flux(m,n,2)=tempU(2);Flux(m,n,3)=tempU(3);Flux(m,n,4)=tempU(4);
    end
end
end

