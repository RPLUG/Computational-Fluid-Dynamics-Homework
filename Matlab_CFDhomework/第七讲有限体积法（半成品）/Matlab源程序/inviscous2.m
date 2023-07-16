%% 物理量重构
%这里使用三阶迎风插值，以及守恒性变量进行重构，需要计算ewsn四边上的两侧值
%左右上下方位与存储矩阵的左右上下一致，并与ewsn对应，这里直接用ewsn代替
%将虚网格量和内部网格量组合成大矩阵，便于矩阵化计算
U=ptoc2(rho,u,v,p);
rho_ew=[rho_dummyW,rho,rho_dummyE];rho_ns=[rho_dummyN;rho;rho_dummyS];
%ew为79*243，ns为83*239
u_ew=[u_dummyW,u,u_dummyE];u_ns=[u_dummyN;u;u_dummyS];
v_ew=[v_dummyW,v,v_dummyE];v_ns=[v_dummyN;v;v_dummyS];
p_ew=[p_dummyW,p,p_dummyE];p_ns=[p_dummyN;p;p_dummyS];
U_ew=ptoc2(rho_ew,u_ew,v_ew,p_ew);U_ns=ptoc2(rho_ns,u_ns,v_ns,p_ns);
%以上，得到了重构四边量用的U
U_E=(-U_ew(:,4:end,:)+5*U_ew(:,3:end-1,:)+2*U_ew(:,2:end-2,:))/6;F_E=U_E;G_E=U_E;Flux_E=U_E;UtE=U_E;
U_W=(-U_ew(:,1:end-3,:)+5*U_ew(:,2:end-2,:)+2*U_ew(:,3:end-1,:))/6;F_W=U_W;G_W=U_W;Flux_W=U_W;UeW=U_W;
U_S=(-U_ns(4:end,:,:)+5*U_ns(3:end-1,:,:)+2*U_ns(2:end-2,:,:))/6;F_S=U_S;G_S=U_S;Flux_S=U_S;UtS=U_S;
U_N=(-U_ns(1:end-3,:,:)+5*U_ns(2:end-2,:,:)+2*U_ns(3:end-1,:,:))/6;F_N=U_N;G_N=U_N;Flux_N=U_N;UtN=U_N;
%以上，重构得到了四边两侧的U值
%% 对重构后的物理量，使用roe方法得到近似黎曼解，参照任玉新《计算流体力学基础》p166以及Toro的书计算
%% ****************************************东西向********************************************************
%% 首先计算各边单位外法向nx，ny
nxE=ones(M-1,N);nyE=ones(M-1,N);
nxW=ones(M-1,N);nyW=ones(M-1,N);
%需要注意的是，之前的法向量与内网格相对应，79*239，这在计算旋转roe平均值的时候需要额外扩充成80*240
for m=1:M-1
    for n=1:N-1
    %基于公式ut=u*nx+v*ny,vt=-u*nx+v*ny
   nxE(m,n+1)=real(nE(m,n));
   nyE(m,n+1)=imag(nE(m,n));
   nxW(m,n)=real(nW(m,n));
   nyW(m,n)=imag(nW(m,n));
    end
end
%额外处理，即e侧外法向，第一列应为w侧外法向第一列的相反数,w侧外法向，最后一列应为e侧外法向最后一列的相反数
for m=1:M-1
    nxE(m,1)=-real(nW(m,1));
    nyE(m,1)=-imag(nW(m,1));
    nxW(m,end)=-real(nE(m,end));
    nyW(m,end)=-imag(nE(m,end));
end
%% 将原始变量变为局部坐标系下的变量，解拟一维问题
utE=U_E(:,:,2).*nxE+U_E(:,:,3).*nyE;vtE=-U_E(:,:,2).*nxE+U_E(:,:,3).*nyE;
utW=U_W(:,:,2).*nxE+U_W(:,:,3).*nyE;vtW=-U_W(:,:,2).*nxE+U_W(:,:,3).*nyE;
rhotW=U_W(:,:,1);etW=U_W(:,:,4)./U_W(:,:,1);%这里的e为质量能量密度，关系为E=rho*e
ptW=(gamma - 1) * (etW - 0.5 * (utW .*utW + vtW.*vtW)).*rhotW;
rhotE=U_E(:,:,1);etE=U_E(:,:,4)./U_E(:,:,1);
ptE=(gamma - 1) * (etE - 0.5 * (utE .*utE + vtE.*vtE)).*rhotE;
UtE=ptoc(rhotE,utE,vtE,ptE);
UtW=ptoc(rhotW,utE,vtE,ptW);
FluxEW=RoeFlux(UtE,UtW,M-1,N,nxE,nxW,nyE,nyW);%东西向，公共边通量，79*240
Flux_E=FluxEW;Flux_W=FluxEW;

Flux_E(:,:,2)=Flux_E(:,:,2).*nxE-Flux_E(:,:,2).*nyE;
Flux_E(:,:,3)=Flux_E(:,:,3).*nyE-Flux_E(:,:,3).*nxE;

Flux_W(:,:,2)=Flux_W(:,:,2).*nxE-Flux_W(:,:,2).*nyE;
Flux_W(:,:,3)=Flux_W(:,:,3).*nyE+Flux_W(:,:,3).*nxE;
 Flux_E(:,1,:)=[];Flux_W(:,end,:)=[];
%% **************************************南北向******************************************************
%% 首先计算各边单位外法向nx，ny
nxN=ones(M,N-1);nyN=ones(M,N-1);
nxS=ones(M,N-1);nyS=ones(M,N-1);
%需要注意的是，之前的法向量与内网格相对应，79*239，这在计算旋转roe平均值的时候需要额外扩充成80*240
for m=1:M-1
    for n=1:N-1
    %基于公式ut=u*nx+v*ny,vt=-u*nx+v*ny
   nxN(m,n)=real(nN(m,n));
   nyN(m,n)=imag(nN(m,n));
   nxS(m+1,n)=real(nS(m,n));
   nyS(m+1,n)=imag(nS(m,n));
    end
end
%额外处理，即s侧外法向，第一行应为n侧外法向第1行的相反数,n侧外法向，最后一行应为s侧外法向最后一行的相反数
for n=1:N-1
    nxN(end,n)=-real(nS(end,n));
    nyN(end,n)=-imag(nS(end,n));
    nxS(1,n)=-real(nN(1,n));
    nyS(1,n)=-imag(nN(1,n));
end
%% 将原始变量变为局部坐标系下的变量，解拟一维问题
utS=U_S(:,:,2).*nxS+U_S(:,:,3).*nyS;vtS=-U_S(:,:,2).*nxS+U_S(:,:,3).*nyS;
utN=U_N(:,:,2).*nxS+U_N(:,:,3).*nyS;vtN=-U_N(:,:,2).*nxS+U_N(:,:,3).*nyS;
rhotN=U_N(:,:,1);etN=U_N(:,:,4)./U_N(:,:,1);
ptN=(gamma - 1) * (etN - 0.5 * (utN .*utN + vtN.*vtN));
rhotS=U_S(:,:,1);etS=U_S(:,:,4)./U_S(:,:,1);
ptS=(gamma - 1) * (etS - 0.5 * (utS .*utS + vtS.*vtS));
UtS=ptoc(rhotS,utS,vtS,ptS);
UtN=ptoc(rhotN,utS,vtS,ptN);
FluxSN=RoeFlux(UtS,UtN,M,N-1,nxS,nxN,nyS,nyN);%东西向，公共边通量，79*240
Flux_S=FluxSN;Flux_N=FluxSN;
Flux_S(:,:,2)=Flux_S(:,:,2).*nxS-Flux_S(:,:,2).*nyS;
Flux_S(:,:,3)=Flux_S(:,:,3).*nyS+Flux_S(:,:,3).*nxS;
Flux_N(:,:,2)=Flux_N(:,:,2).*nxS-Flux_N(:,:,2).*nyS;
Flux_N(:,:,3)=Flux_N(:,:,3).*nyS+Flux_N(:,:,3).*nxS;
 Flux_S(1,:,:)=[];Flux_N(end,:,:)=[];
