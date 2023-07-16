%% �������ع�
%����ʹ������ӭ���ֵ���Լ��غ��Ա��������ع�����Ҫ����ewsn�ı��ϵ�����ֵ
%�������·�λ��洢�������������һ�£�����ewsn��Ӧ������ֱ����ewsn����
%�������������ڲ���������ϳɴ���󣬱��ھ��󻯼���
U=ptoc2(rho,u,v,p);
rho_ew=[rho_dummyW,rho,rho_dummyE];rho_ns=[rho_dummyN;rho;rho_dummyS];
%ewΪ79*243��nsΪ83*239
u_ew=[u_dummyW,u,u_dummyE];u_ns=[u_dummyN;u;u_dummyS];
v_ew=[v_dummyW,v,v_dummyE];v_ns=[v_dummyN;v;v_dummyS];
p_ew=[p_dummyW,p,p_dummyE];p_ns=[p_dummyN;p;p_dummyS];
U_ew=ptoc2(rho_ew,u_ew,v_ew,p_ew);U_ns=ptoc2(rho_ns,u_ns,v_ns,p_ns);
%���ϣ��õ����ع��ı����õ�U
U_E=(-U_ew(:,4:end,:)+5*U_ew(:,3:end-1,:)+2*U_ew(:,2:end-2,:))/6;F_E=U_E;G_E=U_E;Flux_E=U_E;UtE=U_E;
U_W=(-U_ew(:,1:end-3,:)+5*U_ew(:,2:end-2,:)+2*U_ew(:,3:end-1,:))/6;F_W=U_W;G_W=U_W;Flux_W=U_W;UeW=U_W;
U_S=(-U_ns(4:end,:,:)+5*U_ns(3:end-1,:,:)+2*U_ns(2:end-2,:,:))/6;F_S=U_S;G_S=U_S;Flux_S=U_S;UtS=U_S;
U_N=(-U_ns(1:end-3,:,:)+5*U_ns(2:end-2,:,:)+2*U_ns(3:end-1,:,:))/6;F_N=U_N;G_N=U_N;Flux_N=U_N;UtN=U_N;
%���ϣ��ع��õ����ı������Uֵ
%% ���ع������������ʹ��roe�����õ����������⣬���������¡�����������ѧ������p166�Լ�Toro�������
%% ****************************************������********************************************************
%% ���ȼ�����ߵ�λ�ⷨ��nx��ny
nxE=ones(M-1,N);nyE=ones(M-1,N);
nxW=ones(M-1,N);nyW=ones(M-1,N);
%��Ҫע����ǣ�֮ǰ�ķ����������������Ӧ��79*239�����ڼ�����תroeƽ��ֵ��ʱ����Ҫ���������80*240
for m=1:M-1
    for n=1:N-1
    %���ڹ�ʽut=u*nx+v*ny,vt=-u*nx+v*ny
   nxE(m,n+1)=real(nE(m,n));
   nyE(m,n+1)=imag(nE(m,n));
   nxW(m,n)=real(nW(m,n));
   nyW(m,n)=imag(nW(m,n));
    end
end
%���⴦����e���ⷨ�򣬵�һ��ӦΪw���ⷨ���һ�е��෴��,w���ⷨ�����һ��ӦΪe���ⷨ�����һ�е��෴��
for m=1:M-1
    nxE(m,1)=-real(nW(m,1));
    nyE(m,1)=-imag(nW(m,1));
    nxW(m,end)=-real(nE(m,end));
    nyW(m,end)=-imag(nE(m,end));
end
%% ��ԭʼ������Ϊ�ֲ�����ϵ�µı���������һά����
utE=U_E(:,:,2).*nxE+U_E(:,:,3).*nyE;vtE=-U_E(:,:,2).*nxE+U_E(:,:,3).*nyE;
utW=U_W(:,:,2).*nxE+U_W(:,:,3).*nyE;vtW=-U_W(:,:,2).*nxE+U_W(:,:,3).*nyE;
rhotW=U_W(:,:,1);etW=U_W(:,:,4)./U_W(:,:,1);%�����eΪ���������ܶȣ���ϵΪE=rho*e
ptW=(gamma - 1) * (etW - 0.5 * (utW .*utW + vtW.*vtW)).*rhotW;
rhotE=U_E(:,:,1);etE=U_E(:,:,4)./U_E(:,:,1);
ptE=(gamma - 1) * (etE - 0.5 * (utE .*utE + vtE.*vtE)).*rhotE;
UtE=ptoc(rhotE,utE,vtE,ptE);
UtW=ptoc(rhotW,utE,vtE,ptW);
FluxEW=RoeFlux(UtE,UtW,M-1,N,nxE,nxW,nyE,nyW);%�����򣬹�����ͨ����79*240
Flux_E=FluxEW;Flux_W=FluxEW;

Flux_E(:,:,2)=Flux_E(:,:,2).*nxE-Flux_E(:,:,2).*nyE;
Flux_E(:,:,3)=Flux_E(:,:,3).*nyE-Flux_E(:,:,3).*nxE;

Flux_W(:,:,2)=Flux_W(:,:,2).*nxE-Flux_W(:,:,2).*nyE;
Flux_W(:,:,3)=Flux_W(:,:,3).*nyE+Flux_W(:,:,3).*nxE;
 Flux_E(:,1,:)=[];Flux_W(:,end,:)=[];
%% **************************************�ϱ���******************************************************
%% ���ȼ�����ߵ�λ�ⷨ��nx��ny
nxN=ones(M,N-1);nyN=ones(M,N-1);
nxS=ones(M,N-1);nyS=ones(M,N-1);
%��Ҫע����ǣ�֮ǰ�ķ����������������Ӧ��79*239�����ڼ�����תroeƽ��ֵ��ʱ����Ҫ���������80*240
for m=1:M-1
    for n=1:N-1
    %���ڹ�ʽut=u*nx+v*ny,vt=-u*nx+v*ny
   nxN(m,n)=real(nN(m,n));
   nyN(m,n)=imag(nN(m,n));
   nxS(m+1,n)=real(nS(m,n));
   nyS(m+1,n)=imag(nS(m,n));
    end
end
%���⴦����s���ⷨ�򣬵�һ��ӦΪn���ⷨ���1�е��෴��,n���ⷨ�����һ��ӦΪs���ⷨ�����һ�е��෴��
for n=1:N-1
    nxN(end,n)=-real(nS(end,n));
    nyN(end,n)=-imag(nS(end,n));
    nxS(1,n)=-real(nN(1,n));
    nyS(1,n)=-imag(nN(1,n));
end
%% ��ԭʼ������Ϊ�ֲ�����ϵ�µı���������һά����
utS=U_S(:,:,2).*nxS+U_S(:,:,3).*nyS;vtS=-U_S(:,:,2).*nxS+U_S(:,:,3).*nyS;
utN=U_N(:,:,2).*nxS+U_N(:,:,3).*nyS;vtN=-U_N(:,:,2).*nxS+U_N(:,:,3).*nyS;
rhotN=U_N(:,:,1);etN=U_N(:,:,4)./U_N(:,:,1);
ptN=(gamma - 1) * (etN - 0.5 * (utN .*utN + vtN.*vtN));
rhotS=U_S(:,:,1);etS=U_S(:,:,4)./U_S(:,:,1);
ptS=(gamma - 1) * (etS - 0.5 * (utS .*utS + vtS.*vtS));
UtS=ptoc(rhotS,utS,vtS,ptS);
UtN=ptoc(rhotN,utS,vtS,ptN);
FluxSN=RoeFlux(UtS,UtN,M,N-1,nxS,nxN,nyS,nyN);%�����򣬹�����ͨ����79*240
Flux_S=FluxSN;Flux_N=FluxSN;
Flux_S(:,:,2)=Flux_S(:,:,2).*nxS-Flux_S(:,:,2).*nyS;
Flux_S(:,:,3)=Flux_S(:,:,3).*nyS+Flux_S(:,:,3).*nxS;
Flux_N(:,:,2)=Flux_N(:,:,2).*nxS-Flux_N(:,:,2).*nyS;
Flux_N(:,:,3)=Flux_N(:,:,3).*nyS+Flux_N(:,:,3).*nxS;
 Flux_S(1,:,:)=[];Flux_N(end,:,:)=[];
