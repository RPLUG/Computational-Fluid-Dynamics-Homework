%为计算粘性通量，需要计算u,v,T三个量各自在x和y向上的导数，使用雅可比变换到计算空间ε-η中
%为此，需要计算x,y对伊普西隆和伊塔的导数，以上均与各边对应，这里用下标EWNS表示
%与矩阵存储的形式对应，epsilon正方向为矩阵列号增大方向(E)，eta正方向为矩阵行号增大的方向(S)
%% 以下是每一步需要更新的物理量，故放在这里定义
vel=sqrt(u.^2+v.^2);%无量纲速度大小
C=110.4;
miu=(T./T_inf).^1.5.*(1+C)./(T+C);%计算无量纲粘度
a=sqrt(gamma*R*(T*T_inf));%声速
Ma=vel*vel_inf./a;%马赫数
Cp=1./((gamma-1)*Ma.^2);%气体无量纲热容
Pr=0.8;%气体普朗特数，认为是个常数
Re=10000;%流动雷诺数
%% ****************  E边  ************************************
%计算x，y对epsilon，eta的导数，E边共有79*239条，存储其的矩阵规模也为此
%对E边而言，计算上述导数需要N、S、E三侧的虚网格以及角点，这里把角点的值直接做0处理，简化计算
%则存储计算信息，需要的矩阵规模为81*240，这里先构造一个整体的矩阵，再用通用的向量化编程方式计算导数
%定义计算以及存储用矩阵
xcal=makecalmatrixE(M,N,real(xyCell),real(xyCell_dummyN),real(xyCell_dummyS),real(xyCell_dummyE));
ycal=makecalmatrixE(M,N,imag(xyCell),imag(xyCell_dummyN),imag(xyCell_dummyS),imag(xyCell_dummyE));
ucal=makecalmatrixE(M,N,u,u_dummyN(2,:),u_dummyS(1,:),u_dummyE(:,1));
vcal=makecalmatrixE(M,N,v,v_dummyN(2,:),v_dummyS(1,:),v_dummyE(:,1));
Tcal=makecalmatrixE(M,N,T,T_dummyN,T_dummyS,T_dummyE);
%计算各物理量对epsion和eta的导数
xepE=diffepEW(xcal);xetE=diffetEW(xcal);
yepE=diffepEW(ycal);yetE=diffetEW(ycal);
uepE=diffepEW(ucal);uetE=diffetEW(ucal);
vepE=diffepEW(vcal);vetE=diffetEW(vcal);
TepE=diffepEW(Tcal);TetE=diffetEW(Tcal);
JE=1./(xepE.*yetE-yepE.*xetE);%雅可比矩阵值
epxE=JE.*yetE;epyE=-JE.*xetE;etxE=-JE.*yepE;etyE=JE.*xepE;
%计算各物理量对x，y的导数
uxE=diffx(uepE,uetE,etxE,epxE);
vxE=diffx(vepE,vetE,etxE,epxE);
TxE=diffx(TepE,TetE,etxE,epxE);
uyE=diffy(uepE,uetE,etyE,epyE);
vyE=diffy(vepE,vetE,etyE,epyE);
TyE=diffy(TepE,TetE,etyE,epyE);
%计算tau
tau11E=2/3*miu.*(2*uxE-vyE);tau22E=2/3*miu.*(2*vyE-uxE);tau12E=miu.*(uyE+vxE);
%计算E边粘性通量，定义F1E，F2E以及FluxVE,大小均为79*239*4
F1E=zeros(M-1,N-1,4);F2E=zeros(M-1,N-1,4);
F1E(:,:,2)=tau11E;F2E(:,:,2)=tau12E;
F1E(:,:,3)=tau12E;F2E(:,:,3)=tau22E;
F1E(:,:,4)=Cp.*miu/Pr/Re.*TxE+u.*tau11E+v.*tau12E;
F2E(:,:,4)=Cp.*miu/Pr/Re.*TyE+u.*tau12E+v.*tau22E;
FluxVE(:,:,1)=F1E(:,:,1).*real(nE)+F2E(:,:,1).*imag(nE);
FluxVE(:,:,2)=F1E(:,:,2).*real(nE)+F2E(:,:,2).*imag(nE);
FluxVE(:,:,3)=F1E(:,:,3).*real(nE)+F2E(:,:,3).*imag(nE);
FluxVE(:,:,4)=F1E(:,:,4).*real(nE)+F2E(:,:,4).*imag(nE);
%% ************************* W边 ********************************
%计算x，y对epsilon，eta的导数，W边共有79*239条，存储其的矩阵规模也为此
%对W边而言，计算上述导数需要N、S、W三侧的虚网格以及角点，这里把角点的值直接做0处理，简化计算
%则存储计算信息，需要的矩阵规模为81*240，这里先构造一个整体的矩阵，再用通用的向量化编程方式计算导数
%定义计算以及存储用矩阵
xcal=makecalmatrixW(M,N,real(xyCell),real(xyCell_dummyN),real(xyCell_dummyS),real(xyCell_dummyW));
ycal=makecalmatrixW(M,N,imag(xyCell),imag(xyCell_dummyN),imag(xyCell_dummyS),imag(xyCell_dummyW));
ucal=makecalmatrixW(M,N,u,u_dummyN(2,:),u_dummyS(1,:),u_dummyW(:,2));
vcal=makecalmatrixW(M,N,v,v_dummyN(2,:),v_dummyS(1,:),v_dummyW(:,2));
Tcal=makecalmatrixW(M,N,T,T_dummyN,T_dummyS,T_dummyW);
%计算各物理量对epsion和eta的导数
xepW=diffepEW(xcal);xetW=diffetEW(xcal);
yepW=diffepEW(ycal);yetW=diffetEW(ycal);
uepW=diffepEW(ucal);uetW=diffetEW(ucal);
vepW=diffepEW(vcal);vetW=diffetEW(vcal);
TepW=diffepEW(Tcal);TetW=diffetEW(Tcal);
JW=1./(xepW.*yetW-yepW.*xetW);%雅可比矩阵值
epxW=JW.*yetW;epyW=-JW.*xetW;etxW=-JW.*yepW;etyW=JW.*xepW;
%计算各物理量对x，y的导数
uxW=diffx(uepW,uetW,etxW,epxW);
vxW=diffx(vepW,vetW,etxW,epxW);
TxW=diffx(TepW,TetW,etxW,epxW);
uyW=diffy(uepW,uetW,etyW,epyW);
vyW=diffy(vepW,vetW,etyW,epyW);
TyW=diffy(TepW,TetW,etyW,epyW);
%计算tau
tau11W=2/3*miu.*(2*uxW-vyW);tau22W=2/3*miu.*(2*vyW-uxW);tau12W=miu.*(uyW+vxW);
%计算W边粘性通量，定义F1W，F2W以及FluxVW,大小均为79*239*4
F1W=zeros(M-1,N-1,4);F2W=zeros(M-1,N-1,4);
F1W(:,:,2)=tau11W;F2W(:,:,2)=tau12W;
F1W(:,:,3)=tau12W;F2W(:,:,3)=tau22W;
F1W(:,:,4)=Cp.*miu/Pr/Re.*TxW+u.*tau11W+v.*tau12W;
F2W(:,:,4)=Cp.*miu/Pr/Re.*TyW+u.*tau12W+v.*tau22W;
FluxVW(:,:,1)=F1W(:,:,1).*real(nW)+F2W(:,:,1).*imag(nW);
FluxVW(:,:,2)=F1W(:,:,2).*real(nW)+F2W(:,:,2).*imag(nW);
FluxVW(:,:,3)=F1W(:,:,3).*real(nW)+F2W(:,:,3).*imag(nW);
FluxVW(:,:,4)=F1W(:,:,4).*real(nW)+F2W(:,:,4).*imag(nW);
%% ************************* N边 ********************************
%计算x，y对epsilon，eta的导数，N边共有79*239条，存储其的矩阵规模也为此
%对N边而言，计算上述导数需W、E、N三侧的虚网格以及角点，这里把角点的值直接做0处理，简化计算
%则存储计算信息，需要的矩阵规模为79*241，这里先构造一个整体的矩阵，再用通用的向量化编程方式计算导数
%定义计算以及存储用矩阵
xcal=makecalmatrixN(M,N,real(xyCell),real(xyCell_dummyN),real(xyCell_dummyW),real(xyCell_dummyE));
ycal=makecalmatrixN(M,N,imag(xyCell),imag(xyCell_dummyN),imag(xyCell_dummyW),imag(xyCell_dummyE));
ucal=makecalmatrixN(M,N,u,u_dummyN(2,:),u_dummyW(:,2),u_dummyE(:,1));
vcal=makecalmatrixN(M,N,v,v_dummyN(2,:),v_dummyW(:,2),v_dummyE(:,1));
Tcal=makecalmatrixN(M,N,T,T_dummyN,T_dummyW,T_dummyE);
%计算各物理量对epsion和eta的导数
xepN=diffepNS(xcal);xetN=diffetNS(xcal);
yepN=diffepNS(ycal);yetN=diffetNS(ycal);
uepN=diffepNS(ucal);uetN=diffetNS(ucal);
vepN=diffepNS(vcal);vetN=diffetNS(vcal);
TepN=diffepNS(Tcal);TetN=diffetNS(Tcal);
JN=1./(xepN.*yetN-yepN.*xetN);%雅可比矩阵值
epxN=JN.*yetN;epyN=-JN.*xetN;etxN=-JN.*yepN;etyN=JN.*xepN;
%计算各物理量对x，y的导数
uxN=diffx(uepN,uetN,etxN,epxN);
vxN=diffx(vepN,vetN,etxN,epxN);
TxN=diffx(TepN,TetN,etxN,epxN);
uyN=diffy(uepN,uetN,etyN,epyN);
vyN=diffy(vepN,vetN,etyN,epyN);
TyN=diffy(TepN,TetN,etyN,epyN);
%计算tau
tau11N=2/3*miu.*(2*uxN-vyN);tau22N=2/3*miu.*(2*vyN-uxN);tau12N=miu.*(uyN+vxN);
%计算N边粘性通量，定义F1N，F2N以及FluxVN,大小均为79*239*4
F1N=zeros(M-1,N-1,4);F2N=zeros(M-1,N-1,4);
F1N(:,:,2)=tau11N;F2N(:,:,2)=tau12N;
F1N(:,:,3)=tau12N;F2N(:,:,3)=tau22N;
F1N(:,:,4)=Cp.*miu/Pr/Re.*TxN+u.*tau11N+v.*tau12N;
F2N(:,:,4)=Cp.*miu/Pr/Re.*TyN+u.*tau12N+v.*tau22N;
FluxVN(:,:,1)=F1N(:,:,1).*real(nN)+F2N(:,:,1).*imag(nN);
FluxVN(:,:,2)=F1N(:,:,2).*real(nN)+F2N(:,:,2).*imag(nN);
FluxVN(:,:,3)=F1N(:,:,3).*real(nN)+F2N(:,:,3).*imag(nN);
FluxVN(:,:,4)=F1N(:,:,4).*real(nN)+F2N(:,:,4).*imag(nN);
%% ************************* S边 ********************************
%计算x，y对epsilon，eta的导数，S边共有79*239条，存储其的矩阵规模也为此
%对N边而言，计算上述导数需W、E、S三侧的虚网格以及角点，这里把角点的值直接做0处理，简化计算
%则存储计算信息，需要的矩阵规模为79*241，这里先构造一个整体的矩阵，再用通用的向量化编程方式计算导数
%定义计算以及存储用矩阵
xcal=makecalmatrixS(M,N,real(xyCell),real(xyCell_dummyS),real(xyCell_dummyW),real(xyCell_dummyE));
ycal=makecalmatrixS(M,N,imag(xyCell),imag(xyCell_dummyS),imag(xyCell_dummyW),imag(xyCell_dummyE));
ucal=makecalmatrixS(M,N,u,u_dummyS(2,:),u_dummyW(:,2),u_dummyE(:,1));
vcal=makecalmatrixS(M,N,v,v_dummyS(2,:),v_dummyW(:,2),v_dummyE(:,1));
Tcal=makecalmatrixS(M,N,T,T_dummyS,T_dummyW,T_dummyE);
%计算各物理量对epsion和eta的导数
xepS=diffepNS(xcal);xetS=diffetNS(xcal);
yepS=diffepNS(ycal);yetS=diffetNS(ycal);
uepS=diffepNS(ucal);uetS=diffetNS(ucal);
vepS=diffepNS(vcal);vetS=diffetNS(vcal);
TepS=diffepNS(Tcal);TetS=diffetNS(Tcal);
JS=1./(xepS.*yetS-yepS.*xetS);%雅可比矩阵值
epxS=JS.*yetS;epyS=-JS.*xetS;etxS=-JS.*yepS;etyS=JS.*xepS;
%计算各物理量对x，y的导数
uxS=diffx(uepS,uetS,etxS,epxS);
vxS=diffx(vepS,vetS,etxS,epxS);
TxS=diffx(TepS,TetS,etxS,epxS);
uyS=diffy(uepS,uetS,etyS,epyS);
vyS=diffy(vepS,vetS,etyS,epyS);
TyS=diffy(TepS,TetS,etyS,epyS);
%计算tau
tau11S=2/3*miu.*(2*uxS-vyS);tau22S=2/3*miu.*(2*vyS-uxS);tau12S=miu.*(uyS+vxS);
%计算S边粘性通量，定义F1S，F2S以及FluxVS,大小均为79*239*4
F1S=zeros(M-1,N-1,4);F2S=zeros(M-1,N-1,4);
F1S(:,:,2)=tau11S;F2S(:,:,2)=tau12S;
F1S(:,:,3)=tau12S;F2S(:,:,3)=tau22S;
F1S(:,:,4)=Cp.*miu/Pr/Re.*TxS+u.*tau11S+v.*tau12S;
F2S(:,:,4)=Cp.*miu/Pr/Re.*TyS+u.*tau12S+v.*tau22S;
FluxVS(:,:,1)=F1S(:,:,1).*real(nS)+F2S(:,:,1).*imag(nS);
FluxVS(:,:,2)=F1S(:,:,2).*real(nS)+F2S(:,:,2).*imag(nS);
FluxVS(:,:,3)=F1S(:,:,3).*real(nS)+F2S(:,:,3).*imag(nS);
FluxVS(:,:,4)=F1S(:,:,4).*real(nS)+F2S(:,:,4).*imag(nS);