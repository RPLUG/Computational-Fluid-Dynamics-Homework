Vxy=importdata('mesh.txt');%输入节点坐标存入向量Vxy
M=80;%行数
N=240;%列数
%排成80*240矩阵
for m=1:M
Rowhead=1+N*(m-1);
Rowtail=N+N*(m-1);
xy(m,:)=Vxy(Rowhead:Rowtail,1)+1i*Vxy(Rowhead:Rowtail,2);
end
clear Rowhead;clear Rowtail;clear Vxy;
%创建网格中心点，即四边形中心，通过分割出的两个三角形面积加权平均得到
for m=1:M-1
    for n=1:N-1
        S1=-(0.5*((-real(xy(m,n))+real(xy(m+1,n)))*(imag(xy(m,n+1))-imag(xy(m,n))))-...
            ((-imag(xy(m,n))+imag(xy(m+1,n)))*(real(xy(m,n+1))-real(xy(m,n)))));
        S2= 0.5*((-real(xy(m+1,n+1))+real(xy(m+1,n)))*(imag(xy(m,n+1))-imag(xy(m+1,n+1))))-...
            ((-imag(xy(m+1,n+1))+imag(xy(m+1,n)))*(real(xy(m,n+1))-real(xy(m+1,n+1))));
        c1=1/3*(xy(m,n)+xy(m+1,n)+xy(m,n+1));
        c2=1/3*(xy(m+1,n+1)+xy(m+1,n)+xy(m,n+1));
        xyCell(m,n)=(S1*c1+S2*c2)/(S1+S2);
        Scell(m,n)=S1+S2;
    end
end
%最终得到的网格节点坐标矩阵为xy(80,240)，网格中心点坐标矩阵为xyCell(79,239)
%下标相同时，在矩阵空间里，xy对应xyCell的“左上方”
%EdgeEWNS以及VnESWN与各个网格四边对应，也为(79,239)
%创建与网格边对应的单位法向量，为了确保是“外法向”，边向量由网格边界点循环相减，归一化后得到
for m=1:M-1
    for n=1:N-1
EdgeE(m,n)=xy(m,n+1)-xy(m+1,n+1);EdgeN(m,n)=xy(m,n)-xy(m,n+1);
EdgeW(m,n)=xy(m+1,n)-xy(m,n);EdgeS(m,n)=xy(m+1,n+1)-xy(m+1,n);


   nE(m,n)=-(imag(EdgeE(m,n))-real(EdgeE(m,n))*1i)/(abs(EdgeE(m,n)));
   nS(m,n)=-(imag(EdgeS(m,n))-real(EdgeS(m,n))*1i)/(abs(EdgeS(m,n)));
   nW(m,n)=-(imag(EdgeW(m,n))-real(EdgeW(m,n))*1i)/(abs(EdgeW(m,n)));
   nN(m,n)=-(imag(EdgeN(m,n))-real(EdgeN(m,n))*1i)/(abs(EdgeN(m,n)));
    end
end
%虚网格变量命名方式统一为X_dummy*，*为方向
%这里给出nwes方向的坐标位置虚网格，用于粘性项的相关离散公式
%计算方式为靠近虚网格的第一层格点中心坐标减去靠近虚网格的第二层格点中心坐标，直接将得到的向量加在第一层格点中心坐标上
xyCell_dummyN=2*xyCell(1,:)-xyCell(2,:);
xyCell_dummyS=2*xyCell(end,:)-xyCell(end-1,:);
xyCell_dummyW=2*xyCell(:,1)-xyCell(:,2);
xyCell_dummyE=2*xyCell(:,end)-xyCell(:,end-1);


