function cal = makecalmatrixW(M,N,valcell,dummyN,dummyS,dummyW)
%这个函数用来构造计算矩阵
%输入参数分别为内网格，N、S、W三侧虚网格的对应量
cal=zeros(M+1,N);
cal(2:end-1,2:end)=valcell;
cal(1,1:end-1)=dummyN;
cal(end,1:end-1)=dummyS;
cal(2:end-1,1)=dummyW;
end

