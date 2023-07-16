function cal = makecalmatrixN(M,N,valcell,dummyN,dummyW,dummyE)
%这个函数用来构造计算矩阵
%输入参数分别为内网格，EWN三侧虚网格的对应量
cal=zeros(M,N+1);
cal(2:end,2:end-1)=valcell;
cal(1,2:end-1)=dummyN;
cal(2:end,1)=dummyW;
cal(2:end,end)=dummyE;
end

