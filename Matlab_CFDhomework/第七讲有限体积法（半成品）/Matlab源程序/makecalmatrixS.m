function cal = makecalmatrixS(M,N,valcell,dummyS,dummyW,dummyE)
%���������������������
%��������ֱ�Ϊ������EWN����������Ķ�Ӧ��
cal=zeros(M,N+1);
cal(1:end-1,2:end-1)=valcell;
cal(end,2:end-1)=dummyS;
cal(2:end,1)=dummyW;
cal(2:end,end)=dummyE;
end

