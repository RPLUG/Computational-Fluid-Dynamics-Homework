function cal = makecalmatrixE(M,N,valcell,dummyN,dummyS,dummyE)
%���������������������
%��������ֱ�Ϊ������N��S��E����������Ķ�Ӧ��
cal=zeros(M+1,N);
cal(2:end-1,1:end-1)=valcell;
cal(1,1:end-1)=dummyN;
cal(end,1:end-1)=dummyS;
cal(2:end-1,end)=dummyE;
end

