function cal = makecalmatrixW(M,N,valcell,dummyN,dummyS,dummyW)
%���������������������
%��������ֱ�Ϊ������N��S��W����������Ķ�Ӧ��
cal=zeros(M+1,N);
cal(2:end-1,2:end)=valcell;
cal(1,1:end-1)=dummyN;
cal(end,1:end-1)=dummyS;
cal(2:end-1,1)=dummyW;
end

