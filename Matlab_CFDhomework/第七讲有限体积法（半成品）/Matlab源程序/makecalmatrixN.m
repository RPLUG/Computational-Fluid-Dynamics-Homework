function cal = makecalmatrixN(M,N,valcell,dummyN,dummyW,dummyE)
%���������������������
%��������ֱ�Ϊ������EWN����������Ķ�Ӧ��
cal=zeros(M,N+1);
cal(2:end,2:end-1)=valcell;
cal(1,2:end-1)=dummyN;
cal(2:end,1)=dummyW;
cal(2:end,end)=dummyE;
end

