function valeta= diffetEW(cal)
%����������ڼ���W��E���ϱ�����eta�ĵ���
valeta=0.25*(cal(3:end,1:end-1)+cal(3:end,2:end)-cal(1:end-2,1:end-1)-cal(1:end-2,2:end));
end
