function valeps= diffepNS(cal)
%����������ڼ���N��S���ϱ�����epsilon�ĵ���
valeps=0.25*(cal(2:end,3:end)+cal(1:end-1,3:end)-cal(2:end,1:end-2)-cal(1:end-1,1:end-2));
end

