function valeps= diffepNS(cal)
%这个函数用于计算N和S边上变量对epsilon的导数
valeps=0.25*(cal(2:end,3:end)+cal(1:end-1,3:end)-cal(2:end,1:end-2)-cal(1:end-1,1:end-2));
end

