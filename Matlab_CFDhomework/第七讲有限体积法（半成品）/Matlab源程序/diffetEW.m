function valeta= diffetEW(cal)
%这个函数用于计算W和E边上变量对eta的导数
valeta=0.25*(cal(3:end,1:end-1)+cal(3:end,2:end)-cal(1:end-2,1:end-1)-cal(1:end-2,2:end));
end
