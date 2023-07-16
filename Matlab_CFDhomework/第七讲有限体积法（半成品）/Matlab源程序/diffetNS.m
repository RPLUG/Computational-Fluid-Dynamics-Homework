function valeta= diffetNS(cal)
%这个函数用于计算N和S边上变量对eta的导数
valeta=cal(2:end,2:end-1)-cal(1:end-1,2:end-1);
end
