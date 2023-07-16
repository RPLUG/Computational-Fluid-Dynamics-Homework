function valeps= diffepEW(cal)
%这个函数用于计算W和E边上变量对epsilon的导数
valeps=cal(2:end-1,2:end)-cal(2:end-1,1:end-1);
end

