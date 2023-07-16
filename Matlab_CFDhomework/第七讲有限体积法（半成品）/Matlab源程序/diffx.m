function valx = diffx(valep,valet,etx,epx)
%这个函数计算物理量对x的导数
valx=valep.*epx+valet.*etx;
end

