function xyOut = elliptic_mesh(xy) 
    [m, n] = size(xy);%输出矩阵的尺寸，即网格格点数记录 
    convCrit = 1e-5;%收敛判据 
    maxit = 100000;%最大循环次数 
    px = 2:n-1; e = 3:n; w = 1:n-2; 
    py = 2:m-1; n = 3:m; s = 1:m-2; 
    for i = 1:maxit
        %各个点赋初值
        xyE  = xy(py, e); xyN  = xy(n, px); 
        xyW  = xy(py, w); xyS  = xy(s, px); 
        xyNE = xy(n, e ); xySE = xy(s, e ); 
        xyNW = xy(n, w ); xySW = xy(s, w ); 
        %计算导数
        xy_xi = (xyE - xyW)/2; 
        xy_et = (xyN - xyS)/2; 
        xy_xiet = (xyNE - xySE - xyNW + xySW)/4; 
        %计算参数
        a = abs(xy_et).^2; 
        b = real(xy_xi).*real(xy_et) + imag(xy_xi).*imag(xy_et);  
        c = abs(xy_xi).^2; 
        %迭代式
        xyP = (a.*(xyE + xyW) + c.*(xyN + xyS) - 2*b.*(xy_xiet))./(2*(a + c)); 
        % 判断收敛
         err = max(max(abs(xyP - xy(py, px)))); 
         if (err < convCrit) 
             fprintf('----- 网格在循环%i次后收敛\n', i); break; 
         end 
        % 更新网格 
        xy(py, px) = xyP; 
    end 
    xyOut = xy; 
end