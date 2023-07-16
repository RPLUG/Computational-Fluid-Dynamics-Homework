function xyOut = elliptic_mesh(xy) 
    [m, n] = size(xy);%�������ĳߴ磬������������¼ 
    convCrit = 1e-5;%�����о� 
    maxit = 100000;%���ѭ������ 
    px = 2:n-1; e = 3:n; w = 1:n-2; 
    py = 2:m-1; n = 3:m; s = 1:m-2; 
    for i = 1:maxit
        %�����㸳��ֵ
        xyE  = xy(py, e); xyN  = xy(n, px); 
        xyW  = xy(py, w); xyS  = xy(s, px); 
        xyNE = xy(n, e ); xySE = xy(s, e ); 
        xyNW = xy(n, w ); xySW = xy(s, w ); 
        %���㵼��
        xy_xi = (xyE - xyW)/2; 
        xy_et = (xyN - xyS)/2; 
        xy_xiet = (xyNE - xySE - xyNW + xySW)/4; 
        %�������
        a = abs(xy_et).^2; 
        b = real(xy_xi).*real(xy_et) + imag(xy_xi).*imag(xy_et);  
        c = abs(xy_xi).^2; 
        %����ʽ
        xyP = (a.*(xyE + xyW) + c.*(xyN + xyS) - 2*b.*(xy_xiet))./(2*(a + c)); 
        % �ж�����
         err = max(max(abs(xyP - xy(py, px)))); 
         if (err < convCrit) 
             fprintf('----- ������ѭ��%i�κ�����\n', i); break; 
         end 
        % �������� 
        xy(py, px) = xyP; 
    end 
    xyOut = xy; 
end