function [output] = ki(a4,dx,alpha)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
output =a4*dx/10*sin(3*alpha)+(1/2*a4*dx+4/3)*sin(alpha)-(2/5*a4*dx+1/6)*sin(2*alpha);

end

