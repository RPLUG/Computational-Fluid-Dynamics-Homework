function [output] = kr(a4,dx,alpha)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
output =-a4*dx/10*cos(3*alpha)-3/2*a4*dx*cos(alpha)+3/5*a4*dx*cos(2*alpha)+a4*dx;

end


