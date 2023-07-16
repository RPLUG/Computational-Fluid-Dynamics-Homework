function [output] = kr(a4,dx,alpha)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
output =-a4*dx/10*cos(3*alpha)-3/2*a4*dx*cos(alpha)+3/5*a4*dx*cos(2*alpha)+a4*dx;

end


