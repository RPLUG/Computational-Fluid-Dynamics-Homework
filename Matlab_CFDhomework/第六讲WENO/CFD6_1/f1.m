function [f] = f1(p,pi,rhoi)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
ga=1.4;
ci = sqrt(ga * pi / rhoi);
	if (p > pi)
		f = (p - pi) ./ (rhoi * ci * ((ga + 1) / (2 * ga) * (p ./ pi) + (ga - 1) / (2 * ga)).^0.5);
	else
		f = (2 * ci) / (ga - 1) *((p ./ pi).^((ga - 1)./ (2 * ga)) - 1);
end

