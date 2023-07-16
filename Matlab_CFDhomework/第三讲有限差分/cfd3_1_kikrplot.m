clear
clf
dx=2*pi/20;
for i=0.8:0.4:3.2
    alpha=[0:0.1:3];
    k_i=ki(i,dx,alpha);
    k_r=kr(i,dx,alpha);
    lab=['a4=',num2str(i)];
    figure(1);
    plot(alpha,k_i,'DisplayName',lab);hold on;
    figure(2);
    plot(alpha,k_r,'DisplayName',lab);hold on;
    error=norm(k_i-alpha,2);
    disp(i);
    disp(error);
end
figure(1);
plot(alpha,alpha,'DisplayName','�ο���');hold on;
title('��-ki��ϵͼ');xlabel('��');ylabel('ki');
legend('show');
figure(2);
plot(alpha,alpha,'DisplayName','�ο���');hold on;
title('��-kr��ϵͼ');xlabel('��');ylabel('kr');
legend('show');

