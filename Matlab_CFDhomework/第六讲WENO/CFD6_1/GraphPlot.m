%��ͼ
[rho_numerical,u_numerical,p_numerical] = conservation_to_physics(U);
EulerExact;
 plot(x,rho,'m','DisplayName','��ȷ���ܶ�');hold on;
 plot(x,u,'k','DisplayName','��ȷ���ٶ�');hold on;
 plot(x,p,'r','DisplayName','��ȷ��ѹǿ');hold on; 
 plot(x,rho_numerical,'m+','DisplayName','��ֵ���ܶ�');hold on;
 plot(x,u_numerical,'k+','DisplayName','��ֵ���ٶ�');hold on;
 plot(x,p_numerical,'r+','DisplayName','��ֵ��ѹǿ');hold on;
 legend('show');
 title('��ֵ���뾫ȷ��Ա�ͼ');xlabel('x');ylabel('��������ֵ');

