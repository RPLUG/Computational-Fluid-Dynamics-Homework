%% �����������ĳ�ʼ��
%�Գ������������Ϊ�����ٻ��Ĳα�ֵ��IC��rho��u��p��T��Ϊ1������Ϊ0��v��Ϊ0
rho=ones(m,n);u=ones(m,n);v=zeros(m,n);T=ones(m,n);p=ones(m,n);U=ptoc(rho,u,v,p);
T_inf=226.5;%����Զ���¶ȣ���λΪK���������¶Ȳα�ֵ
p_inf=101325;%����Զ��ѹǿ����λΪPa��������ѹ���α�ֵ
R=8.1345;%����Ħ����������λJ/(K��mol)
vel_inf=6*sqrt(gamma*R*(T_inf));%�����ٶȣ���λΪm/s���������ٶȲα�ֵ

%% ����ԭʼ��������������
rho_dummyN=ones(2,n);u_dummyN=ones(2,n);v_dummyN=ones(2,n);p_dummyN=ones(2,n);
rho_dummyS=ones(2,n);u_dummyS=ones(2,n);v_dummyS=ones(2,n);p_dummyS=ones(2,n);
rho_dummyE=ones(m,2);u_dummyE=ones(m,2);v_dummyE=ones(m,2);p_dummyE=ones(m,2);
rho_dummyW=ones(m,2);u_dummyW=ones(m,2);v_dummyW=ones(m,2);p_dummyW=ones(m,2);

