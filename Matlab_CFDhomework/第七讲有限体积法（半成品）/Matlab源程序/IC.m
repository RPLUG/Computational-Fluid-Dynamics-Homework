%% 各个物理量的初始化
%以超音速入口条件为无量纲化的参比值，IC中rho，u，p，T设为1，攻角为0，v设为0
rho=ones(m,n);u=ones(m,n);v=zeros(m,n);T=ones(m,n);p=ones(m,n);U=ptoc(rho,u,v,p);
T_inf=226.5;%无穷远处温度，单位为K，无量纲温度参比值
p_inf=101325;%无穷远处压强，单位为Pa，无量纲压力参比值
R=8.1345;%气体摩尔常数，单位J/(K·mol)
vel_inf=6*sqrt(gamma*R*(T_inf));%来流速度，单位为m/s，无量纲速度参比值

%% 定义原始物理量的虚网格
rho_dummyN=ones(2,n);u_dummyN=ones(2,n);v_dummyN=ones(2,n);p_dummyN=ones(2,n);
rho_dummyS=ones(2,n);u_dummyS=ones(2,n);v_dummyS=ones(2,n);p_dummyS=ones(2,n);
rho_dummyE=ones(m,2);u_dummyE=ones(m,2);v_dummyE=ones(m,2);p_dummyE=ones(m,2);
rho_dummyW=ones(m,2);u_dummyW=ones(m,2);v_dummyW=ones(m,2);p_dummyW=ones(m,2);

