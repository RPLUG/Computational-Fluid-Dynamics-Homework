function dF = FluxDiffRoe(U,dx)
%Roe格式的离散
global gamma;
global N;
%% step 1 计算U_R和U_L
%首先计算左右迎风插值，得到N+1个点
U1 = U(:,1:end-4);%j-2
U2 = U(:,2:end-3);%j-1
U3 = U(:,3:end-2);%j
U4 = U(:,4:end-1);%j+1
U5 = U(:,5:end);%j+2
%weno插值
C1=0.1;C2=0.6;C3=0.3;kesai=1e-6;
%正通量，对应U_L
ISp1=0.25*(U1-4*U2+3*U3).^2+13/12*(U1-2*U2+U3).^2;
ISp2=0.25*(U2-U4).^2+13/12*(U2-2*U3+U4).^2;
ISp3=0.25*(3*U3-4*U4+U5).^2+13/12*(U3-2*U4+U5).^2;
alphap1=C1./(kesai+ISp1).^2;
alphap2=C2./(kesai+ISp2).^2;
alphap3=C3./(kesai+ISp3).^2;
omegap1=alphap1./(alphap1+alphap2+alphap3);
omegap2=alphap2./(alphap1+alphap2+alphap3);
omegap3=alphap3./(alphap1+alphap2+alphap3);
fp1=1/3*U1-7/6*U2+11/6*U3;
fp2=-1/6*U2+5/6*U3+1/3*U4;
fp3=1/3*U3+5/6*U4-1/6*U5;
U_L=fp1.*omegap1+fp2.*omegap2+fp3.*omegap3;
C1=0.1;C2=0.6;C3=0.3;kesai=1e-6;
%负通量，对应U_L
ISn1=0.25*(U5-4*U4+3*U3).^2+13/12*(U5-2*U4+U3).^2;
ISn2=0.25*(U4-U2).^2+13/12*(U4-2*U3+U2).^2;
ISn3=0.25*(3*U3-4*U2+U1).^2+13/12*(U3-2*U2+U1).^2;
alphan1=C1./(kesai+ISn1).^2;
alphan2=C2./(kesai+ISn2).^2;
alphan3=C3./(kesai+ISn3).^2;
omegan1=alphan1./(alphan1+alphan2+alphan3);
omegan2=alphan2./(alphan1+alphan2+alphan3);
omegan3=alphan3./(alphan1+alphan2+alphan3);
fn1=1/3*U5-7/6*U4+11/6*U3;
fn2=-1/6*U4+5/6*U3+1/3*U2;
fn3=1/3*U3+5/6*U2-1/6*U1;
U_R=fn1.*omegan1+fn2.*omegan2+fn3.*omegan3;
U_L(:,end)=[];U_R(:,1)=[];
%% step 2 计算Roe平均
%然后计算Roe平均的值
rho_L = U_L(1,:);
u_L = U_L(2,:)./U_L(1,:);
e_L = U_L(3,:)./U_L(1,:);
p_L = (gamma-1)*(e_L - 0.5*u_L.^2).*rho_L;
H_L = e_L + p_L./rho_L;
rho_R = U_R(1,:);
u_R = U_R(2,:)./U_R(1,:);
e_R = U_R(3,:)./U_R(1,:);
p_R = (gamma-1)*(e_R - 0.5*u_R.^2).*rho_R;
H_R = e_R + p_R./rho_R;
rho_roe = ((sqrt(rho_L)+sqrt(rho_R))/2).^2;
u_roe = (sqrt(rho_L).*u_L + sqrt(rho_R).*u_R)./(sqrt(rho_L)+sqrt(rho_R));
H_roe = (sqrt(rho_L).*H_L + sqrt(rho_R).*H_R)./(sqrt(rho_L)+sqrt(rho_R));
p_roe = (gamma-1)/gamma*(rho_roe.*H_roe-0.5*rho_roe.*u_roe.^2);
c_roe = sqrt((gamma-1)*(H_roe-0.5*u_roe.^2));
U_roe = physics_to_conservation(rho_roe,u_roe,p_roe);
%% step 3 4 计算特征分解
%计算特征分解矩阵
Lambda = zeros(3,3,N+1);
S = zeros(3,3,N+1);
S_inv = zeros(3,3,N+1);
A_abs = zeros(3,3,N+1);
F = zeros(3,N+1);
Lambda(1,1,:) = u_roe;
Lambda(2,2,:) = u_roe - c_roe;
Lambda(3,3,:) = u_roe + c_roe;
S(1,1,:) = 1;
S(1,2,:) = 1;
S(1,3,:) = 1;
S(2,1,:) = u_roe;
S(2,2,:) = u_roe - c_roe;
S(2,3,:) = u_roe + c_roe;
S(3,1,:) = 0.5*u_roe.^2;
S(3,2,:) = 0.5*u_roe.^2 + gamma/(gamma-1)*p_roe./rho_roe - u_roe.*c_roe;
S(3,3,:) = 0.5*u_roe.^2 + gamma/(gamma-1)*p_roe./rho_roe + u_roe.*c_roe;
S_inv(1,1,:) = 1 - (gamma-1)*u_roe.^2./(2*c_roe.^2);
S_inv(1,2,:) = (gamma-1)*u_roe./(c_roe.^2);
S_inv(1,3,:) = -(gamma-1)./c_roe.^2;
S_inv(2,1,:) = (gamma-1)*u_roe.^2./(4*c_roe.^2) + 0.5*u_roe./c_roe;
S_inv(2,2,:) = - (gamma-1)*u_roe./(2*c_roe.^2) - 1./(2*c_roe);
S_inv(2,3,:) = (gamma-1)./(2*c_roe.^2);
S_inv(3,1,:) = (gamma-1)*u_roe.^2./(4*c_roe.^2) - 0.5*u_roe./c_roe;
S_inv(3,2,:) = - (gamma-1)*u_roe./(2*c_roe.^2) + 1./(2*c_roe);
S_inv(3,3,:) = (gamma-1)./(2*c_roe.^2);
for i = 1:N+1
   A_abs(:,:,i) = S(:,:,i)*abs(Lambda(:,:,i))*S_inv(:,:,i);
end
%% step 5 计算半点处重构函数通量
F_L = [rho_L.*u_L; rho_L.*u_L.^2 + p_L; u_L.*(rho_L.*e_L+p_L)];
F_R = [rho_R.*u_R; rho_R.*u_R.^2 + p_R; u_R.*(rho_R.*e_R+p_R)];

for i = 1:N+1
    F(:,i) = 0.5*(F_R(:,i)+F_L(:,i)) - 0.5*A_abs(:,:,i)*(U_R(:,i)-U_L(:,i));
end

%% 最终结果
dF = (F(:,2:end)-F(:,1:end-1))/dx;
end

