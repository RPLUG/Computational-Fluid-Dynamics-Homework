%得到一个三维矩阵存储守恒量 ，这里能量密度用体积能量密度表示，与李新亮老师ppt以及toro书上一致
function U = ptoc2(rho,u,v,p)
global gamma;

%给出守恒量
U(:,:,1) = rho;
U(:,:,2) = rho.*u;
U(:,:,3) = rho.*v;
U(:,:,4) =p/(gamma-1) + 0.5*rho.*(v.^2+u.^2);

end