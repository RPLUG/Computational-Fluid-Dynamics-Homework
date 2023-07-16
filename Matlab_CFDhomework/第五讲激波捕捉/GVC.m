function [U_L,U_R] = GVC(j,U)
if(abs((U(j)-U(j-1)))<abs((U(j+1)-U(j))))
    U_L=(3*U(j)-U(j-1))/2;
    U_R=(3*U(j+1)-U(j+2))/2;
else
    U_L=(U(j+1)+U(j))/2;
    U_R=(U(j+1)+U(j))/2;
end

