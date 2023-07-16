%边界条件处理
function U = BC(U,n)
%U是守恒量，n代表左右各延拓的长度
%输运边界
for i = 1:n
    U = [U(:,1) U U(:,end)];
end 
end
