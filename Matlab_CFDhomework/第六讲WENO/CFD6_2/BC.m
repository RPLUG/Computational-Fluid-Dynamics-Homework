%�߽���������
function U = BC(U,n)
%U���غ�����n�������Ҹ����صĳ���
%���˱߽�
for i = 1:n
    U = [U(:,1) U U(:,end)];
end 
end
