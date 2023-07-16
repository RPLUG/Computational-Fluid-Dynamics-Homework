Vxy=importdata('mesh.txt');%����ڵ������������Vxy
M=80;%����
N=240;%����
%�ų�80*240����
for m=1:M
Rowhead=1+N*(m-1);
Rowtail=N+N*(m-1);
xy(m,:)=Vxy(Rowhead:Rowtail,1)+1i*Vxy(Rowhead:Rowtail,2);
end
clear Rowhead;clear Rowtail;clear Vxy;
%�����������ĵ㣬���ı������ģ�ͨ���ָ�������������������Ȩƽ���õ�
for m=1:M-1
    for n=1:N-1
        S1=-(0.5*((-real(xy(m,n))+real(xy(m+1,n)))*(imag(xy(m,n+1))-imag(xy(m,n))))-...
            ((-imag(xy(m,n))+imag(xy(m+1,n)))*(real(xy(m,n+1))-real(xy(m,n)))));
        S2= 0.5*((-real(xy(m+1,n+1))+real(xy(m+1,n)))*(imag(xy(m,n+1))-imag(xy(m+1,n+1))))-...
            ((-imag(xy(m+1,n+1))+imag(xy(m+1,n)))*(real(xy(m,n+1))-real(xy(m+1,n+1))));
        c1=1/3*(xy(m,n)+xy(m+1,n)+xy(m,n+1));
        c2=1/3*(xy(m+1,n+1)+xy(m+1,n)+xy(m,n+1));
        xyCell(m,n)=(S1*c1+S2*c2)/(S1+S2);
        Scell(m,n)=S1+S2;
    end
end
%���յõ�������ڵ��������Ϊxy(80,240)���������ĵ��������ΪxyCell(79,239)
%�±���ͬʱ���ھ���ռ��xy��ӦxyCell�ġ����Ϸ���
%EdgeEWNS�Լ�VnESWN����������ı߶�Ӧ��ҲΪ(79,239)
%����������߶�Ӧ�ĵ�λ��������Ϊ��ȷ���ǡ��ⷨ�򡱣�������������߽��ѭ���������һ����õ�
for m=1:M-1
    for n=1:N-1
EdgeE(m,n)=xy(m,n+1)-xy(m+1,n+1);EdgeN(m,n)=xy(m,n)-xy(m,n+1);
EdgeW(m,n)=xy(m+1,n)-xy(m,n);EdgeS(m,n)=xy(m+1,n+1)-xy(m+1,n);


   nE(m,n)=-(imag(EdgeE(m,n))-real(EdgeE(m,n))*1i)/(abs(EdgeE(m,n)));
   nS(m,n)=-(imag(EdgeS(m,n))-real(EdgeS(m,n))*1i)/(abs(EdgeS(m,n)));
   nW(m,n)=-(imag(EdgeW(m,n))-real(EdgeW(m,n))*1i)/(abs(EdgeW(m,n)));
   nN(m,n)=-(imag(EdgeN(m,n))-real(EdgeN(m,n))*1i)/(abs(EdgeN(m,n)));
    end
end
%���������������ʽͳһΪX_dummy*��*Ϊ����
%�������nwes���������λ������������ճ����������ɢ��ʽ
%���㷽ʽΪ����������ĵ�һ�������������ȥ����������ĵڶ������������ֱ꣬�ӽ��õ����������ڵ�һ��������������
xyCell_dummyN=2*xyCell(1,:)-xyCell(2,:);
xyCell_dummyS=2*xyCell(end,:)-xyCell(end-1,:);
xyCell_dummyW=2*xyCell(:,1)-xyCell(:,2);
xyCell_dummyE=2*xyCell(:,end)-xyCell(:,end-1);


