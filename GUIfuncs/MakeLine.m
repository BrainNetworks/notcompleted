function [connect_sub,connect_neighb_sub]=MakeLine(P,Q)

%P��P�㣬��1*3�������������±�����[row,column,slice]
%Q��Q�㣬��1*3�������������±�����[row,column,slice]

%connect_sub��P���Q��֮�������connect���±�������Ϊn*3�Ķ�ά���飬nΪconnect����������
%connect_neighb_sub��P���Q��֮�������connect����������ǰ��������connect_neighb���±�������Ϊn*3�Ķ�ά���飬nΪconnect_neighb����������

d0=0.5;
D=norm(P-Q);
d=(0:d0:D)';
d=[d;D];

R1=P(1);
C1=P(2);
S1=P(3);
R2=Q(1);
C2=Q(2);
S2=Q(3);

Rn=d/D*(R2-R1)+R1;
Cn=d/D*(C2-C1)+C1;
Sn=d/D*(S2-S1)+S1;
connect_sub=[Rn,Cn,Sn];
% connect_sub=round(connect_sub);
% connect=false(nx,ny,nz);
% connect_ind=sub2ind(size(connect),connect_sub(:,1),connect_sub(:,2),connect_sub(:,3));
% connect(connect_ind)=1;

connect_neighb_sub=[connect_sub(:,1)-1,connect_sub(:,2),connect_sub(:,3);
                    connect_sub(:,1)+1,connect_sub(:,2),connect_sub(:,3);
                    connect_sub(:,1),connect_sub(:,2)-1,connect_sub(:,3);
                    connect_sub(:,1),connect_sub(:,2)+1,connect_sub(:,3);
                    connect_sub(:,1),connect_sub(:,2),connect_sub(:,3)-1;
                    connect_sub(:,1),connect_sub(:,2),connect_sub(:,3)+1;
                    connect_sub(:,1),connect_sub(:,2),connect_sub(:,3)];
% connect_neighb=false(nx,ny,nz);
% connect_neighb_ind=sub2ind(size(connect_neighb),connect_neighb_sub(:,1),connect_neighb_sub(:,2),connect_neighb_sub(:,3));
% connect_neighb(connect_neighb_ind)=1;

end
