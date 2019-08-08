function [connect,connect_sub,connect_neighb,connect_neighb_sub,VolumeOut]=ConnectTwoConC(ConCj,ConCi,Pcloser,Qcloser)

%ConCj��ĳһ��ά�����ĳһ��ͨ����j��ͨ���Ǿ����б�ȷ����������ConCi���ڵ缫����ͨ����
%ConCi��ĳһ��ά�����ĳһ��ͨ����i,ͨ����ĳһ�缫��ĳһ�缫��һ����
%Pcloser��ConCj���������ص��У�����ConCi�������˵�Ciendpoints_sub����ֱ������ĵ��1*3�±�����
%Qcloser��ConCi�������˵�Ciendpoints_sub�У�����ConCj�Ͻ���һ������1*3���±�����

%connect��ConCi��ConCj֮��ĵ��������ӣ�����ά������
%connect_sub��connect���±�������Ϊn*3�Ķ�ά���飬nΪconnect����������
%connect_neighb��connect����������ǰ��������
%connect_neighb_sub��connect_neighb���±�������Ϊn*3�Ķ�ά���飬nΪconnect_neighb����������
%VolumeOut��VolumeOut=logical(ConCi+ConCj+connect)����ΪConCi��ConCj�����ӽ��

ConCi=logical(ConCi);
ConCj=logical(ConCj);

[nx,ny,nz]=size(ConCj);

d0=0.5;
D=norm(Pcloser-Qcloser);
d=(0:d0:D)';
d=[d;D];

R1=Pcloser(1);
C1=Pcloser(2);
S1=Pcloser(3);
R2=Qcloser(1);
C2=Qcloser(2);
S2=Qcloser(3);

Rn=d/D*(R2-R1)+R1;
Cn=d/D*(C2-C1)+C1;
Sn=d/D*(S2-S1)+S1;
connect_sub=round([Rn,Cn,Sn]);
connect=false(nx,ny,nz);
connect_ind=sub2ind(size(connect),connect_sub(:,1),connect_sub(:,2),connect_sub(:,3));
connect(connect_ind)=1;

connect_neighb_sub=[connect_sub(:,1)-1,connect_sub(:,2),connect_sub(:,3);
                    connect_sub(:,1)+1,connect_sub(:,2),connect_sub(:,3);
                    connect_sub(:,1),connect_sub(:,2)-1,connect_sub(:,3);
                    connect_sub(:,1),connect_sub(:,2)+1,connect_sub(:,3);
                    connect_sub(:,1),connect_sub(:,2),connect_sub(:,3)-1;
                    connect_sub(:,1),connect_sub(:,2),connect_sub(:,3)+1;
                    connect_sub(:,1),connect_sub(:,2),connect_sub(:,3)];

connect_neighb=false(nx,ny,nz);
connect_neighb_ind=sub2ind(size(connect_neighb),connect_neighb_sub(:,1),connect_neighb_sub(:,2),connect_neighb_sub(:,3));
connect_neighb(connect_neighb_ind)=1;

VolumeOut=logical(ConCi+ConCj+connect);
% VolumeOut=logical(ConCi+ConCj+connect_neighb);

end