function [connect,connect_sub,connect_neighb,connect_neighb_sub,VolumeOut]=ConnectTwoConC(ConCj,ConCi,Pcloser,Qcloser)

%ConCj：某一三维数组的某一连通分量j，通常是经过判别确定了其属于ConCi所在电极的连通分量
%ConCi：某一三维数组的某一连通分量i,通常是某一电极或某一电极的一部分
%Pcloser：ConCj的所有像素点中，距离ConCi的两个端点Ciendpoints_sub所在直线最近的点的1*3下标索引
%Qcloser：ConCi的两个端点Ciendpoints_sub中，距离ConCj较近的一个，是1*3的下标索引

%connect：ConCi和ConCj之间的单点线连接，是三维体数据
%connect_sub：connect的下标索引，为n*3的二维数组，n为connect的像素总数
%connect_neighb：connect的上下左右前后六邻域
%connect_neighb_sub：connect_neighb的下标索引，为n*3的二维数组，n为connect_neighb的像素总数
%VolumeOut：VolumeOut=logical(ConCi+ConCj+connect)，即为ConCi和ConCj的连接结果

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