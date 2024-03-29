function [connect,connect_sub,connect_neighb,connect_neighb_sub,VolumeOut]=ConnectTwoConC3(ConCi,ConCj,type)

%ConCi：属于某一电极的第一个连通分量，是三维数组
%ConCj：属于某一电极的第二个连通分量，是三维数组
%type：连通分量的连接方式，type=1时VolumeOut=logical(ConCi+ConCj+connect_neighb);type=2时VolumeOut=logical(ConCi+ConCj+connect)

%connect：ConCi和ConCj之间的单点线连接，是三维体数据
%connect_sub：connect的下标索引，为n*3的二维数组，n为connect的像素总数
%connect_neighb：connect的上下左右前后六邻域
%connect_neighb_sub：connect_neighb的下标索引，为n*3的二维数组，n为connect_neighb的像素总数
%VolumeOut：VolumeOut=logical(ConCi+ConCj+connect_neighb)，即为ConCi和ConCj的连接结果

if nargin<3
    type=1;
end

ConCi=logical(ConCi);
ConCj=logical(ConCj);

[nx,ny,nz]=size(ConCi);
VolumeSize=[nx,ny,nz];

ConCj_Centroid=regionprops3(ConCj,'Centroid');
ConCj_Centroid=ConCj_Centroid{1,1};
ConCj_Centrsub=ConCj_Centroid(:,[2,1,3]);

ConCi_ind=find(ConCi);
[c,r,s]=ind2sub(size(ConCi),ConCi_ind);
ConCi_sub=[c,r,s];
D=pdist2(ConCi_sub,ConCj_Centrsub,'euclidean');
Dmin=min(D);
Label=find(D==Dmin,1);
ConCi_closer=ConCi_sub(Label,:);

[connect,connect_sub,connect_neighb,connect_neighb_sub]=MakeLine2(ConCi_closer,ConCj_Centrsub,VolumeSize);
if type==2
    VolumeOut=logical(ConCi+ConCj+connect);
else
    VolumeOut=logical(ConCi+ConCj+connect_neighb);
end

end