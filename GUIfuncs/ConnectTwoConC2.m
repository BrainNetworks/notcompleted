function [connect,connect_sub,connect_neighb,connect_neighb_sub,VolumeOut]=ConnectTwoConC2(ConCi,ConCj,type)

%ConCi������ĳһ�缫�ĵ�һ����ͨ����������ά����
%ConCj������ĳһ�缫�ĵڶ�����ͨ����������ά����
%type����ͨ���������ӷ�ʽ��type=1ʱVolumeOut=logical(ConCi+ConCj+connect_neighb);type=2ʱVolumeOut=logical(ConCi+ConCj+connect)

%connect��ConCi��ConCj֮��ĵ��������ӣ�����ά������
%connect_sub��connect���±�������Ϊn*3�Ķ�ά���飬nΪconnect����������
%connect_neighb��connect����������ǰ��������
%connect_neighb_sub��connect_neighb���±�������Ϊn*3�Ķ�ά���飬nΪconnect_neighb����������
%VolumeOut��VolumeOut=logical(ConCi+ConCj+connect_neighb)����ΪConCi��ConCj�����ӽ��

if nargin<3
    type=1;
end

ConCi=logical(ConCi);
ConCj=logical(ConCj);

[nx,ny,nz]=size(ConCi);
VolumeSize=[nx,ny,nz];

ConCi_Centroid=regionprops3(ConCi,'Centroid');
ConCj_Centroid=regionprops3(ConCj,'Centroid');
ConCi_Centroid=ConCi_Centroid{1,1};
ConCj_Centroid=ConCj_Centroid{1,1};
ConCi_Centrsub=ConCi_Centroid(:,[2,1,3]);
ConCj_Centrsub=ConCj_Centroid(:,[2,1,3]);

[connect,connect_sub,connect_neighb,connect_neighb_sub]=MakeLine2(ConCi_Centrsub,ConCj_Centrsub,VolumeSize);
if type==2
    VolumeOut=logical(ConCi+ConCj+connect);
else
    VolumeOut=logical(ConCi+ConCj+connect_neighb);
end

end