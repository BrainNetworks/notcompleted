function [Cipca,Cipca_sub,Cipca_neighb,Cipca_neighb_sub,...,
    Ciendpoints,Ciendpoints_sub]=ConcEndpsDetection2(ConCi,p0,range)

%ConCi��ĳһ��ά�����ĳһ��ͨ����
%p0����һ���ɷַ�������ֱ�߾����ĵ���±�������Ҳ�����ɷֿռ��ԭ�㣬Ϊ1*3������
%range������ConCi�ľ�����p0�ĵ�һ���ɷֵ���ֱ�߳���ConCi���˵ĳ̶�

%Cipca��ConCi�ĵ�һ���ɷֵ����ߣ������˳����缫��ͨ������PCAneighb3�е�range��������ά������
%Cipca_sub��Cipca���±���������n*3�Ķ�ά���飬n��Cipca����������
%Cipca_neighb��e12Cpca��ǰ������4���򣨲��������£�������ά������
%Cipca_neighb_sub��Cipca_neighb���±���������n*3�Ķ�ά���飬n��Cipca_neighb����������
%Ciendpoints��ConCi�������˵㣬����ά������
%Ciendpoints_sub��Ciendpoints���±���������2*3�Ķ�ά����

index=find(ConCi);
[r,c,s]=ind2sub(size(ConCi),index);
sub=[r,c,s];

if nargin<3
    range=3;
end

if nargin<2
    p0=mean(sub,1);
end

ConCi=logical(ConCi);
[Cipca_sub,Cipca_neighb_sub]=PCAneighb3(ConCi,p0,range);

Cipca_ind=sub2ind(size(ConCi),Cipca_sub(:,1),Cipca_sub(:,2),Cipca_sub(:,3));
Cipca=false(size(ConCi));
Cipca(Cipca_ind)=1; 

Cipca_neighb_ind=sub2ind(size(ConCi),Cipca_neighb_sub(:,1),Cipca_neighb_sub(:,2),Cipca_neighb_sub(:,3));
Cipca_neighb=false(size(ConCi));
Cipca_neighb(Cipca_neighb_ind)=1;

% cutCipcan=ConCi.*logical(Cipca);
cutCipcan=ConCi.*logical(Cipca_neighb);

cutCipcan_ind=find(cutCipcan);
[x,y,z]=ind2sub(size(cutCipcan),cutCipcan_ind);
Ciendpoints_sub=zeros(2,3);

cutCipcan_sub=[x,y,z];
cutCipcan_coeff=pca(cutCipcan_sub);
cutCipcan_score=(cutCipcan_sub-repmat(p0,size(cutCipcan_sub,1),1))*cutCipcan_coeff;
cutCipcan_score1_min=min(cutCipcan_score(:,1));
cutCipcan_score1_max=max(cutCipcan_score(:,1));

pca_start=[cutCipcan_score1_min,0,0];
pca_end=[cutCipcan_score1_max,0,0];

endps1_sub=pca_start/cutCipcan_coeff+repmat(p0,size(pca_start,1),1);
endps2_sub=pca_end/cutCipcan_coeff+repmat(p0,size(pca_end,1),1);

Ciendpoints_sub(1,:)=round(endps1_sub);
Ciendpoints_sub(2,:)=round(endps2_sub);

Ciendpoints_ind=sub2ind(size(cutCipcan),Ciendpoints_sub(:,1),Ciendpoints_sub(:,2),Ciendpoints_sub(:,3));
Ciendpoints=false(size(cutCipcan));
Ciendpoints(Ciendpoints_ind)=1; 

end