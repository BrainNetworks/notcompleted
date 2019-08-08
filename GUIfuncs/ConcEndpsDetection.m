function [ConCi_gmtm,ConCi_gmt_mean,Cipca,Cipca_sub,Cipca_neighb,Cipca_neighb_sub,...,
    Ciendpoints,Ciendpoints_sub]=ConcEndpsDetection(ConCi,orig,ratio,range)

%ConCi��ĳһ��ά�����ĳһ��ͨ����
%orig��δ������ֵ����ĺ��лҶ�ֵ��Ϣ�ĳ�ʼͼ����electrode1
%ratio����ֵ��Ϊ0~1��С����ConCi�ĵ�һ���ɷֵ���ֱ�߾���ConCi�����������лҶ�ֵ�����ֵratio���ĵ㼯��������ֵ����p0
%range������ConCi�ľ�����p0�ĵ�һ���ɷֵ���ֱ�߳���ConCi���˵ĳ̶�

%ConCi_gmtm����p0������ά������
%ConCi_gmt_mean����p0���±���������1*3�Ķ�ά����
%Cipca��ConCi�ĵ�һ���ɷֵ����ߣ������˳����缫��ͨ������PCAneighb3�е�range��������ά������
%Cipca_sub��Cipca���±���������n*3�Ķ�ά���飬n��Cipca����������
%Cipca_neighb��e12Cpca��ǰ������4���򣨲��������£�������ά������
%Cipca_neighb_sub��Cipca_neighb���±���������n*3�Ķ�ά���飬n��Cipca_neighb����������
%Ciendpoints��ConCi�������˵㣬����ά������
%Ciendpoints_sub��Ciendpoints���±���������2*3�Ķ�ά����

if nargin<4
    range=6;
end

if nargin<3
    ratio=0.9;
end

ConCi=logical(ConCi);
[nx,ny,nz]=size(ConCi);
ConCi_gray=orig.*ConCi;
ConCi_gray_max=max(max(max(ConCi_gray)));
% ConCi_gm=ConCi_gray;
% ConCi_gm(ConCi_gm~=ConCi_gray_max)=0;

ConCi_gm_thre=ConCi_gray_max*ratio;
ConCi_gmt=ConCi_gray;
ConCi_gmt(ConCi_gmt<ConCi_gm_thre)=0;

ConCi_gmt_ind=find(ConCi_gmt);
[x1,y1,z1]=ind2sub(size(ConCi_gmt),ConCi_gmt_ind);
ConCi_gmt_sub=[x1,y1,z1];
ConCi_gmt_mean=round(mean(ConCi_gmt_sub,1));
ConCi_gmtm=zeros(nx,ny,nz);
ConCi_gmtm(ConCi_gmt_mean(1),ConCi_gmt_mean(2),ConCi_gmt_mean(3))=1;

[Cipca_sub,Cipca_neighb_sub]=PCAneighb3(ConCi_gray,ConCi_gmt_mean,range);

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

% dx=max(x)-min(x);
% dy=max(y)-min(y);
% dz=max(z)-min(z);
% d=[dx,dy,dz];
% if max(d)==dx
%     [xnew,xL]=sort(x,1,'ascend');
%     ynew=y(xL);znew=z(xL);
% else
%     if max(d)==dy
%         [ynew,yL]=sort(y,1,'ascend');
%         xnew=x(yL);znew=z(yL);
%     else
%         [znew,zL]=sort(z,1,'ascend');
%         xnew=x(zL);ynew=y(zL);
%     end
% end
% cutCipca_sub=[xnew,ynew,znew];
% Ciendpoints_sub(1,:)=round(cutCipca_sub(1,:));
% Ciendpoints_sub(2,:)=round(cutCipca_sub(end,:));

cutCipcan_sub=[x,y,z];
cutCipcan_coeff=pca(cutCipcan_sub);
cutCipcan_score=(cutCipcan_sub-repmat(ConCi_gmt_mean,size(cutCipcan_sub,1),1))*cutCipcan_coeff;
cutCipcan_score1_min=min(cutCipcan_score(:,1));
cutCipcan_score1_max=max(cutCipcan_score(:,1));

pca_start=[cutCipcan_score1_min,0,0];
pca_end=[cutCipcan_score1_max,0,0];

endps1_sub=pca_start/cutCipcan_coeff+repmat(ConCi_gmt_mean,size(pca_start,1),1);
endps2_sub=pca_end/cutCipcan_coeff+repmat(ConCi_gmt_mean,size(pca_end,1),1);

Ciendpoints_sub(1,:)=round(endps1_sub);
Ciendpoints_sub(2,:)=round(endps2_sub);

Ciendpoints_ind=sub2ind(size(cutCipcan),Ciendpoints_sub(:,1),Ciendpoints_sub(:,2),Ciendpoints_sub(:,3));
Ciendpoints=false(size(cutCipcan));
Ciendpoints(Ciendpoints_ind)=1; 

end