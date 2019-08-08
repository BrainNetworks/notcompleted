function [ConCi_gmtm,ConCi_gmt_mean,Cipca,Cipca_sub,Cipca_neighb,Cipca_neighb_sub,...,
    Ciendpoints,Ciendpoints_sub]=ConcEndpsDetection(ConCi,orig,ratio,range)

%ConCi：某一三维数组的某一连通分量
%orig：未经过阈值处理的含有灰度值信息的初始图像，如electrode1
%ratio：是值域为0~1的小数，ConCi的第一主成分单点直线经过ConCi的所有像素中灰度值达最大值ratio倍的点集的索引均值中心p0
%range：输入ConCi的经过点p0的第一主成分单点直线超出ConCi两端的程度

%ConCi_gmtm：点p0，是三维体数据
%ConCi_gmt_mean：点p0的下标索引，是1*3的二维数组
%Cipca：ConCi的第一主成分单点线，其两端超出电极连通分量（PCAneighb3中的range），是三维体数据
%Cipca_sub：Cipca的下标索引，是n*3的二维数组，n是Cipca的像素总数
%Cipca_neighb：e12Cpca的前后左右4领域（不包括上下），是三维体数据
%Cipca_neighb_sub：Cipca_neighb的下标索引，是n*3的二维数组，n是Cipca_neighb的像素总数
%Ciendpoints：ConCi的两个端点，是三维体数据
%Ciendpoints_sub：Ciendpoints的下标索引，是2*3的二维数组

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