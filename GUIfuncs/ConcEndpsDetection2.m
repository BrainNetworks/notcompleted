function [Cipca,Cipca_sub,Cipca_neighb,Cipca_neighb_sub,...,
    Ciendpoints,Ciendpoints_sub]=ConcEndpsDetection2(ConCi,p0,range)

%ConCi：某一三维数组的某一连通分量
%p0：第一主成分方向所在直线经过的点的下标索引，也是主成分空间的原点，为1*3的向量
%range：输入ConCi的经过点p0的第一主成分单点直线超出ConCi两端的程度

%Cipca：ConCi的第一主成分单点线，其两端超出电极连通分量（PCAneighb3中的range），是三维体数据
%Cipca_sub：Cipca的下标索引，是n*3的二维数组，n是Cipca的像素总数
%Cipca_neighb：e12Cpca的前后左右4领域（不包括上下），是三维体数据
%Cipca_neighb_sub：Cipca_neighb的下标索引，是n*3的二维数组，n是Cipca_neighb的像素总数
%Ciendpoints：ConCi的两个端点，是三维体数据
%Ciendpoints_sub：Ciendpoints的下标索引，是2*3的二维数组

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