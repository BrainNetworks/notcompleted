function [pca_sub,pca_neighb_sub]=PCAneighb3(econc,p0,range)

%econc：某一三维数组的某一连通分量
%p0：第一主成分方向所在直线经过的点的下标索引，也是主成分空间的原点，为1*3的向量
%range：输入econc的经过点p0的第一主成分单点直线超出econc两端的程度

%pca_sub：输入econc的经过点p0的第一主成分单点直线的下标索引
%pca_neighb_sub：输入econc的经过点p0的第一主成分单点直线及其前后左右4邻域的下标索引

econc=logical(econc);
index=find(econc);
[x,y,z]=ind2sub(size(econc),index);
econc_sub=[x,y,z];

if nargin<3
    range=6;
end

if nargin<2
    p0=mean(econc_sub,1);
end

if size(econc_sub,1)>1 %连通分量多于一个像素
    
    econc_coeff=pca(econc_sub);
    econc_score=(econc_sub-repmat(p0,size(econc_sub,1),1))*econc_coeff;
    econc_score1_min=min(econc_score(:,1));
    econc_score1_max=max(econc_score(:,1));
    pca_tmin=econc_score1_min-range;
    pca_tmax=econc_score1_max+range;

    deltat=0.5;
    pca_t=pca_tmin:deltat:pca_tmax;
    pca_xt=econc_coeff(1,1)*pca_t+p0(1);
    pca_yt=econc_coeff(2,1)*pca_t+p0(2);
    pca_zt=econc_coeff(3,1)*pca_t+p0(3);

    pca_sub=[pca_xt',pca_yt',pca_zt']; %若pca_sub不取整，则pca_sub是小数，不是整数
    pca_sub=round(pca_sub); %四舍五入取整
    pca_neighb_sub=[pca_sub(:,1)-1,pca_sub(:,2),pca_sub(:,3);
                    pca_sub(:,1)+1,pca_sub(:,2),pca_sub(:,3);
                    pca_sub(:,1),pca_sub(:,2)-1,pca_sub(:,3);
                    pca_sub(:,1),pca_sub(:,2)+1,pca_sub(:,3);
                    pca_sub(:,1),pca_sub(:,2),pca_sub(:,3)];
    pca_neighb_sub=unique(pca_neighb_sub,'rows'); %若pca_sub不取整，则pca_neighb_sub是小数，不是整数
    
end

end