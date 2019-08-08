function [pca_sub,pca_neighb_sub]=PCAneighb2(econc,p0,range)

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
    
    [econc_coeff,econc_score,econc_latent,econc_tsquared]=pca(econc_sub);
    econc_latent_percent=100*econc_latent/sum(econc_latent);%将latent总和统一为100，便于观察贡献率

    econc_xmin=min(econc_sub(:,1));
    econc_xmax=max(econc_sub(:,1));
    econc_ymin=min(econc_sub(:,2));
    econc_ymax=max(econc_sub(:,2));    
    econc_zmin=min(econc_sub(:,3));
    econc_zmax=max(econc_sub(:,3));
    econc_dx=econc_xmax-econc_xmin;
    econc_dy=econc_ymax-econc_ymin;
    econc_dz=econc_zmax-econc_zmin;
    econc_d=[econc_dx,econc_dy,econc_dz];

    if max(econc_d)==econc_dx
        pca_t1=((econc_xmin-range)-p0(1))/econc_coeff(1,1);
        pca_t2=((econc_xmax+range)-p0(1))/econc_coeff(1,1);
    end
    if max(econc_d)==econc_dy
        pca_t1=((econc_ymin-range)-p0(2))/econc_coeff(2,1);
        pca_t2=((econc_ymax+range)-p0(2))/econc_coeff(2,1);
    end
    if max(econc_d)==econc_dz
        pca_t1=((econc_zmin-range)-p0(3))/econc_coeff(3,1);
        pca_t2=((econc_zmax+range)-p0(3))/econc_coeff(3,1);
    end
    pca_tmin=min(pca_t1,pca_t2);
    pca_tmax=max(pca_t1,pca_t2);

    pca_t=pca_tmin:pca_tmax;
    pca_xt=econc_coeff(1,1)*pca_t+p0(1);
    pca_yt=econc_coeff(2,1)*pca_t+p0(2);
    pca_zt=econc_coeff(3,1)*pca_t+p0(3);

    pca_sub=[pca_xt',pca_yt',pca_zt']; %若pca_sub不取整，则pca_sub是小数，不是整数
    pca_sub=round(pca_sub); %四舍五入取整
    pca_neighb_sub=[];
    el=size(pca_sub,1);
    for n=1:el
        pcan_sub=pca_sub(n,:);
        pcan_neighb_sub=[pcan_sub(1)-1,pcan_sub(2),pcan_sub(3);
                          pcan_sub(1)+1,pcan_sub(2),pcan_sub(3);
                          pcan_sub(1),pcan_sub(2)-1,pcan_sub(3);
                          pcan_sub(1),pcan_sub(2)+1,pcan_sub(3);
                          pcan_sub(1),pcan_sub(2),pcan_sub(3)];
        pca_neighb_sub=[pca_neighb_sub;pcan_neighb_sub];
    end
    pca_neighb_sub=unique(pca_neighb_sub,'rows'); %若pca_sub不取整，则pca_neighb_sub是小数，不是整数
    
end

end