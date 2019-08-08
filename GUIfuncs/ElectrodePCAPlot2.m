function ElectrodePCAPlot2(econc,p0,color,patt)

%econc：某一三维数组的某一连通分量
%p0：第一主成分方向所在直线经过的点的下标索引，也是主成分空间的原点，为1*3的向量
%color：画图时，构成连通分量的点的颜色
%patt：画图时，构成连通分量的点的样式

index=find(econc);
[x,y,z]=ind2sub(size(econc),index);
econc_sub=[x,y,z];

if nargin<4
    patt='-';
end

if nargin<3
    color='b';
end

if nargin<2
    p0=mean(econc_sub,1);
end

if size(econc_sub,1)>1 %连通分量多于一个像素
        
    range1=-6;
    range2=6;

%--------------------------------------------------------------------------
%此小节的错误在于用pca计算得到的econc_score是建立在对econc_sub进行中心化的基础
%上的，相当于默认了主成分空间以mean(econc_sub,1)为原点，即默认第一主成分方向所
%在直线经过的点是mean(econc_sub,1)，但实际上我们想要的原点是输入p0，p0不一定是
%mean(econc_sub,1)。

%     [econc_coeff,econc_score,econc_latent,econc_tsquared]=pca(econc_sub);
%     econc_score_min1=min(econc_score(:,1));
%     econc_score_max1=max(econc_score(:,1));
%     
%     pca_tmin=econc_score_min1+range1;
%     pca_tmax=econc_score_max1+range2;
%--------------------------------------------------------------------------
%此小节才是正确的
    econc_coeff=pca(econc_sub);
    econc_score=(econc_sub-repmat(p0,size(econc_sub,1),1))*econc_coeff;
    econc_score1_min=min(econc_score(:,1));
    econc_score1_max=max(econc_score(:,1));
    pca_tmin=econc_score1_min+range1;
    pca_tmax=econc_score1_max+range2;
%--------------------------------------------------------------------------

    pca_t=pca_tmin:pca_tmax;
    pca_xt=econc_coeff(1,1)*pca_t+p0(1);
    pca_yt=econc_coeff(2,1)*pca_t+p0(2);
    pca_zt=econc_coeff(3,1)*pca_t+p0(3);

    %将索引值xt,yt,zt转换为坐标值X,Y,Z
    pca_X=pca_yt;
    pca_Y=pca_xt;
    pca_Z=pca_zt;

%     plot3(pca_X,pca_Y,pca_Z,patt,'color',color);
    plot3(pca_X,pca_Y,pca_Z,patt,'color',color,'linewidth',2);
    % Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
    Azimuth=41;Elevation=5;view(Azimuth,Elevation);
    [nx,ny,nz]=size(econc);
    axis equal;axis([1 ny 1 nx 1 nz]);
    % axis([150 400 100 270 1 160]);
    xlabel('X');ylabel('Y');zlabel('Z');
    hold on;
    set(gcf,'Position',get(0,'ScreenSize'));
    
end

end