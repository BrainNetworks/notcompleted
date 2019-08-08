function ElectrodePCAPlot(econc,color,patt,p0)

%econc：某一三维数组的某一连通分量
%color：画图时，构成连通分量的点的颜色
%patt：画图时，构成连通分量的点的样式
%p0：第一主成分方向所在直线经过的点的下标索引，也是主成分空间的原点，为1*3的向量

index=find(econc);
[x,y,z]=ind2sub(size(econc),index);
econc_sub=[x,y,z];

if nargin<4
    p0=mean(econc_sub,1);
end

if nargin<3
    patt='-';
end

if nargin<2
    color='r';
end

if size(econc_sub,1)>1 %连通分量多于一个像素
    [econc_coeff,econc_score,econc_latent,econc_tsquared]=pca(econc_sub);
    econc_latent_percent=100*econc_latent/sum(econc_latent);%将latent总和统一为100，便于观察贡献率

    % if econc_latent_percent(1)>90 
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
        range1=-6;
        range2=6;

        %pca_xt=econc_coeff(1,1)*[econc_tmin,econc_tmax]+econc_center(1)=[econc_xmin+range1,econc_xmax+range2]
        %pca_yt=econc_coeff(2,1)*[econc_tmin,econc_tmax]+econc_center(2)=[econc_ymin+range1,econc_ymax+range2]
        %pca_zt=econc_coeff(3,1)*[econc_tmin,econc_tmax]+econc_center(3)=[econc_zmin+range1,econc_zmax+range2]
        %当econc_coeff(1,1)为负时，econc_xmax+range2=econc_coeff(1,1)*econc_tmin+econc_center(1)，
        %也就是说，此时实际上econc_t2才是最小值，而econc_t1是最大值
        if max(econc_d)==econc_dx
            pca_t1=((econc_xmin+range1)-p0(1))/econc_coeff(1,1);
            pca_t2=((econc_xmax+range2)-p0(1))/econc_coeff(1,1);
        end
        if max(econc_d)==econc_dy
            pca_t1=((econc_ymin+range1)-p0(2))/econc_coeff(2,1);
            pca_t2=((econc_ymax+range2)-p0(2))/econc_coeff(2,1);
        end
        if max(econc_d)==econc_dz
            pca_t1=((econc_zmin+range1)-p0(3))/econc_coeff(3,1);
            pca_t2=((econc_zmax+range2)-p0(3))/econc_coeff(3,1);
        end
        pca_tmin=min(pca_t1,pca_t2);
        pca_tmax=max(pca_t1,pca_t2);

        pca_t=pca_tmin:pca_tmax;
        pca_xt=econc_coeff(1,1)*pca_t+p0(1);
        pca_yt=econc_coeff(2,1)*pca_t+p0(2);
        pca_zt=econc_coeff(3,1)*pca_t+p0(3);

        %将索引值xt,yt,zt转换为坐标值X,Y,Z
        pca_X=pca_yt;
        pca_Y=pca_xt;
        pca_Z=pca_zt;

        plot3(pca_X,pca_Y,pca_Z,patt,'color',color);
        Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
        [nx,ny,nz]=size(econc);
        axis equal;axis([1 ny 1 nx 1 nz]);
        % axis([150 400 100 270 1 160]);
        xlabel('X');ylabel('Y');zlabel('Z');
        hold on;
        set(gcf,'Position',get(0,'ScreenSize'));
    % end
end

end