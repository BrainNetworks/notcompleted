function ElectrodePCAPlot2(econc,p0,color,patt)

%econc��ĳһ��ά�����ĳһ��ͨ����
%p0����һ���ɷַ�������ֱ�߾����ĵ���±�������Ҳ�����ɷֿռ��ԭ�㣬Ϊ1*3������
%color����ͼʱ��������ͨ�����ĵ����ɫ
%patt����ͼʱ��������ͨ�����ĵ����ʽ

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

if size(econc_sub,1)>1 %��ͨ��������һ������
        
    range1=-6;
    range2=6;

%--------------------------------------------------------------------------
%��С�ڵĴ���������pca����õ���econc_score�ǽ����ڶ�econc_sub�������Ļ��Ļ���
%�ϵģ��൱��Ĭ�������ɷֿռ���mean(econc_sub,1)Ϊԭ�㣬��Ĭ�ϵ�һ���ɷַ�����
%��ֱ�߾����ĵ���mean(econc_sub,1)����ʵ����������Ҫ��ԭ��������p0��p0��һ����
%mean(econc_sub,1)��

%     [econc_coeff,econc_score,econc_latent,econc_tsquared]=pca(econc_sub);
%     econc_score_min1=min(econc_score(:,1));
%     econc_score_max1=max(econc_score(:,1));
%     
%     pca_tmin=econc_score_min1+range1;
%     pca_tmax=econc_score_max1+range2;
%--------------------------------------------------------------------------
%��С�ڲ�����ȷ��
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

    %������ֵxt,yt,ztת��Ϊ����ֵX,Y,Z
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