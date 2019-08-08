function [pca_sub,pca_neighb_sub]=PCAneighb(econc,p0)

%econc��ĳһ��ά�����ĳһ��ͨ����
%p0����һ���ɷַ�������ֱ�߾����ĵ���±�������Ҳ�����ɷֿռ��ԭ�㣬Ϊ1*3������

%pca_sub������econc�ľ�����p0�ĵ�һ���ɷֵ���ֱ�ߵ��±�����
%pca_neighb_sub������econc�ľ�����p0�ĵ�һ���ɷֵ���ֱ�߼���ǰ������4������±�����

econc=logical(econc);
index=find(econc);
[x,y,z]=ind2sub(size(econc),index);
econc_sub=[x,y,z];

if nargin<2
    p0=mean(econc_sub,1);
end

if size(econc_sub,1)>1 %��ͨ��������һ������
    
    [econc_coeff,econc_score,econc_latent,econc_tsquared]=pca(econc_sub);
    econc_latent_percent=100*econc_latent/sum(econc_latent);%��latent�ܺ�ͳһΪ100�����ڹ۲칱����

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

    pca_sub=[pca_xt',pca_yt',pca_zt']; %��pca_sub��ȡ������pca_sub��С������������
    pca_sub=round(pca_sub); %��������ȡ��
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
    pca_neighb_sub=unique(pca_neighb_sub,'rows'); %��pca_sub��ȡ������pca_neighb_sub��С������������
    
    %������ֵxt,yt,ztת��Ϊ����ֵX,Y,Z
    pca_X=pca_yt;
    pca_Y=pca_xt;
    pca_Z=pca_zt;

    PlotElectrodes(pca_neighb_sub,'b','.');hold on;
    plot3(pca_X,pca_Y,pca_Z,'-+','color','r');hold on;
    Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
    [nx,ny,nz]=size(econc);
    axis equal;axis([1 nx 1 ny 1 nz]);
    % axis([150 400 100 270 1 160]);
    xlabel('X');ylabel('Y');zlabel('Z');
    set(gcf,'Position',get(0,'ScreenSize'));
    
end

end