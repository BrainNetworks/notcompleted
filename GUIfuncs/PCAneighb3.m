function [pca_sub,pca_neighb_sub]=PCAneighb3(econc,p0,range)

%econc��ĳһ��ά�����ĳһ��ͨ����
%p0����һ���ɷַ�������ֱ�߾����ĵ���±�������Ҳ�����ɷֿռ��ԭ�㣬Ϊ1*3������
%range������econc�ľ�����p0�ĵ�һ���ɷֵ���ֱ�߳���econc���˵ĳ̶�

%pca_sub������econc�ľ�����p0�ĵ�һ���ɷֵ���ֱ�ߵ��±�����
%pca_neighb_sub������econc�ľ�����p0�ĵ�һ���ɷֵ���ֱ�߼���ǰ������4������±�����

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

if size(econc_sub,1)>1 %��ͨ��������һ������
    
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

    pca_sub=[pca_xt',pca_yt',pca_zt']; %��pca_sub��ȡ������pca_sub��С������������
    pca_sub=round(pca_sub); %��������ȡ��
    pca_neighb_sub=[pca_sub(:,1)-1,pca_sub(:,2),pca_sub(:,3);
                    pca_sub(:,1)+1,pca_sub(:,2),pca_sub(:,3);
                    pca_sub(:,1),pca_sub(:,2)-1,pca_sub(:,3);
                    pca_sub(:,1),pca_sub(:,2)+1,pca_sub(:,3);
                    pca_sub(:,1),pca_sub(:,2),pca_sub(:,3)];
    pca_neighb_sub=unique(pca_neighb_sub,'rows'); %��pca_sub��ȡ������pca_neighb_sub��С������������
    
end

end