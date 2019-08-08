function [ConCi_coeff1,ConCi_explained,ConCj_coeff1,ConCj_explained,pcadiff,iscolinear,judge]=ConCJudge2(ConCj,ConCi,D_ConCj2Qcloser,p0)

%ConCj��ĳһ��ά�����ĳһ��ͨ����j��ͨ���Ǵ��б����ͨ���������ж�ConCj�Ƿ�����ConCi���ڵĵ缫
%ConCi��ĳһ��ά�����ĳһ��ͨ����i,ͨ����ĳһ�缫��ĳһ�缫��һ����
%D_ConCj2Qcloser����ͨ����ConCj��ÿһ�����ص㵽ConCi�Ľ��˵�Qcloser�ľ��룬��һ��n*1�����飬n��ConCj����������
%p0����һ���ɷַ�������ֱ�߾����ĵ���±�������Ҳ�����ɷֿռ��ԭ�㣬Ϊ1*3��������������ȱʡʱp0ΪConCi�������±�������ֵ���ĵ�

%ConCi_coeff����ͨ����ConCi����pca����õ������coeff��Ϊ3��3��
%ConCi_explained����ͨ����ConCi����pca����õ������explained��Ϊ3��1��
%ConCj_coeff����ͨ����ConCj����pca����õ������coeff��Ϊ3��3��
%ConCj_explained����ͨ����ConCj����pca����õ������explained��Ϊ3��1��
%pcadiff��ConCi_coeff1(ConCi_coeff�ĵ�һ��)��ConCj_coeff1(ConCj_coeff�ĵ�һ��)֮���ŷ�Ͼ��룬���ں���ConCi��ConCj���������ɷַ���Ĳ�����
%iscolinear��iscolinear=1���������ConCi�ĵ�һ���ɷ�����ֱ�߷����ϣ�ConCj��ȫλ��ConCi��ǰ�����ߺ�
%judge��judge=1�����ConCj����ConCi���ڵĵ缫

ConCi=logical(ConCi);
ConCj=logical(ConCj);

ConCi_ind=find(ConCi);
[x1,y1,z1]=ind2sub(size(ConCi),ConCi_ind);
ConCi_sub=[x1,y1,z1];

if nargin<4
    p0=mean(ConCi_sub,1);
end

ConCj_ind=find(ConCj);
[x2,y2,z2]=ind2sub(size(ConCj),ConCj_ind);
ConCj_sub=[x2,y2,z2];

judge=0;

if size(ConCj_sub,1)>1

    [ConCi_coeff,~,~,~,ConCi_explained,~]=pca(ConCi_sub);
    [ConCj_coeff,~,~,~,ConCj_explained,~]=pca(ConCj_sub);
    ConCi_coeff1=ConCi_coeff(:,1);
    ConCj_coeff1=ConCj_coeff(:,1);
    pcadiff=norm(ConCi_coeff1-ConCj_coeff1);

    ConCi_score=(ConCi_sub-repmat(p0,size(ConCi_sub,1),1))*ConCi_coeff;
    ConCj_score=(ConCj_sub-repmat(p0,size(ConCj_sub,1),1))*ConCi_coeff; %ע��˴���ConCi_coeff������ConCj_coeff������

    ConCi_score1_min=min(ConCi_score(:,1));
    ConCi_score1_max=max(ConCi_score(:,1));
    ConCj_score1=ConCj_score(:,1);

    N1=sum(ConCj_score1>ConCi_score1_max);
    N2=sum(ConCj_score1<ConCi_score1_min);
    iscolinear=0; %colinear:���ߵ�
    if (N1/length(ConCj_score1)>0.9)||(N2/length(ConCj_score1)>0.9)
        iscolinear=1;
    end

    D_ConCj2Qcloser_min=min(D_ConCj2Qcloser);

    if (pcadiff<0.12)&&(ConCi_explained(1)>94)&&(ConCj_explained(1)>94)&&(D_ConCj2Qcloser_min<18)&&(iscolinear==1)
        judge=1;
    end

end

end