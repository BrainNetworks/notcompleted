function [Pcloser,Qcloser,D_P2Q1Q2,D_Cc2Q1Q2,D_ConCj2Qcloser,iscolinear,judge]=ConCJudge1(ConCj,ConCi,Ciendpoints_sub,p0)

%ConCj��ĳһ��ά�����ĳһ��ͨ����j��ͨ���Ǵ��б����ͨ���������ж�ConCj�Ƿ�����ConCi���ڵĵ缫
%ConCi��ĳһ��ά�����ĳһ��ͨ����i,ͨ����ĳһ�缫��ĳһ�缫��һ����
%Ciendpoints_sub��ConCi�������˵���±���������2*3�Ķ�ά���飬������ConcEndpsDetection��������
%p0����һ���ɷַ�������ֱ�߾����ĵ���±�������Ҳ�����ɷֿռ��ԭ�㣬Ϊ1*3��������������ȱʡʱp0ΪConCi�������±�������ֵ���ĵ�

%Pcloser��ConCj���������ص��У�����ConCi�������˵�Ciendpoints_sub����ֱ������ĵ��1*3�±�����
%Qcloser��ConCi�������˵�Ciendpoints_sub�У�����ConCj�Ͻ���һ������1*3���±�����
%D_P2Q1Q2����ͨ����ConCj��ÿһ�����ص㵽ConCi�������˵�Ciendpoints_sub����ֱ�ߵľ��룬��һ��n*1�����飬n��ConCj����������
%D_Cc2Q1Q2����ͨ����ConCj�����ģ������±�������ֵ����ConCi�������˵�Ciendpoints_sub����ֱ�ߵľ���
%D_ConCj2Qcloser����ͨ����ConCj��ÿһ�����ص㵽ConCi�Ľ��˵�Qcloser�ľ��룬��һ��n*1�����飬n��ConCj����������
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

    R=size(ConCj_sub,1);
    ConCjcenter_sub=round(mean(ConCj_sub,1));

    ConCi_coeff=pca(ConCi_sub);
    ConCj_coeff=pca(ConCj_sub);
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

    Q1=Ciendpoints_sub(1,:);
    Q2=Ciendpoints_sub(2,:);
    D1=norm(Q1-ConCjcenter_sub); %����������2-��������ŷʽ��������������������������ŷʽ����
    D2=norm(Q2-ConCjcenter_sub);

    if D1<D2
        Qcloser=Q1;
    else
        Qcloser=Q2;
    end

    D_P2Q1Q2=zeros(R,1);
    for i=1:R
        P=ConCj_sub(i,:);
        D_P2Q1Q2(i)=norm(cross(Q2-Q1,P-Q1))/norm(Q2-Q1); %�����P��Q1Q2����ֱ�ߵľ���
    end
    D_P2Q1Q2_min=min(D_P2Q1Q2);
    Label=find(D_P2Q1Q2==D_P2Q1Q2_min,1);
    Pcloser=ConCj_sub(Label,:);
    D_Cc2Q1Q2=norm(cross(Q2-Q1,ConCjcenter_sub-Q1))/norm(Q2-Q1);

    %D=pdist2(X,Y)��D(i,j) corresponds to the pairwise distance between observation i in X and observation j in Y
    D_ConCj2Qcloser=pdist2(ConCj_sub,Qcloser,'euclidean');
    D_ConCj2Qcloser_min=min(D_ConCj2Qcloser);

    if (D_P2Q1Q2_min<2)&&(D_Cc2Q1Q2<6)&&(D_ConCj2Qcloser_min<15)&&(iscolinear==1)
        judge=1;
    end

end

end