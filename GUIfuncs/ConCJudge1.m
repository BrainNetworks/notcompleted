function [Pcloser,Qcloser,D_P2Q1Q2,D_Cc2Q1Q2,D_ConCj2Qcloser,iscolinear,judge]=ConCJudge1(ConCj,ConCi,Ciendpoints_sub,p0)

%ConCj：某一三维数组的某一连通分量j，通常是待判别的连通分量，即判断ConCj是否属于ConCi所在的电极
%ConCi：某一三维数组的某一连通分量i,通常是某一电极或某一电极的一部分
%Ciendpoints_sub：ConCi的两个端点的下标索引，是2*3的二维数组，即函数ConcEndpsDetection的输出结果
%p0：第一主成分方向所在直线经过的点的下标索引，也是主成分空间的原点，为1*3的向量，该输入缺省时p0为ConCi的像素下标索引均值中心点

%Pcloser：ConCj的所有像素点中，距离ConCi的两个端点Ciendpoints_sub所在直线最近的点的1*3下标索引
%Qcloser：ConCi的两个端点Ciendpoints_sub中，距离ConCj较近的一个，是1*3的下标索引
%D_P2Q1Q2：连通分量ConCj中每一个像素点到ConCi的两个端点Ciendpoints_sub所在直线的距离，是一个n*1的数组，n是ConCj的像素总数
%D_Cc2Q1Q2：连通分量ConCj的中心（像素下标索引均值）到ConCi的两个端点Ciendpoints_sub所在直线的距离
%D_ConCj2Qcloser：连通分量ConCj中每一个像素点到ConCi的近端点Qcloser的距离，是一个n*1的数组，n是ConCj的像素总数
%iscolinear：iscolinear=1则表明，在ConCi的第一主成分所在直线方向上，ConCj完全位于ConCi的前方或者后方
%judge：judge=1则表明ConCj属于ConCi所在的电极

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
    ConCj_score=(ConCj_sub-repmat(p0,size(ConCj_sub,1),1))*ConCi_coeff; %注意此处是ConCi_coeff而不是ConCj_coeff！！！

    ConCi_score1_min=min(ConCi_score(:,1));
    ConCi_score1_max=max(ConCi_score(:,1));
    ConCj_score1=ConCj_score(:,1);

    N1=sum(ConCj_score1>ConCi_score1_max);
    N2=sum(ConCj_score1<ConCi_score1_min);
    iscolinear=0; %colinear:共线的
    if (N1/length(ConCj_score1)>0.9)||(N2/length(ConCj_score1)>0.9)
        iscolinear=1;
    end

    Q1=Ciendpoints_sub(1,:);
    Q2=Ciendpoints_sub(2,:);
    D1=norm(Q1-ConCjcenter_sub); %计算向量的2-范数，即欧式范数，即等于组成向量的两点的欧式距离
    D2=norm(Q2-ConCjcenter_sub);

    if D1<D2
        Qcloser=Q1;
    else
        Qcloser=Q2;
    end

    D_P2Q1Q2=zeros(R,1);
    for i=1:R
        P=ConCj_sub(i,:);
        D_P2Q1Q2(i)=norm(cross(Q2-Q1,P-Q1))/norm(Q2-Q1); %计算点P到Q1Q2所在直线的距离
    end
    D_P2Q1Q2_min=min(D_P2Q1Q2);
    Label=find(D_P2Q1Q2==D_P2Q1Q2_min,1);
    Pcloser=ConCj_sub(Label,:);
    D_Cc2Q1Q2=norm(cross(Q2-Q1,ConCjcenter_sub-Q1))/norm(Q2-Q1);

    %D=pdist2(X,Y)：D(i,j) corresponds to the pairwise distance between observation i in X and observation j in Y
    D_ConCj2Qcloser=pdist2(ConCj_sub,Qcloser,'euclidean');
    D_ConCj2Qcloser_min=min(D_ConCj2Qcloser);

    if (D_P2Q1Q2_min<2)&&(D_Cc2Q1Q2<6)&&(D_ConCj2Qcloser_min<15)&&(iscolinear==1)
        judge=1;
    end

end

end