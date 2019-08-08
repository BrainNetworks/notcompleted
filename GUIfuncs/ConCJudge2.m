function [ConCi_coeff1,ConCi_explained,ConCj_coeff1,ConCj_explained,pcadiff,iscolinear,judge]=ConCJudge2(ConCj,ConCi,D_ConCj2Qcloser,p0)

%ConCj：某一三维数组的某一连通分量j，通常是待判别的连通分量，即判断ConCj是否属于ConCi所在的电极
%ConCi：某一三维数组的某一连通分量i,通常是某一电极或某一电极的一部分
%D_ConCj2Qcloser：连通分量ConCj中每一个像素点到ConCi的近端点Qcloser的距离，是一个n*1的数组，n是ConCj的像素总数
%p0：第一主成分方向所在直线经过的点的下标索引，也是主成分空间的原点，为1*3的向量，该输入缺省时p0为ConCi的像素下标索引均值中心点

%ConCi_coeff：连通分量ConCi进行pca运算得到的输出coeff，为3行3列
%ConCi_explained：连通分量ConCi进行pca运算得到的输出explained，为3行1列
%ConCj_coeff：连通分量ConCj进行pca运算得到的输出coeff，为3行3列
%ConCj_explained：连通分量ConCj进行pca运算得到的输出explained，为3行1列
%pcadiff：ConCi_coeff1(ConCi_coeff的第一列)和ConCj_coeff1(ConCj_coeff的第一列)之间的欧氏距离，用于衡量ConCi和ConCj的三个主成分方向的差异性
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

    [ConCi_coeff,~,~,~,ConCi_explained,~]=pca(ConCi_sub);
    [ConCj_coeff,~,~,~,ConCj_explained,~]=pca(ConCj_sub);
    ConCi_coeff1=ConCi_coeff(:,1);
    ConCj_coeff1=ConCj_coeff(:,1);
    pcadiff=norm(ConCi_coeff1-ConCj_coeff1);

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

    D_ConCj2Qcloser_min=min(D_ConCj2Qcloser);

    if (pcadiff<0.12)&&(ConCi_explained(1)>94)&&(ConCj_explained(1)>94)&&(D_ConCj2Qcloser_min<18)&&(iscolinear==1)
        judge=1;
    end

end

end