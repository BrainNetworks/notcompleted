 function [coeff,score,latent,tsquared,explained,mu]=PCACalculation5(IM,conn)

%IM：输入的三维数组
%conn：连通性

%coeff：为n行1列的元胞数组，行号与连通分量编号一致，n为输入IM的连通分量总数，每一个元胞数组为一个连通分量的pca输出coeff
%score：为n行1列的元胞数组，行号与连通分量编号一致，n为输入IM的连通分量总数，每一个元胞数组为一个连通分量的pca输出score
%latent：为n行1列的元胞数组，行号与连通分量编号一致，n为输入IM的连通分量总数，每一个元胞数组为一个连通分量的pca输出latent
%tsquared：为n行1列的元胞数组，行号与连通分量编号一致，n为输入IM的连通分量总数，每一个元胞数组为一个连通分量的pca输出tsquared
%explained：为n行1列的元胞数组，行号与连通分量编号一致，n为输入IM的连通分量总数，每一个元胞数组为一个连通分量的pca输出explained
%mu：为n行1列的元胞数组，行号与连通分量编号一致，n为输入IM的连通分量总数，每一个元胞数组为一个连通分量的pca输出mu

IM=logical(IM);
[IM_L,IM_n]=bwlabeln(IM,conn);
coeff=cell(IM_n,1);
score=cell(IM_n,1);
latent=cell(IM_n,1);
tsquared=cell(IM_n,1);
explained=cell(IM_n,1);
mu=cell(IM_n,1);

for i=1:IM_n
    econc=IM_L;
    econc(econc~=i)=0;
    index=find(econc);
    [x,y,z]=ind2sub(size(econc),index);
    econc_sub=[x,y,z];
    %连通分量少于2个像素即少于2个观测值时，无法计算pca，等于2和3个像素时可以计算但coeff维度只有3*2，latent只有2个值
    if size(econc_sub,1)<2
        coeff{i}=NaN;
        score{i}=NaN;
        latent{i}=NaN;
        tsquared{i}=NaN;
        explained{i}=NaN;
        mu{i}=NaN;
    else
        [econc_coeff,econc_score,econc_latent,econc_tsquared,econc_explained,econc_mu]=pca(econc_sub);
%       econc_latent_percent=100*econc_latent/sum(econc_latent);%将latent总和统一为100，便于观察贡献率
        %注意：econc_latent_percent与econc_explained相等
        coeff{i}=econc_coeff;
        score{i}=econc_score;
        latent{i}=econc_latent;
        tsquared{i}=econc_tsquared;
        explained{i}=econc_explained;
        mu{i}=econc_mu;
    end
end

end