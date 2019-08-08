 function [coeff,score,latent,tsquared,explained,mu]=PCACalculation4(IM)

%IM：输入的三维数组，仅含有一个连通分量，或将所有连通分量视为一个整体
%conn：连通性

%coeff：为输入连通分量的pca输出coeff
%score：为输入连通分量的pca输出score
%latent：为输入连通分量的pca输出latent
%tsquared：为输入连通分量的pca输出tsquared
%explained：为输入连通分量的pca输出explained
%mu：为输入连通分量的pca输出mu

IM=logical(IM);
econc=IM;
index=find(econc);
[x,y,z]=ind2sub(size(econc),index);
econc_sub=[x,y,z];
%连通分量少于2个像素即少于2个观测值时，无法计算pca，等于2和3个像素时可以计算但coeff维度只有3*2，latent只有2个值
if size(econc_sub,1)<2
    coeff=NaN;
    score=NaN;
    latent=NaN;
    tsquared=NaN;
    explained=NaN;
    mu=NaN;
else
    [econc_coeff,econc_score,econc_latent,econc_tsquared,econc_explained,econc_mu]=pca(econc_sub);
%     econc_latent_percent=100*econc_latent/sum(econc_latent);%将latent总和统一为100，便于观察贡献率
    %注意：econc_latent_percent与econc_explained相等
    coeff=econc_coeff;
    score=econc_score;
    latent=econc_latent;
    tsquared=econc_tsquared;
    explained=econc_explained;
    mu=econc_mu;
end

end