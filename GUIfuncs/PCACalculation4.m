 function [coeff,score,latent,tsquared,explained,mu]=PCACalculation4(IM)

%IM���������ά���飬������һ����ͨ��������������ͨ������Ϊһ������
%conn����ͨ��

%coeff��Ϊ������ͨ������pca���coeff
%score��Ϊ������ͨ������pca���score
%latent��Ϊ������ͨ������pca���latent
%tsquared��Ϊ������ͨ������pca���tsquared
%explained��Ϊ������ͨ������pca���explained
%mu��Ϊ������ͨ������pca���mu

IM=logical(IM);
econc=IM;
index=find(econc);
[x,y,z]=ind2sub(size(econc),index);
econc_sub=[x,y,z];
%��ͨ��������2�����ؼ�����2���۲�ֵʱ���޷�����pca������2��3������ʱ���Լ��㵫coeffά��ֻ��3*2��latentֻ��2��ֵ
if size(econc_sub,1)<2
    coeff=NaN;
    score=NaN;
    latent=NaN;
    tsquared=NaN;
    explained=NaN;
    mu=NaN;
else
    [econc_coeff,econc_score,econc_latent,econc_tsquared,econc_explained,econc_mu]=pca(econc_sub);
%     econc_latent_percent=100*econc_latent/sum(econc_latent);%��latent�ܺ�ͳһΪ100�����ڹ۲칱����
    %ע�⣺econc_latent_percent��econc_explained���
    coeff=econc_coeff;
    score=econc_score;
    latent=econc_latent;
    tsquared=econc_tsquared;
    explained=econc_explained;
    mu=econc_mu;
end

end