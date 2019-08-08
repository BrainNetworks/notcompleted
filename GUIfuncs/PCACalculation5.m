 function [coeff,score,latent,tsquared,explained,mu]=PCACalculation5(IM,conn)

%IM���������ά����
%conn����ͨ��

%coeff��Ϊn��1�е�Ԫ�����飬�к�����ͨ�������һ�£�nΪ����IM����ͨ����������ÿһ��Ԫ������Ϊһ����ͨ������pca���coeff
%score��Ϊn��1�е�Ԫ�����飬�к�����ͨ�������һ�£�nΪ����IM����ͨ����������ÿһ��Ԫ������Ϊһ����ͨ������pca���score
%latent��Ϊn��1�е�Ԫ�����飬�к�����ͨ�������һ�£�nΪ����IM����ͨ����������ÿһ��Ԫ������Ϊһ����ͨ������pca���latent
%tsquared��Ϊn��1�е�Ԫ�����飬�к�����ͨ�������һ�£�nΪ����IM����ͨ����������ÿһ��Ԫ������Ϊһ����ͨ������pca���tsquared
%explained��Ϊn��1�е�Ԫ�����飬�к�����ͨ�������һ�£�nΪ����IM����ͨ����������ÿһ��Ԫ������Ϊһ����ͨ������pca���explained
%mu��Ϊn��1�е�Ԫ�����飬�к�����ͨ�������һ�£�nΪ����IM����ͨ����������ÿһ��Ԫ������Ϊһ����ͨ������pca���mu

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
    %��ͨ��������2�����ؼ�����2���۲�ֵʱ���޷�����pca������2��3������ʱ���Լ��㵫coeffά��ֻ��3*2��latentֻ��2��ֵ
    if size(econc_sub,1)<2
        coeff{i}=NaN;
        score{i}=NaN;
        latent{i}=NaN;
        tsquared{i}=NaN;
        explained{i}=NaN;
        mu{i}=NaN;
    else
        [econc_coeff,econc_score,econc_latent,econc_tsquared,econc_explained,econc_mu]=pca(econc_sub);
%       econc_latent_percent=100*econc_latent/sum(econc_latent);%��latent�ܺ�ͳһΪ100�����ڹ۲칱����
        %ע�⣺econc_latent_percent��econc_explained���
        coeff{i}=econc_coeff;
        score{i}=econc_score;
        latent{i}=econc_latent;
        tsquared{i}=econc_tsquared;
        explained{i}=econc_explained;
        mu{i}=econc_mu;
    end
end

end