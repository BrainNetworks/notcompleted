function [IML_PCAcoeff,IML_PCAlatentp]=PCACalculation3(IM_L,IM_n)

%IM_L：连通分量已经编好了号的三维数组
%IM_n：连通分量的最大编号，有可能大于max(max(max(IM_L)))

%IML_PCAcoeff：为n行3列，行号与连通分量编号一致，每一行为一个连通分量的的第一主成分的特征向量（系数）
%IML_PCAlatentp：为n行3列，行号与连通分量编号一致，每一行为一个连通分量的的第一、第二、第三主成分的贡献率百分比

IML_PCAcoeff=zeros(IM_n,3);
IML_PCAlatentp=zeros(IM_n,3);

for i=1:IM_n
    econc=IM_L;
    econc(econc~=i)=0;
    index=find(econc);
    [x,y,z]=ind2sub(size(econc),index);
    econc_sub=[x,y,z];
    %连通分量少于2个像素即少于2个观测值时，无法计算pca，等于2和3个像素时可以计算但coeff维度只有3*2，latent只有2*1
    if size(econc_sub,1)<2
        IML_PCAcoeff(i,:)=[0,0,0];
        IML_PCAlatentp(i,:)=[0,0,0];
    else
        [econc_coeff,~,econc_latent,~]=pca(econc_sub);
        econc_latent_percent=100*econc_latent/sum(econc_latent);%将latent总和统一为100，便于观察贡献率
        if length(econc_latent)==2
            econc_latent_percent=[econc_latent_percent;0];
        end
        IML_PCAcoeff(i,:)=econc_coeff(:,1)';
        IML_PCAlatentp(i,:)=econc_latent_percent';
    end
end

end