function IM_PCA=PCACalculation(IM,conn)

%IM：输入的三维数组
%conn：连通性

%IM_PCA：为n行3列，行号与连通分量编号一致，每一行为一个连通分量的的第一主成分的特征向量（系数）

IM=logical(IM);
[IM_L,IM_n]=bwlabeln(IM,conn);
IM_PCA=zeros(IM_n,3);

for i=1:IM_n
    econc=IM_L;
    econc(econc~=i)=0;
    index=find(econc);
    [x,y,z]=ind2sub(size(econc),index);
    econc_sub=[x,y,z];
    if size(econc_sub,1)<2 %连通分量少于2个像素即少于2个观测值时，无法计算PCA
        IM_PCA(i,:)=[NaN,NaN,NaN];
    else
        [econc_coeff,~,~,~]=pca(econc_sub);
        IM_PCA(i,:)=econc_coeff(:,1)';
    end
end

end