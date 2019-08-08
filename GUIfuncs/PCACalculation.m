function IM_PCA=PCACalculation(IM,conn)

%IM���������ά����
%conn����ͨ��

%IM_PCA��Ϊn��3�У��к�����ͨ�������һ�£�ÿһ��Ϊһ����ͨ�����ĵĵ�һ���ɷֵ�����������ϵ����

IM=logical(IM);
[IM_L,IM_n]=bwlabeln(IM,conn);
IM_PCA=zeros(IM_n,3);

for i=1:IM_n
    econc=IM_L;
    econc(econc~=i)=0;
    index=find(econc);
    [x,y,z]=ind2sub(size(econc),index);
    econc_sub=[x,y,z];
    if size(econc_sub,1)<2 %��ͨ��������2�����ؼ�����2���۲�ֵʱ���޷�����PCA
        IM_PCA(i,:)=[NaN,NaN,NaN];
    else
        [econc_coeff,~,~,~]=pca(econc_sub);
        IM_PCA(i,:)=econc_coeff(:,1)';
    end
end

end