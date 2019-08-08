function [IM_PCA_coeff,IM_PCA_latentp]=PCACalculation2(IM,conn)

%IM���������ά����
%conn����ͨ��

%IM_PCA_coeff��Ϊn��3�У��к�����ͨ�������һ�£�ÿһ��Ϊһ����ͨ�����ĵ�һ���ɷֵ�����������ϵ����
%IM_PCA_latentp��Ϊn��3�У��к�����ͨ�������һ�£�ÿһ��Ϊһ����ͨ�����ĵ�һ���ڶ����������ɷֵĹ����ʰٷֱ�

IM=logical(IM);
[IM_L,IM_n]=bwlabeln(IM,conn);
IM_PCA_coeff=zeros(IM_n,3);
IM_PCA_latentp=zeros(IM_n,3);

for i=1:IM_n
    econc=IM_L;
    econc(econc~=i)=0;
    index=find(econc);
    [x,y,z]=ind2sub(size(econc),index);
    econc_sub=[x,y,z];
    %��ͨ��������2�����ؼ�����2���۲�ֵʱ���޷�����pca������2��3������ʱ���Լ��㵫coeffά��ֻ��3*2��latentֻ��2*1
    if size(econc_sub,1)<2
        IM_PCA_coeff(i,:)=[0,0,0];
        IM_PCA_latentp(i,:)=[0,0,0];
    else
        [econc_coeff,~,econc_latent,~]=pca(econc_sub);
        econc_latent_percent=100*econc_latent/sum(econc_latent);%��latent�ܺ�ͳһΪ100�����ڹ۲칱����
        if length(econc_latent)==2
            econc_latent_percent=[econc_latent_percent;0];
        end
        IM_PCA_coeff(i,:)=econc_coeff(:,1)';
        IM_PCA_latentp(i,:)=econc_latent_percent';
    end
end

end