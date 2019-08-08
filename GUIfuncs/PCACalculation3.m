function [IML_PCAcoeff,IML_PCAlatentp]=PCACalculation3(IM_L,IM_n)

%IM_L����ͨ�����Ѿ�����˺ŵ���ά����
%IM_n����ͨ����������ţ��п��ܴ���max(max(max(IM_L)))

%IML_PCAcoeff��Ϊn��3�У��к�����ͨ�������һ�£�ÿһ��Ϊһ����ͨ�����ĵĵ�һ���ɷֵ�����������ϵ����
%IML_PCAlatentp��Ϊn��3�У��к�����ͨ�������һ�£�ÿһ��Ϊһ����ͨ�����ĵĵ�һ���ڶ����������ɷֵĹ����ʰٷֱ�

IML_PCAcoeff=zeros(IM_n,3);
IML_PCAlatentp=zeros(IM_n,3);

for i=1:IM_n
    econc=IM_L;
    econc(econc~=i)=0;
    index=find(econc);
    [x,y,z]=ind2sub(size(econc),index);
    econc_sub=[x,y,z];
    %��ͨ��������2�����ؼ�����2���۲�ֵʱ���޷�����pca������2��3������ʱ���Լ��㵫coeffά��ֻ��3*2��latentֻ��2*1
    if size(econc_sub,1)<2
        IML_PCAcoeff(i,:)=[0,0,0];
        IML_PCAlatentp(i,:)=[0,0,0];
    else
        [econc_coeff,~,econc_latent,~]=pca(econc_sub);
        econc_latent_percent=100*econc_latent/sum(econc_latent);%��latent�ܺ�ͳһΪ100�����ڹ۲칱����
        if length(econc_latent)==2
            econc_latent_percent=[econc_latent_percent;0];
        end
        IML_PCAcoeff(i,:)=econc_coeff(:,1)';
        IML_PCAlatentp(i,:)=econc_latent_percent';
    end
end

end