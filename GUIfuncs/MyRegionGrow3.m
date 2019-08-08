function [result,judge,t,latentp1]=MyRegionGrow3(origin,seed0,elecenters_n,elecenters_sub,th,ratio,T0,Te)

%origin��������ֵ����ĳ�ʼͼ����electrode1
%seed0���������������ĳ�ʼ���ӣ���ĳһ­����ͨ����
%elecenters_n���缫���elecenters����ͨ��������
%elecenters_sub���缫����������ͨ��������ά�±���������һ��n��1�е�Ԫ������
%th����ֵ�½��̶ȣ�th=2��ʾ�������ε���ֵ���2
%ratio����ʾ­����ͨ����seed0��­�Ǳ߽罵��ֵ����ֱ������λ��­�Ǳ߽��ϵĵ缫�������������ý���������������ռ�õ缫�������������������ratio��
%T0����ֵ��ֵ
%Te����ֵ��ֵ

%result����1�����seed0��ĳһ�缫��­�ڲ��֣���result�����������õ��������缫���
%        ��2�����seed0������ĳһ�缫������­�ڸ��ţ���result=seed0
%judge�����seed0����ĳһ�缫����judge=1�����seed0������ĳһ�缫����judge=0
%t��������������ֵ��ֵ
%latentp1����seed0���������õ��Ľ������pca�����õ��ĵ�һ���ɷֹ����ʵİٷֱ�ֵ

if nargin<7
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<6
    ratio=0.05;
end

if nargin<5
    th=2;
end

T=(T0-th):(-th):Te;
if T(end)~=Te
    T=[T,Te];
end
Tlen=length(T);
flag1=0;
flag3=0;
seed=logical(seed0);

for m=1:Tlen

    seed_last=seed;
    
    seed_ind=find(seed);
    [x,y,z]=ind2sub(size(seed),seed_ind);
    seed_sub=[x,y,z];
    seed_neighb_sub=[x-1,y,z;
                     x+1,y,z;
                     x,y-1,z;
                     x,y+1,z;
                     x,y,z-1;
                     x,y,z+1];
    seed_neighb_sub=unique(seed_neighb_sub,'rows');
    seed_neighb_sub=setdiff(seed_neighb_sub,seed_sub,'rows');
    seed_neighb_ind=sub2ind(size(seed),seed_neighb_sub(:,1),seed_neighb_sub(:,2),seed_neighb_sub(:,3));
    
    seed2=false(size(seed));
    seed2(seed_ind)=1;
    seed2(seed_neighb_ind)=1; 
    t=T(m);
    for l=1:length(seed_neighb_ind)
        if origin(seed_neighb_ind(l))<t
            seed2(seed_neighb_ind(l))=0;
        end
    end
    seed=logical(seed2);
    seed_ind2=find(seed);
    [x2,y2,z2]=ind2sub(size(seed),seed_ind2);
    seed_sub2=[x2,y2,z2];
    
    ist=0;
    for j=1:elecenters_n
        IS=intersect(seed_sub2,elecenters_sub{j},'rows');
        if ~isempty(IS)
            ist=ist+1;
        end        
    end
    
    if ist>1
        flag3=1;
        break;
    end
    
    for j=1:elecenters_n
        IS=intersect(seed_sub2,elecenters_sub{j},'rows');
        if size(IS,1)/size(elecenters_sub{j},1)>ratio
            flag1=1;
            break;
        end
    end
    
    if flag1==1
        break;
    end
    
    seed_pn=length(find(seed));
    if (flag1==0)&&(seed_pn>8000)
        break;
    end
    
end

if size(seed_sub2,1)<2
    latentp1=NaN;
else
    [~,~,seed_latent,~]=pca(seed_sub2);
    seed_latent_percent=100*seed_latent/sum(seed_latent);%��latent�ܺ�ͳһΪ100�����ڹ۲칱����
    latentp1=seed_latent_percent(1);
end
    
if flag1==1
    judge=1;
    result=logical(seed);
else 
    judge=0;
    result=logical(seed0);
end

if flag3==1
    judge=1;
    result=logical(seed_last); %seed_last������­�Ǳ߽�û�н�����������������ڵ缫����޲�ʱ����
end

end