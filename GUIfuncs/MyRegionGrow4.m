function result=MyRegionGrow4(origin,seed0,th,T0,Te)

%origin��������ֵ����ĳ�ʼͼ����electrode1
%seed0���������������ĳ�ʼ����
%th����ֵ�½��̶ȣ�th=2��ʾ�������ε���ֵ���2
%T0����ֵ��ֵ
%Te����ֵ��ֵ

%result�������������

if nargin<4
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<3
    th=2;
end

T=(T0-th):(-th):Te;
if T(end)~=Te
    T=[T,Te];
end
Tlen=length(T);
seed=logical(seed0);

for m=1:Tlen

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
    
end

result=seed;

end
