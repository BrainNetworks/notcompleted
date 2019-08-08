function [result,judge,t,latentp1]=MyRegionGrow3(origin,seed0,elecenters_n,elecenters_sub,th,ratio,T0,Te)

%origin：用于阈值处理的初始图像，如electrode1
%seed0：用于区域生长的初始种子，如某一颅内连通分量
%elecenters_n：电极入点elecenters的连通分量总数
%elecenters_sub：电极入点的所有连通分量的三维下标索引，是一个n行1列的元胞数组
%th：阈值下降程度，th=2表示相邻两次的阈值相差2
%ratio：表示颅内连通分量seed0向颅骨边界降阈值生长直至其与位于颅骨边界上的电极入点产生交集，该交集所含像素总数占该电极入点的像素总数必须大于ratio倍
%T0：阈值初值
%Te：阈值终值

%result：（1）如果seed0是某一电极的颅内部分，则result是区域生长得到的完整电极结果
%        （2）如果seed0不属于某一电极，而是颅内干扰，则result=seed0
%judge：如果seed0属于某一电极，则judge=1；如果seed0不属于某一电极，则judge=0
%t：区域生长的阈值终值
%latentp1：对seed0区域生长得到的结果进行pca分析得到的第一主成分贡献率的百分比值

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
    seed_latent_percent=100*seed_latent/sum(seed_latent);%将latent总和统一为100，便于观察贡献率
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
    result=logical(seed_last); %seed_last可能与颅骨边界没有交集，这种情况可以在电极入点修补时纠正
end

end