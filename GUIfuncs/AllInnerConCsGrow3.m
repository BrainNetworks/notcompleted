function [CranialConCs,InnerConCs,GrowProcess,ImageOut]=AllInnerConCsGrow3(origin,basis,elecenters,th,ratio,T0,Te,fignum)

%origin：用于阈值处理的初始图像，如electrode1
%basis：用于区域生长的输入三维数组
%elecenters：某一三维数组（通常是basis）与颅骨边界的交集部分，称为电极入点
%th：阈值下降程度，th=2表示相邻两次的阈值相差2
%ratio：表示当颅内连通分量向颅骨边界降阈值生长直至其与位于颅骨边界上的电极入点产生交集，且交集所含像素总数占该电极入点的像素总数的
%       ratio倍时，该颅内连通分量并入颅骨边界上的电极入点作为一个新的电极入点；当两个或两个以上颅内连通分量在生长的过程中连接在了
%       一起时，把这些颅内连通分量合并作为一个新的颅内生长种子点，直至该新的颅内生长种子点也与颅骨边界上的电极入点产生交集。
%T0：阈值初值
%Te：阈值终值
%fignum：所画的第一张图片的序号

%CranialConCs：basis的完全位于颅骨边界上或与颅骨边界存在交集的连通分量，称为颅骨连通分量
%InnerConCs：basis的完全位于颅内的与颅骨边界不存在交集的连通分量，称为颅内连通分量
%GrowProcess：记录区域生长过程，是N*4的元胞数组，每一行分别是当前更新的电极入点、当前更新的颅内生长种子、该次区域生长的阈值初值、该次区域生长的阈值终值
%ImageOut：basis区域生长后得到的近乎完整的电极结果，绝大部分电极的颅内连通分量生长至边界

if nargin<8
    fignum=1;
end

if nargin<6
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<4
    th=2;
end

% [basisConC_Sub,elecentersConC_Sub,basis_L,elecenters_L,elecenters_newL,elecenters_recL,basis_n,elecenters_n]=CountElectrodes(basis,elecenters,26);
% basis_recL=(1:1:basis_n)';
% InnerConCs_recL=setdiff(basis_recL,elecenters_recL);
% 
% InnerConCs=basis_L;
% for i=1:length(elecenters_recL)
%     L=elecenters_recL(i);
%     InnerConCs(InnerConCs==L)=0;
% end

basis=logical(basis);
elecenters=logical(elecenters);
[nx,ny,nz]=size(basis);
InnerConCs=false(nx,ny,nz);

conn=26;
[basis_L,basis_n]=bwlabeln(basis,conn);

for i=1:basis_n
    bconci=basis_L;
    bconci(bconci~=i)=0;
    IS=logical(bconci.*elecenters);
    if sum(IS==1)==0
        InnerConCs=logical(InnerConCs+bconci);
    end
end
CranialConCs=logical(basis-InnerConCs);
figure(fignum);PlotElectrodes(CranialConCs,'b');
figure(fignum+1);PlotElectrodes(InnerConCs,'b');
fignum=fignum+2;

seed=InnerConCs;
crconc=CranialConCs;

[~,crconc_n,crconc_sub]=FindIndexOfAllConnectedComponents(crconc,conn);
time=4;
T1=T0-th;
T2=T1-time*th;
N=0;

while T2>=Te
    
    % 注：result=MyRegionGrow4(origin,seed0,th,T1,T2)，若单次生长则T1和T2相等
    result=MyRegionGrow4(origin,seed,th,T1,T2);
    seed2=result;
    [seed2_L,seed2_n]=bwlabeln(seed2,conn);
    
    for i=1:seed2_n
        seed2conci=seed2_L;
        seed2conci(seed2conci~=i)=0;
        seed2conci_ind=find(seed2conci);
        [x,y,z]=ind2sub(size(seed2conci),seed2conci_ind);
        seed2conci_sub=[x,y,z];
        for j=1:crconc_n
            IS=intersect(seed2conci_sub,crconc_sub{j},'rows');
            if size(IS,1)/size(crconc_sub{j},1)>ratio
                N=N+1;
                seed2(seed2conci_ind)=0;
                crconc=logical(crconc+seed2conci);
                [~,crconc_n,crconc_sub]=FindIndexOfAllConnectedComponents(crconc,conn);
                figure(fignum);PlotElectrodes(crconc,'b');
                figure(fignum+1);PlotElectrodes(seed2,'b');
                fignum=fignum+2;
                GrowProcess{N,1}=crconc;
                GrowProcess{N,2}=seed2;
                GrowProcess{N,3}=T1;
                GrowProcess{N,4}=T2;
                break;
            end
        end
    end
    seed=seed2;
    [~,seed_n]=bwlabeln(seed,conn);
    if seed_n==0
        break;
    end
    
    T1=T2-th;
    T2=T1-time*th;
    if T2<Te
        T2=Te;
    end
    
end

ImageOut=logical(crconc+basis);

end