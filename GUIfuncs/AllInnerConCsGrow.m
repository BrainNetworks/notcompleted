function [ImageOut,InCs_stats,InCs_grow]=AllInnerConCsGrow(origin,basis,basis_S,elecenters_n,elecenters_sub,th,ratio,T0,Te,fignum)

%origin：用于阈值处理的初始图像，如electrode1
%basis：用于区域生长的输入三维数组
%basis_S：basis的所有连通分量的统计信息
%elecenters_n：电极入点elecenters的连通分量总数
%elecenters_sub：电极入点的所有连通分量的三维下标索引，是一个n行1列的元胞数组
%th：阈值下降程度，th=2表示相邻两次的阈值相差2
%ratio：表示颅内连通分量向颅骨边界降阈值生长直至其与位于颅骨边界上的电极入点产生交集，且交集所含像素总数占该电极入点的像素总数的ratio倍
%T0：阈值初值
%Te：阈值终值
%fignum：所画的第一张图片的序号

%ImageOut：basis区域生长后得到的近乎完整的电极结果，所有电极的颅内连通分量生长至边界，同时颅内残余干扰被去除
%InCs_stats：basis的所有颅内连通分量的统计信息，是n*4的矩阵，n为basis的颅内连通分量总数
%InCs_stats(:,1)是basis的颅内连通分量编号
%InCs_stats(:,2)是对basis的每个颅内连通分量是否属于电极的判据，值为1表明“是”，值为0表明“否”
%InCs_stats(:,3)是basis的每个颅内连通分量区域生长的阈值终值
%InCs_stats(:,4)是对basis的颅内连通分量区域生长得到的结果进行pca分析得到的第一主成分贡献率的百分比值
%InCs_grow：是一个4维数组，InCs_grow(:,:,:,i)对应第i个颅内连通分量InnerConCi区域生长得到的结果，即函数MyRegionGrow的输出result
%（1）如果第i个颅内连通分量InnerConCi属于某一电极，则InCs_grow(:,:,:,i)是InnerConCi区域生长得到的完整电极结果
%（2）如果第i个颅内连通分量InnerConCi不属于某一电极，而是颅内干扰，则InCs_grow(:,:,:,i)=InnerConCi

if nargin<10
    fignum=1;
end

if nargin<8
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<7
    ratio=0.05;
end

if nargin<6
    th=2;
end

InnerConC_L=[];
conn=26;
basis=logical(basis);
[basis_L,basis_n]=bwlabeln(basis,conn);
for j=1:basis_n
    if basis_S(j,2)==0
        InnerConC_L=[InnerConC_L;j];
    end
end

ICnum=length(InnerConC_L);
[nx,ny,nz]=size(basis);
InCs_grow=false(nx,ny,nz,ICnum);
InCs_judge=zeros(ICnum,1);
InCs_grow_t=zeros(ICnum,1);
InCs_grow_latentp1=zeros(ICnum,1);
InCs_stats=zeros(ICnum,4);
ImageOut=basis;

h=waitbar(0,'Region Growing...'); %创建进度条

for i=1:ICnum
    
    InnerConCi=basis_L;
    InnerConCi(InnerConCi~=InnerConC_L(i))=0;
    InnerConCi=logical(InnerConCi);
    InnerConCi_comp=basis_L;
    InnerConCi_comp(InnerConCi_comp==InnerConC_L(i))=0;
    InnerConCi_comp=logical(InnerConCi_comp);
    figure(fignum+2*(i-1));
    PlotElectrodes(InnerConCi_comp,'b');hold on;
    PlotElectrodes(InnerConCi,'r');hold on;

    [InCi_grow,InCi_judge,InCi_grow_t,InCi_grow_latentp1]=MyRegionGrow3(origin,InnerConCi,elecenters_n,elecenters_sub,th,ratio,T0,Te);
    figure(fignum+2*(i-1)+1);PlotElectrodes(InCi_grow,'b');
    InCs_grow(:,:,:,i)=InCi_grow;
    InCs_judge(i)=InCi_judge;
    InCs_grow_t(i)=InCi_grow_t;
    InCs_grow_latentp1(i)=InCi_grow_latentp1;
    if InCi_judge==1
        ImageOut=logical(ImageOut+InCi_grow);
    else
        InCi_grow_ind=find(InCi_grow);
        ImageOut(InCi_grow_ind)=0;
    end
    
    waitbar(i/ICnum,h); %修改进度显示
            
end

close(h);

InCs_stats(:,1)=InnerConC_L;
InCs_stats(:,2)=InCs_judge;
InCs_stats(:,3)=InCs_grow_t;
InCs_stats(:,4)=InCs_grow_latentp1;

end