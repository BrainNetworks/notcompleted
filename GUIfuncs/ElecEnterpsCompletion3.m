function [InnerConC,label1,label2,ImageOut,CompStats,CompConCs]=ElecEnterpsCompletion3(ImageIn,origin,mask,ratio2,ratio3,th,T0,Te,fignum,conn)

%ImageIn：输入的电极数据，通常为所有颅内连通分量区域生长后的结果，即函数AllInnerConCsGrow的输出ImageOut
%origin：用于阈值处理的初始图像，如regCT,electrode1（常用）,electrode1B=regCT.*brainmask等
%mask：不包含颅骨边界的颅内区域掩模，如brainmask3，用于寻找完全位于颅骨边界上的连通分量
%ratio2：表示输入ImageIn的颅骨边界上的小连通分量向颅内电极降阈值生长直至其与颅内电极产生交集，且交集所含像素总数占该电极的像素总数必须
%        大于ratio2倍，一旦满足这一条件时则立刻保存此时的区域生长结果seed，假设命名为seed_ratio2，seed_ratio2就是电极入点的修补结果
%ratio3：表示输入ImageIn的颅骨边界上的小连通分量向颅内电极降阈值生长直至其与颅内电极产生交集，且交集所含像素总数占该电极的像素总数必须
%        大于ratio3倍，满足这一条件则认为该颅骨边界上的小连通分量属于电极的一部分，即认为其是电极入点断开的连通分量之一
%th：阈值下降程度，th=2表示相邻两次的阈值相差2
%T0：阈值初值
%Te：阈值终值
%fignum：所画的第一张图片的序号
%conn：连通性

%InnerConC：将输入ImageIn的完全位于颅骨边界上的连通分量去除后剩下的部分
%label1：输入ImageIn中完全位于颅骨边界上的连通分量的编号，是列向量
%label2：编号在label1中的像素总数少于500并且其被origin所包含的部分的像素数也小于100的连通分量的编号，是列向量，label2中的连通分量才可用于入点修补
%ImageOut：将输入ImageIn的颅骨边界上属于电极的但却与电极断开的小连通分量连接到电极上，同时去除颅骨边界上不属于电极的残余颅骨干扰后得到的结果
%CompStats：输入ImageIn的颅骨边界上的小连通分量（可能属于电极，也可能是残余颅骨干扰）seed0的统计信息，是一个n*4的矩阵，n是小连通分量的数目
%CompStats(:,1)是每个颅骨边界上的小连通分量seed0在输入ImageIn中的连通分量编号
%CompStats(:,2)是对输入ImageIn的每个颅骨边界上的小连通分量seed0是否属于电极的判据，值为1表明“是”，值为0表明“否”
%Compstats(:,3)是对输入ImageIn的每个颅骨边界上的小连通分量seed0是进行分类，Compstats(:,3)值相同的所有小连通分量被包含在编号等于该值的小连通分量的生长结果中，若值为0表示“不是电极”
%CompStats(:,4)是输入ImageIn的每个颅骨边界上的小连通分量seed0区域生长的阈值终值
%CompConCs：输入ImageIn的颅骨边界上的小连通分量seed0的修补结果，是一个4维数组
%（1）如果第i个颅骨边界上的小连通分量seed0属于某一颅内电极，则CompConCs(:,:,:,i)等于向颅内区域生长直至大于ratio2时得到的结果seed_ratio2
%（2）如果第i个颅骨边界上的小连通分量seed0不属于某一颅内电极，而是残余颅骨干扰，则CompConCs(:,:,:,i)=seed0

if nargin<10
    conn=26;
end

if nargin<9
    fignum=1;
end

if nargin<7
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<6
    th=2;
end

if nargin<5
    ratio3=0.003;
end

if nargin<4
    ratio2=0.0003;
end

ImageIn=logical(ImageIn);
[ImageIn_L,ImageIn_n]=bwlabeln(ImageIn,conn);
%ImageIn_pn：输入ImageIn的所有连通分量的像素总数
% ImageIn_pn=zeros(ImageIn_n,1);
% for i=1:ImageIn_n
%     ImageIn_pn(i)=sum(sum(sum(ImageIn_L==i)));
% end
% ImageIn_pns=sort(ImageIn_pn,'ascend');
% ImageIn_pns_mean=mean(ImageIn_pns(2:(end-1)));

label1=[];
label2=[];
for i=1:ImageIn_n
    IConCi=ImageIn_L;
    IConCi(ImageIn_L~=i)=0;
    IConCi=logical(IConCi);
    IS1=IConCi.*logical(mask);
    if isempty(find(IS1))
        label1=[label1;i];
        IConCiog=logical(IConCi.*origin); %注意此处IConCi.*origin，origin=electrode1B=regCT.*brainmask
        IConCiog_pn=length(find(IConCiog));
        IConCi_pn=length(find(IConCi));
        if (IConCiog_pn<200)&&(IConCi_pn<600)
            label2=[label2;i];
        end
    end
end
ll1=length(label1);
ll2=length(label2);

CompStats=zeros(ll2,4);
CompStats(:,1)=label2;
[nx,ny,nz]=size(ImageIn);
CompConCs=false(nx,ny,nz,ll2);

label1ConC=false(nx,ny,nz);
for i=1:ll1
    label1ConC(ImageIn_L==label1(i))=1;
end
label1ConC=logical(label1ConC);

InnerConC=logical(ImageIn-label1ConC);
[~,InnerConC_n,InnerConC_sub]=FindIndexOfAllConnectedComponents(InnerConC,conn);

label2ConC=false(nx,ny,nz);
for i=1:ll2
    label2ConC(ImageIn_L==label2(i))=1;
end
label2ConC=logical(label2ConC);

figure(fignum);PlotElectrodes(label1ConC,'b');
fignum=fignum+1;
figure(fignum);PlotElectrodes(label2ConC,'b');
fignum=fignum+1;
figure(fignum);PlotElectrodes(InnerConC,'b');
fignum=fignum+1;

T=(T0-th):(-th):Te;
if T(end)~=Te
    T=[T,Te];
end
Tlen=length(T);

ImageOut=ImageIn;
flag3=[];

h=waitbar(0,'Electrodes enter points Completing...'); %创建进度条

for i=1:ll2
    
    k=label2(i);
    
    if isempty(intersect(k,flag3))

        flag=0;
        flag2=0;
        seed0=ImageIn_L;
        seed0(seed0~=k)=0;
        seed0=logical(seed0);
        seed=logical(seed0.*origin); %注意此处seed0.*origin
        figure(fignum);PlotElectrodes(seed,'b');
        fignum=fignum+1;

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
            seed_ind2=find(seed);
            [x2,y2,z2]=ind2sub(size(seed),seed_ind2);
            seed_sub2=[x2,y2,z2];
            
            for j=1:InnerConC_n
                IS2=intersect(seed_sub2,InnerConC_sub{j},'rows');
                if (flag2==0)&&(size(IS2,1)/size(InnerConC_sub{j},1)>ratio2)
                    seed_grow=seed;
                    flag2=1;
                    figure(fignum);PlotElectrodes(seed_grow,'b');
                    fignum=fignum+1;
                end
                if size(IS2,1)/size(InnerConC_sub{j},1)>ratio3
                    flag=1;
                    flag3=[flag3,k];
                    CompStats(i,2)=1;
                    CompStats(i,3)=k;
                    CompConCs(:,:,:,i)=seed_grow;
                    ImageOut=logical(ImageOut+seed_grow);
                    figure(fignum);PlotElectrodes(seed,'b');
                    fignum=fignum+1;
                    
                    for ne=(i+1):ll2
                        ConCj=ImageIn_L;
                        ConCj(ConCj~=label2(ne))=0;
                        ConCj_ind=find(ConCj);
                        IS3=intersect(seed_ind2,ConCj_ind);
                        if length(IS3)/length(ConCj_ind)>0.8
                            CompStats(ne,2)=1;
                            CompStats(ne,3)=k;
                            CompStats(ne,4)=t;
                            CompConCs(:,:,:,ne)=seed_grow;
                            flag3=[flag3,label2(ne)];
                        end
                    end
                    
                    break;
                end
            end

            if flag==1
                break;
            end
            
            seed_pn=length(find(seed));
            if (flag==0)&&(seed_pn>1300)
                CompConCs(:,:,:,i)=seed0;
                seed0_ind=find(seed0);
                ImageOut(seed0_ind)=0;
                break;
            end

        end
        CompStats(i,4)=t;

    else
        i=i+1;
    end
    
    waitbar(i/ll2,h); %修改进度显示

end

close(h);

CompConCs=logical(CompConCs);
ImageOut=logical(ImageOut);

end