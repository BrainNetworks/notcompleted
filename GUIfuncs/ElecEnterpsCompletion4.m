function [InnerConC,label1,label2,ImageOut,CompStats]=ElecEnterpsCompletion4(ImageIn,origin,mask,ratio,th,T0,Te,fignum,conn)

%ImageIn：输入的电极数据，通常为所有颅内连通分量区域生长后的结果，即函数AllInnerConCsGrow的输出ImageOut
%origin：用于阈值处理的初始图像，如regCT,electrode1（常用）,electrode1B=regCT.*brainmask等
%mask：不包含颅骨边界的颅内区域掩模，如brainmask3，用于寻找完全位于颅骨边界上的连通分量
%ratio：表示输入ImageIn的颅骨边界上的小连通分量向颅内电极降阈值生长直至其与颅内电极产生交集，且交集所含像素总数占该电极的像素总数必须
%       大于ratio倍，满足这一条件则认为该颅骨边界上的小连通分量属于电极的一部分，即初步认为其是电极入点断开的连通分量之一
%th：阈值下降程度，th=2表示相邻两次的阈值相差2
%T0：阈值初值
%Te：阈值终值
%fignum：所画的第一张图片的序号
%conn：连通性

%InnerConC：将输入ImageIn的完全位于颅骨边界上的连通分量去除后剩下的部分
%label1：输入ImageIn中完全位于颅骨边界上的连通分量的编号，是列向量
%label2：编号在label1中的像素总数少于1500的连通分量的编号，是列向量，label2中的连通分量才可用于入点修补
%ImageOut：将输入ImageIn的颅骨边界上属于电极的但却与电极断开的小连通分量连接到电极上，同时去除颅骨边界上不属于电极的残余颅骨干扰后得到的结果
%CompStats：输入ImageIn的颅骨边界上的小连通分量（可能属于电极，也可能是残余颅骨干扰）seed0的统计信息，是一个n*4的矩阵，n是小连通分量的数目
%CompStats(:,1)是每个颅骨边界上的小连通分量seed0在输入ImageIn中的连通分量编号
%CompStats(:,2)是对输入ImageIn的每个颅骨边界上的小连通分量seed0是否属于电极的判据，值为1表明“是”，值为0表明“否”
%Compstats(:,3)是对输入ImageIn的每个颅骨边界上的小连通分量seed0是进行分类，Compstats(:,3)值相同的所有小连通分量被包含在编号等于该值的小连通分量的生长结果中，若值为0表示“不是电极”
%CompStats(:,4)是输入ImageIn的每个颅骨边界上的小连通分量seed0区域生长的阈值终值

if nargin<9
    conn=26;
end

if nargin<8
    fignum=1;
end

if nargin<6
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<5
    th=2;
end

if nargin<4
    ratio=0.001;
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
        IConCi_pn=length(find(IConCi));
%         IConCiog=logical(IConCi.*origin); %注意当origin=electrode1B=regCT.*brainmask，此处有IConCi.*origin
%         IConCiog_pn=length(find(IConCiog));
%         if (IConCiog_pn<400)&&(IConCi_pn<1500)
        if IConCi_pn<1500
            label2=[label2;i];
        end
    end
end
ll1=length(label1);
ll2=length(label2);
label2_renew=label2;

CompStats=zeros(ll2,4);
CompStats(:,1)=label2;
[nx,ny,nz]=size(ImageIn);

label1ConC=false(nx,ny,nz);
for i=1:ll1
    label1ConC(ImageIn_L==label1(i))=1;
end
label1ConC=logical(label1ConC);

InnerConC=logical(ImageIn-label1ConC);
[InnerConC_L,InnerConC_n,InnerConC_sub]=FindIndexOfAllConnectedComponents(InnerConC,conn);

InnerConC_allind=find(InnerConC);
[xi,yi,zi]=ind2sub(size(InnerConC),InnerConC_allind);
InnerConC_allsub=[xi,yi,zi];
            
label2ConC=false(nx,ny,nz);
for i=1:ll2
    label2ConC(ImageIn_L==label2(i))=1;
end
label2ConC=logical(label2ConC);

% figure(fignum);PlotElectrodes(label1ConC,'b');
% fignum=fignum+1;
% figure(fignum);PlotElectrodes(label2ConC,'b');
% fignum=fignum+1;
% figure(fignum);PlotElectrodes(InnerConC,'b');
% fignum=fignum+1;

T=(T0-th):(-th):Te;
if T(end)~=Te
    T=[T,Te];
end
Tlen=length(T);

ImageOut=ImageIn;
flag2=[];

% h=waitbar(0,'Electrode enter point repair...'); %创建进度条

for i=1:ll2
    
    k=label2(i);
    
    if isempty(intersect(k,flag2))

        flag=0;
        seed0=ImageIn_L;
        seed0(seed0~=k)=0;
        seed0=logical(seed0);
%         figure(fignum);PlotElectrodes(seed0,'b');
%         fignum=fignum+1;
        
        seed0_ind=find(seed0);
        [x0,y0,z0]=ind2sub(size(seed0),seed0_ind);
        seed0_sub=[x0,y0,z0];
        
        D=pdist2(seed0_sub,InnerConC_allsub,'euclidean');
        Dmin=min(min(D));
        if Dmin>=40
            ImageOut(seed0_ind)=0;
            flag2=[flag2,k];
        else

%             seed=logical(seed0.*origin); %注意当origin=electrode1B=regCT.*brainmask，此处seed0.*origin
            seed=logical(seed0); %注意当origin=electrode1=regCT.*brainmask2时，此处IConCi不能再点乘origin

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
                    if size(IS2,1)/size(InnerConC_sub{j},1)>ratio
                        InnerConCj=InnerConC_L;
                        InnerConCj(InnerConC_L~=j)=0;
                        InnerConCj=logical(InnerConCj);

                        [~,~,connect_neighb,~,~]=ConnectTwoConC3(InnerConCj,seed0);
                        ICcnS=logical(InnerConCj+connect_neighb+seed0);
%                         [~,~,~,~,ICcnS_explained,~]=PCACalculation4(ICcnS);

                        InnerConCjPAL=regionprops3(InnerConCj,'PrincipalAxisLength');
                        InnerConCj_PAL=InnerConCjPAL.PrincipalAxisLength;
                        ICcnSPAL=regionprops3(ICcnS,'PrincipalAxisLength');
                        ICcnS_PAL=ICcnSPAL.PrincipalAxisLength;
                        PALdiff=ICcnS_PAL(2)-InnerConCj_PAL(2);

%                         if (ICcnS_explained(1)>90)||(PALdiff<8)
                        if PALdiff<8
                            flag=1;
                            CompStats(i,2)=1;
                            CompStats(i,3)=k;
                            ImageOut=logical(ImageOut+connect_neighb);
%                             figure(fignum);PlotElectrodes(logical(InnerConCj+seed0),'b');hold on;
%                             figure(fignum);PlotElectrodes(logical(connect_neighb),'r');hold on;
%                             fignum=fignum+1;
%                             figure(fignum);PlotElectrodes(ImageOut,'b');
%                             fignum=fignum+1;
                        else
                            label2_renew=setdiff(label2_renew,k);
                        end
                        flag2=[flag2,k];

                        for ne=(i+1):ll2
                            ConCj=ImageIn_L;
                            ConCj(ConCj~=label2(ne))=0;
                            ConCj_ind=find(ConCj);
                            IS3=intersect(seed_ind2,ConCj_ind);
                            if length(IS3)/length(ConCj_ind)>0.8
                                CompStats(ne,2)=1;
                                CompStats(ne,3)=k;
                                CompStats(ne,4)=t;
                                [~,~,connect2_neighb,~,~]=ConnectTwoConC3(InnerConCj,ConCj);
                                ImageOut=logical(ImageOut+connect2_neighb);
                                flag2=[flag2,label2(ne)];

                                [~,~,connect2_neighb,~,~]=ConnectTwoConC3(InnerConCj,ConCj);
                                ICcnS2=logical(InnerConCj+connect2_neighb+ConCj);

                                ICcnS2PAL=regionprops3(ICcnS2,'PrincipalAxisLength');
                                ICcnS2_PAL=ICcnS2PAL.PrincipalAxisLength;
                                PALdiff2=ICcnS2_PAL(2)-InnerConCj_PAL(2);

                                if PALdiff2<8
                                    CompStats(ne,2)=1;
                                    CompStats(ne,3)=k;
                                    CompStats(ne,4)=t;
                                    ImageOut=logical(ImageOut+connect2_neighb);
                                    flag2=[flag2,label2(ne)];
                                end
                            end
                        end

                        break;
                    end
                end

                if flag==1
                    break;
                end

                seed_pn=length(find(seed));
                if (flag==0)&&(seed_pn>6000)
                    label2_renew=setdiff(label2_renew,k);
                    break;
                end

            end
            CompStats(i,4)=t;
            
        end

    else
        i=i+1;
    end
    
%     waitbar(i/ll2,h); %修改进度显示

end

% close(h);

label2=label2_renew;
ImageOut=logical(ImageOut);

end