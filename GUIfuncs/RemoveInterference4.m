function [VolumeIn_explained,VIintc_explained,VolumeIn_PAL,VIintc_PAL,...,
    VIintc_s1ratio,VIc_s1ratio,VolumeOut]=RemoveInterference4(VolumeIn,Mask,conn)

%VolumeIn：输入的三维体数据
%Mask：掩模，也是三维数组，维度与VolumeIn相同，如brainmask,brainmask3等
%conn：连通性

%VolumeIn_explained：为n*3的二维数组，n是VolumeIn的连通分量总数，每一行为VolumeIn的一个连通分量的pca输出explained
%VIintc_explained：为n*3的二维数组，n是VolumeIn的连通分量总数，每一行为VolumeIn的一个连通分量的颅内部分VIintc的pca输出explained
%VolumeIn_PAL：为n*3的二维数组，n是VolumeIn的连通分量总数，每一行为VolumeIn的一个连通分量在其三个主成分空间坐标轴的投影长度
%VIintc_PAL：为n*3的二维数组，n是VolumeIn的连通分量总数，每一行为VolumeIn的一个连通分量的颅内部分VIintc在其三个主成分空间坐标轴的投影长度
%设VolumeIn的每一个连通分量VI的颅内部分为VIintc，颅骨边界部分为VIc，VIintc在VI的第一主成分坐标轴的投影为L1，VIc在VI的第一主成分坐标轴的投影为L2，设L1和L2的交集为L1L2IS，则：
%VIintc_s1ratio：为n*1的列向量，n是VolumeIn的连通分量总数，每一行为VolumeIn的一个连通分量的L1L2IS长度与L1长度的比值，值域为[0,1]
%VIc_s1ratio：为n*1的列向量，n是VolumeIn的连通分量总数，每一行为VolumeIn的一个连通分量的L1L2IS长度与L2长度的比值，值域为[0,1]
%VolumeOut：将在RemoveInterference3中未能去除的与颅骨边界存在交集并且其颅内部分像素总数不少于30个的残余干扰去除后得到的输出三维体数据

if nargin<3
    conn=26;
end

VolumeIn=logical(VolumeIn);
Mask=logical(Mask);

[VolumeIn_L,VolumeIn_n]=bwlabeln(VolumeIn,conn);
% VolumeInConC_ID=(1:VolumeIn_n)';

VolumeIn_explained=zeros(VolumeIn_n,3);
VIintc_explained=zeros(VolumeIn_n,3);
VolumeIn_PAL=zeros(VolumeIn_n,3);
VIintc_PAL=zeros(VolumeIn_n,3);
VIintc_s1ratio=zeros(VolumeIn_n,1);
VIc_s1ratio=zeros(VolumeIn_n,1);

for i=1:VolumeIn_n
    
    VIconci=VolumeIn_L;
    VIconci(VIconci~=i)=0;
    VIconci=logical(VIconci);
    VIciIntc=logical(VIconci.*Mask);
    VIciC=logical(VIconci.*(~Mask));
    
    [VIciIntc_L,VIciIntc_n]=bwlabeln(VIciIntc,conn);
    if VIciIntc_n>1
        VIciIntc_PN=zeros(VIciIntc_n,1);
        for j=1:VIciIntc_n
            VIciIntcj=VIciIntc_L;
            VIciIntcj(VIciIntc_L~=j)=0;
            VIciIntcj=logical(VIciIntcj);
            VIciIntc_PN(j)=length(find(VIciIntcj));
        end
        L=find(VIciIntc_PN==max(VIciIntc_PN));
        VIciIntc(VIciIntc_L~=L)=0;
    end
    
    [VIciC_L,VIciC_n]=bwlabeln(VIciC,conn);
    if VIciC_n>1
        VIciC_PN=zeros(VIciC_n,1);
        for j=1:VIciC_n
            VIciCj=VIciC_L;
            VIciCj(VIciC_L~=j)=0;
            VIciCj=logical(VIciCj);
            VIciC_PN(j)=length(find(VIciCj));
        end
        L=find(VIciC_PN==max(VIciC_PN));
        VIciC(VIciC_L~=L)=0;
    end
    
    VIconci_ind=find(VIconci);
    [x1,y1,z1]=ind2sub(size(VIconci),VIconci_ind);
    VIconci_sub=[x1,y1,z1];
    VIconci_ctrsub=mean(VIconci_sub,1);
    
    VIciIntc_ind=find(VIciIntc);
    [x2,y2,z2]=ind2sub(size(VIciIntc),VIciIntc_ind);
    VIciIntc_sub=[x2,y2,z2];

    VIciC_ind=find(VIciC);
    [x2,y2,z2]=ind2sub(size(VIciC),VIciC_ind);
    VIciC_sub=[x2,y2,z2];
    
    [VIconci_coeff,~,~,~,VIconci_explained,~]=pca(VIconci_sub);
    [~,~,~,~,VIciIntc_explained,~]=pca(VIciIntc_sub);
    VolumeIn_explained(i,:)=VIconci_explained;
    VIintc_explained(i,:)=VIciIntc_explained;
    
    VIconci_PAL=regionprops3(VIconci,'PrincipalAxisLength');
    VolumeIn_PAL(i,:)=VIconci_PAL.PrincipalAxisLength;
    VIciIntc_PAL=regionprops3(VIciIntc,'PrincipalAxisLength');
    VIintc_PAL(i,:)=VIciIntc_PAL.PrincipalAxisLength;
    
    VIciIntc_score=(VIciIntc_sub-repmat(VIconci_ctrsub,size(VIciIntc_sub,1),1))*VIconci_coeff;
    VIciC_score=(VIciC_sub-repmat(VIconci_ctrsub,size(VIciC_sub,1),1))*VIconci_coeff;
    VIciIntc_s1min=min(VIciIntc_score(:,1));
    VIciIntc_s1max=max(VIciIntc_score(:,1));
    VIciC_s1min=min(VIciC_score(:,1));
    VIciC_s1max=max(VIciC_score(:,1));
    VIciIntc_s1range=round((VIciIntc_s1min*100)):1:round((VIciIntc_s1max*100));
    VIciC_s1range=round((VIciC_s1min*100)):1:round((VIciC_s1max*100));
    VIciIC_s1IS=intersect(VIciIntc_s1range,VIciC_s1range);
    VIciIntc_s1ratio=length(VIciIC_s1IS)/length(VIciIntc_s1range);
    VIciC_s1ratio=length(VIciIC_s1IS)/length(VIciC_s1range);
    
    VIintc_s1ratio(i)=VIciIntc_s1ratio;
    VIc_s1ratio(i)=VIciC_s1ratio;

end

VolumeOut=VolumeIn;
judge1=0;
judge2=0;
judge3=0;
for i=1:VolumeIn_n
    if (VolumeIn_explained(i,1)<85)||(VIintc_explained(i,1)<85)
        VolumeOut(VolumeIn_L==i)=0;
        judge1=judge1+1;
    else
        if (VolumeIn_PAL(i,1)<32)||(VIintc_PAL(i,1)<25)
            VolumeOut(VolumeIn_L==i)=0;
            judge2=judge2+1;
        else
%             if (VIintc_s1ratio(i)>0.35)||(VIc_s1ratio(i)>0.85)
            if VIintc_s1ratio(i)>0.35
                VolumeOut(VolumeIn_L==i)=0;
                judge3=judge3+1;
            end
        end
    end
end
VolumeOut=logical(VolumeOut);

end