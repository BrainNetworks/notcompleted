function [ces1,ces2,judge]=ConnectionJudge(conc,maskedge,origin,mask,T0,conn,fignum)

%conc：输入的怀疑由两根电极连在一起的某一连通分量，是三维数组
%maskedge：颅骨边界，通常是brainmask_edge2
%origin：未经阈值处理的初始图像，此处通常是未去除颅外干扰的regCT
%mask：大脑掩模（通常是brainmask2）
%T0：阈值初值，也就是对origin进行统一阈值处理时所使用的阈值
%conn：连通性
%fignum：所画的第一张图片的序号

%ces1：输入连通分量conc与颅骨边界maskedge的第一个交集，也是第一根电极的入点，是三维体数据
%ces2：输入连通分量conc与颅骨边界maskedge的第二个交集，也是第二根电极的入点，是三维体数据
%judge：judge=1表明输入连通分量conc确实是由两根电极连在一起

if nargin<7
    fignum=1;
end

if nargin<6
    conn=26;
end

judge=0;

conc=logical(conc);
maskedge=logical(maskedge);

ces=logical(conc.*maskedge); %ces:concs_elecenters
%注意：此处不能用ces=logical(conc.*elecenters)来计算ces，因为此时计算得到的结果
%是初始时的电极入点，没有经过区域生长（或闭操作）、入点修正和最终的去除干扰，这样
%一来得到的ces很可能不止两个连通分量，从而对判断造成影响。
ces=bwareaopen(ces,3,conn); %去除由于maskedge的边缘不连续不光滑性而产生的零散小像素点

[ces_L,ces_n]=bwlabeln(ces,conn);

if ces_n<2
    ces1=[];
    ces2=[];
end

if ces_n==2
    
    figure(fignum);PlotElectrodes(ces,'b');hold on;

    ces1=ces_L;ces1(ces1~=1)=0;ces1=logical(ces1);
    ces2=ces_L;ces2(ces2~=2)=0;ces2=logical(ces2);

    [conc_coeff,~,~,~,conc_explained,~]=PCACalculation4(conc);
    conc_coeff1=conc_coeff(:,1);
 
    if conc_explained(1)<60
        judge=1;
    end
    
    if (conc_explained(1)>=60)&&(conc_explained(1)<85)

        ces1_Centroid=regionprops3(ces1,'Centroid');
        ces2_Centroid=regionprops3(ces2,'Centroid');
        ces1_Centroid=ces1_Centroid{1,1};
        ces2_Centroid=ces2_Centroid{1,1};
        
%         ces1_Centrsub=round(ces1_Centroid(:,[2,1,3]));
%         ces2_Centrsub=round(ces2_Centroid(:,[2,1,3]));
        ces1_Centrsub=ces1_Centroid(:,[2,1,3]); %不对ces1_Centrsub取整，以求后续得到更加精确的ces12_coeff1
        ces2_Centrsub=ces2_Centroid(:,[2,1,3]);
        figure(fignum);plot3(ces1_Centrsub(2),ces1_Centrsub(1),ces1_Centrsub(3),'o','color','r');hold on;
        figure(fignum);plot3(ces2_Centrsub(2),ces2_Centrsub(1),ces2_Centrsub(3),'o','color','r');hold on;
        fignum=fignum+1;

        ces12=ces1_Centrsub-ces2_Centrsub;
        ces12_d=norm(ces12);
        ces12_coeff1=(ces12./ces12_d)';
        diff=norm(conc_coeff1-ces12_coeff1); %值域在0~sqrt(6)=2.4495之间
        
        if diff>1.2
            judge=1;
        end
        
    end
    
    if judge==0
        fignum=fignum+1;
        originT0=origin;
        originT0(origin>T0)=T0;
        originT0(origin<T0)=0;
        originT0=logical(originT0);
        figure(fignum);PlotElectrodes(originT0);
        fignum=fignum+1;
        [originT0_L,~]=bwlabeln(originT0);
        ces1_ind=find(ces1);
        ces2_ind=find(ces2);
        ces1_Lv=originT0_L(ces1_ind(round(end/2)));
        ces2_Lv=originT0_L(ces2_ind(round(end/2)));
        originT0_ces12=originT0_L;
        originT0_ces12(originT0_L==ces2_Lv)=ces1_Lv;
        originT0_ces12(originT0_ces12~=ces1_Lv)=0;
        originT0_ces12=logical(originT0_ces12);
        figure(fignum);PlotElectrodes(originT0_ces12);
        fignum=fignum+1;
        orTc_outmask=logical(originT0_ces12.*(~mask));
        orTc_outmask=bwareaopen(orTc_outmask,3,conn); %去除由于maskedge的边缘不连续不光滑性而产生的零散小像素点
        figure(fignum);PlotElectrodes(orTc_outmask);
        [~,orTcom_n]=bwlabeln(orTc_outmask);
        if orTcom_n==2
            judge=1;
        end
    end
    
end

if ces_n>2
    
    judge=1;

    figure(fignum);PlotElectrodes(ces,'b');
    fignum=fignum+1;

    originT0=origin;
    originT0(origin>T0)=T0;
    originT0(origin<T0)=0;
    originT0=logical(originT0);
    figure(fignum);PlotElectrodes(originT0);
    fignum=fignum+1;
    [originT0_L,~]=bwlabeln(originT0);
    
    i_rec=[];
    for i=1:ces_n
        cesi=ces_L;cesi(cesi~=i)=0;cesi=logical(cesi);
        cesi_ind=find(cesi);
        cesi_Lv=originT0_L(cesi_ind(round(end/2)));
        originT0_cesi=originT0_L;
        originT0_cesi(originT0_L~=cesi_Lv)=0;
        originT0_cesi=logical(originT0_cesi);
        figure(fignum);PlotElectrodes(originT0_cesi);
        fignum=fignum+1;
        orTci_outmask=logical(originT0_cesi.*(~mask));
        figure(fignum);PlotElectrodes(orTci_outmask);
        if sum(sum(sum(orTci_outmask~=0)))~=0
            i_rec=[i_rec,i];
        end
    end 
    ces1=ces_L;ces1(ces1~=i_rec(1))=0;ces1=logical(ces1);
    ces2=ces_L;ces2(ces2~=i_rec(2))=0;ces2=logical(ces2);

end

end