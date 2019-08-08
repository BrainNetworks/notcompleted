function [ces1,ces2,judge]=ConnectionJudge3(conc,maskedge,origin,mask,T0,conn,fignum)

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
ces=bwareaopen(ces,2,conn); %去除由于maskedge的边缘不连续不光滑性而产生的零散小像素点

[~,ces_n]=bwlabeln(ces,conn);

if ces_n<2
    ces1=[];
    ces2=[];
end

if ces_n>=2
    
    figure(fignum);PlotElectrodes(ces,'b');hold on;
    fignum=fignum+1;

    originT0=origin;
    originT0(origin>T0)=T0;
    originT0(origin<T0)=0;
    originT0=logical(originT0);
    figure(fignum);PlotElectrodes(originT0,'b');
    fignum=fignum+1;
    
    [originT0_L,originT0_n]=bwlabeln(originT0);
    orgT0Lces=ces.*originT0_L;
    
    existL1=[];
    for i=1:originT0_n
        if ~isempty(find(orgT0Lces==i))
            existL1=[existL1;i];
        end
    end
    N1=length(existL1);
    
    orgTc=zeros(size(originT0_L));
    for j=1:N1
        orgTc(originT0_L==existL1(j))=existL1(j);
    end
    figure(fignum);PlotElectrodes(orgTc,'b');
    fignum=fignum+1;
        
    orgTc_outmask=orgTc.*logical(~mask);
    figure(fignum);PlotElectrodes(orgTc_outmask,'b');
    fignum=fignum+1;

    existL2=[];
    for i=1:originT0_n
        if ~isempty(find(orgTc_outmask==i))
            existL2=[existL2;i];
        end
    end
    N2=length(existL2);

    if N2<2
        ces1=[];
        ces2=[];
    else
        judge=1;
        lv1=existL2(1);orgTc1=originT0;orgTc1(originT0_L~=lv1)=0;orgTc1=logical(orgTc1);
        lv2=existL2(2);orgTc2=originT0;orgTc2(originT0_L~=lv2)=0;orgTc2=logical(orgTc2);
        
        orgTc1_outskull=logical(orgTc1.*(~mask));
        orgTc2_outskull=logical(orgTc2.*(~mask));
        [~,orgTc1os_n]=bwlabeln(orgTc1_outskull,conn);
        [~,orgTc2os_n]=bwlabeln(orgTc2_outskull,conn);

%         orgTc1os_ind=find(orgTc1_outskull);
%         [c1,r1,s1]=ind2sub(size(orgTc1_outskull),orgTc1os_ind);
%         orgTc1os_sub=[c1,r1,s1];
%         orgTc2os_ind=find(orgTc2_outskull);
%         [c2,r2,s2]=ind2sub(size(orgTc2_outskull),orgTc2os_ind);
%         orgTc2os_sub=[c2,r2,s2];

        orgTc1ces=logical(orgTc1.*ces);
        [orgTc1ces_L,orgTc1ces_n]=bwlabeln(orgTc1ces,conn);
        for j=1:orgTc1ces_n
            orgTc1ces_Lj=orgTc1ces;
            orgTc1ces_Lj(orgTc1ces_L~=j)=0;
            orgTc1ces_Lj=logical(orgTc1ces_Lj);
            oc1j=logical(orgTc1ces_Lj+orgTc1_outskull);
            [~,oc1j_n]=bwlabeln(oc1j,conn);
            if oc1j_n==orgTc1os_n
                ces1=orgTc1ces_Lj;
            end
        end

        orgTc2ces=logical(orgTc2.*ces);
        [orgTc2ces_L,orgTc2ces_n]=bwlabeln(orgTc2ces,conn);
        for j=1:orgTc2ces_n
            orgTc2ces_Lj=orgTc2ces;
            orgTc2ces_Lj(orgTc2ces_L~=j)=0;
            orgTc2ces_Lj=logical(orgTc2ces_Lj);
            oc2j=logical(orgTc2ces_Lj+orgTc2_outskull);
            [~,oc2j_n]=bwlabeln(oc2j,conn);
            if oc2j_n==orgTc2os_n
                ces2=orgTc2ces_Lj;
            end
        end
        
        figure(fignum);PlotElectrodes(ces1,'b');hold on;
        figure(fignum);PlotElectrodes(ces2,'k');hold on;

    end
    
end

end