function [ces_all,cesn,judge]=ConnectionJudge4(conc,maskedge,origin,mask,T0,conn,fignum)

%conc：输入的怀疑由两或两根以上的电极连在一起的某一连通分量，是三维数组
%maskedge：颅骨边界，通常是brainmask_edge2
%origin：未经阈值处理的初始图像，此处通常是未去除颅外干扰的regCT
%mask：大脑掩模（通常是brainmask2）
%T0：本函数中对origin进行统一阈值处理时所使用的阈值
%conn：连通性
%fignum：所画的第一张图片的序号

%ces_all：输入连通分量conc的全部电极入点，是四维数组，ces_all(:,:,:,i)是第i个电极入点
%cesn:电极入点总数，是一个标量数值
%judge：judge=1表明输入连通分量conc确实是由两或两根以上的电极连在一起

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
    ces_all=[];
    cesn=1;
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
    
    [originT0_L,originT0_n]=bwlabeln(originT0,conn);
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
    [orgTc_outmask_L,orgTc_outmask_n]=bwlabeln(orgTc_outmask,conn);
    for j=1:orgTc_outmask_n
        if length(find(orgTc_outmask(orgTc_outmask_L==j)))<3
            orgTc_outmask(orgTc_outmask_L==j)=0;
        end
    end
    figure(fignum);PlotElectrodes(orgTc_outmask,'b');
    fignum=fignum+1;

    existL2=[];
    for i=1:originT0_n
        if ~isempty(find(orgTc_outmask==i))
            existL2=[existL2;i];
        end
    end
    N2=length(existL2);
    cesn=N2;

    if N2<2
        ces_all=[];
    else
        [nx,ny,nz]=size(conc);
        ces_all=zeros(nx,ny,nz,N2);
        judge=1;
        for i=1:N2
            lvi=existL2(i);orgTci=originT0;orgTci(originT0_L~=lvi)=0;orgTci=logical(orgTci);
            orgTci_outskull=logical(orgTci.*(~mask));
%             figure;PlotElectrodes(orgTci_outskull,'b');
%         end
            [~,orgTcios_n]=bwlabeln(orgTci_outskull,conn);

            orgTcices=logical(orgTci.*ces);
%             figure;PlotElectrodes(orgTcices,'b');
%         end
            [orgTcices_L,orgTcices_n]=bwlabeln(orgTcices,conn);
            for j=1:orgTcices_n
                orgTcices_Lj=orgTcices;
                orgTcices_Lj(orgTcices_L~=j)=0;
                orgTcices_Lj=logical(orgTcices_Lj);
                ocij=logical(orgTcices_Lj+orgTci_outskull);
%                 figure;PlotElectrodes(ocij,'b');
                [~,ocij_n]=bwlabeln(ocij,conn);
                if ocij_n==1 %注意不是ocij_n==orgTcices_n
                    cesi=orgTcices_Lj;
                end
            end

            figure(fignum);PlotElectrodes(cesi,'b');
            fignum=fignum+1;
            ces_all(:,:,:,i)=cesi;
        end
        
    end
    
end

end