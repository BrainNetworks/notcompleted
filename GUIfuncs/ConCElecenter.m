function ces_all=ConCElecenter(conc,maskedge,origin,mask,T0,conn,fignum)

%conc：某一电极的某一连通分量，是三维数组
%maskedge：颅骨边界，通常是brainmask_edge2
%origin：未经阈值处理的初始图像，此处通常是未去除颅外干扰的regCT
%mask：大脑掩模（通常是brainmask2）
%T0：本函数中对origin进行统一阈值处理时所使用的阈值
%conn：连通性
%fignum：所画的第一张图片的序号

%ces_all：输入连通分量conc的电极入点，在本函数中输入conc应该均只有一个电极入点

if nargin<7
    fignum=1;
end

if nargin<6
    conn=26;
end

conc_ecp=[];
conc=logical(conc);
maskedge=logical(maskedge);
[nx,ny,nz]=size(conc);

ces=logical(conc.*maskedge); %ces:concs_elecenters
ces=bwareaopen(ces,2,conn); %去除由于maskedge的边缘不连续不光滑性而产生的零散小像素点
[ces_L,ces_n]=bwlabeln(ces,conn);

if ces_n>1
    
%     figure(fignum);PlotElectrodes(ces,'b');hold on;
%     fignum=fignum+1;

    originT0=origin;
    originT0(origin>T0)=T0;
    originT0(origin<T0)=0;
    originT0=logical(originT0);
    [originT0_L,~]=bwlabeln(originT0,conn);
%     figure(fignum);PlotElectrodes(originT0,'b');
%     fignum=fignum+1;
    
    Num_all=[];
    for i=1:ces_n
        cesi=ces_L;
        cesi(ces_L~=i)=0;
        cesi=logical(cesi);
        cesiorgL=cesi.*originT0_L;
        cesiorgL_sub=find(cesiorgL);
        Num=unique(cesiorgL(cesiorgL_sub));
        NN=length(Num);
        if NN==1
            Num_all=[Num_all;Num];
        end
        if NN>1
            L=zeros(NN,1);
            for j=1:NN
                cesiorgLj=cesiorgL;
                cesiorgLj(cesiorgL~=Num(j))=0;
                L(j)=length(find(cesiorgLj));
            end
            LmaxL=find(L==max(L));
            Num_all=[Num_all;Num(LmaxL)];
        end
    end
    
    N1=length(Num_all);
    if N1>0
        orgTc=zeros(nx,ny,nz);
        for i=1:N1
            orgTc(originT0_L==Num_all(i))=Num_all(i);
        end
%         figure(fignum);PlotElectrodes(orgTc,'b');
%         fignum=fignum+1;

        orgTc_outmask=orgTc.*logical(~mask);
        orgTc_outmask_sub=find(orgTc_outmask);
        Num2_all=unique(orgTc_outmask(orgTc_outmask_sub));
        
        N2=length(Num2_all);
        if N2>0
            orgom=false(nx,ny,nz);
            for i=1:N2
                orgomi=orgTc_outmask;
                orgomi(orgTc_outmask~=Num2_all(i))=0;
                [orgomi_L,orgomi_n]=bwlabeln(orgomi,conn);
                if orgomi_n==1
                    orgom=logical(orgom+orgomi);
                end
                if orgomi_n>1
                    L2=zeros(orgomi_n,1);
                    for j=1:orgomi_n
                        orgomiLj=orgomi_L;
                        orgomiLj(orgomi_L~=j)=0;
                        L2(j)=length(find(orgomiLj));                
                    end            
                    L2maxL=find(L2==max(L2));
                    orgomiLj=orgomi_L;
                    orgomiLj(orgomi_L~=L2maxL)=0;
                    orgomiLj=logical(orgomiLj);
                    orgom=logical(orgom+orgomiLj);
                end
            end
            
            [orgom_L,orgom_n]=bwlabeln(orgom,conn);
%             figure(fignum);PlotElectrodes(orgom,'b');
%             fignum=fignum+1;
            
            if orgom_n>0
                cesn=orgom_n;
                ces_all=zeros(nx,ny,nz,cesn);
                for i=1:orgom_n
                    orgomi=orgom;
                    orgomi(orgom_L~=i)=0;
                    orgomi=logical(orgomi);
                    for j=1:ces_n
                        cesj=ces_L;
                        cesj(ces_L~=j)=0;
                        cesj=logical(cesj);
                        orgomicj=logical(orgomi+cesj);
                        [~,orgomicj_n]=bwlabeln(orgomicj,conn);
                        if orgomicj_n==1
                            ces_all(:,:,:,i)=cesj;
%                             figure(fignum);PlotElectrodes(cesj,'b');
%                             fignum=fignum+1;
                            break;
                        end
                    end
                end
            end
        end
    end    
    
end

end