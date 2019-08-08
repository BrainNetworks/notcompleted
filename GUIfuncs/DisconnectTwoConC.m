function [E1concs,E2concs]=DisconnectTwoConC(conc,maskedge,ces1,ces2,origelec,conn,fignum)

%conc：经ConnectionJudge判定是由两根电极连在一起的某一连通分量，是三维数组
%maskedge：颅骨边界，通常是brainmask_edge2
%ces1：输入连通分量conc与颅骨边界maskedge的第一个交集，也是第一根电极的入点，是三维体数据
%ces2：输入连通分量conc与颅骨边界maskedge的第二个交集，也是第二根电极的入点，是三维体数据
%origelec：regCT去除颅外干扰和阈值处理后得到的初步电极结果，通常是electrode2,electrode2B等
%conn：连通性
%fignum：所画的第一张图片的序号

%E1concs：E1concs等于ces1加上origelec中属于第一根电极的所有连通分量，是三维体数据
%E2concs：E2concs等于ces2加上origelec中属于第二根电极的所有连通分量，是三维体数据

if nargin<7
    fignum=1;
end

if nargin<6
    conn=26;
end

conc=logical(conc);
maskedge=logical(maskedge);
origelec=logical(origelec);

%蓝色：电极1
%黑色：电极2
%绿色：未归属的连通分量
%红色：电极的中点线
oics=logical(conc.*(origelec-origelec.*maskedge)); %oics:origin_intracranial_concs
% oics=bwareaopen(oics,2,conn); %去除由于maskedge的边缘不连续不光滑性而产生的零散小像素点
% figure(fignum);PlotElectrodes(ces1,'b');hold on; %ces:concs_elecenters
% figure(fignum);PlotElectrodes(ces2,'k');hold on;
% figure(fignum);PlotElectrodes(oics,'g');hold on;

ces1_Centroid=regionprops3(ces1,'Centroid');
ces2_Centroid=regionprops3(ces2,'Centroid');
ces1_Centroid=ces1_Centroid{1,1};
ces2_Centroid=ces2_Centroid{1,1};
ces1_Centrsub=ces1_Centroid(:,[2,1,3]);
ces2_Centrsub=ces2_Centroid(:,[2,1,3]);
% figure(fignum);plot3(ces1_Centrsub(2),ces1_Centrsub(1),ces1_Centrsub(3),'o','color','r');hold on;
% figure(fignum);plot3(ces2_Centrsub(2),ces2_Centrsub(1),ces2_Centrsub(3),'o','color','r');hold on;

[oics_L,oics_n]=bwlabeln(oics,conn);
oics_ind=find(oics);
[oics_r,oics_c,oics_s]=ind2sub(size(oics),oics_ind);
oics_sub=[oics_r,oics_c,oics_s];

D_oics2ces1Ct=pdist2(oics_sub,ces1_Centrsub,'euclidean');
dmin1=min(D_oics2ces1Ct);
DM=25;
if dmin1<DM
    L1=(D_oics2ces1Ct==dmin1);
    Lv1=oics_L(oics_ind(L1));
    oic_closer1=oics_L;oic_closer1(oics_L~=Lv1)=0;oic_closer1=logical(oic_closer1);
else
    ratio1_rec=zeros(oics_n,1);
    oics_concs=cell(oics_n,1);
    conc_ind=find(conc);
    [conc_r,conc_c,conc_s]=ind2sub(size(conc),conc_ind);
    conc_sub=[conc_r,conc_c,conc_s];
    for i=1:oics_n
        oicsi=oics_L;oicsi(oics_L~=i)=0;oicsi=logical(oicsi);
        oics_concs{i}=oicsi;
        oicsi_Centroid=regionprops3(oicsi,'Centroid');
        oicsi_Centroid=oicsi_Centroid{1,1};
        oicsi_Centrsub=oicsi_Centroid(:,[2,1,3]);
        [linei_sub,~]=MakeLine(ces1_Centrsub,oicsi_Centrsub);
        linei_sub=round(linei_sub);
        IS1=intersect(linei_sub,conc_sub,'rows');
        ratio1=size(IS1,1)/size(linei_sub,1);
        ratio1_rec(i)=ratio1;
    end
    ratio1max=max(ratio1_rec);
    L1=(ratio1_rec==ratio1max);
    oic_closer1=oics_concs{L1};
end
oicc1_Centroid=regionprops3(oic_closer1,'Centroid');
oicc1_Centroid=oicc1_Centroid{1,1};
oicc1_Centrsub=oicc1_Centroid(:,[2,1,3]);
[line1_sub,~]=MakeLine(ces1_Centrsub,oicc1_Centrsub);

D_oics2ces2Ct=pdist2(oics_sub,ces2_Centrsub,'euclidean');
dmin2=min(D_oics2ces2Ct);
if dmin2<DM
    L2=(D_oics2ces2Ct==dmin2);
    Lv2=oics_L(oics_ind(L2));
    oic_closer2=oics_L;oic_closer2(oics_L~=Lv2)=0;oic_closer2=logical(oic_closer2);
else
    ratio2_rec=zeros(oics_n,1);
    if dmin1<DM
        conc_ind=find(conc);
        [conc_r,conc_c,conc_s]=ind2sub(size(conc),conc_ind);
        conc_sub=[conc_r,conc_c,conc_s];
        oics_concs=cell(oics_n,1);
        for i=1:oics_n
            oicsi=oics_L;oicsi(oics_L~=i)=0;oicsi=logical(oicsi);
            oics_concs{i}=oicsi;
        end
    end
    for i=1:oics_n
        oicsi=oics_concs{i};
        oicsi_Centroid=regionprops3(oicsi,'Centroid');
        oicsi_Centroid=oicsi_Centroid{1,1};
        oicsi_Centrsub=oicsi_Centroid(:,[2,1,3]);
        [linei_sub,~]=MakeLine(ces2_Centrsub,oicsi_Centrsub);
        linei_sub=round(linei_sub);
        IS2=intersect(linei_sub,conc_sub,'rows');
        ratio2=size(IS2,1)/size(linei_sub,1);
        ratio2_rec(i)=ratio2;
    end
    ratio2max=max(ratio2_rec);
    L2=find(ratio2_rec==ratio2max);
    oic_closer2=oics_L;oic_closer2(oics_L~=L2)=0;oic_closer2=logical(oic_closer2);    
end
oicc2_Centroid=regionprops3(oic_closer2,'Centroid');
oicc2_Centroid=oicc2_Centroid{1,1};
oicc2_Centrsub=oicc2_Centroid(:,[2,1,3]);
[line2_sub,~]=MakeLine(ces2_Centrsub,oicc2_Centrsub);

[nx,ny,nz]=size(conc);XYZWorldLimits=[1,ny,1,nx,1,nz];
% % figure(fignum+1);PlotElectrodes(ces1,'b');hold on;
% % figure(fignum+1);PlotElectrodes(oic_closer1,'b');hold on;
% % figure(fignum+1);plot3(ces1_Centrsub(2),ces1_Centrsub(1),ces1_Centrsub(3),'o','color','r');hold on;
% % figure(fignum+1);plot3(oicc1_Centrsub(2),oicc1_Centrsub(1),oicc1_Centrsub(3),'o','color','r');hold on;
% % figure(fignum+1);PlotElectrodes2(line1_sub,XYZWorldLimits,'r','-.');hold on;
% % 
% figure(fignum+1);PlotElectrodes(ces2,'k');hold on;
% figure(fignum+1);PlotElectrodes(oic_closer2,'k');hold on;
% figure(fignum+1);plot3(ces2_Centrsub(2),ces2_Centrsub(1),ces2_Centrsub(3),'o','color','r');hold on;
% figure(fignum+1);plot3(oicc2_Centrsub(2),oicc2_Centrsub(1),oicc2_Centrsub(3),'o','color','r');hold on;
% figure(fignum+1);PlotElectrodes2(line2_sub,XYZWorldLimits,'r','-.');hold on;

oics2=logical(oics-oic_closer1-oic_closer2);
% figure(fignum+1);PlotElectrodes(oics2,'g');hold on;
[oics2_L,oics2_n]=bwlabeln(oics2,conn);

if oics2_n==0
    E1concs=logical(ces1+oic_closer1);
    E2concs=logical(ces2+oic_closer2);
else
    
    oics2_stats=zeros(oics2_n,5);
    oics2_stats(:,1)=(1:1:oics2_n)';
    oics2_concs=cell(oics2_n,1);
    for i=1:oics2_n
        oics2i=oics2_L;oics2i(oics2_L~=i)=0;oics2i=logical(oics2i);
        oics2_concs{i}=oics2i;
        oics2i_Centroid=regionprops3(oics2i,'Centroid');
        oics2i_Centroid=oics2i_Centroid{1,1};
        oics2i_Centrsub=oics2i_Centroid(:,[2,1,3]);
        oics2_stats(i,2:4)=oics2i_Centrsub;
    end

    E1concs=logical(ces1+oic_closer1);
    E2concs=logical(ces2+oic_closer2);
    P1=ces1_Centrsub;
    Q1=oicc1_Centrsub;
    P2=ces2_Centrsub;
    Q2=oicc2_Centrsub;

    Q1_renew=oicc1_Centrsub;
    Q2_renew=oicc2_Centrsub;
    CitoP1Q1_rec=zeros(oics2_n);
    CitoP2Q2_rec=zeros(oics2_n);
    L_rec=[];
    for j=1:oics2_n
        CitoP1Q1=zeros(oics2_n,1);
        CitoP2Q2=zeros(oics2_n,1);
        for i=1:oics2_n
            if isempty(intersect(L_rec,i))
                Ci=oics2_stats(i,2:4);
                CitoP1Q1(i)=norm(cross(P1-Q1,Ci-Q1))/norm(P1-Q1);
                CitoP2Q2(i)=norm(cross(P2-Q2,Ci-Q2))/norm(P2-Q2);
            end
        end
        CitoP1Q1_rec(:,j)=CitoP1Q1;
        CitoP2Q2_rec(:,j)=CitoP2Q2;
        diff=abs(CitoP1Q1-CitoP2Q2);
        diffmax=max(diff);
        L=find(diff==diffmax);
        L_rec=[L_rec;L];
        if CitoP1Q1(L)<CitoP2Q2(L)
            E1concs=logical(E1concs+oics2_concs{L});
            oics2_stats(L,5)=1;
            Q1_renew=[Q1_renew;oics2_stats(L,2:4)];
            Q1=mean(Q1_renew,1);
            [line1_sub,~]=MakeLine(P1,Q1);
        else
            E2concs=logical(E2concs+oics2_concs{L});
            oics2_stats(L,5)=2;
            Q2_renew=[Q2_renew;oics2_stats(L,2:4)];
            Q2=mean(Q2_renew,1);
            [line2_sub,~]=MakeLine(P2,Q2);
        end
        oics_leave=logical(oics+ces1+ces2-E1concs-E2concs);
%         figure(fignum+1+j);PlotElectrodes(E1concs,'b');hold on;
%         figure(fignum+1+j);plot3(P1(2),P1(1),P1(3),'o','color','r');hold on;
%         figure(fignum+1+j);plot3(Q1(2),Q1(1),Q1(3),'o','color','r');hold on;
%         figure(fignum+1+j);PlotElectrodes2(line1_sub,XYZWorldLimits,'r','-.');hold on;
%         figure(fignum+1+j);PlotElectrodes(E2concs,'k');hold on;
%         figure(fignum+1+j);plot3(P2(2),P2(1),P2(3),'o','color','r');hold on;
%         figure(fignum+1+j);plot3(Q2(2),Q2(1),Q2(3),'o','color','r');hold on;
%         figure(fignum+1+j);PlotElectrodes2(line2_sub,XYZWorldLimits,'r','-.');hold on;
%         figure(fignum+1+j);PlotElectrodes(oics_leave,'g');hold on;
    end

end

% for i=fignum:1:(fignum+1+oics2_n)
%     figure(i);
%     Azimuth=30;Elevation=13;view(Azimuth,Elevation);
% end

% filename1='SubjectName_DiscennectElectrodes.gif'; %生成文件名
% for frame=1:1:(oics2_n+2)
%     figure(fignum+frame-1);
%     Azimuth=34;Elevation=10;view(Azimuth,Elevation);
%     f(frame)=getframe(fignum+frame-1); %捕获坐标区作为影片帧
%     imind=frame2im(f(frame)); %返回与影片帧关联的图像句柄
%     [imind,cm]=rgb2ind(imind,256); %将RGB图像转换为索引图像
%     if frame==1 %如果为第一帧
%         imwrite(imind,cm,filename1,'gif','Loopcount',inf,'DelayTime',0.5); %每隔0.5s写入一帧
%     else %如果不是第一帧
%         imwrite(imind,cm,filename1,'gif','WriteMode','append','DelayTime',0.5); %每隔0.5s写入一帧
%     end
% end

end