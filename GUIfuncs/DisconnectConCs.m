function Eallconcs=DisconnectConCs(conc,maskedge,ces_all,origelec,conn,fignum)

%conc：经ConnectionJudge判定是由三或三根以上的电极连在一起的某一连通分量，是三维数组
%maskedge：颅骨边界，通常是brainmask_edge2
%ces_all：输入连通分量conc与颅骨边界maskedge的全部交集，也即全部电极入点，是四维数组，ces_all(:,:,:,i)是第i个电极入点cesi
%origelec：regCT去除颅外干扰和阈值处理后得到的初步电极结果，通常是electrode2,electrode2B等
%conn：连通性
%fignum：所画的第一张图片的序号

%Eallconcs：是一个N*1的元胞数组，N是电极总数，每一个元胞数组Eallconcs{i}等于cesi加上origelec中属于第i电极的所有连通分量，Eallconcs{i}是三维体数据

if nargin<7
    fignum=1;
end

if nargin<6
    conn=26;
end

conc=logical(conc);
maskedge=logical(maskedge);
origelec=logical(origelec);

%parula色系：电极
%黑色：未归属的连通分量
%红色：电极的中点线
oics=logical(conc.*(origelec-origelec.*maskedge)); %oics:origin_intracranial_concs
oics=bwareaopen(oics,2,conn); %去除由于maskedge的边缘不连续不光滑性而产生的零散小像素点

N=size(ces_all,4);
parulac=parula;
cl=size(parulac,1);
dc=round(cl/(N-1));
ces_allcolor=zeros(N,3);
ces_allcolor(1,:)=parulac(1,:);
ces_allcolor(N,:)=parulac(end,:);
for i=2:(N-1)
    ces_allcolor(i,:)=parulac(dc*(i-1),:);
end

for i=1:N
    cesi=ces_all(:,:,:,i);
    cesi_color=ces_allcolor(i,:);
%     figure(fignum);PlotElectrodes(cesi,cesi_color);hold on;
end
% figure(fignum);PlotElectrodes(oics,'k');hold on;

ces_Centrsub=zeros(N,3);
for i=1:N
    cesi=ces_all(:,:,:,i);
    cesi_Centroid=regionprops3(cesi,'Centroid');
    cesi_Centroid=cesi_Centroid{1,1};
    cesi_Centrsub=cesi_Centroid(:,[2,1,3]);
    ces_Centrsub(i,:)=cesi_Centrsub;
%     figure(fignum);plot3(cesi_Centrsub(2),cesi_Centrsub(1),cesi_Centrsub(3),'o','color','r');hold on;
end

[oics_L,oics_n]=bwlabeln(oics,conn);
oics_ind=find(oics);
[oics_r,oics_c,oics_s]=ind2sub(size(oics),oics_ind);
oics_sub=[oics_r,oics_c,oics_s];

[nx,ny,nz]=size(conc);XYZWorldLimits=[1,ny,1,nx,1,nz];
oicc_Centrsub=zeros(N,3);
oic_closer=cell(N,1);
line_sub=cell(N,1);
for i=1:N
    cesi_Centrsub=ces_Centrsub(i,:);
    D_oics2cesiCt=pdist2(oics_sub,cesi_Centrsub,'euclidean');
    dmini=min(D_oics2cesiCt);
    DM=25;
    if dmini<DM
        Li=(D_oics2cesiCt==dmini);
        Lvi=oics_L(oics_ind(Li));
        oic_closeri=oics_L;oic_closeri(oics_L~=Lvi)=0;oic_closeri=logical(oic_closeri);
    else
        ratioi_rec=zeros(oics_n,1);
        oics_concs=cell(oics_n,1);
        conc_ind=find(conc);
        [conc_r,conc_c,conc_s]=ind2sub(size(conc),conc_ind);
        conc_sub=[conc_r,conc_c,conc_s];
        for k=1:oics_n
            oicsk=oics_L;oicsk(oics_L~=k)=0;oicsk=logical(oicsk);
            oics_concs{k}=oicsk;
            oicsk_Centroid=regionprops3(oicsk,'Centroid');
            oicsk_Centroid=oicsk_Centroid{1,1};
            oicsk_Centrsub=oicsk_Centroid(:,[2,1,3]);
            [linek_sub,~]=MakeLine(cesi_Centrsub,oicsk_Centrsub);
            linek_sub=round(linek_sub);
            ISi=intersect(linek_sub,conc_sub,'rows');
            ratioi=size(ISi,1)/size(linek_sub,1);
            ratioi_rec(k)=ratioi;
        end
        ratioimax=max(ratioi_rec);
        Li=(ratioi_rec==ratioimax);
        oic_closeri=oics_concs{Li};
    end
    oic_closer{i}=oic_closeri;
    oicci_Centroid=regionprops3(oic_closeri,'Centroid');
    oicci_Centroid=oicci_Centroid{1,1};
    oicci_Centrsub=oicci_Centroid(:,[2,1,3]);
    oicc_Centrsub(i,:)=oicci_Centrsub;
    [linei_sub,~]=MakeLine(cesi_Centrsub,oicci_Centrsub);
    line_sub{i}=linei_sub;

    cesi=ces_all(:,:,:,i);
    cesi_color=ces_allcolor(i,:);
%     figure(fignum+1);PlotElectrodes(cesi,cesi_color);hold on;
%     figure(fignum+1);PlotElectrodes(oic_closeri,cesi_color);hold on;
%     figure(fignum+1);plot3(cesi_Centrsub(2),cesi_Centrsub(1),cesi_Centrsub(3),'o','color','r');hold on;
%     figure(fignum+1);plot3(oicci_Centrsub(2),oicci_Centrsub(1),oicci_Centrsub(3),'o','color','r');hold on;
%     figure(fignum+1);PlotElectrodes2(linei_sub,XYZWorldLimits,'r','-.');hold on;
end

oics2=oics;
for i=1:N
    oic_closeri=oic_closer{i};
    oics2=logical(oics2-oic_closeri);
end
% figure(fignum+1);PlotElectrodes(oics2,'k');hold on;
[oics2_L,oics2_n]=bwlabeln(oics2,conn);

Eallconcs=cell(N,1);
for i=1:N
    cesi=ces_all(:,:,:,i);
    oic_closeri=oic_closer{i};
    Eiconcs=logical(cesi+oic_closeri);
    Eallconcs{i}=Eiconcs;
end

if oics2_n>0
    
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

    Qi_renew=cell(N,1);
    for i=1:N
        Qi_renew{i}=oicc_Centrsub(i,:);
    end

    CktoPiQi_rec=zeros(oics2_n,N,oics2_n);
    L_rec=[];
    oics_leave=oics2;
    for j=1:oics2_n

        for i=1:N
            Pi=ces_Centrsub(i,:);
            Qi=oicc_Centrsub(i,:);

            CktoPiQi=zeros(oics2_n,1);
            for k=1:oics2_n
                if isempty(intersect(L_rec,k))
                    Ck=oics2_stats(k,2:4);
                    CktoPiQi(k)=norm(cross(Pi-Qi,Ck-Qi))/norm(Pi-Qi);
                else
                    CktoPiQi(k)=NaN;
                end
            end
            CktoPiQi_rec(:,i,j)=CktoPiQi;
        end

        Dj=CktoPiQi_rec(:,:,j);
        Djmin=min(min(Dj));
        Djmin_ind=find(Dj==Djmin,1);
        [Djmin_r,Djmin_c]=ind2sub(size(Dj),Djmin_ind);
        L_rec=[L_rec;Djmin_r];

        Ejconcs=logical(Eallconcs{Djmin_c}+oics2_concs{Djmin_r});
        Eallconcs{Djmin_c}=Ejconcs;
        oics2_stats(Djmin_r,5)=Djmin_c;
        Qi_renew{Djmin_c}=[Qi_renew{Djmin_c};oics2_stats(Djmin_r,2:4)];
        Qi=mean(Qi_renew{Djmin_c},1);
        oicc_Centrsub(Djmin_c,:)=Qi;
        oics_leave=logical(oics_leave-oics2_concs{Djmin_r});
        
        for i=1:N
            cesi_color=ces_allcolor(i,:);
            Pi=ces_Centrsub(i,:);
            Qi=oicc_Centrsub(i,:);
            [linei_sub,~]=MakeLine(Pi,Qi);
%             figure(fignum+1+j);PlotElectrodes(Eallconcs{i},cesi_color);hold on;
%             figure(fignum+1+j);plot3(Pi(2),Pi(1),Pi(3),'o','color','r');hold on;
%             figure(fignum+1+j);plot3(Qi(2),Qi(1),Qi(3),'o','color','r');hold on;
%             figure(fignum+1+j);PlotElectrodes2(linei_sub,XYZWorldLimits,'r','-.');hold on;
        end
%         figure(fignum+1+j);PlotElectrodes(oics_leave,'k');hold on;
    end
    
end

% for i=fignum:1:(fignum+1+oics2_n)
%     figure(i);
%     Azimuth=34;Elevation=10;view(Azimuth,Elevation);
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

% filename1='SubjectName_DiscennectElectrodes.gif'; %文件名
% for frame=1:1:n
%     figure(frame);
%     Azimuth=34;Elevation=10;view(Azimuth,Elevation); %调整视角
%     f(frame)=getframe(frame); %捕获坐标区作为影片帧
%     imind=frame2im(f(frame)); %返回与影片帧关联的图像句柄
%     [imind,cm]=rgb2ind(imind,256); %将RGB图像转换为索引图像
%     if frame==1 %如果为第一帧
%         imwrite(imind,cm,filename1,'gif','Loopcount',inf,'DelayTime',0.5); %每隔0.5s写入一帧
%     else %如果不是第一帧
%         imwrite(imind,cm,filename1,'gif','WriteMode','append','DelayTime',0.5); %每隔0.5s写入一帧
%     end
% end

end