function [e17Cendpoints,e17Cendpoints_judge,e17Centerps,e17Ctargetps,e17Cepsj_ras,e17ElesLength,e17EiContacts_RASC,...
    e17EiContacts_sub,e17EiContacts_cell,e17Contacts]=ContactSeg3(T0C,T1_IOP,T1_IPP1,regCT,brainmask2,brainmask3,...
    brainmask_edge2,electrode12,electrode17,M_T1_recon,RAS_T1_recon,f)

waitbar(0.75,f,'Contact modeling...');

conn=26;
[nx,ny,nz]=size(electrode17);
e17Cendpoints=zeros(nx,ny,nz);
[electrode17_L,electrode17_n]=bwlabeln(electrode17,conn);
e17Cendpoints_sub=zeros(2,3,electrode17_n);
for i=1:electrode17_n
    
    e17ConCi=electrode17_L;
    e17ConCi(e17ConCi~=i)=0;
    e17ConCi=logical(e17ConCi);

    [~,~,~,~,e17Ciendpoints,e17Ciendpoints_sub]=ConcEndpsDetection2(e17ConCi);
    e17Cendpoints=logical(e17Cendpoints+e17Ciendpoints);
    e17Cendpoints_sub(:,:,i)=e17Ciendpoints_sub;
end
% figure(105);PlotElectrodes(electrode17,'b','.');hold on;
% figure(105);PlotElectrodes(e17Cendpoints,'r','o');hold on;
[~,e17Cendpoints_n]=bwlabeln(e17Cendpoints,26)

e17Cendps_sked=logical((~brainmask3).*e17Cendpoints); %sked：skull endpoints
% figure(108);PlotElectrodes(e17Cendps_sked,'r','o');
[~,e17Cendps_sked_n]=bwlabeln(e17Cendps_sked,26) %有些电极可能穿透整个大脑，从而其入点和靶点均在颅骨边界上

e17Cendpoints_judge=zeros(2,4,electrode17_n);
e17Cendpoints_judge(:,1:3,:)=e17Cendpoints_sub;
%e17Cendpoints_judge(i,4,j)=1则表明点e17Cendpoints_judge(i,1:3,j)是电极入点
%e17Cendpoints_judge(i,4,j)=0则表明点e17Cendpoints_judge(i,1:3,j)是电极靶点
e17Cendps_sked_ind=find(e17Cendps_sked);
[x1,y1,z1]=ind2sub(size(e17Cendps_sked),e17Cendps_sked_ind);
e17Cendps_sked_sub=[x1,y1,z1];
fignum=109;

for i=1:electrode17_n
    
    for j=1:2
        IS=intersect(e17Cendpoints_sub(j,:,i),e17Cendps_sked_sub,'rows');
        if size(IS,1)==1
            e17Cendpoints_judge(j,4,i)=1;
        end
    end
    
    IS2=intersect(e17Cendpoints_sub(:,:,i),e17Cendps_sked_sub,'rows');
    if size(IS2,1)==2
        e17ConCi=electrode17_L;
        e17ConCi(e17ConCi~=i)=0;
        e17ConCi=logical(e17ConCi);
        e17ConCi_ind=find(e17ConCi);
%         figure(fignum);PlotElectrodes(e17ConCi,'b','.');fignum=fignum+1;

        [electrode12_L,electrode12_n]=bwlabeln(electrode12,conn);
        e17Cie12Lv=unique(electrode12_L(e17ConCi_ind));
        e12ConCLv=electrode12;
        e12ConCLv(electrode12_L~=e17Cie12Lv)=0;
        e12ConCLv=logical(e12ConCLv);
%         figure(fignum);PlotElectrodes(e12ConCLv,'b','.');fignum=fignum+1;

        conc=e12ConCLv;maskedge=brainmask_edge2;origin=regCT;mask=brainmask2;T0=T0C;
        ces_all=ConCElecenter(conc,maskedge,origin,mask,T0,conn,fignum);
        e12ConCLv_ces=ces_all;
        
        e12ConCLv_ces_ind=find(e12ConCLv_ces);
        [x2,y2,z2]=ind2sub(size(e12ConCLv_ces),e12ConCLv_ces_ind);
        e12ConCLv_ces_sub=[x2,y2,z2];
        e12ConCLv_ces_ctrsub=[mean(x2),mean(y2),mean(z2)];
        e12ConCLv_ces_ctrsub=round(e12ConCLv_ces_ctrsub);
        D1=norm(e17Cendpoints_sub(1,:,i)-e12ConCLv_ces_ctrsub); %计算向量的2-范数，即欧式范数，即等于组成向量的两点的欧式距离
        D2=norm(e17Cendpoints_sub(2,:,i)-e12ConCLv_ces_ctrsub);
        if D1<D2
            e17Cendpoints_judge(1,4,i)=1;
            e17Cendpoints_judge(2,4,i)=0;
        else
            e17Cendpoints_judge(1,4,i)=0;
            e17Cendpoints_judge(2,4,i)=1;
        end
    end

end

e17Centerps=false(nx,ny,nz);
e17Ctargetps=false(nx,ny,nz);
for i=1:electrode17_n
    for j=1:2
        if e17Cendpoints_judge(j,4,i)==1
            e17Centerps(e17Cendpoints_judge(j,1,i),e17Cendpoints_judge(j,2,i),e17Cendpoints_judge(j,3,i))=1;
        else
            e17Ctargetps(e17Cendpoints_judge(j,1,i),e17Cendpoints_judge(j,2,i),e17Cendpoints_judge(j,3,i))=1;            
        end
    end
end
% figure(116);PlotElectrodes(e17Centerps,'b','o');
% [~,e17Centerps_n]=bwlabeln(e17Centerps,26)
% figure(117);PlotElectrodes(e17Ctargetps,'b','o');
% [~,e17Ctargetps_n]=bwlabeln(e17Ctargetps,26)

figure(118);
PlotElectrodes(electrode17,'b','.');hold on;
PlotElectrodes(e17Centerps,'r','o');hold on;
PlotElectrodes(e17Ctargetps,'k','o');hold on;

waitbar(0.8,f,'Contact modeling...');

if isempty(RAS_T1_recon)
    PixelSpacing=[0.6;0.6];SliceThickness=0.6;coordinate_judge=1;
    [M_e17Ceps,RAS_e17Ceps]=transform2RAS(e17Centerps,PixelSpacing,SliceThickness,T1_IOP,T1_IPP1,coordinate_judge);
else
    M_e17Ceps=M_T1_recon;
    RAS_e17Ceps=RAS_T1_recon;
end
        
electrode17_n=size(e17Cendpoints_judge,3);
e17Cepsj_ras=e17Cendpoints_judge;
e17ElesLength=zeros(electrode17_n,1);
for k=1:electrode17_n
    for j=1:2
        %将索引转换为坐标，因为输入RAS_e17Ceps是坐标集，而不是下标索引集，所以需要进行转换以与其对应
        X=e17Cendpoints_judge(j,2,k);
        Y=e17Cendpoints_judge(j,1,k);
        Z=e17Cendpoints_judge(j,3,k);
        R=RAS_e17Ceps(1,X,Y,Z);
        A=RAS_e17Ceps(2,X,Y,Z);
        S=RAS_e17Ceps(3,X,Y,Z);
        e17Cepsj_ras(j,1:3,k)=[R,A,S];
    end
    e17ElesLength(k)=norm(e17Cepsj_ras(1,1:3,k)-e17Cepsj_ras(2,1:3,k));
end

d0=1+1.5+1; %相邻两个电极触点中心之间的距离
e17EiContacts_RASC=cell(electrode17_n,1);
e17EiContacts_sub=cell(electrode17_n,1);
e17EiContacts_cell=cell(electrode17_n,1);
[nx1,ny1,nz1]=size(e17Centerps);
e17Contacts=false(nx1,ny1,nz1);
for k=1:electrode17_n
    for j=1:2
        if e17Cepsj_ras(j,4,k)==0 %电极靶点
            r1=e17Cepsj_ras(j,1,k);
            a1=e17Cepsj_ras(j,2,k);
            s1=e17Cepsj_ras(j,3,k);
        else %电极入点
            r2=e17Cepsj_ras(j,1,k);
            a2=e17Cepsj_ras(j,2,k);
            s2=e17Cepsj_ras(j,3,k);            
        end
    end
    D=e17ElesLength(k);
    if D>=90 %两组触点
        d1=(1:d0:(1+5*d0))';
        d2=((1+5*d0+12):d0:(1+5*d0+12+9*d0))';
        d=[d1;d2];
    end
    if D<90
        d=(1:d0:D)';
        if length(d)>18 %电极触点最多只有18个
            d=d(1:18);
        end
    end
    rn=d/D*(r2-r1)+r1;
    an=d/D*(a2-a1)+a1;
    sn=d/D*(s2-s1)+s1;
    rasn=[rn,an,sn];
    e17EiContacts_RASC{k}=rasn;
    O=ones(1,length(d));
    rasn1=[rasn';O];
    xyz1=M_e17Ceps\rasn1;
    xyz=round((xyz1(1:3,:))');
    sub=xyz(:,[2,1,3]); %将坐标值xyz转换为索引值sub
    e17EiContacts_sub{k}=sub;
    e17EiContacts=false(nx1,ny1,nz1);
    for i=1:length(d)
        e17EiContacts(sub(i,1),sub(i,2),sub(i,3))=1;
    end
    e17EiContacts_cell{k}=e17EiContacts;
    e17Contacts=e17Contacts+e17EiContacts;
end

for k=1:electrode17_n
    EiCsub=e17EiContacts_sub{k};
    figure(118);scatter3(EiCsub(:,2),EiCsub(:,1),EiCsub(:,3),25,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r');hold on;
end
legend({'电极连通分量','电极入点','电极靶点','电极触点'},'FontSize',14,'Location','northwest','NumColumns',3);

end