function [electrode3_n,electrode3_sub,electrode8,electrode9,e9InnerConC,electrode10,electrode11,T0B,T0C,electrode2B,...
    electrode12,electrode13,electrode14,electrode15,electrode16,electrode17]=ContactSeg2(ElecNum,T0,Te,regCT,brainmask,...
    brainmask2,brainmask3,brainmask_edge2,electrode1,electrode3,electrode5,electrode5_S,f)
%%
waitbar(0.2,f,'Threshold-reduction region growing...');

conn=26;
[~,electrode3_n,electrode3_sub]=FindIndexOfAllConnectedComponents(electrode3,conn);
origin=electrode1;basis=electrode5;basis_S=electrode5_S;elecenters_n=electrode3_n;elecenters_sub=electrode3_sub;
th=2;ratio=0.01;fignum=82;
[ImageOut,InCs_stats,InCs_grow]=AllInnerConCsGrow4(origin,basis,basis_S,elecenters_n,elecenters_sub,th,ratio,T0,Te,fignum);
electrode8=ImageOut;
e5InCs_stats=InCs_stats;
e5InCs_grow=InCs_grow;

figure(93);PlotElectrodes(electrode5,'b');
figure(94);PlotElectrodes(electrode8,'b');
[~,electrode8_n]=bwlabeln(electrode8,conn)
% ShowAllConnectedComponents(electrode8,conn,95);
% [electrode8_PCAcoeff,electrode8_PCAlatentp]=PCACalculation2(electrode8,conn);

waitbar(0.32,f,'Threshold-reduction region growing...');

VolumeIn=electrode8;Mask=brainmask;conn=26; %注意的Mask用的是brainmask，而不是brainmask3
%原来的RemoveInterference2去除的是完全位于颅骨边界上的连通分量（像素总数大小不限），
%最新的RemoveInterference2去除的是完全位于颅骨边界上并且像素总数少于20个的连通分量。
VolumeOut=RemoveInterference2(VolumeIn,Mask,conn);
electrode9=VolumeOut; 
figure(95);PlotElectrodes(electrode9,'b');
[~,electrode9_n]=bwlabeln(electrode9,conn)
%%
waitbar(0.35,f,'Electrode enter point repair...');

ImageIn=electrode9;
origin=electrode1; %注意此处用回electrode1=regCT.*brainmask2，而不是electrode1B=regCT.*brainmask
mask=brainmask3; %注意此处不是用brainmask或者brainmask2，而是用brainmask3=logical(brainmask2-brainmask_edge2)
ratio=0.001;th=2;fignum=96;conn=26;
[InnerConC,label1,label2,ImageOut,CompStats]=ElecEnterpsCompletion4(ImageIn,origin,mask,ratio,th,T0,Te,fignum,conn);
e9InnerConC=InnerConC;
e9label1=label1;
e9label2=label2;
electrode10=ImageOut;
e10CompStats=CompStats;
figure(97);PlotElectrodes(electrode10,'b');
[~,electrode10_n]=bwlabeln(electrode10,conn)

e9label3=setdiff(e9label1,e9label2);
[electrode9_L,~]=bwlabeln(electrode9,conn);
[electrode10_L,~]=bwlabeln(electrode10,conn);

waitbar(0.4,f,'Electrode enter point repair...');

electrode11=electrode10;

for i=1:length(e9label3)
    
    e9conci=electrode9_L;
    e9conci(electrode9_L~=e9label3(i))=0;
    e9conci=logical(e9conci);
%     figure(96+2*i);PlotElectrodes(e9conci,'r');hold on;
    e9conci_ind=find(e9conci);
    [x1,y1,z1]=ind2sub(size(e9conci),e9conci_ind);
    e9conci_sub=[x1,y1,z1];
    e9conci_ctrsub=mean(e9conci_sub,1);
%     figure(96+2*i);plot3(e9conci_ctrsub(2),e9conci_ctrsub(1),e9conci_ctrsub(3),'o','color','g');hold on;

    e9IC_ind=find(e9InnerConC);
    [x2,y2,z2]=ind2sub(size(e9InnerConC),e9IC_ind);
    e9IC_sub=[x2,y2,z2];
    
    D=pdist2(e9IC_sub,e9conci_ctrsub,'euclidean');
    Dmin=min(D);
    %注意D是e9InnerConC（而不是electrode10）中的连通分量与点e9conci_ctrsub的距离，两者是不同的，因为electrode10中
    %含有e9label3的连通分量，同时还包含了e9label2及其生长结果，而e9InnerConC中都不含有，这就导致了D和Dmin的不同
    
    if Dmin<40
        
        Label=find(D==Dmin,1);
        closer=e9IC_sub(Label,:);

        e10LV=electrode10_L(closer(1),closer(2),closer(3));
        e10conci=electrode10_L;
        e10conci(electrode10_L~=e10LV)=0;
        e10conci=logical(e10conci);
%         figure(96+2*i);PlotElectrodes(e10conci,'b');hold on;
%         figure(96+2*i);plot3(closer(2),closer(1),closer(3),'o','color','g');hold on;
        
        [~,~,~,~,e10conci_explained,~]=PCACalculation4(e10conci);
        
        if e10conci_explained(1)>90

            [e10Cipca,~,~,~,~,e10Ciendpoints_sub]=ConcEndpsDetection2(e10conci);
%             figure(96+2*i);PlotElectrodes(e10Cipca,'k');hold on;
%             figure(96+2*i);plot3(e10Ciendpoints_sub(1,2),e10Ciendpoints_sub(1,1),e10Ciendpoints_sub(1,3),'o','color','k');hold on;
%             figure(96+2*i);plot3(e10Ciendpoints_sub(2,2),e10Ciendpoints_sub(2,1),e10Ciendpoints_sub(2,3),'o','color','k');hold on;

            [Pcloser,Qcloser,D_P2Q1Q2,D_Cc2Q1Q2,D_ConCj2Qcloser,e10conci_coeff,e10conci_explained,...,
                e9conci_coeff,e9conci_explained,pcadiff,iscolinear,judge]=ConCJudge(e9conci,e10conci,e10Ciendpoints_sub);

            if judge==1
                [connect,connect_sub,~,~,e9e10conci]=ConnectTwoConC(e9conci,e10conci,Pcloser,Qcloser);
%                 figure(96+2*i);PlotElectrodes(connect,'c');hold on;
%                 figure(97+2*i);PlotElectrodes(e9e10conci,'b');
                electrode11=logical(electrode11+e9e10conci);
            else
                electrode11(e9conci_ind)=0;
            end
            
        else
            
            [e9Cipca,~,~,~,~,e9Ciendpoints_sub]=ConcEndpsDetection2(e9conci);
%             figure(96+2*i);PlotElectrodes(e9Cipca,'k');hold on;
%             figure(96+2*i);plot3(e9Ciendpoints_sub(1,2),e9Ciendpoints_sub(1,1),e9Ciendpoints_sub(1,3),'o','color','k');hold on;
%             figure(96+2*i);plot3(e9Ciendpoints_sub(2,2),e9Ciendpoints_sub(2,1),e9Ciendpoints_sub(2,3),'o','color','k');hold on;            
            P1=e9Ciendpoints_sub(1,:);
            Q1=e9Ciendpoints_sub(2,:);
            
            e10conci_ind=find(e10conci);
            [x2,y2,z2]=ind2sub(size(e10conci),e10conci_ind);
            e10conci_sub=[x2,y2,z2];
            L=size(e10conci_sub,1);
            e10CtoPiQ1=zeros(L,1);
            for k=1:L
                e10Ck=e10conci_sub(k,:);
                e10CtoPiQ1(k)=norm(cross(P1-Q1,e10Ck-Q1))/norm(P1-Q1);
            end
            ratio=sum(e10CtoPiQ1<10)/L;
            
            if ratio>0.25
                [~,~,connect_neighb,~,e9e10conci]=ConnectTwoConC3(e10conci,e9conci,1); %注意此处e10conci和e9conci的顺序不能弄反
%                 figure(96+2*i);PlotElectrodes(connect_neighb,'c');hold on;
%                 figure(97+2*i);PlotElectrodes(e9e10conci,'b');
                electrode11=logical(electrode11+e9e10conci);
            else
                electrode11(e9conci_ind)=0;
            end
            
        end
        
        
    else
        electrode11(e9conci_ind)=0;
    end
    
end
conn=26;
[~,electrode11_n]=bwlabeln(electrode11,conn)
figure(98);PlotElectrodes(electrode11,'b');
%%
waitbar(0.5,f,'Interconnected Electrodes Determination and Separation...');

T0B=round(max(max(max(electrode1)))*0.8);
electrode2B=electrode1;
electrode2B(electrode1>T0B)=T0B;
electrode2B(electrode1<T0B)=0;
electrode2B=logical(electrode2B);
% [~,electrode2B_n]=bwlabeln(electrode2B,conn)
% figure(99);PlotElectrodes(electrode2B);

[electrode11_L,~]=bwlabeln(electrode11,conn);
electrode12=electrode11;
e12CL=[]; %记录由多根电极连在一起的连通分量的编号
T0C=round(max(max(max(electrode1)))*0.75);

maskedge=brainmask_edge2;origin=regCT;mask=brainmask2;conn=26;fignum=100;
origelec=electrode2B; %注意此处的origelec用的是electrode2B而不是electrode2
e11cesn_rec=zeros(electrode11_n,1);
e11judge_rec=zeros(electrode11_n,1);
for i=1:electrode11_n
    
    e11conci=electrode11_L;
    e11conci(electrode11_L~=i)=0;
    e11conci=logical(e11conci);
    conc=e11conci;
    [ces_all,cesn,judge]=ConnectionJudge5(conc,maskedge,origin,mask,T0C,conn,fignum); %注意此处是T0C而不是T0
    e11cesn_rec(i)=cesn;
    e11judge_rec(i)=judge;

    if judge==1
        
        e12CL=[e12CL,i];
        
        if cesn==2
            fignum=fignum+20; %ConnectionJudge5可能产生较多图片
            e11conci_ind=find(e11conci);
            electrode12(e11conci_ind)=0;
            ces1=ces_all(:,:,:,1);ces2=ces_all(:,:,:,2);
            
            conc=e11conci;
            [nx,ny,nz]=size(conc);
            [E1concs,E2concs]=DisconnectTwoConC(conc,maskedge,ces1,ces2,origelec,conn,fignum);
%             fignum=fignum+40; %DisconnectTwoConC可能产生较多图片
%             figure(fignum);PlotElectrodes(E1concs,'b');hold on;
%             figure(fignum);PlotElectrodes(E2concs,'r');hold on;
%             fignum=fignum+1;

            E1=ConnectConCsToElectrode(E1concs,conn,1);
            E2=ConnectConCsToElectrode(E2concs,conn,1);
            EC=logical(E1+E2);
%             [~,EC_n]=bwlabeln(EC,conn);
%             if EC_n<2
%                 E1=ConnectConCsToElectrode(E1concs,conn,2);
%                 E2=ConnectConCsToElectrode(E2concs,conn,2);
%             end
            [~,EC_n]=bwlabeln(EC,conn);
            if EC_n<2
                [~,~,~,~,~,E1endpoints_sub]=ConcEndpsDetection2(E1);
                P1=E1endpoints_sub(1,:);
                Q1=E1endpoints_sub(2,:);
                [~,~,~,~,~,E2endpoints_sub]=ConcEndpsDetection2(E2);
                P2=E2endpoints_sub(1,:);
                Q2=E2endpoints_sub(2,:);
                VolumeSize=[nx,ny,nz];
                [P1Q1line,~,~,~]=MakeLine2(P1,Q1,VolumeSize);
                [P2Q2line,~,~,~]=MakeLine2(P2,Q2,VolumeSize);
                E1=P1Q1line;
                E2=P2Q2line;
            end
            
%             figure(fignum);PlotElectrodes(E1,'b');hold on;
%             figure(fignum);PlotElectrodes(E2,'r');hold on;
%             fignum=fignum+1;

            electrode12=logical(electrode12+E1);
            electrode12=logical(electrode12+E2);
        end
        
        if cesn>2
            fignum=fignum+20; %ConnectionJudge5可能产生较多图片
            e11conci_ind=find(e11conci);
            electrode12(e11conci_ind)=0;
            conc=e11conci;
            [nx,ny,nz]=size(conc);
            Eallconcs=DisconnectConCs(conc,maskedge,ces_all,origelec,conn,fignum);
            fignum=fignum+40; %DisconnectConCs可能产生较多图片
            
            EC=false(nx,ny,nz);
            Ej_rec=false(nx,ny,nz,cesn);
            for j=1:cesn
%                 figure(fignum);PlotElectrodes(Eallconcs{j},'b');hold on;
                Ej=ConnectConCsToElectrode(Eallconcs{j},conn,1);
                EC=logical(EC+Ej);
                Ej_rec(:,:,:,j)=Ej;
            end
%             fignum=fignum+1;
%             [~,EC_n]=bwlabeln(EC,conn);
%             if EC_n<cesn
%                 EC=false(size(conc));
%                 for j=1:cesn
%                     Ej=ConnectConCsToElectrode(Eallconcs{j},conn,2);
%                     EC=logical(EC+Ej);
%                 end
%             end
            [~,EC_n]=bwlabeln(EC,conn);
            if EC_n<cesn
                EC=false(nx,ny,nz);
                for j=1:cesn
                    Ej=Ej_rec(:,:,:,j);
                    [~,~,~,~,~,Ejendpoints_sub]=ConcEndpsDetection2(Ej);
                    Pj=Ejendpoints_sub(1,:);
                    Qj=Ejendpoints_sub(2,:);
                    VolumeSize=size(conc);
                    [PjQjline,~,~,~]=MakeLine2(Pj,Qj,VolumeSize);
                    Ej=PjQjline;
                    EC=logical(EC+Ej);
                end
            end
            
%             figure(fignum);PlotElectrodes(EC,'b');
%             fignum=fignum+1;
            electrode12=logical(electrode12+EC);
        end
        
    end
    
end

[~,electrode12_n]=bwlabeln(electrode12,conn)
figure(99);PlotElectrodes(electrode12,'b');
%%
waitbar(0.65,f,'CranioCerebral Interference Removing...');

[electrode12_L,~]=bwlabeln(electrode12,conn);
electrode13=false(size(electrode12));
mm=0;
for i=1:electrode12_n
    e12conci=electrode12_L;
    e12conci(e12conci~=i)=0;
    e12conci=logical(e12conci);
    e12ciIntc=logical(e12conci.*brainmask);
    [e12ciIntc_L,e12ciIntc_n]=bwlabeln(e12ciIntc,conn);
    if e12ciIntc_n>1
        mm=mm+1;
        e12ciIntc_PN=zeros(e12ciIntc_n,1);
        for j=1:e12ciIntc_n
            e12ciIntcj=e12ciIntc_L;
            e12ciIntcj(e12ciIntc_L~=j)=0;
            e12ciIntcj=logical(e12ciIntcj);
            e12ciIntc_PN(j)=length(find(e12ciIntcj));
        end
        L=find(e12ciIntc_PN==max(e12ciIntc_PN));
        e12ciIntc(e12ciIntc_L~=L)=0;
    end
    electrode13=logical(electrode13+e12ciIntc);
end
[~,electrode13_n]=bwlabeln(electrode13,conn)
figure(100);PlotElectrodes(electrode13,'b');

waitbar(0.68,f,'CranioCerebral Interference Removing...');

if electrode13_n>ElecNum
    VolumeIn=electrode13;Mask=brainmask3;conn=26;
    [VolumeIn_PN,VolumeOut]=RemoveInterference3(VolumeIn,Mask,conn);
    electrode13_PN=VolumeIn_PN;
    electrode14=VolumeOut;
    save electrode13_PN 'electrode13_PN';
else
    electrode14=electrode13;
end
[~,electrode14_n]=bwlabeln(electrode14,conn)
figure(101);PlotElectrodes(electrode14,'b');

waitbar(0.7,f,'CranioCerebral Interference Removing...');

if electrode14_n>ElecNum
    VolumeIn=electrode14;Mask=brainmask3;conn=26;
    [VolumeIn_explained,VIintc_explained,VolumeIn_PAL,VIintc_PAL,...,
        VIintc_s1ratio,VIc_s1ratio,VolumeOut]=RemoveInterference4(VolumeIn,Mask,conn);
    e14_explained=VolumeIn_explained;
    e14VIintc_explained=VIintc_explained;
    e14_PAL=VolumeIn_PAL;
    e14VIintc_PAL=VIintc_PAL;
    e14VIintc_s1ratio=VIintc_s1ratio;
    e14VIc_s1ratio=VIc_s1ratio;
    electrode15=VolumeOut;
%     save e14_explained 'e14_explained'
%     save e14VIintc_explained 'e14VIintc_explained'
%     save e14_PAL 'e14_PAL'
%     save e14VIintc_PAL 'e14VIintc_PAL'
%     save e14VIintc_s1ratio 'e14VIintc_s1ratio'
%     save e14VIc_s1ratio 'e14VIc_s1ratio'
else
    electrode15=electrode14;
end
[~,electrode15_n]=bwlabeln(electrode15,conn)
figure(102);PlotElectrodes(electrode15,'b');

[electrode15_L,~]=bwlabeln(electrode15,conn);
e15_PAL=zeros(electrode15_n,3);
for i=1:electrode15_n
    e15conci=electrode15_L;
    e15conci(electrode15_L~=i)=0;
    e15conci=logical(e15conci);
    e15ciPAL=regionprops3(e15conci,'PrincipalAxisLength');
    e15_PAL(i,:)=e15ciPAL.PrincipalAxisLength;
end

waitbar(0.72,f,'CranioCerebral Interference Removing...');

electrode16=electrode15;fignum=103;conn=26;
for i=1:electrode15_n
    if e15_PAL(i,2)>7.5 %该值可能还需要调整
        e15conci=electrode15_L;
        e15conci(electrode15_L~=i)=0;
        e15conci=logical(e15conci);
%         figure(fignum);PlotElectrodes(e15conci,'b');
%         fignum=fignum+1;
        
        electrode16=logical(electrode16-e15conci);
        e15ciorg=logical(e15conci.*electrode5);
%         e15ciorg=bwareaopen(e15ciorg,2,conn); %去除由于electrode5的边缘不连续不光滑性而产生的零散小像素点
%         figure(fignum);PlotElectrodes(e15ciorg,'b');
%         fignum=fignum+1;
        
        e15cinew=ConnectConCsToElectrode(e15ciorg,conn,1);
%         figure(fignum);PlotElectrodes(e15cinew,'b');
%         fignum=fignum+1;
        
        electrode16=logical(electrode16+e15cinew);
    end
end

[~,electrode16_n]=bwlabeln(electrode16,conn)
figure(103);PlotElectrodes(electrode16,'b');

[electrode16_L,~]=bwlabeln(electrode16,conn);
if electrode16_n>ElecNum
    electrode17=false(size(electrode16));
    e16_PAL=zeros(electrode16_n,3);
    for i=1:electrode16_n
        e16conci=electrode16_L;
        e16conci(electrode16_L~=i)=0;
        e16conci=logical(e16conci);
        e16ciPAL=regionprops3(e16conci,'PrincipalAxisLength');
        e16_PAL(i,:)=e16ciPAL.PrincipalAxisLength;
    end
%     save e16_PAL 'e16_PAL';
    e16_PAL1=e16_PAL(:,1);
    [e16_PAL1s,e16_PAL1sL]=sort(e16_PAL1,'descend');
    e17L=e16_PAL1sL(1:ElecNum);
    for j=1:ElecNum
        electrode17(electrode16_L==e17L(j))=1;
    end
else
    electrode17=electrode16;
end
[~,electrode17_n]=bwlabeln(electrode17,conn)
figure(104);PlotElectrodes(electrode17,'b');

end