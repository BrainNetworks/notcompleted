function [brainmask,brainmask2,brainmask3,brainmask_edge,brainmask_edge2,brainmask_edge3,T0,Te,...
    electrode1,electrode1B,electrode2,electrode3,electrode4,imclose_judge,imclose_r,...
    VIclose_nrec,electrode5,electrode5_S]=ContactSeg1(ElecNum,regCT,mask2regCT,f)

waitbar(0.02,f,'Extracranial tissue stripping...');

brainmask=logical(mask2regCT);
[r1,c1,p1]=size(brainmask);

brainmask_edge=brainmask;
for i=1:p1
    brainmask_edge(:,:,i)=edge(brainmask(:,:,i),'Canny');
end

nhood=ones(20,20); %创建一个二维结构元素
brainmask_edge2=brainmask_edge;
for i=1:p1
    brainmask_edge2(:,:,i)=imdilate(brainmask_edge(:,:,i),nhood);
end

brainmask_edge3=logical(brainmask_edge2.*(~brainmask));
brainmask2=logical(brainmask+brainmask_edge3);
brainmask3=logical(brainmask2-brainmask_edge2);

electrode1=regCT.*brainmask2;
electrode1B=regCT.*brainmask;

waitbar(0.05,f,'Threshold processing...');

conn=26;
T0=round(max(max(max(electrode1)))*0.7); %阈值初值
Te=round(max(max(max(electrode1)))*0.01); %阈值终值
electrode2=electrode1;
electrode2(electrode1>T0)=T0;
electrode2(electrode1<T0)=0;
electrode2=logical(electrode2);
[~,electrode2_n]=bwlabeln(electrode2,conn)
figure(57);PlotElectrodes(electrode2);

electrode3=electrode2.*brainmask_edge2;
[~,electrode3_n]=bwlabeln(electrode3,conn)
% figure(58);PlotElectrodes(electrode3);
% electrode3_T0=electrode3;

electrode4=electrode2.*logical(brainmask2-brainmask_edge2);
[~,electrode4_n]=bwlabeln(electrode4,conn)
% figure(59);PlotElectrodes(electrode4);

waitbar(0.08,f,'Morphological closing determination...');

VolumeIn=electrode4;conn=26;fignum=60;
[VIclose,imclose_judge,imclose_r,VIclose_nrec]=ImcloseJudge(VolumeIn,ElecNum,conn,fignum);

if imclose_judge==1
    SE=strel('sphere',imclose_r); %Create a 3-D sphere-shaped structuring element with a radius of r.
    electrode5=imclose(electrode2,SE);
    maskedge=brainmask_edge2;origin=regCT;mask=brainmask2;conn=26;fignum=60;
    [electrode5_L,electrode5_n]=bwlabeln(electrode5,conn);
    e5cesn_rec=zeros(electrode5_n,1);
    for i=1:electrode5_n
        e5conci=electrode5_L;
        e5conci(electrode5_L~=i)=0;
        e5conci=logical(e5conci);
        conc=e5conci;
        [~,cesn,~]=ConnectionJudge5(conc,maskedge,origin,mask,T0,conn,fignum);
        e5cesn_rec(i)=cesn;
        if cesn>=3
            imclose_judge=0;
            electrode5=electrode2;
        end
    end
else
    electrode5=electrode2;
end

[~,electrode5_n]=bwlabeln(electrode5,conn)
% figure(60);PlotElectrodes(electrode5);

E1=electrode5;E2=electrode3;conn=26;
[~,~,~,~,E1_S,E2_newL,E2_recL]=InnerConCsDetection(E1,E2,conn);
electrode5_S=E1_S;
% electrode3_newL=E2_newL;
% electrode3_recL=E2_recL;
% electrode3_checkL=unique(electrode3_recL); %去除重复元素并升序排序
% electrode3_judgeL=(length(electrode3_recL)~=length(electrode3_checkL));

end