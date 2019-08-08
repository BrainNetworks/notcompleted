function [VIclose,imclose_judge,imclose_r,VIclose_nrec]=ImcloseJudge(VolumeIn,ElecNum,conn,fignum)

%VolumeIn���������ά�����ݣ�ͨ����ȥ����­�Ǳ߽��ĳ����缫��ȡ���electrode4
%ElecNum���û��ߵ�­�ڵ缫��������ֵ��Ҫ��ҽ����������
%conn����ͨ��

%VIclose������VolumeIn�ղ����Ľ��
%imclose_judge��judge=1�������Ҫ�Ե缫��ͨ���ǻ�δȥ��­�Ǳ߽��electrode2��ִ�бղ�����judge=0ʱ��������ܽ��бղ���
%imclose_r����judge=1ʱ��r�ǶԵ缫��ͨ���ǻ�δȥ��­�Ǳ߽��electrode2��ִ�бղ���ʱʹ�õ����νṹԪ����Ѱ뾶����judge=0ʱr������
%VIclose_nrec������VolumeInÿ�ν��бղ�����Ľ��VIclose����ͨ������ĿVIclose_n�ļ�¼����һ��������

if nargin<4
    fignum=1;
end

if nargin<3
    conn=26;
end

VolumeIn=logical(VolumeIn);

imclose_judge=0;
r_all=[3,6,8,9,10];
rN=length(r_all);
VIclose_nrec=zeros(rN,1);

for i=1:rN
    imclose_r=r_all(i);
    SE=strel('sphere',imclose_r); %Create a 3-D sphere-shaped structuring element with a radius of r.
    VIclose=imclose(VolumeIn,SE);
    [VIclose_L,VIclose_n]=bwlabeln(VIclose,conn);
    VIclose_nrec(i)=VIclose_n;
    if VIclose_n<ElecNum
        break;
    end
    if (VIclose_n>=ElecNum)&&(VIclose_n<=round(ElecNum*1.3))
        imclose_judge=1;
        if VIclose_n==ElecNum
            break;
        else
            L=find(VIclose_nrec==min(VIclose_nrec),1);
            imclose_r=r_all(L);
        end
    end
end

if (VIclose_n>round(ElecNum*1.3))&&(VIclose_n<=round(ElecNum*2))
    for i=1:VIclose_n
        VIcConCi=VIclose_L;
        VIcConCi(VIclose_L~=i)=0;
        VIcConCi=logical(VIcConCi);
        VIC=logical(VIcConCi.*VolumeIn);
        [VIC_L,VIC_n]=bwlabeln(VIC,conn);
        if (VIC_n>=5)&&(VIC_n<=18)
            PixelNum=zeros(VIC_n,1);
            CenterSub=zeros(VIC_n,3);
%             figure(fignum);PlotElectrodes(VIC,'b');hold on;
            for j=1:VIC_n
                VICj=VIC_L;
                VICj(VIC_L~=j)=0;
                VICj=logical(VICj);
                VICj_PixelNum=length(find(VICj));
                PixelNum(j)=VICj_PixelNum;
                VICj_Centroid=regionprops3(VICj,'Centroid');
                VICj_Centroid=VICj_Centroid{1,1};
                VICj_Centrsub=VICj_Centroid(:,[2,1,3]);
                CenterSub(j,:)=VICj_Centrsub;
            end
            R=CenterSub(:,1);
            C=CenterSub(:,2);
            S=CenterSub(:,3);
            Rdiff=max(R)-min(R);
            Cdiff=max(C)-min(C);
            Sdiff=max(S)-min(S);
            diff=[Rdiff,Cdiff,Sdiff];
            if max(diff)==Rdiff
                [~,RL]=sort(R,'ascend');
                PixelNum=PixelNum(RL,:);
                CenterSub=CenterSub(RL,:);
            end
            if max(diff)==Cdiff
                [~,CL]=sort(C,'ascend');
                PixelNum=PixelNum(CL,:);
                CenterSub=CenterSub(CL,:);
            end
            if max(diff)==Sdiff
                [~,SL]=sort(S,'ascend');
                PixelNum=PixelNum(SL,:);
                CenterSub=CenterSub(SL,:);
            end
            PixelNum2=PixelNum(2:(end-1),:);
            CenterSub2=CenterSub(2:(end-1),:);
%             figure(fignum);plot3(CenterSub2(:,2),CenterSub2(:,1),CenterSub2(:,3),'o','color','r');hold on;
%             fignum=fignum+1;            
            PQD=zeros(VIC_n-3,1);
            for k=1:(VIC_n-3)
                P=CenterSub2(k,:);
                Q=CenterSub2(k+1,:);
                PQDk=norm(P-Q);
                PQD(k)=PQDk;
            end
            PixelNum_std=std(PixelNum2./max(PixelNum2));
            PQD_std=std(PQD./max(PQD));
            if (PixelNum_std<0.3)||(PQD_std<0.2)
                imclose_judge=1;
                break;
            end
        end
    end
end

% if (VIclose_n>round(ElecNum*1.3))&&(VIclose_n<=round(ElecNum*2))
%     CN_rec=[];
%     PNstd_rec=[];
%     PQDstd_rec=[];
%     for i=1:VIclose_n
%         VIcConCi=VIclose_L;
%         VIcConCi(VIclose_L~=i)=0;
%         VIcConCi=logical(VIcConCi);
%         VIC=logical(VIcConCi.*VolumeIn);
%         [VIC_L,VIC_n]=bwlabeln(VIC,conn);
%         if (VIC_n>=5)&&(VIC_n<=18)
%             CN_rec=[CN_rec;VIC_n];
%             PixelNum=zeros(VIC_n,1);
%             CenterSub=zeros(VIC_n,3);
%             figure(fignum);PlotElectrodes(VIC,'b');hold on;
%             for j=1:VIC_n
%                 VICj=VIC_L;
%                 VICj(VIC_L~=j)=0;
%                 VICj=logical(VICj);
%                 VICj_PixelNum=length(find(VICj));
%                 PixelNum(j)=VICj_PixelNum;
%                 VICj_Centroid=regionprops3(VICj,'Centroid');
%                 VICj_Centroid=VICj_Centroid{1,1};
%                 VICj_Centrsub=VICj_Centroid(:,[2,1,3]);
%                 CenterSub(j,:)=VICj_Centrsub;
%             end
%             R=CenterSub(:,1);
%             C=CenterSub(:,2);
%             S=CenterSub(:,3);
%             Rdiff=max(R)-min(R);
%             Cdiff=max(C)-min(C);
%             Sdiff=max(S)-min(S);
%             diff=[Rdiff,Cdiff,Sdiff];
%             if max(diff)==Rdiff
%                 [~,RL]=sort(R,'ascend');
%                 PixelNum=PixelNum(RL,:);
%                 CenterSub=CenterSub(RL,:);
%             end
%             if max(diff)==Cdiff
%                 [~,CL]=sort(C,'ascend');
%                 PixelNum=PixelNum(CL,:);
%                 CenterSub=CenterSub(CL,:);
%             end
%             if max(diff)==Sdiff
%                 [~,SL]=sort(S,'ascend');
%                 PixelNum=PixelNum(SL,:);
%                 CenterSub=CenterSub(SL,:);
%             end
%             PixelNum2=PixelNum(2:(end-1),:);
%             CenterSub2=CenterSub(2:(end-1),:);
%             figure(fignum);plot3(CenterSub2(:,2),CenterSub2(:,1),CenterSub2(:,3),'o','color','r');hold on;
%             fignum=fignum+1;            
%             PQD=zeros(VIC_n-3,1);
%             for k=1:(VIC_n-3)
%                 P=CenterSub2(k,:);
%                 Q=CenterSub2(k+1,:);
%                 PQDk=norm(P-Q);
%                 PQD(k)=PQDk;
%             end
%             PixelNum_std=std(PixelNum2./max(PixelNum2));
%             PQD_std=std(PQD./max(PQD));
%             PNstd_rec=[PNstd_rec;PixelNum_std];
%             PQDstd_rec=[PQDstd_rec;PQD_std];
%         end
%     end
%     stats=[CN_rec,PNstd_rec,PQDstd_rec];
%     if (sum(stats(:,2)<0.3)~=0)||(sum(stats(:,3)<0.2)~=0)
%         judge=1;
%     end
% end

end