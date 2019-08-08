function [ces1,ces2,judge]=ConnectionJudge(conc,maskedge,origin,mask,T0,conn,fignum)

%conc������Ļ����������缫����һ���ĳһ��ͨ����������ά����
%maskedge��­�Ǳ߽磬ͨ����brainmask_edge2
%origin��δ����ֵ����ĳ�ʼͼ�񣬴˴�ͨ����δȥ��­����ŵ�regCT
%mask��������ģ��ͨ����brainmask2��
%T0����ֵ��ֵ��Ҳ���Ƕ�origin����ͳһ��ֵ����ʱ��ʹ�õ���ֵ
%conn����ͨ��
%fignum�������ĵ�һ��ͼƬ�����

%ces1��������ͨ����conc��­�Ǳ߽�maskedge�ĵ�һ��������Ҳ�ǵ�һ���缫����㣬����ά������
%ces2��������ͨ����conc��­�Ǳ߽�maskedge�ĵڶ���������Ҳ�ǵڶ����缫����㣬����ά������
%judge��judge=1����������ͨ����concȷʵ���������缫����һ��

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
%ע�⣺�˴�������ces=logical(conc.*elecenters)������ces����Ϊ��ʱ����õ��Ľ��
%�ǳ�ʼʱ�ĵ缫��㣬û�о���������������ղ�������������������յ�ȥ�����ţ�����
%һ���õ���ces�ܿ��ܲ�ֹ������ͨ�������Ӷ����ж����Ӱ�졣
ces=bwareaopen(ces,3,conn); %ȥ������maskedge�ı�Ե���������⻬�Զ���������ɢС���ص�

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
        ces1_Centrsub=ces1_Centroid(:,[2,1,3]); %����ces1_Centrsubȡ������������õ����Ӿ�ȷ��ces12_coeff1
        ces2_Centrsub=ces2_Centroid(:,[2,1,3]);
        figure(fignum);plot3(ces1_Centrsub(2),ces1_Centrsub(1),ces1_Centrsub(3),'o','color','r');hold on;
        figure(fignum);plot3(ces2_Centrsub(2),ces2_Centrsub(1),ces2_Centrsub(3),'o','color','r');hold on;
        fignum=fignum+1;

        ces12=ces1_Centrsub-ces2_Centrsub;
        ces12_d=norm(ces12);
        ces12_coeff1=(ces12./ces12_d)';
        diff=norm(conc_coeff1-ces12_coeff1); %ֵ����0~sqrt(6)=2.4495֮��
        
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
        orTc_outmask=bwareaopen(orTc_outmask,3,conn); %ȥ������maskedge�ı�Ե���������⻬�Զ���������ɢС���ص�
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