function [ces1,ces2,judge]=ConnectionJudge3(conc,maskedge,origin,mask,T0,conn,fignum)

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
ces=bwareaopen(ces,2,conn); %ȥ������maskedge�ı�Ե���������⻬�Զ���������ɢС���ص�

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