function [ces_all,cesn,judge]=ConnectionJudge4(conc,maskedge,origin,mask,T0,conn,fignum)

%conc������Ļ����������������ϵĵ缫����һ���ĳһ��ͨ����������ά����
%maskedge��­�Ǳ߽磬ͨ����brainmask_edge2
%origin��δ����ֵ����ĳ�ʼͼ�񣬴˴�ͨ����δȥ��­����ŵ�regCT
%mask��������ģ��ͨ����brainmask2��
%T0���������ж�origin����ͳһ��ֵ����ʱ��ʹ�õ���ֵ
%conn����ͨ��
%fignum�������ĵ�һ��ͼƬ�����

%ces_all��������ͨ����conc��ȫ���缫��㣬����ά���飬ces_all(:,:,:,i)�ǵ�i���缫���
%cesn:�缫�����������һ��������ֵ
%judge��judge=1����������ͨ����concȷʵ���������������ϵĵ缫����һ��

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
                if ocij_n==1 %ע�ⲻ��ocij_n==orgTcices_n
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