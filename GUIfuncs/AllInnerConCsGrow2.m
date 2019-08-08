function [ImageOut,InCs_stats,InCs_grow]=AllInnerConCsGrow2(origin,basis,basis_S,elecenters_n,elecenters_sub,th,ratio,T0,Te,fignum)

%origin��������ֵ����ĳ�ʼͼ����electrode1
%basis����������������������ά����
%basis_S��basis��������ͨ������ͳ����Ϣ
%elecenters_n���缫���elecenters����ͨ��������
%elecenters_sub���缫����������ͨ��������ά�±���������һ��n��1�е�Ԫ������
%th����ֵ�½��̶ȣ�th=2��ʾ�������ε���ֵ���2
%ratio����ʾ­����ͨ������­�Ǳ߽罵��ֵ����ֱ������λ��­�Ǳ߽��ϵĵ缫�������������ҽ���������������ռ�õ缫��������������ratio��
%T0����ֵ��ֵ
%Te����ֵ��ֵ
%fignum�������ĵ�һ��ͼƬ�����

%ImageOut��basis����������õ��Ľ��������ĵ缫��������е缫��­����ͨ�����������߽磬ͬʱ­�ڲ�����ű�ȥ��
%InCs_stats��basis������­����ͨ������ͳ����Ϣ����n*5�ľ���nΪbasis��­����ͨ��������
%InCs_stats(:,1)��basis��­����ͨ�������
%InCs_stats(:,2)�Ƕ�basis��ÿ��­����ͨ�����Ƿ����ڵ缫���оݣ�ֵΪ1�������ǡ���ֵΪ0��������
%InCs_stats(:,3)�Ƕ�basis��­����ͨ�������з��࣬InCs_stats(:,3)ֵ��ͬ��������ͨ�����������ڱ�ŵ��ڸ�ֵ����ͨ��������������У���ֵΪ0��ʾ�����ǵ缫��
%InCs_stats(:,4)��basis��ÿ��­����ͨ����������������ֵ��ֵ
%InCs_stats(:,5)�Ƕ�basis��­����ͨ�������������õ��Ľ������pca�����õ��ĵ�һ���ɷֹ����ʵİٷֱ�ֵ
%InCs_grow����һ��4ά���飬InCs_grow(:,:,:,i)��Ӧ��i��­����ͨ����InnerConCi���������õ��Ľ����������MyRegionGrow�����result
%��1�������i��­����ͨ����InnerConCi����ĳһ�缫����InCs_grow(:,:,:,i)��InnerConCi���������õ��������缫���
%��2�������i��­����ͨ����InnerConCi������ĳһ�缫������­�ڸ��ţ���InCs_grow(:,:,:,i)=InnerConCi

if nargin<10
    fignum=1;
end

if nargin<8
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<7
    ratio=0.05;
end

if nargin<6
    th=2;
end

InnerConC_L=[];
conn=26;
basis=logical(basis);
[basis_L,basis_n]=bwlabeln(basis,conn);
for j=1:basis_n
    if basis_S(j,2)==0 %­����ͨ����ָ��­�Ǳ߽���ȫû�н�������ͨ����
        InnerConC_L=[InnerConC_L;j];
    end
end

ICnum=length(InnerConC_L);
[nx,ny,nz]=size(basis);
InCs_grow=false(nx,ny,nz,ICnum);
InCs_judge=zeros(ICnum,1);
InnerConC_EL=zeros(ICnum,1);
InCs_grow_t=zeros(ICnum,1);
InCs_grow_latentp1=zeros(ICnum,1);
InCs_stats=zeros(ICnum,5);
ImageOut=basis;
flag=[];

h=waitbar(0,'Region Growing...'); %����������

for i=1:ICnum
    
    if isempty(intersect(InnerConC_L(i),flag))
        InnerConCi=basis_L;
        InnerConCi(InnerConCi~=InnerConC_L(i))=0;
        InnerConCi=logical(InnerConCi);
        InnerConCi_comp=basis_L;
        InnerConCi_comp(InnerConCi_comp==InnerConC_L(i))=0;
        InnerConCi_comp=logical(InnerConCi_comp);
        figure(fignum);
        PlotElectrodes(InnerConCi_comp,'b');hold on;
        PlotElectrodes(InnerConCi,'r');hold on;
        fignum=fignum+1;

        [InCi_grow,InCi_judge,InCi_grow_t,InCi_grow_latentp1]=MyRegionGrow3(origin,InnerConCi,elecenters_n,elecenters_sub,th,ratio,T0,Te);
        figure(fignum);PlotElectrodes(InCi_grow,'b');
        fignum=fignum+1;

        InCs_grow(:,:,:,i)=InCi_grow;
        InCs_judge(i)=InCi_judge;
        InCs_grow_t(i)=InCi_grow_t;
        InCs_grow_latentp1(i)=InCi_grow_latentp1;
        InnerConC_EL(i)=InnerConC_L(i);

        InCi_grow_ind=find(InCi_grow);
        if InCi_judge==1
            ImageOut=logical(ImageOut+InCi_grow);
        else
            ImageOut(InCi_grow_ind)=0;
        end
        flag=[flag,InnerConC_L(i)];
        
        for j=(i+1):ICnum
            InnerConCj=basis_L;
            InnerConCj(InnerConCj~=InnerConC_L(j))=0;
            InnerConCj_ind=find(InnerConCj);
            IS=intersect(InCi_grow_ind,InnerConCj_ind);
            if length(IS)/length(InnerConCj_ind)>0.8
                InCs_grow(:,:,:,j)=InCi_grow;
                InCs_judge(j)=InCi_judge;
                InCs_grow_t(j)=InCi_grow_t;
                InCs_grow_latentp1(j)=InCi_grow_latentp1;
                InnerConC_EL(j)=InnerConC_L(i);
                flag=[flag,InnerConC_L(j)];
            end
        end
    else
        i=i+1;
        waitbar(i/ICnum,h); %�޸Ľ�����ʾ
    end
    
    waitbar(i/ICnum,h); %�޸Ľ�����ʾ
            
end

close(h);

InCs_stats(:,1)=InnerConC_L;
InCs_stats(:,2)=InCs_judge;
InCs_stats(:,3)=InnerConC_EL;
InCs_stats(:,4)=InCs_grow_t;
InCs_stats(:,5)=InCs_grow_latentp1;

end