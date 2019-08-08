function [CranialConCs,InnerConCs,GrowProcess,ImageOut]=AllInnerConCsGrow3(origin,basis,elecenters,th,ratio,T0,Te,fignum)

%origin��������ֵ����ĳ�ʼͼ����electrode1
%basis����������������������ά����
%elecenters��ĳһ��ά���飨ͨ����basis����­�Ǳ߽�Ľ������֣���Ϊ�缫���
%th����ֵ�½��̶ȣ�th=2��ʾ�������ε���ֵ���2
%ratio����ʾ��­����ͨ������­�Ǳ߽罵��ֵ����ֱ������λ��­�Ǳ߽��ϵĵ缫�������������ҽ���������������ռ�õ缫��������������
%       ratio��ʱ����­����ͨ��������­�Ǳ߽��ϵĵ缫�����Ϊһ���µĵ缫��㣻����������������­����ͨ�����������Ĺ�������������
%       һ��ʱ������Щ­����ͨ�����ϲ���Ϊһ���µ�­���������ӵ㣬ֱ�����µ�­���������ӵ�Ҳ��­�Ǳ߽��ϵĵ缫������������
%T0����ֵ��ֵ
%Te����ֵ��ֵ
%fignum�������ĵ�һ��ͼƬ�����

%CranialConCs��basis����ȫλ��­�Ǳ߽��ϻ���­�Ǳ߽���ڽ�������ͨ��������Ϊ­����ͨ����
%InnerConCs��basis����ȫλ��­�ڵ���­�Ǳ߽粻���ڽ�������ͨ��������Ϊ­����ͨ����
%GrowProcess����¼�����������̣���N*4��Ԫ�����飬ÿһ�зֱ��ǵ�ǰ���µĵ缫��㡢��ǰ���µ�­���������ӡ��ô�������������ֵ��ֵ���ô�������������ֵ��ֵ
%ImageOut��basis����������õ��Ľ��������ĵ缫��������󲿷ֵ缫��­����ͨ�����������߽�

if nargin<8
    fignum=1;
end

if nargin<6
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<4
    th=2;
end

% [basisConC_Sub,elecentersConC_Sub,basis_L,elecenters_L,elecenters_newL,elecenters_recL,basis_n,elecenters_n]=CountElectrodes(basis,elecenters,26);
% basis_recL=(1:1:basis_n)';
% InnerConCs_recL=setdiff(basis_recL,elecenters_recL);
% 
% InnerConCs=basis_L;
% for i=1:length(elecenters_recL)
%     L=elecenters_recL(i);
%     InnerConCs(InnerConCs==L)=0;
% end

basis=logical(basis);
elecenters=logical(elecenters);
[nx,ny,nz]=size(basis);
InnerConCs=false(nx,ny,nz);

conn=26;
[basis_L,basis_n]=bwlabeln(basis,conn);

for i=1:basis_n
    bconci=basis_L;
    bconci(bconci~=i)=0;
    IS=logical(bconci.*elecenters);
    if sum(IS==1)==0
        InnerConCs=logical(InnerConCs+bconci);
    end
end
CranialConCs=logical(basis-InnerConCs);
figure(fignum);PlotElectrodes(CranialConCs,'b');
figure(fignum+1);PlotElectrodes(InnerConCs,'b');
fignum=fignum+2;

seed=InnerConCs;
crconc=CranialConCs;

[~,crconc_n,crconc_sub]=FindIndexOfAllConnectedComponents(crconc,conn);
time=4;
T1=T0-th;
T2=T1-time*th;
N=0;

while T2>=Te
    
    % ע��result=MyRegionGrow4(origin,seed0,th,T1,T2)��������������T1��T2���
    result=MyRegionGrow4(origin,seed,th,T1,T2);
    seed2=result;
    [seed2_L,seed2_n]=bwlabeln(seed2,conn);
    
    for i=1:seed2_n
        seed2conci=seed2_L;
        seed2conci(seed2conci~=i)=0;
        seed2conci_ind=find(seed2conci);
        [x,y,z]=ind2sub(size(seed2conci),seed2conci_ind);
        seed2conci_sub=[x,y,z];
        for j=1:crconc_n
            IS=intersect(seed2conci_sub,crconc_sub{j},'rows');
            if size(IS,1)/size(crconc_sub{j},1)>ratio
                N=N+1;
                seed2(seed2conci_ind)=0;
                crconc=logical(crconc+seed2conci);
                [~,crconc_n,crconc_sub]=FindIndexOfAllConnectedComponents(crconc,conn);
                figure(fignum);PlotElectrodes(crconc,'b');
                figure(fignum+1);PlotElectrodes(seed2,'b');
                fignum=fignum+2;
                GrowProcess{N,1}=crconc;
                GrowProcess{N,2}=seed2;
                GrowProcess{N,3}=T1;
                GrowProcess{N,4}=T2;
                break;
            end
        end
    end
    seed=seed2;
    [~,seed_n]=bwlabeln(seed,conn);
    if seed_n==0
        break;
    end
    
    T1=T2-th;
    T2=T1-time*th;
    if T2<Te
        T2=Te;
    end
    
end

ImageOut=logical(crconc+basis);

end