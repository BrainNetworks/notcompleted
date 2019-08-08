function [InnerConC,label1,label2,ImageOut,CompStats,CompConCs]=ElecEnterpsCompletion3(ImageIn,origin,mask,ratio2,ratio3,th,T0,Te,fignum,conn)

%ImageIn������ĵ缫���ݣ�ͨ��Ϊ����­����ͨ��������������Ľ����������AllInnerConCsGrow�����ImageOut
%origin��������ֵ����ĳ�ʼͼ����regCT,electrode1�����ã�,electrode1B=regCT.*brainmask��
%mask��������­�Ǳ߽��­��������ģ����brainmask3������Ѱ����ȫλ��­�Ǳ߽��ϵ���ͨ����
%ratio2����ʾ����ImageIn��­�Ǳ߽��ϵ�С��ͨ������­�ڵ缫����ֵ����ֱ������­�ڵ缫�����������ҽ���������������ռ�õ缫��������������
%        ����ratio2����һ��������һ����ʱ�����̱����ʱ�������������seed����������Ϊseed_ratio2��seed_ratio2���ǵ缫�����޲����
%ratio3����ʾ����ImageIn��­�Ǳ߽��ϵ�С��ͨ������­�ڵ缫����ֵ����ֱ������­�ڵ缫�����������ҽ���������������ռ�õ缫��������������
%        ����ratio3����������һ��������Ϊ��­�Ǳ߽��ϵ�С��ͨ�������ڵ缫��һ���֣�����Ϊ���ǵ缫���Ͽ�����ͨ����֮һ
%th����ֵ�½��̶ȣ�th=2��ʾ�������ε���ֵ���2
%T0����ֵ��ֵ
%Te����ֵ��ֵ
%fignum�������ĵ�һ��ͼƬ�����
%conn����ͨ��

%InnerConC��������ImageIn����ȫλ��­�Ǳ߽��ϵ���ͨ����ȥ����ʣ�µĲ���
%label1������ImageIn����ȫλ��­�Ǳ߽��ϵ���ͨ�����ı�ţ���������
%label2�������label1�е�������������500�����䱻origin�������Ĳ��ֵ�������ҲС��100����ͨ�����ı�ţ�����������label2�е���ͨ�����ſ���������޲�
%ImageOut��������ImageIn��­�Ǳ߽������ڵ缫�ĵ�ȴ��缫�Ͽ���С��ͨ�������ӵ��缫�ϣ�ͬʱȥ��­�Ǳ߽��ϲ����ڵ缫�Ĳ���­�Ǹ��ź�õ��Ľ��
%CompStats������ImageIn��­�Ǳ߽��ϵ�С��ͨ�������������ڵ缫��Ҳ�����ǲ���­�Ǹ��ţ�seed0��ͳ����Ϣ����һ��n*4�ľ���n��С��ͨ��������Ŀ
%CompStats(:,1)��ÿ��­�Ǳ߽��ϵ�С��ͨ����seed0������ImageIn�е���ͨ�������
%CompStats(:,2)�Ƕ�����ImageIn��ÿ��­�Ǳ߽��ϵ�С��ͨ����seed0�Ƿ����ڵ缫���оݣ�ֵΪ1�������ǡ���ֵΪ0��������
%Compstats(:,3)�Ƕ�����ImageIn��ÿ��­�Ǳ߽��ϵ�С��ͨ����seed0�ǽ��з��࣬Compstats(:,3)ֵ��ͬ������С��ͨ�����������ڱ�ŵ��ڸ�ֵ��С��ͨ��������������У���ֵΪ0��ʾ�����ǵ缫��
%CompStats(:,4)������ImageIn��ÿ��­�Ǳ߽��ϵ�С��ͨ����seed0������������ֵ��ֵ
%CompConCs������ImageIn��­�Ǳ߽��ϵ�С��ͨ����seed0���޲��������һ��4ά����
%��1�������i��­�Ǳ߽��ϵ�С��ͨ����seed0����ĳһ­�ڵ缫����CompConCs(:,:,:,i)������­����������ֱ������ratio2ʱ�õ��Ľ��seed_ratio2
%��2�������i��­�Ǳ߽��ϵ�С��ͨ����seed0������ĳһ­�ڵ缫�����ǲ���­�Ǹ��ţ���CompConCs(:,:,:,i)=seed0

if nargin<10
    conn=26;
end

if nargin<9
    fignum=1;
end

if nargin<7
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<6
    th=2;
end

if nargin<5
    ratio3=0.003;
end

if nargin<4
    ratio2=0.0003;
end

ImageIn=logical(ImageIn);
[ImageIn_L,ImageIn_n]=bwlabeln(ImageIn,conn);
%ImageIn_pn������ImageIn��������ͨ��������������
% ImageIn_pn=zeros(ImageIn_n,1);
% for i=1:ImageIn_n
%     ImageIn_pn(i)=sum(sum(sum(ImageIn_L==i)));
% end
% ImageIn_pns=sort(ImageIn_pn,'ascend');
% ImageIn_pns_mean=mean(ImageIn_pns(2:(end-1)));

label1=[];
label2=[];
for i=1:ImageIn_n
    IConCi=ImageIn_L;
    IConCi(ImageIn_L~=i)=0;
    IConCi=logical(IConCi);
    IS1=IConCi.*logical(mask);
    if isempty(find(IS1))
        label1=[label1;i];
        IConCiog=logical(IConCi.*origin); %ע��˴�IConCi.*origin��origin=electrode1B=regCT.*brainmask
        IConCiog_pn=length(find(IConCiog));
        IConCi_pn=length(find(IConCi));
        if (IConCiog_pn<200)&&(IConCi_pn<600)
            label2=[label2;i];
        end
    end
end
ll1=length(label1);
ll2=length(label2);

CompStats=zeros(ll2,4);
CompStats(:,1)=label2;
[nx,ny,nz]=size(ImageIn);
CompConCs=false(nx,ny,nz,ll2);

label1ConC=false(nx,ny,nz);
for i=1:ll1
    label1ConC(ImageIn_L==label1(i))=1;
end
label1ConC=logical(label1ConC);

InnerConC=logical(ImageIn-label1ConC);
[~,InnerConC_n,InnerConC_sub]=FindIndexOfAllConnectedComponents(InnerConC,conn);

label2ConC=false(nx,ny,nz);
for i=1:ll2
    label2ConC(ImageIn_L==label2(i))=1;
end
label2ConC=logical(label2ConC);

figure(fignum);PlotElectrodes(label1ConC,'b');
fignum=fignum+1;
figure(fignum);PlotElectrodes(label2ConC,'b');
fignum=fignum+1;
figure(fignum);PlotElectrodes(InnerConC,'b');
fignum=fignum+1;

T=(T0-th):(-th):Te;
if T(end)~=Te
    T=[T,Te];
end
Tlen=length(T);

ImageOut=ImageIn;
flag3=[];

h=waitbar(0,'Electrodes enter points Completing...'); %����������

for i=1:ll2
    
    k=label2(i);
    
    if isempty(intersect(k,flag3))

        flag=0;
        flag2=0;
        seed0=ImageIn_L;
        seed0(seed0~=k)=0;
        seed0=logical(seed0);
        seed=logical(seed0.*origin); %ע��˴�seed0.*origin
        figure(fignum);PlotElectrodes(seed,'b');
        fignum=fignum+1;

        for m=1:Tlen

            seed_ind=find(seed);
            [x,y,z]=ind2sub(size(seed),seed_ind);
            seed_sub=[x,y,z];
            seed_neighb_sub=[x-1,y,z;
                             x+1,y,z;
                             x,y-1,z;
                             x,y+1,z;
                             x,y,z-1;
                             x,y,z+1];
            seed_neighb_sub=unique(seed_neighb_sub,'rows');
            seed_neighb_sub=setdiff(seed_neighb_sub,seed_sub,'rows');
            seed_neighb_ind=sub2ind(size(seed),seed_neighb_sub(:,1),seed_neighb_sub(:,2),seed_neighb_sub(:,3));
            
            seed2=false(size(seed));
            seed2(seed_ind)=1;
            seed2(seed_neighb_ind)=1;

            t=T(m);
            for l=1:length(seed_neighb_ind)
                if origin(seed_neighb_ind(l))<t
                    seed2(seed_neighb_ind(l))=0;
                end
            end
            
            seed=logical(seed2);
            seed_ind2=find(seed);
            [x2,y2,z2]=ind2sub(size(seed),seed_ind2);
            seed_sub2=[x2,y2,z2];
            
            for j=1:InnerConC_n
                IS2=intersect(seed_sub2,InnerConC_sub{j},'rows');
                if (flag2==0)&&(size(IS2,1)/size(InnerConC_sub{j},1)>ratio2)
                    seed_grow=seed;
                    flag2=1;
                    figure(fignum);PlotElectrodes(seed_grow,'b');
                    fignum=fignum+1;
                end
                if size(IS2,1)/size(InnerConC_sub{j},1)>ratio3
                    flag=1;
                    flag3=[flag3,k];
                    CompStats(i,2)=1;
                    CompStats(i,3)=k;
                    CompConCs(:,:,:,i)=seed_grow;
                    ImageOut=logical(ImageOut+seed_grow);
                    figure(fignum);PlotElectrodes(seed,'b');
                    fignum=fignum+1;
                    
                    for ne=(i+1):ll2
                        ConCj=ImageIn_L;
                        ConCj(ConCj~=label2(ne))=0;
                        ConCj_ind=find(ConCj);
                        IS3=intersect(seed_ind2,ConCj_ind);
                        if length(IS3)/length(ConCj_ind)>0.8
                            CompStats(ne,2)=1;
                            CompStats(ne,3)=k;
                            CompStats(ne,4)=t;
                            CompConCs(:,:,:,ne)=seed_grow;
                            flag3=[flag3,label2(ne)];
                        end
                    end
                    
                    break;
                end
            end

            if flag==1
                break;
            end
            
            seed_pn=length(find(seed));
            if (flag==0)&&(seed_pn>1300)
                CompConCs(:,:,:,i)=seed0;
                seed0_ind=find(seed0);
                ImageOut(seed0_ind)=0;
                break;
            end

        end
        CompStats(i,4)=t;

    else
        i=i+1;
    end
    
    waitbar(i/ll2,h); %�޸Ľ�����ʾ

end

close(h);

CompConCs=logical(CompConCs);
ImageOut=logical(ImageOut);

end