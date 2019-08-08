function [InnerConC,label1,label2,ImageOut,CompStats]=ElecEnterpsCompletion4(ImageIn,origin,mask,ratio,th,T0,Te,fignum,conn)

%ImageIn������ĵ缫���ݣ�ͨ��Ϊ����­����ͨ��������������Ľ����������AllInnerConCsGrow�����ImageOut
%origin��������ֵ����ĳ�ʼͼ����regCT,electrode1�����ã�,electrode1B=regCT.*brainmask��
%mask��������­�Ǳ߽��­��������ģ����brainmask3������Ѱ����ȫλ��­�Ǳ߽��ϵ���ͨ����
%ratio����ʾ����ImageIn��­�Ǳ߽��ϵ�С��ͨ������­�ڵ缫����ֵ����ֱ������­�ڵ缫�����������ҽ���������������ռ�õ缫��������������
%       ����ratio����������һ��������Ϊ��­�Ǳ߽��ϵ�С��ͨ�������ڵ缫��һ���֣���������Ϊ���ǵ缫���Ͽ�����ͨ����֮һ
%th����ֵ�½��̶ȣ�th=2��ʾ�������ε���ֵ���2
%T0����ֵ��ֵ
%Te����ֵ��ֵ
%fignum�������ĵ�һ��ͼƬ�����
%conn����ͨ��

%InnerConC��������ImageIn����ȫλ��­�Ǳ߽��ϵ���ͨ����ȥ����ʣ�µĲ���
%label1������ImageIn����ȫλ��­�Ǳ߽��ϵ���ͨ�����ı�ţ���������
%label2�������label1�е�������������1500����ͨ�����ı�ţ�����������label2�е���ͨ�����ſ���������޲�
%ImageOut��������ImageIn��­�Ǳ߽������ڵ缫�ĵ�ȴ��缫�Ͽ���С��ͨ�������ӵ��缫�ϣ�ͬʱȥ��­�Ǳ߽��ϲ����ڵ缫�Ĳ���­�Ǹ��ź�õ��Ľ��
%CompStats������ImageIn��­�Ǳ߽��ϵ�С��ͨ�������������ڵ缫��Ҳ�����ǲ���­�Ǹ��ţ�seed0��ͳ����Ϣ����һ��n*4�ľ���n��С��ͨ��������Ŀ
%CompStats(:,1)��ÿ��­�Ǳ߽��ϵ�С��ͨ����seed0������ImageIn�е���ͨ�������
%CompStats(:,2)�Ƕ�����ImageIn��ÿ��­�Ǳ߽��ϵ�С��ͨ����seed0�Ƿ����ڵ缫���оݣ�ֵΪ1�������ǡ���ֵΪ0��������
%Compstats(:,3)�Ƕ�����ImageIn��ÿ��­�Ǳ߽��ϵ�С��ͨ����seed0�ǽ��з��࣬Compstats(:,3)ֵ��ͬ������С��ͨ�����������ڱ�ŵ��ڸ�ֵ��С��ͨ��������������У���ֵΪ0��ʾ�����ǵ缫��
%CompStats(:,4)������ImageIn��ÿ��­�Ǳ߽��ϵ�С��ͨ����seed0������������ֵ��ֵ

if nargin<9
    conn=26;
end

if nargin<8
    fignum=1;
end

if nargin<6
    T0=round(max(max(max(origin)))*0.6);
    Te=round(max(max(max(origin)))*0.01);
end

if nargin<5
    th=2;
end

if nargin<4
    ratio=0.001;
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
        IConCi_pn=length(find(IConCi));
%         IConCiog=logical(IConCi.*origin); %ע�⵱origin=electrode1B=regCT.*brainmask���˴���IConCi.*origin
%         IConCiog_pn=length(find(IConCiog));
%         if (IConCiog_pn<400)&&(IConCi_pn<1500)
        if IConCi_pn<1500
            label2=[label2;i];
        end
    end
end
ll1=length(label1);
ll2=length(label2);
label2_renew=label2;

CompStats=zeros(ll2,4);
CompStats(:,1)=label2;
[nx,ny,nz]=size(ImageIn);

label1ConC=false(nx,ny,nz);
for i=1:ll1
    label1ConC(ImageIn_L==label1(i))=1;
end
label1ConC=logical(label1ConC);

InnerConC=logical(ImageIn-label1ConC);
[InnerConC_L,InnerConC_n,InnerConC_sub]=FindIndexOfAllConnectedComponents(InnerConC,conn);

InnerConC_allind=find(InnerConC);
[xi,yi,zi]=ind2sub(size(InnerConC),InnerConC_allind);
InnerConC_allsub=[xi,yi,zi];
            
label2ConC=false(nx,ny,nz);
for i=1:ll2
    label2ConC(ImageIn_L==label2(i))=1;
end
label2ConC=logical(label2ConC);

% figure(fignum);PlotElectrodes(label1ConC,'b');
% fignum=fignum+1;
% figure(fignum);PlotElectrodes(label2ConC,'b');
% fignum=fignum+1;
% figure(fignum);PlotElectrodes(InnerConC,'b');
% fignum=fignum+1;

T=(T0-th):(-th):Te;
if T(end)~=Te
    T=[T,Te];
end
Tlen=length(T);

ImageOut=ImageIn;
flag2=[];

% h=waitbar(0,'Electrode enter point repair...'); %����������

for i=1:ll2
    
    k=label2(i);
    
    if isempty(intersect(k,flag2))

        flag=0;
        seed0=ImageIn_L;
        seed0(seed0~=k)=0;
        seed0=logical(seed0);
%         figure(fignum);PlotElectrodes(seed0,'b');
%         fignum=fignum+1;
        
        seed0_ind=find(seed0);
        [x0,y0,z0]=ind2sub(size(seed0),seed0_ind);
        seed0_sub=[x0,y0,z0];
        
        D=pdist2(seed0_sub,InnerConC_allsub,'euclidean');
        Dmin=min(min(D));
        if Dmin>=40
            ImageOut(seed0_ind)=0;
            flag2=[flag2,k];
        else

%             seed=logical(seed0.*origin); %ע�⵱origin=electrode1B=regCT.*brainmask���˴�seed0.*origin
            seed=logical(seed0); %ע�⵱origin=electrode1=regCT.*brainmask2ʱ���˴�IConCi�����ٵ��origin

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
                    if size(IS2,1)/size(InnerConC_sub{j},1)>ratio
                        InnerConCj=InnerConC_L;
                        InnerConCj(InnerConC_L~=j)=0;
                        InnerConCj=logical(InnerConCj);

                        [~,~,connect_neighb,~,~]=ConnectTwoConC3(InnerConCj,seed0);
                        ICcnS=logical(InnerConCj+connect_neighb+seed0);
%                         [~,~,~,~,ICcnS_explained,~]=PCACalculation4(ICcnS);

                        InnerConCjPAL=regionprops3(InnerConCj,'PrincipalAxisLength');
                        InnerConCj_PAL=InnerConCjPAL.PrincipalAxisLength;
                        ICcnSPAL=regionprops3(ICcnS,'PrincipalAxisLength');
                        ICcnS_PAL=ICcnSPAL.PrincipalAxisLength;
                        PALdiff=ICcnS_PAL(2)-InnerConCj_PAL(2);

%                         if (ICcnS_explained(1)>90)||(PALdiff<8)
                        if PALdiff<8
                            flag=1;
                            CompStats(i,2)=1;
                            CompStats(i,3)=k;
                            ImageOut=logical(ImageOut+connect_neighb);
%                             figure(fignum);PlotElectrodes(logical(InnerConCj+seed0),'b');hold on;
%                             figure(fignum);PlotElectrodes(logical(connect_neighb),'r');hold on;
%                             fignum=fignum+1;
%                             figure(fignum);PlotElectrodes(ImageOut,'b');
%                             fignum=fignum+1;
                        else
                            label2_renew=setdiff(label2_renew,k);
                        end
                        flag2=[flag2,k];

                        for ne=(i+1):ll2
                            ConCj=ImageIn_L;
                            ConCj(ConCj~=label2(ne))=0;
                            ConCj_ind=find(ConCj);
                            IS3=intersect(seed_ind2,ConCj_ind);
                            if length(IS3)/length(ConCj_ind)>0.8
                                CompStats(ne,2)=1;
                                CompStats(ne,3)=k;
                                CompStats(ne,4)=t;
                                [~,~,connect2_neighb,~,~]=ConnectTwoConC3(InnerConCj,ConCj);
                                ImageOut=logical(ImageOut+connect2_neighb);
                                flag2=[flag2,label2(ne)];

                                [~,~,connect2_neighb,~,~]=ConnectTwoConC3(InnerConCj,ConCj);
                                ICcnS2=logical(InnerConCj+connect2_neighb+ConCj);

                                ICcnS2PAL=regionprops3(ICcnS2,'PrincipalAxisLength');
                                ICcnS2_PAL=ICcnS2PAL.PrincipalAxisLength;
                                PALdiff2=ICcnS2_PAL(2)-InnerConCj_PAL(2);

                                if PALdiff2<8
                                    CompStats(ne,2)=1;
                                    CompStats(ne,3)=k;
                                    CompStats(ne,4)=t;
                                    ImageOut=logical(ImageOut+connect2_neighb);
                                    flag2=[flag2,label2(ne)];
                                end
                            end
                        end

                        break;
                    end
                end

                if flag==1
                    break;
                end

                seed_pn=length(find(seed));
                if (flag==0)&&(seed_pn>6000)
                    label2_renew=setdiff(label2_renew,k);
                    break;
                end

            end
            CompStats(i,4)=t;
            
        end

    else
        i=i+1;
    end
    
%     waitbar(i/ll2,h); %�޸Ľ�����ʾ

end

% close(h);

label2=label2_renew;
ImageOut=logical(ImageOut);

end