function rasc=Scatter3InRASC(Volume,VolumeRAS,color,patt,RASWorldLimits)
%%
%Volume��������ά������
%VolumeRAS��������ά������Volume��Ӧ��RAS����Ⱥ���Ǻ���transform2RAS�����RAS
%color�����ص����ɫ
%patt�����ص����ʽ
%RASWorldLimits��RAS��������ϵ��Χ����һ����Ԫ������[Rmin,Rmax,Amin,Amax,Smin,Smax]

%rasc��������ά������Volume�����ػҶ�ֵ������ras���꣬Ϊn*3�Ķ�ά���飬
%      ��1~3������Ϊr,a,s���꣬nΪVolume�����ػҶ�ֵ����������Ŀ
%%
if nargin<5
    VolumeRAS_rmin=min(min(min(VolumeRAS(1,:,:,:))));
    VolumeRAS_rmax=max(max(max(VolumeRAS(1,:,:,:))));
    VolumeRAS_amin=min(min(min(VolumeRAS(2,:,:,:))));
    VolumeRAS_amax=max(max(max(VolumeRAS(2,:,:,:))));
    VolumeRAS_smin=min(min(min(VolumeRAS(3,:,:,:))));
    VolumeRAS_smax=max(max(max(VolumeRAS(3,:,:,:))));
end

if nargin<4
    patt='.';
end

if nargin<3
    color='b';
end

Volume_ind=find(Volume);
[x,y,z]=ind2sub(size(Volume),Volume_ind);
X=y;Y=x;Z=z; %������ת��Ϊ���꣬��Ϊ����VolumeRAS������transform2RAS������������꼯���������±���������������Ҫ����ת���������Ӧ
L=length(X);
r=zeros(L,1);
a=zeros(L,1);
s=zeros(L,1);
for i=1:L
    r(i)=VolumeRAS(1,X(i),Y(i),Z(i));
    a(i)=VolumeRAS(2,X(i),Y(i),Z(i));
    s(i)=VolumeRAS(3,X(i),Y(i),Z(i));    
end
rasc=[r,a,s];
%% ���ǲ��Եģ���Ϊ�õ���r,a,s�Ѿ�������ֵ�ˣ�����������ֵ
% R=a;A=r;S=s;
% scatter3(R,A,S,patt,color);
% % plot3(R,A,S,patt,'color',color);
%% ����������ȷ��
R=r;A=a;S=s;
% scatter3(R,A,S,patt,color);
plot3(R,A,S,patt,'color',color);
xlabel('R');ylabel('A');zlabel('S');
axis equal
if nargin==5
    axis(RASWorldLimits);
else
    axis([VolumeRAS_rmin VolumeRAS_rmax VolumeRAS_amin VolumeRAS_amax VolumeRAS_smin VolumeRAS_smax]);
end
hold on
% Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
set(gcf,'Position',get(0,'ScreenSize'));
end