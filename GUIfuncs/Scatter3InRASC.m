function rasc=Scatter3InRASC(Volume,VolumeRAS,color,patt,RASWorldLimits)
%%
%Volume：输入三维体数据
%VolumeRAS：输入三维体数据Volume对应的RAS坐标群，是函数transform2RAS的输出RAS
%color：体素点的颜色
%patt：体素点的样式
%RASWorldLimits：RAS世界坐标系范围，是一个六元素数组[Rmin,Rmax,Amin,Amax,Smin,Smax]

%rasc：输入三维体数据Volume中体素灰度值非零点的ras坐标，为n*3的二维数组，
%      第1~3列依次为r,a,s坐标，n为Volume中体素灰度值非零点的总数目
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
X=y;Y=x;Z=z; %将索引转换为坐标，因为输入VolumeRAS（函数transform2RAS的输出）是坐标集，而不是下标索引集，所以需要进行转换以与其对应
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
%% 这是不对的，因为得到的r,a,s已经是坐标值了，而不是索引值
% R=a;A=r;S=s;
% scatter3(R,A,S,patt,color);
% % plot3(R,A,S,patt,'color',color);
%% 这样才是正确的
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