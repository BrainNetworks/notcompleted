function Scatter3InRASC2(rasc,VolumeRAS,color,patt,RASWorldLimits)

%rasc：输入n*3的二维数组，第1、2、3列分别为r、a、s坐标，n为点的总数
%VolumeRAS：输入rasc对应的RAS坐标群，是函数transform2RAS的输出RAS
%color：体素点的颜色
%patt：体素点的样式
%RASWorldLimits：RAS世界坐标系范围，是一个六元素数组[Rmin,Rmax,Amin,Amax,Smin,Smax]

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

R=rasc(:,1);
A=rasc(:,2);
S=rasc(:,3);
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