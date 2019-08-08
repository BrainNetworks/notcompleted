function Scatter3InRASC2(rasc,VolumeRAS,color,patt,RASWorldLimits)

%rasc������n*3�Ķ�ά���飬��1��2��3�зֱ�Ϊr��a��s���꣬nΪ�������
%VolumeRAS������rasc��Ӧ��RAS����Ⱥ���Ǻ���transform2RAS�����RAS
%color�����ص����ɫ
%patt�����ص����ʽ
%RASWorldLimits��RAS��������ϵ��Χ����һ����Ԫ������[Rmin,Rmax,Amin,Amax,Smin,Smax]

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