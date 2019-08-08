function SliceInRASC(Volume,VolumeRAS,slicedim,slicenum,cmap,RASWorldLimits)

%Volume��������ά������
%VolumeRAS��������ά������Volume��ά��Ϊm��n��s����Ӧ��RAS����Ⱥ��ά��Ϊ4��n��m��s�����Ǻ���transform2RAS�����RAS
%slicedim����Ҫ�鿴�Ĳ��ά�ȣ���slicedim=1��ʾ�̶�Volume�ĵ�1ά��
%slicenum����Ҫ�鿴�Ĳ�ı�ţ����ڼ���
%cmap��colormap������ͼ��ɫ��
%RASWorldLimits��RAS��������ϵ��Χ����һ����Ԫ������[Rmin,Rmax,Amin,Amax,Smin,Smax]

[x,y,z]=size(Volume);
nx1=y;
ny1=x;
nz1=z;
% [a,b,c,d]=size(VolumeRAS);
% nx1=b;
% ny1=c;
% nz1=d;

if nargin<6
    VolumeRAS_rmin=min(min(min(VolumeRAS(1,:,:,:))));
    VolumeRAS_rmax=max(max(max(VolumeRAS(1,:,:,:))));
    VolumeRAS_amin=min(min(min(VolumeRAS(2,:,:,:))));
    VolumeRAS_amax=max(max(max(VolumeRAS(2,:,:,:))));
    VolumeRAS_smin=min(min(min(VolumeRAS(3,:,:,:))));
    VolumeRAS_smax=max(max(max(VolumeRAS(3,:,:,:))));
end

if nargin<5
    cmap=gray;
end

if nargin<3
    slicedim=3;
    slicenum=round(nz1/2);
end

if slicedim==1
    Volume_i=1:nx1;
    Volume_j=slicenum; %ע��Volume��VolumeRAS������ά���෴
    Volume_k=1:nz1;
    Volume_r=VolumeRAS(1,Volume_i,Volume_j,Volume_k);
    Volume_a=VolumeRAS(2,Volume_i,Volume_j,Volume_k);
    Volume_s=VolumeRAS(3,Volume_i,Volume_j,Volume_k);
    Volume_r=flipud(squeeze(Volume_r)'); %squeeze������ά����������ߵ�
    Volume_a=flipud(squeeze(Volume_a)');
    Volume_s=flipud(squeeze(Volume_s)');
    slice=flipud(squeeze(Volume(slicenum,:,:))'); %squeeze������ά����������ߵ�
end

if slicedim==2
    Volume_i=slicenum; %ע��Volume��VolumeRAS������ά���෴
    Volume_j=1:ny1;
    Volume_k=1:nz1;
    Volume_r=VolumeRAS(1,Volume_i,Volume_j,Volume_k);
    Volume_a=VolumeRAS(2,Volume_i,Volume_j,Volume_k);
    Volume_s=VolumeRAS(3,Volume_i,Volume_j,Volume_k);
    Volume_r=squeeze(Volume_r); %����squeeze������ά��û����������ߵ�
    Volume_a=squeeze(Volume_a);
    Volume_s=squeeze(Volume_s);
    slice=squeeze(Volume(:,slicenum,:)); %����squeeze������ά��û����������ߵ�
end

if slicedim==3
    Volume_i=1:nx1;
    Volume_j=1:ny1;
    Volume_k=slicenum;
    Volume_r=VolumeRAS(1,Volume_i,Volume_j,Volume_k);
    Volume_a=VolumeRAS(2,Volume_i,Volume_j,Volume_k);
    Volume_s=VolumeRAS(3,Volume_i,Volume_j,Volume_k);
    Volume_r=squeeze(Volume_r)'; %squeeze������ά����������ߵ�
    Volume_a=squeeze(Volume_a)';
    Volume_s=squeeze(Volume_s)';
    slice=squeeze(Volume(:,:,slicenum)); %����squeeze������ά��û����������ߵ�
end

surf(Volume_r,Volume_a,Volume_s,slice);
xlabel('R');ylabel('A');zlabel('S');
axis equal
if nargin==6
    axis(RASWorldLimits);
else
    axis([VolumeRAS_rmin VolumeRAS_rmax VolumeRAS_amin VolumeRAS_amax VolumeRAS_smin VolumeRAS_smax]);
end
hold on
shading interp %ͨ����ÿ�����������ж���ɫͼ���������ɫֵ���в�ֵ���ı�����������е���ɫ
colormap(cmap)
% Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
set(gcf,'Position',get(0,'ScreenSize')); %ȫ��
% set(gcf,'Position',[350 40 700 650]); %��߽磬�±߽磬��ȣ��߶�

end