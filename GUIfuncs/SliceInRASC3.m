function SliceInRASC3(Volume,VolumeRAS,slicedim,slicenum,RASWorldLimits,axes,cmap,EC,FC,FA,FL)

%Volume��������ά������
%VolumeRAS��������ά������Volume��ά��Ϊm��n��s����Ӧ��RAS����Ⱥ��ά��Ϊ4��n��m��s�����Ǻ���transform2RAS�����RAS
%slicedim����Ҫ�鿴�Ĳ��ά�ȣ���slicedim=1��ʾ�̶�Volume�ĵ�1ά��
%slicenum����Ҫ�鿴�Ĳ�ı�ţ����ڼ���
%RASWorldLimits��RAS��������ϵ��Χ����һ����Ԫ������[Rmin,Rmax,Amin,Amax,Smin,Smax]
%axes�����������
%cmap��colormap������ͼ��ɫ��
%EC��EdgeColor
%FC��FaceColor
%FA��FaceAlpha
%FL��FaceLighting

[x,y,z]=size(Volume);
nx1=y;
ny1=x;
nz1=z;
% [a,b,c,d]=size(VolumeRAS);
% nx1=b;
% ny1=c;
% nz1=d;

if nargin<7
    cmap=gray;
    EC='none';
    FC='interp';
    FA=1;
    FL='flat';
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

surf(Volume_r,Volume_a,Volume_s,slice,'EdgeColor',EC,'FaceColor',FC,'FaceAlpha',FA,'FaceLighting',FL);
% s=surf(Volume_r,Volume_a,Volume_s,slice,'EdgeColor',EC,'FaceColor',FC,'FaceAlpha',FA,'FaceLighting',FL);
xlabel('R');ylabel('A');zlabel('S');
axis equal
if nargin>=5
    axis(RASWorldLimits);
else
    axis([VolumeRAS_rmin VolumeRAS_rmax VolumeRAS_amin VolumeRAS_amax VolumeRAS_smin VolumeRAS_smax]);
end
hold on
if length(cmap)>1
    colormap(axes,cmap);
end
% Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
% set(gcf,'Position',get(0,'ScreenSize')); %ȫ��
% set(gcf,'Position',[350 40 700 650]); %��߽磬�±߽磬��ȣ��߶�
% alpha(s,0.5)

end