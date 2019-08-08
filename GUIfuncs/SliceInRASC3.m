function SliceInRASC3(Volume,VolumeRAS,slicedim,slicenum,RASWorldLimits,axes,cmap,EC,FC,FA,FL)

%Volume：输入三维体数据
%VolumeRAS：输入三维体数据Volume（维度为m×n×s）对应的RAS坐标群（维度为4×n×m×s），是函数transform2RAS的输出RAS
%slicedim：所要查看的层的维度，如slicedim=1表示固定Volume的第1维度
%slicenum：所要查看的层的标号，即第几层
%RASWorldLimits：RAS世界坐标系范围，是一个六元素数组[Rmin,Rmax,Amin,Amax,Smin,Smax]
%axes：坐标区句柄
%cmap：colormap，即绘图颜色栏
%EC：EdgeColor
%FC：FaceColor
%FA：FaceAlpha
%FL：FaceLighting

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
    Volume_j=slicenum; %注意Volume与VolumeRAS有两个维度相反
    Volume_k=1:nz1;
    Volume_r=VolumeRAS(1,Volume_i,Volume_j,Volume_k);
    Volume_a=VolumeRAS(2,Volume_i,Volume_j,Volume_k);
    Volume_s=VolumeRAS(3,Volume_i,Volume_j,Volume_k);
    Volume_r=flipud(squeeze(Volume_r)'); %squeeze函数降维后造成索引颠倒
    Volume_a=flipud(squeeze(Volume_a)');
    Volume_s=flipud(squeeze(Volume_s)');
    slice=flipud(squeeze(Volume(slicenum,:,:))'); %squeeze函数降维后造成索引颠倒
end

if slicedim==2
    Volume_i=slicenum; %注意Volume与VolumeRAS有两个维度相反
    Volume_j=1:ny1;
    Volume_k=1:nz1;
    Volume_r=VolumeRAS(1,Volume_i,Volume_j,Volume_k);
    Volume_a=VolumeRAS(2,Volume_i,Volume_j,Volume_k);
    Volume_s=VolumeRAS(3,Volume_i,Volume_j,Volume_k);
    Volume_r=squeeze(Volume_r); %本次squeeze函数降维后没有造成索引颠倒
    Volume_a=squeeze(Volume_a);
    Volume_s=squeeze(Volume_s);
    slice=squeeze(Volume(:,slicenum,:)); %本次squeeze函数降维后没有造成索引颠倒
end

if slicedim==3
    Volume_i=1:nx1;
    Volume_j=1:ny1;
    Volume_k=slicenum;
    Volume_r=VolumeRAS(1,Volume_i,Volume_j,Volume_k);
    Volume_a=VolumeRAS(2,Volume_i,Volume_j,Volume_k);
    Volume_s=VolumeRAS(3,Volume_i,Volume_j,Volume_k);
    Volume_r=squeeze(Volume_r)'; %squeeze函数降维后造成索引颠倒
    Volume_a=squeeze(Volume_a)';
    Volume_s=squeeze(Volume_s)';
    slice=squeeze(Volume(:,:,slicenum)); %本次squeeze函数降维后没有造成索引颠倒
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
% set(gcf,'Position',get(0,'ScreenSize')); %全屏
% set(gcf,'Position',[350 40 700 650]); %左边界，下边界，宽度，高度
% alpha(s,0.5)

end