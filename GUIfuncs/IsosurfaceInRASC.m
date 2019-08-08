function IsosurfaceInRASC(Volume,VolumeRAS)

%Volume：输入三维体数据
%VolumeRAS：输入三维体数据Volume（维度为m×n×s）对应的RAS坐标群（维度为4×n×m×s），是函数transform2RAS的输出RAS

[nx1,ny1,nz1]=size(Volume);

VolumeRAS_rmin=min(min(min(VolumeRAS(1,:,:,:))));
VolumeRAS_rmax=max(max(max(VolumeRAS(1,:,:,:))));
VolumeRAS_amin=min(min(min(VolumeRAS(2,:,:,:))));
VolumeRAS_amax=max(max(max(VolumeRAS(2,:,:,:))));
VolumeRAS_smin=min(min(min(VolumeRAS(3,:,:,:))));
VolumeRAS_smax=max(max(max(VolumeRAS(3,:,:,:))));

Volume_i=1:ny1; %%注意Volume与VolumeRAS有两个维度相反
Volume_j=1:nx1;
Volume_k=1:nz1;
Volume_r=VolumeRAS(1,Volume_i,Volume_j,Volume_k);
Volume_a=VolumeRAS(2,Volume_i,Volume_j,Volume_k);
Volume_s=VolumeRAS(3,Volume_i,Volume_j,Volume_k);
Volume_r=squeeze(Volume_r); %squeeze函数用于降维
Volume_a=squeeze(Volume_a);
Volume_s=squeeze(Volume_s);
Volume_r=permute(Volume_r,[2,1,3]); %squeeze函数降维后造成索引颠倒
Volume_a=permute(Volume_a,[2,1,3]);
Volume_s=permute(Volume_s,[2,1,3]);

Volume_p=patch(isosurface(Volume_r,Volume_a,Volume_s,Volume==0,0));
xlabel('R');ylabel('A');zlabel('S');
axis equal
axis([VolumeRAS_rmin VolumeRAS_rmax VolumeRAS_amin VolumeRAS_amax VolumeRAS_smin VolumeRAS_smax]);
hold on
% Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
set(gcf,'Position',get(0,'ScreenSize'));

end