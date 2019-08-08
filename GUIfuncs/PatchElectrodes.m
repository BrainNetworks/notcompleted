function PatchElectrodes(IM,color)

%IM：输入的三维数组
%color：三维可视化时的表面颜色

if nargin<2
    color='black';
end

p=patch(isosurface(IM==1,0));
p.EdgeColor=color;
Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
[nx,ny,nz]=size(IM);
axis equal;axis([1 ny 1 nx 1 nz]);
% axis([150 400 100 270 1 160]);
xlabel('X');ylabel('Y');zlabel('Z');
set(gcf,'Position',get(0,'ScreenSize'));

end