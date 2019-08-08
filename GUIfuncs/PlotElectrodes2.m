function PlotElectrodes2(IM,XYZWorldLimits,color,patt)

%IM：三维数组的下标索引[rows,columns,slices]，为n行3列的二维数组
%XYZWorldLimits：XYZ世界坐标系范围，是一个六元素数组[Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]
%color：画图时，构成连通分量的点的颜色
%patt：画图时，构成连通分量的点的样式

if nargin<3
    color='b';
end

if nargin<4
    patt='.';
end

x=IM(:,1);y=IM(:,2);z=IM(:,3);
X=y;Y=x;Z=z;
% scatter3(X,Y,Z,patt,color);
plot3(X,Y,Z,patt,'color',color);
axis equal;axis(XYZWorldLimits);
Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
xlabel('X');ylabel('Y');zlabel('Z');
hold on;
set(gcf,'Position',get(0,'ScreenSize'));

end