function PlotElectrodes(IM,color,patt)

%IM：输入的三维数组或该三维数组的下标索引[rows,columns,slices]，为n行3列的二维数组
%color：画图时，构成连通分量的点的颜色
%patt：画图时，构成连通分量的点的样式

if nargin<2
    color='b';
end

if nargin<3
    patt='.';
end

if length(size(IM))==3 %要求此时的输入IM为三维数组
    index=find(IM);
    [x,y,z]=ind2sub(size(IM),index);
    %x是行的索引，y是列的索引，z是层的索引
    %但是对于函数plot3或者scatter的输入（X,Y,Z）
    %（X,Y,Z）分别代表的却是点在X,Y,Z坐标轴上的坐标值
    %而不是点的索引值（x,y,z）
    %又因为X坐标轴上是列索引y，而Y坐标轴上是行索引x
    %从而X=y，Y=x
    X=y;Y=x;Z=z;
%     scatter3(X,Y,Z,patt,color);
    plot3(X,Y,Z,patt,'color',color);
else %要求此时的输入IM为某三维数组的下标索引[rows,columns,slices]，为n行3列的二维数组
    x=IM(:,1);y=IM(:,2);z=IM(:,3);
    X=y;Y=x;Z=z;
%     scatter3(X,Y,Z,patt,color);
    plot3(X,Y,Z,patt,'color',color);
end

if length(size(IM))==3
    Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
    [nx,ny,nz]=size(IM);
    axis equal;axis([1 ny 1 nx 1 nz]);
    xlabel('X');ylabel('Y');zlabel('Z');
    hold on;
    set(gcf,'Position',get(0,'ScreenSize'));
else
%     Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
    Azimuth=41;Elevation=5;view(Azimuth,Elevation);
    xlabel('X');ylabel('Y');zlabel('Z');
    hold on;
    set(gcf,'Position',get(0,'ScreenSize'));
end

end