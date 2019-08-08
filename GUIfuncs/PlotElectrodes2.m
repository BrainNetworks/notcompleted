function PlotElectrodes2(IM,XYZWorldLimits,color,patt)

%IM����ά������±�����[rows,columns,slices]��Ϊn��3�еĶ�ά����
%XYZWorldLimits��XYZ��������ϵ��Χ����һ����Ԫ������[Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]
%color����ͼʱ��������ͨ�����ĵ����ɫ
%patt����ͼʱ��������ͨ�����ĵ����ʽ

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