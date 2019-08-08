function PatchElectrodes2(IM,facecolor,edgecolor)

%IM���������ά����
%facecolor����ά���ӻ�ʱ�ı�����ɫ
%edgecolor����ά���ӻ�ʱ�ı�Ե��ɫ

if nargin<3
    edgecolor='none';
end

if nargin<2
    facecolor=[0.8 0.8 0.8];
end

[r,c,s]=size(IM);
[X,Y,Z]=meshgrid(1:c,1:r,1:s);
V=IM;
isovalue=0;
p=patch(isosurface(X,Y,Z,V,isovalue));

p.FaceColor=facecolor;
p.EdgeColor=edgecolor;
daspect([1 1 1])
% Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
Azimuth=-15;Elevation=3;view(Azimuth,Elevation);
% camlight 
% lighting gouraud

axis equal;axis([1 c 1 r 1 s]);
xlabel('X');ylabel('Y');zlabel('Z');
set(gcf,'Position',get(0,'ScreenSize'));

end