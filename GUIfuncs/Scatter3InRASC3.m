function Scatter3InRASC3(rasc,RASWorldLimits,color,patt)

%rasc：输入n*3的二维数组，第1、2、3列分别为r、a、s坐标，n为点的总数
%RASWorldLimits：RAS世界坐标系范围，是一个六元素数组[Rmin,Rmax,Amin,Amax,Smin,Smax]
%color：体素点的颜色
%patt：体素点的样式

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
axis(RASWorldLimits);
hold on
% Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
set(gcf,'Position',get(0,'ScreenSize'));
% set(gcf,'Position',[350 40 700 650]); %左边界，下边界，宽度，高度

end