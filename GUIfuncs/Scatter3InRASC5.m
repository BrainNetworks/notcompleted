function Scatter3InRASC5(rasc,RASWorldLimits,Area,patt,LW,MEC,MFC,MEA,MFA)

%rasc：输入n*3的二维数组，第1、2、3列分别为r、a、s坐标，n为点的总数
%RASWorldLimits：RAS世界坐标系范围，是一个六元素数组[Rmin,Rmax,Amin,Amax,Smin,Smax]
%color：体素点的颜色
%patt：体素点的样式

if nargin<4
    Area=25;
    patt='.';
    LW=0.5;
    MEC='r';
    MFC='r';
    MEA=1;
    MFA=1;
end

R=rasc(:,1);
A=rasc(:,2);
S=rasc(:,3);

scatter3(R,A,S,Area,patt,'Linewidth',LW,'MarkerEdgeColor',MEC,'MarkerFaceColor',MFC,'MarkerEdgeAlpha',MEA,'MarkerFaceAlpha',MFA);hold on;
xlabel('R');ylabel('A');zlabel('S');
axis equal
axis(RASWorldLimits);
hold on
% Azimuth=-15;Elevation=13;view(Azimuth,Elevation);
% set(gcf,'Position',get(0,'ScreenSize'));
% set(gcf,'Position',[350 40 700 650]); %左边界，下边界，宽度，高度

end