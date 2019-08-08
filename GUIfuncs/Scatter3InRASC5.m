function Scatter3InRASC5(rasc,RASWorldLimits,Area,patt,LW,MEC,MFC,MEA,MFA)

%rasc������n*3�Ķ�ά���飬��1��2��3�зֱ�Ϊr��a��s���꣬nΪ�������
%RASWorldLimits��RAS��������ϵ��Χ����һ����Ԫ������[Rmin,Rmax,Amin,Amax,Smin,Smax]
%color�����ص����ɫ
%patt�����ص����ʽ

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
% set(gcf,'Position',[350 40 700 650]); %��߽磬�±߽磬��ȣ��߶�

end