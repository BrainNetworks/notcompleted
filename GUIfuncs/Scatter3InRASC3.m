function Scatter3InRASC3(rasc,RASWorldLimits,color,patt)

%rasc������n*3�Ķ�ά���飬��1��2��3�зֱ�Ϊr��a��s���꣬nΪ�������
%RASWorldLimits��RAS��������ϵ��Χ����һ����Ԫ������[Rmin,Rmax,Amin,Amax,Smin,Smax]
%color�����ص����ɫ
%patt�����ص����ʽ

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
% set(gcf,'Position',[350 40 700 650]); %��߽磬�±߽磬��ȣ��߶�

end