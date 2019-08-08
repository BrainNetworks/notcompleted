function PlotElectrodes(IM,color,patt)

%IM���������ά��������ά������±�����[rows,columns,slices]��Ϊn��3�еĶ�ά����
%color����ͼʱ��������ͨ�����ĵ����ɫ
%patt����ͼʱ��������ͨ�����ĵ����ʽ

if nargin<2
    color='b';
end

if nargin<3
    patt='.';
end

if length(size(IM))==3 %Ҫ���ʱ������IMΪ��ά����
    index=find(IM);
    [x,y,z]=ind2sub(size(IM),index);
    %x���е�������y���е�������z�ǲ������
    %���Ƕ��ں���plot3����scatter�����루X,Y,Z��
    %��X,Y,Z���ֱ�����ȴ�ǵ���X,Y,Z�������ϵ�����ֵ
    %�����ǵ������ֵ��x,y,z��
    %����ΪX����������������y����Y����������������x
    %�Ӷ�X=y��Y=x
    X=y;Y=x;Z=z;
%     scatter3(X,Y,Z,patt,color);
    plot3(X,Y,Z,patt,'color',color);
else %Ҫ���ʱ������IMΪĳ��ά������±�����[rows,columns,slices]��Ϊn��3�еĶ�ά����
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