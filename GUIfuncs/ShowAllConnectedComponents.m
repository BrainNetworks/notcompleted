function ShowAllConnectedComponents(electrodes,conn,fignum,L,n)

%electrodes���������ά����
%conn����ͨ��
%fignum�������ĵ�һ��ͼƬ�����
%L��electrodes��ԭ��ͨ������ţ������Ǹ���/У׼����±��
%n��electrodes��ԭ��ͨ���������������Ǹ���/У׼����±�ŵ����ֵ

if nargin<2
    conn=26;
end

if nargin<3
    fignum=1;
end

if nargin<5
    [L,n]=bwlabeln(electrodes,conn);
end

electrodes=logical(electrodes);
for i=1:n
    electrode_i=L;
    electrode_i(electrode_i~=i)=0;
    electrode_ni=L;
    electrode_ni(electrode_ni==i)=0;
    figure(fignum+floor((i-1)/4));
    subplot(2,2,(mod(i,4)==0)*4+mod(i,4));
%     PatchElectrodes(electrode_ni,'b');hold on;
%     PatchElectrodes(electrode_i,'r');hold on;
    PlotElectrodes(electrode_ni,'b');hold on;
    PlotElectrodes(electrode_i,'r');hold on;
    title(['��ͨ����',num2str(i)]);
    set(gcf,'Position',get(0,'ScreenSize'));
end

end