function ShowAllConnectedComponents(electrodes,conn,fignum,L,n)

%electrodes：输入的三维数组
%conn：连通性
%fignum：所画的第一张图片的序号
%L：electrodes的原连通分量编号，或者是更改/校准后的新编号
%n：electrodes的原连通分量总数，或者是更改/校准后的新编号的最大值

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
    title(['连通分量',num2str(i)]);
    set(gcf,'Position',get(0,'ScreenSize'));
end

end