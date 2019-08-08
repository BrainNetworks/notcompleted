function PlotPCAForAllConnectedComponents(electrodes,conn,PCApatt,EConCpatt,fignum)

%electrodes：输入的三维数组
%conn：连通性
%PCApatt：画图时，构成PCA直线的点的样式
%EConCpatt：画图时，构成连通分量的点的样式
%fignum：所画的第一张图片的序号

if nargin<2
    conn=26;
end

if nargin<3
    PCApatt='-';
    EConCpatt='.';
end

if nargin<4
    EConCpatt='.';
end

if nargin<5
    fignum=1;
end

electrodes=logical(electrodes);
[L,n]=bwlabeln(electrodes,conn);

for i=1:n
    econc=L;
    econc(econc~=i)=0;
    
    figure(fignum+floor((i-1)/4));
    subplot(2,2,(mod(i,4)==0)*4+mod(i,4));
    PlotElectrodes(econc,'b',EConCpatt);hold on;
    ElectrodePCAPlot2(econc,'r',PCApatt);hold on;
    title(['连通分量',num2str(i)]);
end

end