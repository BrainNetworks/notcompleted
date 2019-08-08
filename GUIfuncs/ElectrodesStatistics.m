function [E1_S,E2_newL,E2_recL,E2_checkL,E2_judgeL]=ElectrodesStatistics(E1,E2,conn)

%E1：经过regCT阈值处理等操作得到的电极初步提取结果
%E2：某一三维数组（通常是E1）与颅骨边界的交集部分，称为电极入点
%以输入E1的连通分量编号为标准，对输入E2的连通分量进行重新编号
%conn：连通性

%E1_S：输入E1的所有连通分量的统计信息
%E1_S(:,1)为连通分量编号
%E1_S(:,2)为连通分量的像素数目
%E1_S(:,3)为连通分量是否完全在颅骨边界上或者与颅骨交界有交集，值1表示是，值0表示否
%E1_S(:,4:6)为连通分量第一主成分的特征向量
%E2_newL：输入E2的新的连通分量编号，是三维数组
%E2_recL：输入E2含有的E1连通分量的编号数值，E2_recL的行号是E2原来的连通分量编号，因此E2_recL可反映出原来E2和E1的编号对应关系，rec即record，是一维列向量
%E2_checkL：是对E2_recL去除重复元素并升序排序的结果
%E2_judgeL：E2_judgeL=1则表明，存在着一个或一个以上的连通分量，它或它们与颅骨各有两个交点

E1=logical(E1);
E2=logical(E2);
[E1_L,E1_n]=bwlabeln(E1,conn);
E1_S=zeros(E1_n,6);
E1_S(:,1)=1:E1_n;

for i=1:E1_n
    E1_S(i,2)=sum(sum(sum(E1_L==i)));
end

[~,~,~,~,E2_newL,E2_recL,~,~]=CountElectrodes(E1,E2,conn);
E2_checkL=unique(E2_recL); %去除重复元素并升序排序
E1_S(E2_checkL,3)=1;
E2_judgeL=(length(E2_recL)~=length(E2_checkL));
%E2_judgeL=1则表明，存在着一个或一个以上的连通分量，它或它们与颅骨各有两个交点

E1_pca=PCACalculation(E1,conn);
E1_S(:,4:6)=E1_pca;

end