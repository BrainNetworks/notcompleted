function [E1_L,E2_L,E1_n,E2_n,E1_S,E2_newL,E2_recL]=InnerConCsDetection(E1,E2,conn)

%E1：经过regCT阈值处理等操作得到的电极初步提取结果
%E2：某一三维数组（通常是E1）与颅骨边界的交集部分，称为电极入点
%conn：连通性

%E1_L：输入E1的初始连通分量编号，是三维数组
%E2_L：输入E2的初始连通分量编号，是三维数组
%E1_n：输入E1的连通分量总数
%E2_n：输入E2的连通分量总数
%E1_S：输入E1的所有连通分量的统计信息
%E1_S(:,1)为连通分量编号
%E1_S(:,2)为连通分量是否完全在颅骨边界上或者与颅骨交界有交集，值1表示是，值0表示否；
%         值为1称为颅骨边界连通分量，值为0则称为颅内连通分量
%E2_newL：输入E2的新的连通分量编号，是三维数组
%E2_recL：输入E2含有的E1连通分量的编号数值，E2_recL的行号是E2原来的连通分量编号，
%         因此E2_recL可反映出原来E2和E1的编号对应关系，rec即record，是一维列向量

if nargin<3
    conn=26;
end

E1=logical(E1);
E2=logical(E2);

[E1_L,E1_n]=bwlabeln(E1,conn);
[E2_L,E2_n]=bwlabeln(E2,conn);
E2_newL=zeros(size(E2_L));
E2_recL=zeros(E2_n,1);
                
for i=1:E2_n
    
    E2ConCi=E2_L;
    E2ConCi(E2_L~=i)=0;
    E2ConCi=logical(E2ConCi);
    
    IS=E2ConCi.*E1_L;
    if sum(sum(sum(IS~=0)))~=0
        IS_ind=find(IS);
        LV=IS(IS_ind(1));
        E2_recL(i)=LV;
        E2_newL(E2ConCi)=LV;
    end
    
end

E1_S=zeros(E1_n,2);
E1_S(:,1)=(1:1:E1_n)';
E2_checkL=unique(E2_recL); %去除重复元素并升序排序
E1_S(E2_checkL,2)=1;

end