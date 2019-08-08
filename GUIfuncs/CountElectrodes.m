function [E1ConC_Sub,E2ConC_Sub,E1_L,E2_L,E2_newL,E2_recL,E1_n,E2_n]=CountElectrodes(E1,E2,conn)

%E1：经过regCT阈值处理等操作得到的电极初步提取结果
%E2：某一三维数组（通常是E1）与颅骨边界的交集部分，称为电极入点
%conn：连通性

%E1ConC_Sub：输入E1的所有连通分量的三维下标索引
%E2ConC_Sub：输入E2的所有连通分量的三维下标索引
%E1_L：输入E1的初始连通分量编号，是三维数组
%E2_L：输入E2的初始连通分量编号，是三维数组
%E2_newL：输入E2的新的连通分量编号，是三维数组
%E2_recL：输入E2含有的E1连通分量的编号数值，E2_recL的行号是E2原来的连通分量编号，因此E2_recL可反映出原来E2和E1的编号对应关系，rec即record，是一维列向量
%E1_n：输入E1的连通分量总数
%E2_n：输入E2的连通分量总数

if nargin<3
    conn=26;
end

E1=logical(E1);
E2=logical(E2);
[E1_L,E1_n,E1ConC_Sub]=FindIndexOfAllConnectedComponents(E1,conn);
[E2_L,E2_n,E2ConC_Sub]=FindIndexOfAllConnectedComponents(E2,conn);
E2_newL=zeros(size(E2_L));
E2_recL=zeros(E2_n,1);
                
for i=1:E2_n
    E2ConCIi=E2ConC_Sub{i};
    for j=1:E1_n
        E1ConCIj=E1ConC_Sub{j};
%         E2ConCIi_xmin=min(E2ConCIi(:,1));
%         E2ConCIi_xmax=max(E2ConCIi(:,1));
%         E2ConCIi_xcenter=floor((E2ConCIi_xmax-E2ConCIi_xmin)/2);
%         E2ConCIi_range=(E2ConCIi_xcenter-3):1:(E2ConCIi_xcenter+3);
%         similarity=intersect(E2ConCIi_range,E1ConCIj(:,1)');
%         if isempty(~similarity)
            IS=intersect(E2ConCIi,E1ConCIj,'rows');
            E2ConCIi_nr=size(E2ConCIi,1);
            ISnr=size(IS,1);
            ratio=ISnr/E2ConCIi_nr;
            if ratio>0.1
                E2_newL_ind=sub2ind(size(E2_newL),E2ConCIi(:,1),E2ConCIi(:,2),E2ConCIi(:,3));
                E2_newL(E2_newL_ind)=j;
                E2_recL(i)=j;
            end
%         end
    end
end

end