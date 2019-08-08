function [IM_L,IM_n,IMConC_Sub]=FindIndexOfAllConnectedComponents(IM,conn)

%IM：输入的三维数组
%conn：连通性

%IM_L：输入IM的连通分量编号
%IM_n：输入IM的连通分量总数
%IMConC_Sub：输入IM的所有连通分量的三维下标索引，是一个n行1列的元胞数组，IMConC_Sub{i}
%              表示第i个连通分量的三维下标索引[rows,columns,slices]（n行3列）

if nargin<2
    conn=26;
end

IM=logical(IM);

[IM_L,IM_n]=bwlabeln(IM,conn);
IMConC_Sub=cell(IM_n,1);

for i=1:IM_n
    IMConC=IM_L;
    IMConC(IMConC~=i)=0;
    index=find(IMConC);
    [x,y,z]=ind2sub(size(IMConC),index);
    sub=[x,y,z];
    IMConC_Sub{i}=sub;
end

end
