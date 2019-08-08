function [IM_L,IM_n,IMConC_Sub]=FindIndexOfAllConnectedComponents(IM,conn)

%IM���������ά����
%conn����ͨ��

%IM_L������IM����ͨ�������
%IM_n������IM����ͨ��������
%IMConC_Sub������IM��������ͨ��������ά�±���������һ��n��1�е�Ԫ�����飬IMConC_Sub{i}
%              ��ʾ��i����ͨ��������ά�±�����[rows,columns,slices]��n��3�У�

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
