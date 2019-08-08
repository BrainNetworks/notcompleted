function [E1_L,E2_L,E1_n,E2_n,E1_S,E2_newL,E2_recL]=InnerConCsDetection(E1,E2,conn)

%E1������regCT��ֵ����Ȳ����õ��ĵ缫������ȡ���
%E2��ĳһ��ά���飨ͨ����E1����­�Ǳ߽�Ľ������֣���Ϊ�缫���
%conn����ͨ��

%E1_L������E1�ĳ�ʼ��ͨ������ţ�����ά����
%E2_L������E2�ĳ�ʼ��ͨ������ţ�����ά����
%E1_n������E1����ͨ��������
%E2_n������E2����ͨ��������
%E1_S������E1��������ͨ������ͳ����Ϣ
%E1_S(:,1)Ϊ��ͨ�������
%E1_S(:,2)Ϊ��ͨ�����Ƿ���ȫ��­�Ǳ߽��ϻ�����­�ǽ����н�����ֵ1��ʾ�ǣ�ֵ0��ʾ��
%         ֵΪ1��Ϊ­�Ǳ߽���ͨ������ֵΪ0���Ϊ­����ͨ����
%E2_newL������E2���µ���ͨ������ţ�����ά����
%E2_recL������E2���е�E1��ͨ�����ı����ֵ��E2_recL���к���E2ԭ������ͨ������ţ�
%         ���E2_recL�ɷ�ӳ��ԭ��E2��E1�ı�Ŷ�Ӧ��ϵ��rec��record����һά������

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
E2_checkL=unique(E2_recL); %ȥ���ظ�Ԫ�ز���������
E1_S(E2_checkL,2)=1;

end