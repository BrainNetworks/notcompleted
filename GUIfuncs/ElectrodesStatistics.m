function [E1_S,E2_newL,E2_recL,E2_checkL,E2_judgeL]=ElectrodesStatistics(E1,E2,conn)

%E1������regCT��ֵ�����Ȳ����õ��ĵ缫������ȡ���
%E2��ĳһ��ά���飨ͨ����E1����­�Ǳ߽�Ľ������֣���Ϊ�缫���
%������E1����ͨ�������Ϊ��׼��������E2����ͨ�����������±��
%conn����ͨ��

%E1_S������E1��������ͨ������ͳ����Ϣ
%E1_S(:,1)Ϊ��ͨ�������
%E1_S(:,2)Ϊ��ͨ������������Ŀ
%E1_S(:,3)Ϊ��ͨ�����Ƿ���ȫ��­�Ǳ߽��ϻ�����­�ǽ����н�����ֵ1��ʾ�ǣ�ֵ0��ʾ��
%E1_S(:,4:6)Ϊ��ͨ������һ���ɷֵ���������
%E2_newL������E2���µ���ͨ������ţ�����ά����
%E2_recL������E2���е�E1��ͨ�����ı����ֵ��E2_recL���к���E2ԭ������ͨ������ţ����E2_recL�ɷ�ӳ��ԭ��E2��E1�ı�Ŷ�Ӧ��ϵ��rec��record����һά������
%E2_checkL���Ƕ�E2_recLȥ���ظ�Ԫ�ز���������Ľ��
%E2_judgeL��E2_judgeL=1�������������һ����һ�����ϵ���ͨ����������������­�Ǹ�����������

E1=logical(E1);
E2=logical(E2);
[E1_L,E1_n]=bwlabeln(E1,conn);
E1_S=zeros(E1_n,6);
E1_S(:,1)=1:E1_n;

for i=1:E1_n
    E1_S(i,2)=sum(sum(sum(E1_L==i)));
end

[~,~,~,~,E2_newL,E2_recL,~,~]=CountElectrodes(E1,E2,conn);
E2_checkL=unique(E2_recL); %ȥ���ظ�Ԫ�ز���������
E1_S(E2_checkL,3)=1;
E2_judgeL=(length(E2_recL)~=length(E2_checkL));
%E2_judgeL=1�������������һ����һ�����ϵ���ͨ����������������­�Ǹ�����������

E1_pca=PCACalculation(E1,conn);
E1_S(:,4:6)=E1_pca;

end