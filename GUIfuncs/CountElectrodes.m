function [E1ConC_Sub,E2ConC_Sub,E1_L,E2_L,E2_newL,E2_recL,E1_n,E2_n]=CountElectrodes(E1,E2,conn)

%E1������regCT��ֵ����Ȳ����õ��ĵ缫������ȡ���
%E2��ĳһ��ά���飨ͨ����E1����­�Ǳ߽�Ľ������֣���Ϊ�缫���
%conn����ͨ��

%E1ConC_Sub������E1��������ͨ��������ά�±�����
%E2ConC_Sub������E2��������ͨ��������ά�±�����
%E1_L������E1�ĳ�ʼ��ͨ������ţ�����ά����
%E2_L������E2�ĳ�ʼ��ͨ������ţ�����ά����
%E2_newL������E2���µ���ͨ������ţ�����ά����
%E2_recL������E2���е�E1��ͨ�����ı����ֵ��E2_recL���к���E2ԭ������ͨ������ţ����E2_recL�ɷ�ӳ��ԭ��E2��E1�ı�Ŷ�Ӧ��ϵ��rec��record����һά������
%E1_n������E1����ͨ��������
%E2_n������E2����ͨ��������

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