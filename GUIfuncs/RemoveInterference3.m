function [VolumeIn_PN,VolumeOut]=RemoveInterference3(VolumeIn,Mask,conn)

%VolumeIn���������ά������
%Mask����ģ��Ҳ����ά���飬ά����VolumeIn��ͬ����brainmask,brainmask3��
%conn����ͨ��

%VolumeIn_PN��Ϊn*4�Ķ�ά���飬��1,2,3,4�зֱ�Ϊ������ά������VolumeIn����ͨ������ţ�����������­����������­�Ǳ߽���������
%VolumeOut��������ȫλ��­�Ǳ߽��ϵĸ��ţ�­�ڲ�������������30���ĸ����Լ���ȫλ��­�ڵĸ���ȥ����õ��������ά������

if nargin<3
    conn=26;
end

VolumeIn=logical(VolumeIn);
Mask=logical(Mask);

[VolumeIn_L,VolumeIn_n]=bwlabeln(VolumeIn,conn);
VIL_intracranial=VolumeIn_L.*Mask; %intracranial��­�ڵ�
VIL_cranial=VolumeIn_L.*(~Mask); %cranial��­���ϵ�

VolumeInConC_ID=(1:VolumeIn_n)';
VolumeInConC_PN=zeros(VolumeIn_n,1);
VILintcConC_PN=zeros(VolumeIn_n,1); %VILintc��VIL_intracranial
VILcConC_PN=zeros(VolumeIn_n,1); %VILc��VIL_cranial
for i=1:VolumeIn_n
    VolumeInConC_PN(i)=sum(sum(sum(VolumeIn_L==i))); %ConCPN��number of pixels
    VILintcConC_PN(i)=sum(sum(sum(VIL_intracranial==i)));
    VILcConC_PN(i)=sum(sum(sum(VIL_cranial==i)));
end

VolumeIn_PN=[VolumeInConC_ID,VolumeInConC_PN,VILintcConC_PN,VILcConC_PN];

VolumeOut=VolumeIn_L;
for i=1:VolumeIn_n
    if VILintcConC_PN(i)==0
        VolumeOut(VolumeOut==i)=0;
    end
    if VILcConC_PN(i)==0
        VolumeOut(VolumeOut==i)=0;
    end
    if (VILintcConC_PN(i)~=0)&&(VILcConC_PN(i)~=0)
        if VILintcConC_PN(i)<30
            VolumeOut(VolumeOut==i)=0;
        end
    end
end

VolumeOut=logical(VolumeOut);

end