function VolumeOut=RemoveInterference2(VolumeIn,Mask,conn)

%VolumeIn���������ά������
%Mask����ģ��Ҳ����ά���飬ά����VolumeIn��ͬ����brainmask,brainmask3��
%conn����ͨ��

%VolumeOut�����������ؾ�λ��­���ϵĸ��ţ�����ȫλ��­�Ǳ߽��ϵ���ͨ������ȥ����õ��������ά������

if nargin<3
    conn=26;
end

VolumeIn=logical(VolumeIn);
Mask=logical(Mask);

[VolumeIn_L,VolumeIn_n]=bwlabeln(VolumeIn,conn);
VIL_intcracranial=VolumeIn_L.*Mask; %intcracranial��­�ڵ�

VolumeInConC_PN=zeros(VolumeIn_n,1);
VILintcConC_PN=zeros(VolumeIn_n,1); %VILintc��VIL_intcracranial
for i=1:VolumeIn_n
    VolumeInConC_PN(i)=sum(sum(sum(VolumeIn_L==i))); %ConCPN��number of pixels
    VILintcConC_PN(i)=sum(sum(sum(VIL_intcracranial==i)));
end

VolumeOut=VolumeIn_L;
for i=1:VolumeIn_n
%     if VILintcConC_PN(i)<2
    if (VILintcConC_PN(i)<2)&&(VolumeInConC_PN(i)<20)
        VolumeOut(VolumeOut==i)=0;
    end
end

VolumeOut=logical(VolumeOut);

end