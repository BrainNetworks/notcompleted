function VolumeOut=RemoveInterference2(VolumeIn,Mask,conn)

%VolumeIn：输入的三维体数据
%Mask：掩模，也是三维数组，维度与VolumeIn相同，如brainmask,brainmask3等
%conn：连通性

%VolumeOut：将所有像素均位于颅骨上的干扰（即完全位于颅骨边界上的连通分量）去除后得到的输出三维体数据

if nargin<3
    conn=26;
end

VolumeIn=logical(VolumeIn);
Mask=logical(Mask);

[VolumeIn_L,VolumeIn_n]=bwlabeln(VolumeIn,conn);
VIL_intcracranial=VolumeIn_L.*Mask; %intcracranial：颅内的

VolumeInConC_PN=zeros(VolumeIn_n,1);
VILintcConC_PN=zeros(VolumeIn_n,1); %VILintc：VIL_intcracranial
for i=1:VolumeIn_n
    VolumeInConC_PN(i)=sum(sum(sum(VolumeIn_L==i))); %ConCPN：number of pixels
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