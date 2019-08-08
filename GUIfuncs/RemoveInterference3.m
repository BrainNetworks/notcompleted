function [VolumeIn_PN,VolumeOut]=RemoveInterference3(VolumeIn,Mask,conn)

%VolumeIn：输入的三维体数据
%Mask：掩模，也是三维数组，维度与VolumeIn相同，如brainmask,brainmask3等
%conn：连通性

%VolumeIn_PN：为n*4的二维数组，第1,2,3,4列分别为输入三维体数据VolumeIn的连通分量编号，像素总数，颅内像素数和颅骨边界上像素数
%VolumeOut：将所完全位于颅骨边界上的干扰，颅内部分像素数少于30个的干扰以及完全位于颅内的干扰去除后得到的输出三维体数据

if nargin<3
    conn=26;
end

VolumeIn=logical(VolumeIn);
Mask=logical(Mask);

[VolumeIn_L,VolumeIn_n]=bwlabeln(VolumeIn,conn);
VIL_intracranial=VolumeIn_L.*Mask; %intracranial：颅内的
VIL_cranial=VolumeIn_L.*(~Mask); %cranial：颅骨上的

VolumeInConC_ID=(1:VolumeIn_n)';
VolumeInConC_PN=zeros(VolumeIn_n,1);
VILintcConC_PN=zeros(VolumeIn_n,1); %VILintc：VIL_intracranial
VILcConC_PN=zeros(VolumeIn_n,1); %VILc：VIL_cranial
for i=1:VolumeIn_n
    VolumeInConC_PN(i)=sum(sum(sum(VolumeIn_L==i))); %ConCPN：number of pixels
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