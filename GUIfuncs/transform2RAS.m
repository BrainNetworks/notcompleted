function [M,RAS]=transform2RAS(data,Spacing,Thickness,IOP,IPP,coordinate_judge)
%%���ں���������
    %dataΪ�������ά�����ݣ���ά��Ϊm��n��s
    %SpacingΪ��������ĳһɨ����������������ؼ�࣬Ϊ������Ԫ�ص�һά����
    %ThicknessΪ��������������ɨ���Ĳ����
    %IOPΪ�������ݵ�ImageOrientationPatient
    %IPPΪ�������ݵĵ�һ��ɨ��ͼƬ��ImagePositionPatient
    %coordinate_judgeΪ����ϵ�жϣ����������任���ĵ�����
    %��1����dataʹ�õ�����������ϵʱ��coordinate_judge=1��
    %��2����dataʹ�õ�����������ϵʱ��coordinate_judge=0��
%%���ں��������
    %MΪ��XYZ����ϵ��RAS����ϵ�任��4*4��ά��������(r;a;s;1)=M*(x;y;z;1)
    %RASΪ��������ȫ��������RAS����ϵ�е�����Ⱥ��Ϊ4ά����ά��Ϊ4��n��m��s(ע�ⲻ��4��m��n��s)
%% ------����������ͼ��XYZ����ϵ�еĳ�ʼ����------ %%
[row_data,column_data,slice_data]=size(data); 
xyz=zeros(4,column_data,row_data,slice_data); %������ά����
for x=1:column_data %xΪ������,�����������ص�X���꣨X�����᷽�������е�������
    xyz(1,x,:,:)=x; 
end
for y=1:row_data %yΪ������,�����������ص�Y���꣨Y�����᷽�������е�������
    xyz(2,:,y,:)=y; 
end
for z=1:slice_data %zΪ������,�����������ص�Z���꣨Z�����᷽�����ǲ��������
    xyz(3,:,:,z)=z; 
end
xyz(4,:,:,:)=1; %���һά��Ϊ1
%% ------����T0����------ %%
%T0�����ǽ�ͼƬ�����Ͻ����ص��е�����XYZ����ϵ��ԭ��(0,0,0)�������ʼΪ(1,1,1)
T0=[1,0,0,-1;0,1,0,-1;0,0,1,-1;0,0,0,1]; 
%% ------����S����------ %%
%XYZ����ϵ�����صļ��Ĭ��Ϊ1����LPH����ϵ�����صļ����spacing��������Ҫ��ͼ��������ŵ���
%Spacing(2)��X�����ᣨ������������ģ�Spacing(1)��Y�����ᣨ�������������
S=[Spacing(2),0,0,0;0,Spacing(1),0,0;0,0,Thickness,0;0,0,0,1];
%% ------����R����------ %%
%R�����ǽ�ͼƬ��XYZ����ϵ�䵽LPH����ϵ
if coordinate_judge==1 %��������ϵ
    Slice=cross(IOP(1:3),IOP(4:6)); %���ò�����ɨ�跽��ķ�������
end
if coordinate_judge==0 %��������ϵ
    Slice=cross(IOP(4:6),IOP(1:3)); %���ò�����ɨ�跽��ķ�������
end
R=[[IOP(1:3);0],[IOP(4:6);0],[Slice;0],[0;0;0;1]];
%% ------����TIFF����------ %%
%R�����ǽ��任��LPH����ϵ���ͼƬ�Ƶ�ԭ��IPP���涨��λ��
TIFF=[1,0,0,IPP(1,1);0,1,0,IPP(1,2);0,0,1,IPP(1,3);0,0,0,1];
%% ------����lph2ras����------ %%
%lph2ras�����ǽ�LPH����ϵת��ΪRAS����ϵ
lph2ras=[-1,0,0,0;0,-1,0,0;0,0,1,0;0,0,0,1];
%% ------����M����------ %%
M=lph2ras*TIFF*R*S*T0;
%% ------ʵ������任------ %%
RAS=zeros(4,column_data,row_data,slice_data); %������ά����
h=waitbar(0,'Coordinate transforming. Please wait...'); %����������
for z=1:slice_data
    for y=1:row_data
        RAS(:,:,y,z)=M*xyz(:,:,y,z); %xΪ��������yΪ��������zΪ������
    end 
    waitbar(z/slice_data,h); %�޸Ľ�����ʾ
end
close(h);
end