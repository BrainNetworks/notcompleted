function [M,RAS]=transform2RAS(data,Spacing,Thickness,IOP,IPP,coordinate_judge)
%%关于函数的输入
    %data为输入的三维体数据，设维度为m×n×s
    %Spacing为输入数据某一扫描层的两个方向的像素间距，为包含两元素的一维向量
    %Thickness为输入数据相邻两扫描层的层间间距
    %IOP为输入数据的ImageOrientationPatient
    %IPP为输入数据的第一张扫描图片的ImagePositionPatient
    %coordinate_judge为坐标系判断（详见《坐标变换》文档）：
    %（1）当data使用的是右手坐标系时，coordinate_judge=1；
    %（2）当data使用的是左手坐标系时，coordinate_judge=0。
%%关于函数的输出
    %M为从XYZ坐标系向RAS坐标系变换的4*4二维矩阵，满足(r;a;s;1)=M*(x;y;z;1)
    %RAS为输入数据全体体素在RAS坐标系中的坐标群，为4维矩阵，维度为4×n×m×s(注意不是4×m×n×s)
%% ------定义输入组图在XYZ坐标系中的初始坐标------ %%
[row_data,column_data,slice_data]=size(data); 
xyz=zeros(4,column_data,row_data,slice_data); %定义四维矩阵
for x=1:column_data %x为列索引,构建所有体素的X坐标（X坐标轴方向上是列的索引）
    xyz(1,x,:,:)=x; 
end
for y=1:row_data %y为行索引,构建所有体素的Y坐标（Y坐标轴方向上是行的索引）
    xyz(2,:,y,:)=y; 
end
for z=1:slice_data %z为层索引,构建所有体素的Z坐标（Z坐标轴方向上是层的索引）
    xyz(3,:,:,z)=z; 
end
xyz(4,:,:,:)=1; %最后一维恒为1
%% ------定义T0矩阵------ %%
%T0作用是将图片的左上角像素的中点移至XYZ坐标系的原点(0,0,0)处，其初始为(1,1,1)
T0=[1,0,0,-1;0,1,0,-1;0,0,1,-1;0,0,0,1]; 
%% ------定义S矩阵------ %%
%XYZ坐标系中体素的间距默认为1，而LPH坐标系中体素的间距是spacing，所以需要对图像进行缩放调整
%Spacing(2)是X坐标轴（列索引）方向的，Spacing(1)是Y坐标轴（行索引）方向的
S=[Spacing(2),0,0,0;0,Spacing(1),0,0;0,0,Thickness,0;0,0,0,1];
%% ------定义R矩阵------ %%
%R作用是将图片从XYZ坐标系变到LPH坐标系
if coordinate_judge==1 %右手坐标系
    Slice=cross(IOP(1:3),IOP(4:6)); %利用叉乘算出扫描方向的方向余弦
end
if coordinate_judge==0 %左手坐标系
    Slice=cross(IOP(4:6),IOP(1:3)); %利用叉乘算出扫描方向的方向余弦
end
R=[[IOP(1:3);0],[IOP(4:6);0],[Slice;0],[0;0;0;1]];
%% ------定义TIFF矩阵------ %%
%R作用是将变换到LPH坐标系后的图片移到原有IPP所规定的位置
TIFF=[1,0,0,IPP(1,1);0,1,0,IPP(1,2);0,0,1,IPP(1,3);0,0,0,1];
%% ------定义lph2ras矩阵------ %%
%lph2ras矩阵是将LPH坐标系转变为RAS坐标系
lph2ras=[-1,0,0,0;0,-1,0,0;0,0,1,0;0,0,0,1];
%% ------定义M矩阵------ %%
M=lph2ras*TIFF*R*S*T0;
%% ------实现坐标变换------ %%
RAS=zeros(4,column_data,row_data,slice_data); %定义四维矩阵
h=waitbar(0,'Coordinate transforming. Please wait...'); %创建进度条
for z=1:slice_data
    for y=1:row_data
        RAS(:,:,y,z)=M*xyz(:,:,y,z); %x为列索引，y为行索引，z为层索引
    end 
    waitbar(z/slice_data,h); %修改进度显示
end
close(h);
end