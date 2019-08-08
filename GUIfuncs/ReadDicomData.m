function [Volume,Volume_IPP,Volume_IOP,Volume_PixelSpacing,Volume_SliceThickness]=ReadDicomData(path)

%path��Ҫ��ȡ�����������ļ��е�·��

%Volume����ֱ�Ӷ�ȡ�õ���Volume_orig����������������������ǰ��ת�ȴ����ʹ�ý����淽��һ�»�����ά������
%Volume_IPP��Volume��ÿһ���Ӧ��IPP����Volume(:,:,i)��IPPΪVolume_IPP(i,:)����������Volume�����˴������ÿ����LPH�������IPPֵ�����غܿ��ܲ���λ�����Ͻ�
%Volume_IOP��׼ȷ����Ӧ����Volume_orig��IOP����Ϊ���������õ���Volume�����淽��һ�£�IOPʵ����Ӧ�ö�����round(IOP)=[0;1;0;0;0;-1]
%Volume_PixelSpacing�����������õ���Volume��PixelSpacing���ܿ�����Volume_orig��PixelSpacing��ͬ
%Volume_SliceThickness�����������õ���Volume��SliceThickness���ܿ�����Volume_orig��SliceThickness��ͬ

Volume_picstr=dir([path,'*.dcm']); %dir�г���ǰ�ļ����е��ļ�������Ϣ
Volume_info1=dicominfo([path,Volume_picstr(1).name]); %������1���ͼƬ��Ϣ

Volume_rows=Volume_info1.Rows;
Volume_columns=Volume_info1.Columns;
Volume_slices=size(Volume_picstr,1); %Volumeɨ�����/֡��

Volume_orig=zeros(Volume_rows,Volume_columns,Volume_slices);
Volume_IPP=zeros(Volume_slices,3); %����IPP��ÿһ��IPPֵ����ͬ��IPPֵ�Ǹò�ͼƬ���Ͻ����ص�LPH���꣬��λ��mm
Volume_IOP=Volume_info1.ImageOrientationPatient; %����IOP������Ƭ���IOPֵ������ͬ��
Volume_PixelSpacing=Volume_info1.PixelSpacing; %����ͬһƬ�����������صļ�࣬����Ƭ���PixelSpacingֵ������ͬ�ģ���λ��mm
Volume_SliceThickness=Volume_info1.SliceThickness; %����������Ƭ��֮��ļ�࣬����Ƭ���Thicknessֵ������ͬ�ģ���λ��mm

for i=1:Volume_slices
    Volume_infoi=dicominfo([path,Volume_picstr(i).name]); %������i���ͼƬ��Ϣ
    Volume_orig(:,:,i)=dicomread([path,Volume_picstr(i).name]); %������i���VolumeͼƬ
    Volume_IPP(i,:)=Volume_infoi.ImagePositionPatient; %����IPP����λ��mm����������Ƭ���IOPֵ֮����ɨ�跽�������һ�����
end

IOPround=round(Volume_IOP);

if isequal(abs(IOPround),[0;1;0;0;0;1])
    
    Volume_L=Volume_IPP(:,1);
    [~,Volume_LL]=sort(Volume_L,'descend');
    Volume_IPP=Volume_IPP(Volume_LL,:);
    Volume=Volume_orig(:,:,Volume_LL);
    
    if isequal(IOPround,[0;1;0;0;0;1])
        Volume=flipud(Volume);
    end

    if isequal(IOPround,[0;-1;0;0;0;-1])
        Volume=fliplr(Volume);
    end
    
    if isequal(IOPround,[0;-1;0;0;0;1])
        Volume=rot90(Volume,2); %��ͬ��Volume=fliplr(flipud(Volume));
    end
    
end

if isequal(abs(IOPround),[1;0;0;0;0;1])
    
    Volume_P=Volume_IPP(:,2);
    [~,Volume_PL]=sort(Volume_P,'ascend');
    Volume_IPP=Volume_IPP(Volume_PL,:);
    Volume=Volume_orig(:,:,Volume_PL);
    
    if isequal(IOPround,[1;0;0;0;0;1])
        Volume=flipud(Volume);
    end

    if isequal(IOPround,[-1;0;0;0;0;-1])
        Volume=fliplr(Volume);
    end
    
    if isequal(IOPround,[-1;0;0;0;0;1])
        Volume=rot90(Volume,2); %��ͬ��Volume=fliplr(flipud(Volume));
    end
    
    Volume=permute(Volume,[1,3,2]);
    Volume=Volume(:,:,size(Volume,3):(-1):1);
    
    Volume_PixelSpacing=[Volume_PixelSpacing(1);Volume_SliceThickness];
    Volume_SliceThickness=Volume_PixelSpacing(2);

end

if isequal(abs(IOPround),[1;0;0;0;1;0])
    
    Volume_H=Volume_IPP(:,3);
    [~,Volume_HL]=sort(Volume_H,'ascend');
    Volume_IPP=Volume_IPP(Volume_HL,:);
    Volume=Volume_orig(:,:,Volume_HL);
    
    if isequal(IOPround,[1;0;0;0;-1;0])
        Volume=flipud(Volume);
    end

    if isequal(IOPround,[-1;0;0;0;1;0])
        Volume=fliplr(Volume);
    end
    
    if isequal(IOPround,[-1;0;0;0;-1;0])
        Volume=rot90(Volume,2); %��ͬ��Volume=fliplr(flipud(Volume));
    end
    
    Volume=permute(Volume,[3,2,1]);
    Volume=flipud(Volume);
    Volume=permute(Volume,[1,3,2]);
    Volume=Volume(:,:,size(Volume,3):(-1):1);

%--------------------------------------------------------------------------
%ע����С����A��
%     Volume_PixelSpacing=[Volume_SliceThickness;Volume_PixelSpacing(2)];
%     Volume_SliceThickness=Volume_PixelSpacing(1);
%     Volume_PixelSpacing=[Volume_PixelSpacing(1);Volume_SliceThickness];
%     Volume_SliceThickness=Volume_PixelSpacing(2);
%--------------------------------------------------------------------------
%ע����С����B�ڣ�B�ڵ�ͬ��A��
    Volume_PixelSpacing=[Volume_SliceThickness;Volume_PixelSpacing(1)];
    Volume_SliceThickness=Volume_PixelSpacing(2);
%--------------------------------------------------------------------------

end

end