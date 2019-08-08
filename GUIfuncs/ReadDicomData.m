function [Volume,Volume_IPP,Volume_IOP,Volume_PixelSpacing,Volume_SliceThickness]=ReadDicomData(path)

%path：要读取的数据所在文件夹的路径

%Volume：对直接读取得到的Volume_orig进行索引交换、上下左右前后翻转等处理后，使得解剖面方向一致化的三维体数据
%Volume_IPP：Volume的每一层对应的IPP，如Volume(:,:,i)的IPP为Volume_IPP(i,:)，但是由于Volume经过了处理，因此每层中LPH坐标等于IPP值的像素很可能不再位于左上角
%Volume_IOP：准确来讲应该是Volume_orig的IOP，因为经过处理后得到的Volume解剖面方向一致，IOP实际上应该都满足round(IOP)=[0;1;0;0;0;-1]
%Volume_PixelSpacing：经过处理后得到的Volume的PixelSpacing，很可能与Volume_orig的PixelSpacing不同
%Volume_SliceThickness：经过处理后得到的Volume的SliceThickness，很可能与Volume_orig的SliceThickness不同

Volume_picstr=dir([path,'*.dcm']); %dir列出当前文件夹中的文件及其信息
Volume_info1=dicominfo([path,Volume_picstr(1).name]); %读出第1层的图片信息

Volume_rows=Volume_info1.Rows;
Volume_columns=Volume_info1.Columns;
Volume_slices=size(Volume_picstr,1); %Volume扫描层数/帧数

Volume_orig=zeros(Volume_rows,Volume_columns,Volume_slices);
Volume_IPP=zeros(Volume_slices,3); %读出IPP，每一层IPP值不相同，IPP值是该层图片左上角像素的LPH坐标，单位是mm
Volume_IOP=Volume_info1.ImageOrientationPatient; %读出IOP，所有片层的IOP值都是相同的
Volume_PixelSpacing=Volume_info1.PixelSpacing; %读出同一片层内相邻像素的间距，所有片层的PixelSpacing值都是相同的，单位是mm
Volume_SliceThickness=Volume_info1.SliceThickness; %读出相邻两片层之间的间距，所有片层的Thickness值都是相同的，单位是mm

for i=1:Volume_slices
    Volume_infoi=dicominfo([path,Volume_picstr(i).name]); %读出第i层的图片信息
    Volume_orig(:,:,i)=dicomread([path,Volume_picstr(i).name]); %读出第i层的Volume图片
    Volume_IPP(i,:)=Volume_infoi.ImagePositionPatient; %读出IPP，单位是mm，相邻两个片层的IOP值之间在扫描方向上相差一个层厚
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
        Volume=rot90(Volume,2); %等同于Volume=fliplr(flipud(Volume));
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
        Volume=rot90(Volume,2); %等同于Volume=fliplr(flipud(Volume));
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
        Volume=rot90(Volume,2); %等同于Volume=fliplr(flipud(Volume));
    end
    
    Volume=permute(Volume,[3,2,1]);
    Volume=flipud(Volume);
    Volume=permute(Volume,[1,3,2]);
    Volume=Volume(:,:,size(Volume,3):(-1):1);

%--------------------------------------------------------------------------
%注：本小节是A节
%     Volume_PixelSpacing=[Volume_SliceThickness;Volume_PixelSpacing(2)];
%     Volume_SliceThickness=Volume_PixelSpacing(1);
%     Volume_PixelSpacing=[Volume_PixelSpacing(1);Volume_SliceThickness];
%     Volume_SliceThickness=Volume_PixelSpacing(2);
%--------------------------------------------------------------------------
%注：本小节是B节，B节等同于A节
    Volume_PixelSpacing=[Volume_SliceThickness;Volume_PixelSpacing(1)];
    Volume_SliceThickness=Volume_PixelSpacing(2);
%--------------------------------------------------------------------------

end

end