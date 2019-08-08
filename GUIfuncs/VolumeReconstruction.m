function Volume_recon=VolumeReconstruction(Volume,Volume_PixelSpacing,Volume_SliceThickness)

%Volume����ֱ�Ӷ�ȡ�õ���Volume_orig����������������������ǰ��ת�ȴ����ʹ�ý����淽��һ�»�����ά������
%Volume_PixelSpacing�����������õ���Volume��PixelSpacing���ܿ�����Volume_orig��PixelSpacing��ͬ
%Volume_SliceThickness�����������õ���Volume��SliceThickness���ܿ�����Volume_orig��SliceThickness��ͬ

%Volume_recon������PixelSpacing��SliceThickness��Volume���в�ֵ�ؽ��õ��Ľ��

[nx,ny,nz]=size(Volume);
spacing=[Volume_PixelSpacing;Volume_SliceThickness];
if (spacing(1)==spacing(2))&&(spacing(2)==spacing(3))
    Volume_recon=CT;
else
    if spacing(1)==spacing(2)
        numrows=nx;
        numcols=ny;
        numplanes=round(nz*spacing(3)/spacing(1));
    else
        if spacing(1)==spacing(3)
            numrows=nx;
            numcols=round(ny*spacing(2)/spacing(1));
            numplanes=nz;
        end
        if spacing(2)==spacing(3)
            numrows=round(nx*spacing(1)/spacing(2));
            numcols=ny;
            numplanes=nz;
        end
    end
    Volume_recon=imresize3(Volume,[numrows numcols numplanes],'nearest');
end

end