function Volume_recon=VolumeReconstruction(Volume,Volume_PixelSpacing,Volume_SliceThickness)

%Volume：对直接读取得到的Volume_orig进行索引交换、上下左右前后翻转等处理后，使得解剖面方向一致化的三维体数据
%Volume_PixelSpacing：经过处理后得到的Volume的PixelSpacing，很可能与Volume_orig的PixelSpacing不同
%Volume_SliceThickness：经过处理后得到的Volume的SliceThickness，很可能与Volume_orig的SliceThickness不同

%Volume_recon：根据PixelSpacing和SliceThickness对Volume进行插值重建得到的结果

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