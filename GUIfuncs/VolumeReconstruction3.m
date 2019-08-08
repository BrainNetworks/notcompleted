function Volume_recon=VolumeReconstruction3(Volume,Volume_PixelSpacing,Volume_SliceThickness)

%Volume：对直接读取得到的Volume_orig进行索引交换、上下左右前后翻转等处理后，使得解剖面方向一致化的三维体数据
%Volume_PixelSpacing：经过处理后得到的Volume的PixelSpacing，很可能与Volume_orig的PixelSpacing不同
%Volume_SliceThickness：经过处理后得到的Volume的SliceThickness，很可能与Volume_orig的SliceThickness不同

%Volume_recon：对输入Volume进行插值重建得到的结果，Volume_recon的PixelSpacing和SliceThickness都是0.6

[nx,ny,nz]=size(Volume);
oldspacing=[Volume_PixelSpacing;Volume_SliceThickness];
newspacing=[0.6;0.6;0.6];

numrows=round(nx*oldspacing(1)/newspacing(1));
numcols=round(ny*oldspacing(2)/newspacing(2));
numplanes=round(nz*oldspacing(3)/newspacing(3));

Volume_recon=imresize3(Volume,[numrows numcols numplanes],'nearest');

end