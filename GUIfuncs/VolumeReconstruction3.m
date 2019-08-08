function Volume_recon=VolumeReconstruction3(Volume,Volume_PixelSpacing,Volume_SliceThickness)

%Volume����ֱ�Ӷ�ȡ�õ���Volume_orig����������������������ǰ��ת�ȴ����ʹ�ý����淽��һ�»�����ά������
%Volume_PixelSpacing�����������õ���Volume��PixelSpacing���ܿ�����Volume_orig��PixelSpacing��ͬ
%Volume_SliceThickness�����������õ���Volume��SliceThickness���ܿ�����Volume_orig��SliceThickness��ͬ

%Volume_recon��������Volume���в�ֵ�ؽ��õ��Ľ����Volume_recon��PixelSpacing��SliceThickness����0.6

[nx,ny,nz]=size(Volume);
oldspacing=[Volume_PixelSpacing;Volume_SliceThickness];
newspacing=[0.6;0.6;0.6];

numrows=round(nx*oldspacing(1)/newspacing(1));
numcols=round(ny*oldspacing(2)/newspacing(2));
numplanes=round(nz*oldspacing(3)/newspacing(3));

Volume_recon=imresize3(Volume,[numrows numcols numplanes],'nearest');

end