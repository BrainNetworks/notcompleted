function [movingRegisteredVolume,optimizer_small,geomtform_small]=Registration(fixedHeader,fixedVolume,...
    movingHeader,movingVolume,ratio,InitialRadius,MaximumIterations,f)

fixedVolume_small=fixedVolume(1:ratio:end,1:ratio:end,1:ratio:end);
movingVolume_small=movingVolume(1:ratio:end,1:ratio:end,1:ratio:end);

[optimizer_small,metric_small]=imregconfig('multimodal');

Rfixed_small=imref3d(size(fixedVolume_small),fixedHeader.PixelSpacing(2).*ratio,...
    fixedHeader.PixelSpacing(1).*ratio,fixedHeader.SliceThickness.*ratio);
Rmoving_small=imref3d(size(movingVolume_small),movingHeader.PixelSpacing(2).*ratio,...
    movingHeader.PixelSpacing(1).*ratio,movingHeader.SliceThickness.*ratio);

Rfixed=imref3d(size(fixedVolume),fixedHeader.PixelSpacing(2),...
    fixedHeader.PixelSpacing(1),fixedHeader.SliceThickness);
Rmoving=imref3d(size(movingVolume),movingHeader.PixelSpacing(2),...
    movingHeader.PixelSpacing(1),movingHeader.SliceThickness);

optimizer_small.InitialRadius=InitialRadius;
optimizer_small.MaximumIterations=MaximumIterations; 
geomtform_small=imregtform(movingVolume_small,Rmoving_small,...
    fixedVolume_small,Rfixed_small,'rigid',optimizer_small,metric_small);
waitbar(0.5,f);

% movingRegisteredVolume_small=imwarp(movingVolume_small,...
%     Rmoving_small,geomtform2_small,'bicubic','OutputView',Rfixed_small);

movingRegisteredVolume=imwarp(movingVolume,...
    Rmoving,geomtform_small,'bicubic','OutputView',Rfixed);
waitbar(0.7,f);

end