function rotateVol=RotateAvg(origVol)

    [~,~,nz5]=size(origVol);
    origVol_ci=origVol;
    for i=1:nz5
        origVol_ci(:,:,i)=fliplr(origVol(:,:,i)); %���ҷ�ת
    end
    origVol_ci2=origVol_ci;
    
    for i=1:nz5
        origVol_ci2(:,:,i)=flipud(origVol_ci(:,:,i)); %���·�ת
    end
    origVol_ci3=permute(origVol_ci2,[3,2,1]); %�ı�����˳��
    
    [~,~,nz6]=size(origVol_ci3);
    rotateVol=origVol_ci3;
    for i=1:nz6
        rotateVol(:,:,i)=flipud(origVol_ci3(:,:,i)); %���·�ת
    end
        
end