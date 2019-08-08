function rotateVol=InvRotateAvg(origVol)

    [~,~,nz5]=size(origVol);
    origVol_ci=origVol;
    for i=1:nz5
        origVol_ci(:,:,i)=flipud(origVol(:,:,i)); %���·�ת
    end
    
    origVol_ci2=permute(origVol_ci,[3,2,1]); %�ı�����˳��
    
    origVol_ci3=origVol_ci2;
    [~,~,nz6]=size(origVol_ci2);
    for i=1:nz6
        origVol_ci3(:,:,i)=flipud(origVol_ci2(:,:,i)); %���·�ת
    end
    
    rotateVol=origVol_ci3;
    for i=1:nz6
        rotateVol(:,:,i)=fliplr(origVol_ci3(:,:,i)); %���ҷ�ת
    end
    
end