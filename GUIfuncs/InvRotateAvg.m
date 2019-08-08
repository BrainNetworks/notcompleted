function rotateVol=InvRotateAvg(origVol)

    [~,~,nz5]=size(origVol);
    origVol_ci=origVol;
    for i=1:nz5
        origVol_ci(:,:,i)=flipud(origVol(:,:,i)); %上下翻转
    end
    
    origVol_ci2=permute(origVol_ci,[3,2,1]); %改变索引顺序
    
    origVol_ci3=origVol_ci2;
    [~,~,nz6]=size(origVol_ci2);
    for i=1:nz6
        origVol_ci3(:,:,i)=flipud(origVol_ci2(:,:,i)); %上下翻转
    end
    
    rotateVol=origVol_ci3;
    for i=1:nz6
        rotateVol(:,:,i)=fliplr(origVol_ci3(:,:,i)); %左右翻转
    end
    
end