function [modalImage] = TNC_ImgStackMode(imageArray,range,modeFlag)

dimsOfImage = size(imageArray{1,1}(:,:));
totalPixels = dimsOfImage(1,1).*dimsOfImage(1,2);

for h = 1: dimsOfImage(1,1)
    for i = 1: dimsOfImage(1,2)
    
        for j=range(1,1):range(1,2)
            B = medfilt2(imageArray{1,i}(:,:));
            zStack(j) = B(h,i);
        end
        
        if modeFlag==1
            pixelMode = mode(zStack);            
        else
            pixelMode = mean(zStack);
        end
        
        modalImage(h,i) = pixelMode;

    end
end
