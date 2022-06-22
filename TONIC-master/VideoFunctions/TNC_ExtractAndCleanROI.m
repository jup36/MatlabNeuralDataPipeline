function [subImage] = TNC_ExtractAndCleanROI(imageArray,roi,dispFlag)
% compress the video data in an image array by first calculating the mean
% image and then measuring frame by frame changes in the image. store an
% array of pixels that change and the magnitude of the change

% Format of the 'roi' structure fields
% roi(1).rowLo
% roi(1).rowHi
% roi(1).colLo
% roi(1).colHi

figure(1); colormap(gray);

for i = 1:size(imageArray,2)

    B = medfilt2(imageArray{1,i}(:,:),[2 2]);
    subImage(i).frame = B(roi(1).rowLo:roi(1).rowHi,roi(1).colLo:roi(1).colHi);
    if dispFlag
        figure(1);
        imagesc(subImage(i).frame,[0 255]);
    end
end
