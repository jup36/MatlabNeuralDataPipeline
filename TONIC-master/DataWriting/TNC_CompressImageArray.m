function [compressedImageStructure] = TNC_CompressImageArray(modalImage,imageArray,deltaThresh)
% compress the video data in an image array by first calculating the mean
% image and then measuring frame by frame changes in the image. store an
% array of pixels that change and the magnitude of the change

compressedImageStructure.modalImage = modalImage;

for i = 2:size(imageArray,2)

    B = medfilt2(imageArray{1,i}(:,:));

    diffImage   = B - modalImage;
    indices     = find(abs(diffImage) > deltaThresh);
    deltas      = diffImage(indices);
    
    compressedImageStructure.deltas(i).values  = int8(deltas);    
    compressedImageStructure.deltas(i).indices = uint32(indices);
    
end
