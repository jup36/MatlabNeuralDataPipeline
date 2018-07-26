function [] = TNC_VideoCompressWrapper(inName,outName,chunkSize,difference)

if ischar(chunkSize)
    chunkSize=str2num(chunkSize);
end

if ischar(difference)
    difference=str2num(difference);
end

tstart = tic;

% get some file properties and put them in the image structure
[seq_info, fid] = TNC_ReadSeqHeader(inName);

numChunks = ceil(seq_info.NumberFrames./chunkSize);

% load the first chunk of frames in the file
seq_image = TNC_ReadSeqImages(seq_info, fid, 1:1:chunkSize);
disp('Loaded initial chunk of frames');
disp(sprintf('Image size is %g [w] by %g [h] pixels',seq_info.Width,seq_info.Height));

% calculate the modal image
[modalImage] = TNC_ImgStackMode(seq_image,[1,50],1);
disp('Calculated modal image for background');
compressedImageStructure.modalImage = modalImage;

[chunkStructure] = TNC_CompressImageArray(modalImage,seq_image,difference); % look for changes of greater than bitDepth bits
compressedImageStructure.deltas = chunkStructure.deltas;
disp('Compressed first chunk of frames');
clear chunkStructure;
        
for i = 2:numChunks

    % load a set of images
    if i == numChunks
        seq_image = TNC_ReadSeqImages(seq_info, fid,((i-1).*chunkSize)+1:1:seq_info.NumberFrames );        
    else
        seq_image = TNC_ReadSeqImages(seq_info, fid,((i-1).*chunkSize)+1:1:(i.*chunkSize) );
    end
    disp(sprintf('Chunk %g of %g ... is loaded',i,numChunks));
    
    % create the compression sequence arrays
    [chunkStructure] = TNC_CompressImageArray(modalImage,seq_image,difference); %look for changes of greater than bitDepth bits
    compressedImageStructure.deltas = [compressedImageStructure.deltas,chunkStructure.deltas];
    disp(sprintf('Chunk %g of %g ... is compressed',i,numChunks));

    clear chunkStructure; % can this stop the bad memory leak?
    disp(sprintf('Chunk %g of %g ... swap data is cleared',i,numChunks));

end

disp(sprintf('Compressed file in %g sec...',toc(tstart)));

disp(sprintf('Saving the compressed file as %s...',outName));
eval(sprintf('save %s compressedImageStructure -v7.3',outName));

disp(sprintf('TNC_VideoCompressWrapper completed in %g sec.',toc(tstart)));