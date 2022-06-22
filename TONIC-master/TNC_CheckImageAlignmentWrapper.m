function [] = TNC_CheckImageAlignmentWrapper(inName,outName)

tstart = tic;

% get some file properties and put them in the image structure
[seq_info, fid] = TNC_ReadSeqHeader(inName);
chunkSize = 2000;

numChunks = ceil(seq_info.NumberFrames./chunkSize);

% load the first chunk of frames in the file
seq_image = TNC_ReadSeqImages(seq_info, fid, 1:1:chunkSize );

disp('Loaded initial chunk of frames');
disp(sprintf('Image size is %g [w] by %g [h] pixels',seq_info.Width,seq_info.Height));

[alignmentMetrics] = TNC_ExtractImageAlignment(seq_image,invert,compressed);

structToSave = alignmentMetrics;
disp('Extracted parameters from initial chunk of frames');

for i = 2:numChunks

    % load a set of images
    if i == numChunks
        seq_image = read_seq_images(seq_info, fid,((i-1).*chunkSize)+1:1:seq_info.NumberFrames );        
    else
        seq_image = read_seq_images(seq_info, fid,((i-1).*chunkSize)+1:1:(i.*chunkSize) );
    end
    
    disp(sprintf('Chunk %g of %g ... is loaded',i,numChunks));

    [alignmentMetrics] = TNC_ExtractImageAlignment(seq_image,invert,compressed);

    % scalars
    structToSave.corrShift = [structToSave.corrShift,alignmentMetrics.corrShift];
    structToSave.corrScore = [structToSave.corrScore,alignmentMetrics.corrScore];

    % matrices
    structToSave.profile25 = [structToSave.profile25;alignmentMetrics.profile25];
    structToSave.profile75 = [structToSave.profile75;alignmentMetrics.profile75];
    
    clear seq_image;

end

disp(sprintf('Extracted all alignment data in %g sec...',toc(tstart)));

disp(sprintf('Saving the analyzed output as %s...',outName));
eval(sprintf('save %s structToSave -v7.3',outName));
