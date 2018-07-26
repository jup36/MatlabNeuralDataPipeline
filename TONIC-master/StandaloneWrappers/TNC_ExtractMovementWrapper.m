function [] = TNC_ExtractMovementWrapper(inName,outName,roi,invert,compressed)
% FUNCTION DETAILS: Wrapper to allow extraction of all movement data from an entire seq file through a series of chunks that are assembled into one data structure.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

tstart = tic;

% get some file properties and put them in the image structure
[seq_info, fid] = TNC_ReadSeqHeader(inName);
chunkSize = 2000;

numChunks = ceil(seq_info.NumberFrames./chunkSize);
% numChunks = 10;

% load the first chunk of frames in the file
seq_image = TNC_ReadSeqImages( seq_info, fid, 1:1:chunkSize );
disp('Loaded initial chunk of frames');
disp(sprintf('Image size is %g [w] by %g [h] pixels',seq_info.Width,seq_info.Height));

figure(1); clf;% show ROI
imagesc(seq_image{1}(:,:),[0 255]); colormap(gray); hold on;
drawnow;

numROIs = length(roi);
for j = 1:numROIs        
    line([roi(j).colLo,roi(j).colHi],[roi(j).rowHi,roi(j).rowHi]);
    line([roi(j).colLo,roi(j).colHi],[roi(j).rowLo,roi(j).rowLo]);
    line([roi(j).colLo,roi(j).colLo],[roi(j).rowLo,roi(j).rowHi]);
    line([roi(j).colHi,roi(j).colHi],[roi(j).rowLo,roi(j).rowHi]);
end
drawnow;
 
[outputStructure] = TNC_ExtractMovement(seq_image,roi,invert,compressed);
structToSave = outputStructure;
disp('Extracted parameters from initial chunk of frames');

for i = 2:numChunks

    % load a set of images
    if i == numChunks
        seq_image = TNC_ReadSeqImages(seq_info, fid,((i-1).*chunkSize)+1:1:seq_info.NumberFrames );        
    else
        seq_image = TNC_ReadSeqImages(seq_info, fid,((i-1).*chunkSize)+1:1:(i.*chunkSize) );
    end
    
    disp(sprintf('Chunk %g of %g ... is loaded',i,numChunks));

    [outputStructure] = TNC_ExtractMovement(seq_image,roi,invert,compressed);

    for indROI = 1:length(roi)

        % scalars
        structToSave.roi(indROI).varR = [structToSave.roi(indROI).varR,outputStructure.roi(indROI).varR];
        structToSave.roi(indROI).varC = [structToSave.roi(indROI).varC,outputStructure.roi(indROI).varC];
        structToSave.roi(indROI).maxVR = [structToSave.roi(indROI).maxVR,outputStructure.roi(indROI).maxVR];
        structToSave.roi(indROI).maxLR = [structToSave.roi(indROI).maxLR,outputStructure.roi(indROI).maxLR];
        structToSave.roi(indROI).maxVC = [structToSave.roi(indROI).maxVC,outputStructure.roi(indROI).maxVC];
        structToSave.roi(indROI).maxLC = [structToSave.roi(indROI).maxLC,outputStructure.roi(indROI).maxLC];

        % matrices
        structToSave.roi(indROI).xcorrR = [structToSave.roi(indROI).xcorrR;outputStructure.roi(indROI).xcorrR];
        structToSave.roi(indROI).xcorrC = [structToSave.roi(indROI).xcorrC;outputStructure.roi(indROI).xcorrC];
        structToSave.roi(indROI).lotCumSumR = [structToSave.roi(indROI).lotCumSumR;outputStructure.roi(indROI).lotCumSumR];
        structToSave.roi(indROI).lotCumSumC = [structToSave.roi(indROI).lotCumSumC;outputStructure.roi(indROI).lotCumSumC];
    
    end
    
    clear seq_image;

end

disp(sprintf('Extracted all motion data in %g sec...',toc(tstart)));

disp(sprintf('Saving the analyzed output as %s...',outName));
eval(sprintf('save %s structToSave -v7.3',outName));
