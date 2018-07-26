function [] = TNC_ROIpicker(imageArray,roi,imgIndArray)
% FUNCTION DETAILS: This function goes through a single channel of filtered recording data and looks for threshold crossings. A second stage then tests these threshold crossings according to a template matching heuristic to try to classify significant events.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% OUTPUT STRUCTURE ELEMENTS

numROIs = length(roi);

% display image
for i = 1:length(imgIndArray)

    % display the current image
    figure(10); clf;
subplot(211);
% imagesc(imageArray{1,i}(:,:),[0 255]); colormap(gray); 
subplot(212);
B = medfilt2(imageArray{1,i}(:,:),[2,2]);
imagesc(B); colormap(gray); 

    hold on;

    % display first roi
    for j = 1:numROIs        
        line([roi(j).colLo,roi(j).colHi],[roi(j).rowHi,roi(j).rowHi]);
        line([roi(j).colLo,roi(j).colHi],[roi(j).rowLo,roi(j).rowLo]);
        line([roi(j).colLo,roi(j).colLo],[roi(j).rowLo,roi(j).rowHi]);
        line([roi(j).colHi,roi(j).colHi],[roi(j).rowLo,roi(j).rowHi]);
    end
    
    drawnow;
    
end