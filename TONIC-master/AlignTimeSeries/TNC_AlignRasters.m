function [response, psthZ] = TNC_AlignRasters(delta,spkStamps,stableTrials,alignStamps,window,rasterFlag,boxcar)
% FUNCTION DETAILS: Using a set of timestamps provided in ALIGNSTAMPS creates a matrix of DATA vectors that span -WINDOW(1,1) to +WINDOW(1,2) around each timestamp.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

alignStamps = round(alignStamps);

if size(alignStamps,1) > size(alignStamps,2)
    alignStamps = alignStamps';
end

count = 1;
numStamps = size(alignStamps,2); 

% overall mean and std rate
if alignStamps(numStamps) > numel(delta)
    stdRate     = std(delta(alignStamps(1):numel(delta)));     
    meanRate     = mean(delta(alignStamps(1):numel(delta)));     
else
    stdRate     = std(delta(alignStamps(1):alignStamps(numStamps))); 
    meanRate     = mean(delta(alignStamps(1):alignStamps(numStamps)));     
end

response.image.aligned = zeros(numStamps, window(1,1)+window(1,2)+1);

for i = 1:numStamps

    currStamp = round(alignStamps(1,i)); % current behavioral timestamp

    % get the list of timestamps within the current window
    if rasterFlag
        validSpks = find(spkStamps>currStamp-window(1,1) & spkStamps<currStamp+window(1,2));
        response.raster.trial(i).ts = spkStamps(validSpks) - currStamp;        
    end
    
    % case 'image' (use all indices by padding)
    if currStamp + window(1,2) < size(delta,2)
        if currStamp-window(1,1) > 0
            response.image.aligned(count,:)   = delta(1,currStamp-window(1,1):currStamp+window(1,2));
            response.image.indicesUsed(count) = i;
            response.image.indices(count)     = currStamp;
        else
%             response.image.aligned(count,:)   = [zeros(1,abs(currStamp-window(1,1)-1)),delta(1,1:currStamp+window(1,2))];
%             response.image.indicesUsed(count) = i;
            count = count-1;
        end
    else
        if currStamp-window(1,1) > 0
            tmp = delta(1,currStamp-window(1,1):length(delta));
            response.image.aligned(count,1:length(tmp))   = tmp;
            response.image.indicesUsed(count) = i;
        end
    end

    count = count+1;
        
end

% write out the response variables...
count = count-1;

% size(response.image.aligned,2)
    
% boxcar is just delta function averaged across trials more sparsely
% sampled across time (e.g. every 6 time points)
if boxcar == 1
    centers = (0:6:size(response.image.aligned,2)-6) + 3;
    boxcarTMP = sum(response.image.aligned,1)./numStamps;
    for i = 1:length(centers)
        response.image.boxcar(i) = sum(boxcarTMP(centers(i)-2:centers(i)+3));
    end   
end

% get mean and std FR across trials
if stableTrials==-1 % in this case, include all trials without excluding any trials in the beginning
    response.image.alignedtime    = -window(1,1):1:window(1,2);
    response.image.psthAVG        = mean(response.image.aligned,1);% meanRate;
    response.image.psthSEM        = std( response.image.aligned,0,1) ./ sqrt(size(response.image.aligned,1)-1);
else                % in this case, exclude initial trials from mean and std calculation
    response.image.alignedtime    = -window(1,1):1:window(1,2);
    response.image.psthAVG        = mean(response.image.aligned(stableTrials:count,:),1);% - meanRate;
    response.image.psthSEM        = std( response.image.aligned(stableTrials:count,:),0,1) ./ sqrt(size(response.image.aligned(stableTrials:count,:),1)-1);
end

% use the entire left window (half a second) as an estimate of the mean
avgRate = mean(response.image.psthAVG(1:round(window(1,1)./2)));
if stdRate>0.0001
    response.image.psthZ         = (response.image.psthAVG-avgRate) ./ stdRate;  
    response.image.psthZe         = response.image.psthSEM ./ stdRate;  
    response.noZ                  = 0;
else % do not calculate z-scores when std is way too low
    disp(['STD=' num2str( mean(response.image.psthAVG(1,1:window(1,1)),2) ) ' | cannot z-score']);
    response.image.psthZ          = response.image.psthAVG;
    response.image.psthZe         = response.image.psthAVG;
    response.noZ                  = 1;
end

psthZ = response.image.psthZ;
response.meanRate = meanRate; % overall mean firing rate per 1 ms needs to be multiplied by 1000 to convert it to Hz
response.stdRate = stdRate;   % overall std firing rate per 1 ms needs to be multiplied by 1000 to convert it to Hz
