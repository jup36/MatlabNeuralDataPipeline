function [ ts_new, binEdges ] = binTimeSeriesByMean( ts_org, binSize_org, binSize_new, varargin )

p = parse_input_binTimeSeries( ts_org, binSize_org, binSize_new, varargin );

% sanity check
assert(binSize_new > binSize_org)
if ~isrow(ts_org)
    ts_org = ts_org'; % ensure it's a row vector 
end

ts_new = []; % initialize binned spike count mat
step = 0; 
steps = []; 
stepSize = floor(binSize_new/binSize_org); 

if strcmp(p.Results.align,'subsequent')
    while step + stepSize <= size(ts_org,2) % this prevents stepping out of the range
        currts_new = ts_org(:,step+1:step+stepSize,:); % current portion of the spike count mat
        steps = [steps,step];   % accumulate steps
        step = step + stepSize; % increment step
        ts_new = [ts_new,nanmean(currts_new,2)]; % accumulate binned spike count mat
    end
    binEdges(1,:) = steps(1,:) + 1; % steps
    binEdges(2,:) = binEdges(1,:) + stepSize -1; % where bin ends for each step
else
    halfStepSize = round(stepSize/2); % take equally from previous and subsequent bins
    edges1 = []; edges2 = []; 
    while step + stepSize <= size(ts_org,2) % this prevents stepping out of the range
        tempRange = max(1,step+1-halfStepSize):min(step+halfStepSize,size(ts_org,2));
        currts_new = ts_org(:,tempRange,:); % current portion of the spike count mat
        steps = [steps,step];   % accumulate steps
        edges1 = [edges1,max(1,step+1-halfStepSize)]; 
        edges2 = [edges2,min(step+halfStepSize,size(ts_org,2))]; 
        ts_new = [ts_new,nanmean(currts_new,2)]; % accumulate binned spike count mat
        step = step + stepSize; % increment step
    end 
    binEdges = [ edges1; edges2 ]; 
end

% normalize the spikeCounts for the size of each bin
ts_new = ts_new./repmat((binEdges(2,:)-binEdges(1,:)+1)/stepSize,[size(ts_new,1), 1, size(ts_new,3)]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% nested helper function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_bin1msSpkCountMat( ts_org, binSize, stepSize, vargs ) 
        %parse input, and extract name-value pairs for the main function 'bin1msSpkCountMat.m'
        
        default_align = 'center'; 
        
        p = inputParser; % create parser object
        
        addRequired(p,'ts_org'); % SpkCountMat to be binned
        addRequired(p,'binSize');  % binSize
        addRequired(p,'stepSize'); % stepSize
        
        addParameter(p,'align', default_align) % either 'center' to take preceding and following bins or 'subsequent' to just take following bins

        parse(p, ts_org, binSize, stepSize, vargs{:})
    end
end