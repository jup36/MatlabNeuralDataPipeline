function [ binMat, binEdges ] = bin1msSpkCountMat( SpkCountMat1ms, binSize, stepSize, varargin )
% bin1msSpkCountMat takes a trial-by-1msBin spike count matrix and generates 
%  a binned spikeCountMat (binMat) using the binSize and stepSize specified by the user.   
%  binEdges contains: the range of each column of the binMat

p = parse_input_bin1msSpkCountMat( SpkCountMat1ms, binSize, stepSize, varargin );
%p = parse_input_bin1msSpkCountMat( SpkCountMat1ms, binSize, stepSize, {'align','center'} );

step = 0;    % initialize step
binMat = []; % initialize binned spike count mat
steps  = []; % initialize accumulated steps

if strcmp(p.Results.align,'subsequent')
    while step + binSize <= size(SpkCountMat1ms,2) % this prevents stepping out of the range
        currBinMat = SpkCountMat1ms(:,step+1:step+binSize,:); % current portion of the spike count mat
        steps = [steps,step];   % accumulate steps
        step = step + stepSize; % increment step
        binMat = [binMat,sum(currBinMat,2)]; % accumulate binned spike count mat
    end
    binEdges(1,:) = steps(1,:) + 1; % steps
    binEdges(2,:) = binEdges(1,:) + binSize -1; % where bin ends for each step
else
    halfBinSize = round(binSize/2); % take equally from previous and subsequent bins
    edges1 = []; edges2 = []; 
    while step + binSize <= size(SpkCountMat1ms,2) % this prevents stepping out of the range
        tempRange = max(1,step+1-halfBinSize):min(step+halfBinSize,size(SpkCountMat1ms,2));
        currBinMat = SpkCountMat1ms(:,tempRange,:); % current portion of the spike count mat
        steps = [steps,step];   % accumulate steps
        edges1 = [edges1,max(1,step+1-halfBinSize)]; 
        edges2 = [edges2,min(step+halfBinSize,size(SpkCountMat1ms,2))]; 
        binMat = [binMat,sum(currBinMat,2)]; % accumulate binned spike count mat
        step = step + stepSize; % increment step
    end 
    binEdges = [ edges1; edges2 ]; 
end

% normalize the spikeCounts for the size of each bin
binMat = binMat./repmat((binEdges(2,:)-binEdges(1,:)+1)/binSize,[size(binMat,1), 1, size(binMat,3)]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% nested helper function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_bin1msSpkCountMat( SpkCountMat1ms, binSize, stepSize, vargs ) 
        %parse input, and extract name-value pairs for the main function 'bin1msSpkCountMat.m'
        
        default_align = 'center'; 
        
        p = inputParser; % create parser object
        
        addRequired(p,'SpkCountMat1ms'); % SpkCountMat to be binned
        addRequired(p,'binSize');  % binSize
        addRequired(p,'stepSize'); % stepSize
        
        addParameter(p,'align', default_align) % either 'center' to take preceding and following bins or 'subsequent' to just take following bins

        parse(p, SpkCountMat1ms, binSize, stepSize, vargs{:})
    end
end