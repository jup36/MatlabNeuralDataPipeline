function [ binMat, binEdges ] = bin1msSpkCountMat( SpkCountMat1ms, binSize, stepSize )
% bin1msSpkCountMat takes a trial-by-1msBin spike count matrix and generates 
%  a binned spikeCountMat (binMat) using the binSize and stepSize specified by the user.   
%  binEdges contains: the range of each column of the binMat

step = 0;    % initialize step
binMat = []; % initialize binned spike count mat
steps  = []; % initialize accumulated steps
while step + binSize <= size(SpkCountMat1ms,2) % this prevents stepping out of the range
    currBinMat = SpkCountMat1ms(:,step+1:step+binSize); % current portion of the spike count mat
    steps = [steps,step];   % accumulate steps
    step = step + stepSize; % increment step
    binMat = [binMat,sum(currBinMat,2)]; % accumulate binned spike count mat
end

binEdges(1,:) = steps(1,:); % steps
binEdges(2,:) = binEdges(1,:) + binSize; % where bin ends for each step 

