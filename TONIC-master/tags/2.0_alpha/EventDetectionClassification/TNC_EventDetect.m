function [events] = TNC_EventDetect(bandPassedData,sampleRate,snrThresh)
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
% events.wfs
% events.inds
% events.sampleRate
% events.channel
% events.snrThresh
% events.winL
% events.winR
% events.numEvents
 
% README for the TONIC SPIKE SORTING SUB-PACKAGE
 
% Starting to use a new design for packages. Top-level directory will
% contain a few subfolders that organize the package. The distributable
% core of the package is a library or functions used to transform, analyze,
% and plot data.
% 
% Within this directory, the useable form of the distribution for the lab
% will contain a subfolder that holds "Analysis Workflows". This directory
% contains 2nd order functions that call upon the TONIC library to perform
% analysis by stringing together library functions into a workflow.
% 
 
% Algorithm for spike sorting invovles a few stages:
% 
% (0) filtering / denoising
%     acausal filtering is applied to the raw data
% (1) detection 
%     find moments at which the recording exceeds a given SNR threshold and also has a rapid recovery to a level below the SNR threshold.
% (2) extraction
%     sort the extracted events to allow largest amp events (compared to other channels) to trigger extraction of data windows from other rec channels.
% (3) alignment
%     does alignment make sense in cases of multichannel?
% (4) ts-quant
%     calculate projections onto freq or timeseries space. other scalar methods implemented as well.
% (5) sp-quant
%     calculate projections of distribution of ts-quant projections in electrode space (i.e. across channels).
% (6) id-ing
%     cluster in sp-quant or ts-quant space according to the watershed like algorithm
% 
 
%% EVENT DETECTION
tstart = tic;
 
    peakTime = 1;
 
    events.sampleRate = sampleRate;
    samples = numel(bandPassedData);
    
% find the rms values of the data to set the SNR threshold
    % standard estimator
    % stdRawData = std(bandPassedData);
    
    % robust estimator
    avgRawData = median(bandPassedData);
    stdRawData = median(abs(bandPassedData-avgRawData)) ./ 0.6745;
 
    if snrThresh < 0
        threshold = -snrThresh;
    else
        threshold = stdRawData.*snrThresh;
    end
    events.ampThresh = threshold;
    disp(['threshold=' num2str(threshold) ' stdRawData=' num2str(stdRawData) ' avgRawData=' num2str(avgRawData) ' min(bandPassedData)=' num2str(min(bandPassedData))]);
%     disp(['SD calculated to be: ' num2str(stdRawData) ' ... giving a threshold of: ' num2str(threshold)]);
 
% display the threshold choice
%     figure(2); clf;
    xValsForThresh = 1:1000:length(bandPassedData);
    yValsForThresh = ones(1,length(xValsForThresh)) .* threshold;
%     plot(1:1:length(bandPassedData),bandPassedData(1,:),'b',xValsForThresh,yValsForThresh,'r-',xValsForThresh,-yValsForThresh,'r-')
%     drawnow;
    
% first stage is to run find on the array of data values
    allSupraThreshIndices = find(bandPassedData<-threshold);

%     disp(sprintf('First pass detected events: %g',length(allSupraThreshIndices)));
%     plot(1:1:length(bandPassedData),bandPassedData(1,:),'b',allSupraThreshIndices,bandPassedData(allSupraThreshIndices),'r.');
%     drawnow;
    
    if numel(allSupraThreshIndices)==0
 
        candidateEvents = [];
        events.inds = [];
        events.ts   = [];
        
    else
 
    % eliminate any crossings without a return below threshold
        candidateEvents(1,1) = allSupraThreshIndices(1,1);
        j=2;
        for i=2:length(allSupraThreshIndices)
            if  allSupraThreshIndices(1,i)<(samples-60) && allSupraThreshIndices(1,i)>(60)
                if allSupraThreshIndices(1,i) ~= allSupraThreshIndices(1,i-1) + 1
                    if peakTime==1
                        % use the peak
                        pkLoc = find(bandPassedData(1,allSupraThreshIndices(1,i):allSupraThreshIndices(1,i)+12) == min(bandPassedData(1,allSupraThreshIndices(1,i):allSupraThreshIndices(1,i)+12)),1);
                        candidateEvents(1,j) = allSupraThreshIndices(1,i)+pkLoc-1; 
                        j=j+1;
                    else
                        candidateEvents(1,j) = allSupraThreshIndices(1,i)-2;
                        j=j+1;
                    end
 
                end
            end  
        end
 
        events.inds = candidateEvents;
        events.ts   = candidateEvents./sampleRate;
 
    end
    
telapsed = toc(tstart);
 
    numberOfEvents = length(candidateEvents);    
 
    disp(sprintf('     EVENT DETECT | Number of Events: %g | Elapsed time (seconds): %g',numberOfEvents,telapsed));
 

