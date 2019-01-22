function [meanFRxCorrUnitIdx] = getGlobalUnitIdx(filePath, fileName, eventName, varargin)
%This function identifies a unit index satisfying the mean FR and the
% pairwise cross-correlation criteria, this index can be used to pick units
% to include in subsequent dimensionality reduction analyses. 

p = parse_input_getGlobalUnitIdx( filePath, fileName, eventName, varargin );
% p = parse_input_getGlobalUnitIdx( filePath, fileName, eventName, {'frHighPass', 2, 'frLowPass', 50} )

unitTimeTrialC = cell(1,length(eventName)); 
%% inspect mean FR of the single units 
for e = 1:length(eventName)
    S = load(fullfile(filePath,fileName),eventName{e}); 
    S = S.(eventName{e}); 
    unitTimeTrialC{1,e} = S.unitTimeTrial; 
end
clearvars e

meanUnitTimeTrialC = cellfun(@(x) nanmean(nanmean(x,2),3).*1000, unitTimeTrialC, 'UniformOutput', false);
unitMeanFR  = nanmean(cell2mat(meanUnitTimeTrialC),2); 

%% inspect pairwise high cross-correlation (remove duplicate/split units)
if sum(strcmpi(eventName,'reach'))==0
    error('reach PSTH has not been provided to identify the high cross-corr neuronal pairs with!')
elseif sum(strcmpi(eventName,'reach'))==1
    reachUnitTimeTrial = unitTimeTrialC{strcmpi(eventName,'reach')};
end

unitTimeAcrossTrials = reshape(reachUnitTimeTrial, [size(reachUnitTimeTrial,1), size(reachUnitTimeTrial,2)*size(reachUnitTimeTrial,3)]); % reshape the unitTimeTrial matrix to get unitTimeAcrossTrials (concatenate across trials)
if p.Results.redefineXcorrUnits
    xcorUnitIdx = getUnitOfHighXcorr( unitTimeAcrossTrials, p.Results.xcorThresholdPer );
    if contains(fileName,'Ctx','IgnoreCase',true)
        save(fullfile(p.Results.filePath,strcat('xcorUnitIdx','Ctx')), 'xcorUnitIdx') % to speed up running again in the future save this idx, and bypass the process of examining all unit pairs
    elseif contains(fileName,'Str','IgnoreCase',true)
        save(fullfile(p.Results.filePath,strcat('xcorUnitIdx','Str')), 'xcorUnitIdx') % to speed up running again in the future save this idx, and bypass the process of examining all unit pairs
    end    
else
    if xor(contains(fileName,'Ctx','IgnoreCase',true), contains(fileName,'Str','IgnoreCase',true))
        if contains(fileName,'Ctx','IgnoreCase',true)
            if exist(fullfile(filePath,'xcorUnitIdxCtx.mat'),'file')==2 % if the file exists
                load(fullfile(filePath,'xcorUnitIdxCtx.mat'),'xcorUnitIdx'); % just load it instead of repeating the xcorr process
            else
                xcorUnitIdx = getUnitOfHighXcorr( unitTimeAcrossTrials, p.Results.xcorThresholdPer );
                save(fullfile(p.Results.filePath,strcat('xcorUnitIdx','Ctx')), 'xcorUnitIdx') % to speed up running again in the future save this idx, and bypass the process of examining all unit pairs
            end
        elseif contains(fileName,'Str','IgnoreCase',true)
            if exist(fullfile(filePath,'xcorUnitIdxStr.mat'),'file')==2 % if the file exists
                load(fullfile(filePath,'xcorUnitIdxStr.mat'),'xcorUnitIdx'); % just load it instead of repeating the xcorr process
            else
                xcorUnitIdx = getUnitOfHighXcorr( unitTimeAcrossTrials, p.Results.xcorThresholdPer );
                save(fullfile(p.Results.filePath,strcat('xcorUnitIdx','Str')), 'xcorUnitIdx') % to speed up running again in the future save this idx, and bypass the process of examining all unit pairs
            end
        end
    else
        error('Check if the input variable fileName contains a proper brain region info!!!')
    end
end

meanFRxCorrUnitIdx = unitMeanFR > p.Results.frHighPass & unitMeanFR < p.Results.frLowPass & xcorUnitIdx; % units to be used for pca

% save the unitIdx
if xor(contains(fileName,'Ctx','IgnoreCase',true),contains(fileName,'Str','IgnoreCase',true))
    if contains(fileName,'Ctx','IgnoreCase',true)
        save(fullfile(p.Results.filePath,strcat('meanFRxCorrUnitIdx','Ctx')), 'meanFRxCorrUnitIdx')
    elseif contains(fileName,'Str','IgnoreCase',true)
        save(fullfile(p.Results.filePath,strcat('meanFRxCorrUnitIdx','Str')), 'meanFRxCorrUnitIdx')
    else
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_getGlobalUnitIdx( filePath, fileName, eventName, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
        
        %parse input, and extract name-value pairs for the main function
        % 'runPCA'
        default_frHighPass  = 1;  % low FR cut (3Hz)
        default_frLowPass   = 60; % high FR cut (60Hz) % low pass might not be necessary
        default_xcorThresholdPer = 0.2; % default cross-correlation spike cooccurrence tolerance (i.e., 20% of the total spike count of the fewer spiking unit of the pair)
        default_redefineXcorrUnits = false; 
        
        p = inputParser; % create parser object
        addRequired(p,'filePath'); % file directory
        addRequired(p,'fileName'); % file name
        addRequired(p,'eventName'); % event name, e.g. reach, tagLaser - the events around which the psths were taken
        
        addParameter(p,'frHighPass', default_frHighPass)
        addParameter(p,'frLowPass', default_frLowPass) % FR high cut might not be necessary
        addParameter(p,'xcorThresholdPer', default_xcorThresholdPer)
        addParameter(p,'redefineXcorrUnits', default_redefineXcorrUnits)
        
        parse(p, filePath, fileName, eventName, vargs{:})
    end

end