function [stimE] = stimEffectReachToGrasp( filePath, fileName, probeDepth, varargin )
%'stimEffectReachToGrasp' examines the effect of laser stim on the activity of individual units. 
%  

% filePath = '/Volumes/Beefcake/Junchol_Data/jayReachToGrasp/M314Sim1GtACr2Junchol/M314_20200427_1000um_g0'; 
cd(filePath)

p = parse_input_stimEffectReachToGrasp( filePath, fileName, probeDepth, varargin );
% p = parse_input_stimEffectReachToGrasp( filePath, 'binSpkCountCTXM314_20200427_1000um', 1000, {} ); % use this when running line-by-line

if strcmp(p.Results.fileName(1,end-3:end),'.mat')
    fileList = dir(p.Results.fileName);
else
    fileList = dir(strcat(p.Results.fileName,'.mat'));
end

S=load(fileList.name);    % load the binSpkCount structure
%S= S.(p.Results.varName); % rename the structure simply as S (get rid of the original varName)

cueOn = find(p.Results.cueWinEdges==p.Results.cueOnTime); % reach-On (start) point within the reachWinEdges
tagStimOn = find(p.Results.tagWinEdges==p.Results.tagOnTime);   % tagStim-On (start) point within the tagWinEdges

gaussianKernel = TNC_CreateGaussian(2*25,2,2*50,1); % get a gaussian kernel for convolution (mean,std,width)

for u = 1:length(S.cueNoLaser.SpkCountMatZ) % increment units
    
    stimE.meanCueNoLaser(u,1) = nanmean(S.cueNoLaser.SpkCountMatZ{u}(1,cueOn:cueOn+p.Results.reachDur));    % mean FR during reach period
    stimE.maxCueNoLaser(u,1)  = nanmax(S.cueNoLaser.SpkCountMatZ{u}(1,cueOn:cueOn+p.Results.reachDur));     % max  FR during reach period
    
    if isfield(S,'laserCue2s')
        stimE.meanLaserCue2s(u,1) = nanmean(S.laserCue2s.SpkCountMatZ{u}(1,cueOn:cueOn+p.Results.reachDur)); % mean FR during stimlaser period
        stimE.maxLaserCue2s(u,1)  = nanmax(S.laserCue2s.SpkCountMatZ{u}(1,cueOn:cueOn+p.Results.reachDur));  % max  FR during stimlaser period
    end
    
    if isfield(S,'laserOnly2s') && isfield(S,'tagLaser1s')
        tagStmSpkCntMat = cell2mat(getSpkCntMatFromSpkTimes( [S.tagLaser1s.SpkTimes{u}; S.laserOnly2s.SpkTimes{u}], S.tagLaser1s.params)); % get the current unit's spikeCountMat (trial-by-1msBin)
    elseif ~isfield(S,'laserOnly2s') && isfield(S,'tagLaser1s')
        tagStmSpkCntMat = cell2mat(getSpkCntMatFromSpkTimes( S.tagLaser1s.SpkTimes{u}, S.tagLaser1s.params)); % get the current unit's spikeCountMat (trial-by-1msBin)
    elseif isfield(S,'laserOnly2s') && ~isfield(S,'tagLaser1s')
        tagStmSpkCntMat = cell2mat(getSpkCntMatFromSpkTimes( S.laserOnly2s.SpkTimes{u}, S.laserOnly2s.params)); % get the current unit's spikeCountMat (trial-by-1msBin)
    end
    
    stimE.sumTagStmOn(u,1)    = sum(sum(tagStmSpkCntMat(:,tagStimOn+50:tagStimOn+p.Results.tagDur))); % total spike counts during tag stim on 
    stimE.sumPreTagStmOn(u,1) = sum(sum(tagStmSpkCntMat(:,tagStimOn-p.Results.tagDur:tagStimOn))); % total spike counts before tag stim on
    
    stimE.tagStmOn{u,1} = sum(tagStmSpkCntMat(:,tagStimOn+50:tagStimOn+p.Results.tagDur),2);
    stimE.preTagStmOn{u,1} = sum(tagStmSpkCntMat(:,tagStimOn-p.Results.tagDur:tagStimOn),2);
    
    stimE.meanTagStmOn(u,1)    = sum(sum(tagStmSpkCntMat(:,tagStimOn+50:tagStimOn+p.Results.tagDur)))/(size(tagStmSpkCntMat,1))/(p.Results.tagDur/1000); % mean spike counts during tag stim on
    stimE.meanPreTagStmOn(u,1) = sum(sum(tagStmSpkCntMat(:,tagStimOn-p.Results.tagDur:tagStimOn)))/(size(tagStmSpkCntMat,1))/(p.Results.tagDur/1000); % mean spike counts before tag stim on
    
    % paired t-test for the tag effect
    [stimE.ttesth(u,1),stimE.ttestp(u,1)] = ttest(stimE.preTagStmOn{u,1},stimE.tagStmOn{u,1},'Alpha',0.05); % paired ttest for preTag vs Tag summed spike counts
    
    if stimE.ttesth(u,1)==1 && (stimE.sumTagStmOn(u,1)>stimE.sumPreTagStmOn(u,1))
        stimE.tagE(u,1) = 1; % put 1 for units with excitatory tagging effect 
        stimE.periTagMeanSpkC{u,1} = mean(bin1msSpkCountMat(tagStmSpkCntMat(:,tagStimOn:tagStimOn+p.Results.tagDur-1),p.Results.binSize,p.Results.binStepSize),1).*(1000/p.Results.binSize); % bin
        stimE.periTagMeanSpkC{u,1} = stimE.periTagMeanSpkC{u,1}-stimE.meanPreTagStmOn(u,1); % rezero the bin spike count mat (in Hz)
        stimE.periTagMeanSpkC{u,1} = conv(stimE.periTagMeanSpkC{u,1},gaussianKernel,'same'); % convolution of the delta function with the Gaussian kernel
        [stimE.tagEPeak(u,1)]  = max(stimE.periTagMeanSpkC{u,1}); % get the maximal effect size
        stimE.tagEHalfPeakT(u,1) =  find(stimE.periTagMeanSpkC{u,1}>stimE.tagEPeak(u,1)/2,1)*p.Results.binStepSize;
        
    elseif stimE.ttesth(u,1)==1 && (stimE.sumTagStmOn(u,1)<stimE.sumPreTagStmOn(u,1))
        stimE.tagE(u,1) = -1; % put -1 for units with inhibitory tagging effect
        stimE.periTagMeanSpkC{u,1} = mean(bin1msSpkCountMat(tagStmSpkCntMat(:,tagStimOn:tagStimOn+p.Results.tagDur-1),p.Results.binSize,p.Results.binStepSize),1).*(1000/p.Results.binSize); % bin
        stimE.periTagMeanSpkC{u,1} = stimE.periTagMeanSpkC{u,1}-stimE.meanPreTagStmOn(u,1); % rezero the bin spike count mat (in Hz)
        stimE.periTagMeanSpkC{u,1} = conv(stimE.periTagMeanSpkC{u,1},gaussianKernel,'same'); % convolution of the delta function with the Gaussian kernel
        [stimE.tagEPeak(u,1)] = min(stimE.periTagMeanSpkC{u,1}); % get the maximal effect size
        stimE.tagEHalfPeakT(u,1) =  find(stimE.periTagMeanSpkC{u,1}<stimE.tagEPeak(u,1)/2,1)*p.Results.binStepSize;
        
    elseif stimE.ttesth(u,1)==0
        stimE.tagE(u,1) = 0; % put 0 for units without effect 
        stimE.tagEPeak(u,1) = NaN;
        stimE.tagEHalfPeakT(u,1) = NaN;
    end
    
    SpkCueNoLaser = cell2mat(getSpkCntMatFromSpkTimes( S.cueNoLaser.SpkTimes{u}, S.cueNoLaser.params )); % get the current unit's spikeCountMat (trial-by-1msBin)
    SpkCueNoLaserBase = SpkCueNoLaser(:,3000:5000); 
    %SpkReach = full(cell2mat(S.reach.SpkCountMat{u}));  % spike count mat during the whole reach period
    MeanFR   = nanmean(SpkCueNoLaserBase(:))*1000/2;   % mean FR during the entire reach period
    
    stimE.FRidx(u,1) = MeanFR > p.Results.lowFRcut; % in case mean FR > 0.5 Hz
end
clearvars u

%% Plot X-Y relationship with gramm
group = zeros(sum(stimE.FRidx),1); % there's only one group currently, but modify this to include multiple groups
if  isfield(S,'laserCue2s')
    clear g
    
    % X-Y plot cueNoLaser vs. cueLaser (triggered by a slight reach)
    g(1,1)= gramm('x', stimE.meanCueNoLaser(stimE.FRidx), 'y', stimE.meanLaserCue2s(stimE.FRidx), 'color', group);
    g(1,1).set_color_options('lightness_range',[70 40],'chroma_range',[60 70],'legend','separate');
    g(1,1).geom_point();
    g(1,1).stat_cornerhist('edges',-10:0.5:10,'aspect',0.6);
    g(1,1).geom_abline();
    g(1,1).set_names('x','cueNoLaser','y','cueLaser');
    g(1,1).set_title('CueNoLaserVsCueLaser z-score');
    
    % X-Y plot Tag Laser ON vs. OFF
    g(1,2)= gramm('x',stimE.sumPreTagStmOn(stimE.FRidx), 'y', stimE.sumTagStmOn(stimE.FRidx), 'color', group);
    g(1,2).geom_point();
    %g(1,3).stat_cornerhist('edges',-50:5:50,'aspect',0.6);
    g(1,2).geom_abline();
    g(1,2).set_names('x','Stim-Off','y','Stim-On');
    g(1,2).set_title('preTagVsTag spike counts');
    
    fig1 = figure('Position',[100 100 1650 550]);
    g.draw();
    
    if isfolder(fullfile(filePath,'Figure')) % if there's Figure folder already
        print(fig1, fullfile(filePath,'Figure',strcat(p.Results.fileName,'_stimE_Summary')), '-dpdf','-bestfit'); % print figure as pdf
    else
        mkdir(fullfile(filePath,'Figure')) % otherwise, make a folder
        print(fig1, fullfile(filePath,'Figure',strcat(p.Results.fileName,'_stimE_Summary')), '-dpdf','-bestfit');
    end
    
    % plot meanLaser-meanReach
    laserVSreach=stimE.meanLaserCue2s-stimE.meanCueNoLaser;
    laserVSreach(:,2)=1:length(laserVSreach);
    laserVSreach=sortrows(laserVSreach,1);
    stimE.laserVSreach = laserVSreach;
end
sites  = cell2mat(S.laserOnly2s.Site); % site IDs
geoms  = cell2mat(S.laserOnly2s.geometry); % electrode x-y positions
depths = geoms(:,2); % depth from the pial surface

% plot sumTagStmOn-sumPreTagStmOn 
tagVSpreTag = stimE.sumTagStmOn - stimE.sumPreTagStmOn; 
tagVSpreTag(:,2)=1:length(tagVSpreTag);
tagVSpreTag =sortrows(tagVSpreTag,1);
stimE.tagVSpreTag = tagVSpreTag;


% plot tag effect
scatterColorSpace = linspace(min(abs(stimE.tagEHalfPeakT)),max(abs(stimE.tagEHalfPeakT)),20); % to color code each scatter by the tagEHalfPeakT
scatterColor = nan(length(stimE.tagEHalfPeakT),1);
for u = 1:length(stimE.tagEHalfPeakT)
    if ~isnan(stimE.tagEHalfPeakT(u))
        scatterColor(u,1) = find(scatterColorSpace>=stimE.tagEHalfPeakT(u),1); % match the color from the scatterColorSpace
    end
end

% generate an X-Y plot; Tag off vs. Tag on activity of each units with color-coding by the latency to the half-max tagging effect
figure; 
hold on;
scatter(stimE.meanPreTagStmOn(stimE.tagE~=0),stimE.meanTagStmOn(stimE.tagE~=0),50,scatterColor(stimE.tagE~=0),'filled') % plot the units with significant tagging effect
scatter(stimE.meanPreTagStmOn(stimE.tagE==0),stimE.meanTagStmOn(stimE.tagE==0),20,'b') % plot the units without significant tagging effect
plot(-5:max([stimE.meanPreTagStmOn;stimE.meanTagStmOn])+5,-5:max([stimE.meanPreTagStmOn;stimE.meanTagStmOn])+5,':r') % just to add a diagonal line
xlim([-5 max([stimE.meanPreTagStmOn;stimE.meanTagStmOn])+5]);
ylim([-5 max([stimE.meanPreTagStmOn;stimE.meanTagStmOn])+5]);
pbaspect([1 1 1])
hold off
print(fullfile(filePath,'Figure',strcat(p.Results.fileName,'_tagEpreTagOffOnXY')), '-dpdf','-bestfit'); % print figure as pdf

% generate an X-Y plot; activity change by tagging vs. depth from the pial surface with color-coding by the latency to the half-max tagging effect
figure; 
hold on; 
scatter(stimE.meanTagStmOn(stimE.tagE~=0)-stimE.meanPreTagStmOn(stimE.tagE~=0),depths(stimE.tagE~=0),50,scatterColor(stimE.tagE~=0),'filled') % plot the units with significant tagging effect
scatter(stimE.meanTagStmOn(stimE.tagE==0)-stimE.meanPreTagStmOn(stimE.tagE==0),depths(stimE.tagE==0),20,'b') % plot the units without significant tagging effect
%xlim([-40 15]);
%ylim([]);
pbaspect([1 1 1])
set(gca,'TickDir','out')
set(gca,'YDir','reverse')
hold off
print(fullfile(filePath,'Figure',strcat(p.Results.fileName,'_tagEbyDepth')), '-dpdf','-bestfit'); % print figure as pdf

%% Save stimE
%stimE.S = S; % do not save the binSpkCount*.mat due to the file size
stimE.p = p; 

saveName = strcat(p.Results.fileName,'_stimE');
save(saveName,'stimE')

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_stimEffectReachToGrasp( filePath, fileName, probeDepth, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
        %parse input, and extract name-value pairs for the main function
        % 'stimEffect.m'
        
        default_cueWinEdges = -5e3+1:5e3; % default reach window
        default_cueOnTime   = 0;   % time at which reach is on
        default_reachDur      = 1000; % default reach duration
        
        default_tagWinEdges = -5e3+1:5e3; % default tag window
        default_tagOnTime   = 0;   % time at which tag stim is on
        default_tagDur      = 500; % default tag stim duration
        default_tagColorAxis = [-5 5]; % default tag color axis
        default_binStepSize = 10;  % 10 ms to be used for spike counts
        default_binSize = 50; % 50ms bin to be used for spike counts
        
        default_lowFRcut  = 0.25;    % to exclude units with low FR
        default_colorAxis = [-5 5]; % color axis for stim effect probe colormap
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addRequired(p,'fileName');
        addRequired(p,'probeDepth');
        
        addParameter(p,'cueWinEdges', default_cueWinEdges)
        addParameter(p,'cueOnTime', default_cueOnTime)
        addParameter(p,'reachDur', default_reachDur)
        addParameter(p,'tagWinEdges', default_tagWinEdges)
        addParameter(p,'tagOnTime', default_tagOnTime)
        addParameter(p,'tagDur', default_tagDur)
        addParameter(p,'lowFRcut', default_lowFRcut)
        addParameter(p,'colorAxis', default_colorAxis)
        addParameter(p,'tagColorAxis', default_tagColorAxis)
        addParameter(p,'binStepSize', default_binStepSize)
        addParameter(p,'binSize', default_binSize)
        
        parse(p, filePath, fileName, probeDepth, vargs{:})
        
    end

end