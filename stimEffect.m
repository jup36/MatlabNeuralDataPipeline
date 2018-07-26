function [stimE] = stimEffect( filePath, fileName, probeDepth, varargin )
%'stimEffect' examines the effect of laser stim on the activity of individual units. 
%   It generates/stores a structure 'stimE' which contains relevant quantities such as meanReach and meanLaser Z-scores. 
%   Also, plots X-Y plots are generated to show the effect of stimulation
%   during behavioral and tagging trials. 
%   Usage example: stimEffect( filePath, 'binSpkCountCTXIT06_040218', 'binSpkCountCTX', {} )
%   Use IndividualUnitPlotGramm.m to plot individual unit PSTHs based on the outcome of 'stimEffect'.  

cd(filePath)

p = parse_input_stimEffect( filePath, fileName, probeDepth, varargin );
% p = parse_input_stimEffect( filePath, fileName, probeDepth, {} ); % use this when running line-by-line

if strcmp(p.Results.fileName(1,end-3:end),'.mat')
    fileList = dir(p.Results.fileName);
else
    fileList = dir(strcat(p.Results.fileName,'.mat'));
end

S=load(fileList.name);    % load the binSpkCount structure
%S= S.(p.Results.varName); % rename the structure simply as S (get rid of the original varName)

reachOn = find(p.Results.reachWinEdges==p.Results.reachOnTime); % reach-On (start) point within the reachWinEdges
tagStimOn = find(p.Results.tagWinEdges==p.Results.tagOnTime);   % tagStim-On (start) point within the tagWinEdges

gaussianKernel = TNC_CreateGaussian(2*25,2,2*50,1); % get a gaussian kernel for convolution (mean,std,width)

for u = 1:length(S.reach.SpkCountMatZ) % increment units
    
    stimE.meanReach(u,1) = nanmean(S.reach.SpkCountMatZ{u}(1,reachOn:reachOn+p.Results.reachDur));    % mean FR during reach period
    stimE.maxReach(u,1)  = nanmax(S.reach.SpkCountMatZ{u}(1,reachOn:reachOn+p.Results.reachDur));     % max  FR during reach period
    
    stimE.meanLaser(u,1) = nanmean(S.stmLaser.SpkCountMatZ{u}(1,reachOn:reachOn+p.Results.reachDur)); % mean FR during stimlaser period
    stimE.maxLaser(u,1)  = nanmax(S.stmLaser.SpkCountMatZ{u}(1,reachOn:reachOn+p.Results.reachDur));  % max  FR during stimlaser period
    
    stimE.meanStmReach(u,1) = nanmean(S.stmReach.SpkCountMatZ{u}(1,reachOn:reachOn+p.Results.reachDur)); % mean FR during reach period with laser on
    stimE.maxStmReach(u,1)  = nanmax(S.stmReach.SpkCountMatZ{u}(1,reachOn:reachOn+p.Results.reachDur));  % max  FR during reach period with laser on
    
    tagStmSpkCntMat   = full(cell2mat(S.tagLaser.SpkCountMat{u})); % spikecountmat around the tag stims
    stimE.sumTagStmOn(u,1)    = sum(sum(tagStmSpkCntMat(:,tagStimOn:tagStimOn+p.Results.tagDur))); % total spike counts during tag stim on 
    stimE.sumPreTagStmOn(u,1) = sum(sum(tagStmSpkCntMat(:,tagStimOn-p.Results.tagDur:tagStimOn))); % total spike counts before tag stim on
    
    stimE.tagStmOn{u,1} = sum(tagStmSpkCntMat(:,tagStimOn:tagStimOn+p.Results.tagDur),2);
    stimE.preTagStmOn{u,1} = sum(tagStmSpkCntMat(:,tagStimOn-p.Results.tagDur:tagStimOn),2);
    
    stimE.meanTagStmOn(u,1)    = sum(sum(tagStmSpkCntMat(:,tagStimOn:tagStimOn+p.Results.tagDur)))/(size(tagStmSpkCntMat,1))/(p.Results.tagDur/1000); % mean spike counts during tag stim on
    stimE.meanPreTagStmOn(u,1) = sum(sum(tagStmSpkCntMat(:,tagStimOn-p.Results.tagDur:tagStimOn)))/(size(tagStmSpkCntMat,1))/(p.Results.tagDur/1000); % mean spike counts before tag stim on
    
    % paired t-test for the tag effect
    [stimE.ttesth(u,1),stimE.ttestp(u,1)] = ttest(stimE.preTagStmOn{u,1},stimE.tagStmOn{u,1},'Alpha',0.01); % paired ttest for preTag vs Tag summed spike counts
    
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
    
    SpkReach = full(cell2mat(S.reach.SpkCountMat{u}));  % spike count mat during the whole reach period
    MeanFR   = nanmean(SpkReach(:))*1000/(length(p.Results.reachWinEdges)/1000);   % mean FR during the entire reach period
    
    stimE.FRidx(u,1) = MeanFR > p.Results.lowFRcut; % in case mean FR > 0.5 Hz
end
clearvars u

%% Plot X-Y relationship with gramm
group = zeros(sum(stimE.FRidx),1); % there's only one group currently, but modify this to include multiple groups

clear g

% X-Y plot Reach vs. stimLaser (triggered by a slight reach)
g(1,1)= gramm('x', stimE.meanReach(stimE.FRidx), 'y', stimE.meanLaser(stimE.FRidx), 'color', group);
g(1,1).set_color_options('lightness_range',[70 40],'chroma_range',[60 70],'legend','separate');
g(1,1).geom_point();
g(1,1).stat_cornerhist('edges',-10:0.5:10,'aspect',0.6);
g(1,1).geom_abline();
g(1,1).set_names('x','Reach','y','Laser');
g(1,1).set_title('ReachVsLaser z-score');

% X-Y plot Reach (laser off) vs. Reach (laser on) (triggered by a slight reach)
g(1,2)= gramm('x', stimE.meanReach(stimE.FRidx), 'y', stimE.meanStmReach(stimE.FRidx), 'color', group);
g(1,2).geom_point();
g(1,2).stat_cornerhist('edges',-10:0.5:10,'aspect',0.6);
g(1,2).geom_abline();
g(1,2).set_names('x','Reach','y','Reach (laser on)');
g(1,2).set_title('ReachVsLaser z-score');

% X-Y plot Laser ON vs. Laser OFF 
g(1,3)= gramm('x',stimE.sumPreTagStmOn(stimE.FRidx), 'y', stimE.sumTagStmOn(stimE.FRidx), 'color', group);
g(1,3).geom_point();
%g(1,3).stat_cornerhist('edges',-50:5:50,'aspect',0.6);
g(1,3).geom_abline();
g(1,3).set_names('x','Stim-Off','y','Stim-On');
g(1,3).set_title('preTagVsTag spike counts');

fig1 = figure('Position',[100 100 1650 550]);
g.draw();
print(fig1, strcat(p.Results.fileName,'_stimE_Summary'), '-dpdf'); % print figure as pdf

sites  = cell2mat(S.reach.Site); % site IDs
geoms  = cell2mat(S.reach.geometry); % electrode x-y positions 
depths = geoms(:,2); % depth from the pial surface 

% plot meanLaser-meanReach
laserVSreach=stimE.meanLaser-stimE.meanReach;
laserVSreach(:,2)=1:length(laserVSreach);
laserVSreach=sortrows(laserVSreach,1);
stimE.laserVSreach = laserVSreach;

fig2 = imecOpt3GeomColorMapZ( sites(laserVSreach(stimE.FRidx,2)), laserVSreach(stimE.FRidx,1), 'rb', true, p.Results.colorAxis, p.Results.probeDepth ); % imecOpt3GeomColorMapZ( imecSites, Z, colorScheme, drawAllSites, colorAxis, varargin )
print(fig2, strcat(p.Results.fileName,'_stimE_siteColorMap'), '-dpdf'); % print figure as pdf

% plot sumTagStmOn-sumPreTagStmOn 
tagVSpreTag = stimE.sumTagStmOn - stimE.sumPreTagStmOn; 
tagVSpreTag(:,2)=1:length(tagVSpreTag);
tagVSpreTag =sortrows(tagVSpreTag,1);
stimE.tagVSpreTag = tagVSpreTag;

fig3 = imecOpt3GeomColorMapZ( sites(tagVSpreTag(stimE.FRidx,2)), tagVSpreTag(stimE.FRidx,1), 'rb', true, p.Results.tagColorAxis, p.Results.probeDepth ); % imecOpt3GeomColorMapZ( imecSites, Z, colorScheme, drawAllSites, colorAxis, varargin )
print(fig3, strcat(p.Results.fileName,'_tagE_siteColorMap'), '-dpdf'); % print figure as pdf

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


%% Save stimE
%stimE.S = S; % do not save the binSpkCount*.mat due to the file size
stimE.p = p; 

saveName = strcat(p.Results.fileName,'_stimE');
save(saveName,'stimE')

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_stimEffect( filePath, fileName, probeDepth, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
        %parse input, and extract name-value pairs for the main function
        % 'stimEffect.m'
        
        default_reachWinEdges = -2e3+1:2e3; % default reach window
        default_reachOnTime   = 0;   % time at which reach is on
        default_reachDur      = 750; % default reach duration
        
        default_tagWinEdges = -5e3+1:5e3; % default tag window
        default_tagOnTime   = 0;   % time at which tag stim is on
        default_tagDur      = 500; % default tag stim duration
        default_tagColorAxis = [-5 5]; % default tag color axis 
        default_binStepSize = 10;  % 10 ms to be used for spike counts
        default_binSize = 50; % 50ms bin to be used for spike counts
        
        default_lowFRcut  = 0.5;    % to exclude units with low FR
        default_colorAxis = [-5 5]; % color axis for stim effect probe colormap
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addRequired(p,'fileName');
        addRequired(p,'probeDepth');
        
        addParameter(p,'reachWinEdges', default_reachWinEdges)
        addParameter(p,'reachOnTime', default_reachOnTime)
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

