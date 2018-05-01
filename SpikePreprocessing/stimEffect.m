function [stimE] = stimEffect( filePath, fileName, probeDepth, varargin )
%'stimEffect' examines the effect of laser stim on the activity of individual units. 
%   It generates/stores a structure 'stimE' which contains relevant quantities such as meanReach and meanLaser Z-scores. 
%   Also, plots X-Y plots are generated to show the effect of stimulation
%   during behavioral and tagging trials. 
%   Usage example: stimEffect( filePath, 'binSpkCountCTXIT06_040218', 'binSpkCountCTX', {} )
%   Use IndividualUnitPlotGramm.m to plot individual unit PSTHs based on the outcome of 'stimEffect'.  

cd(filePath)

p = parse_input_stimEffect( filePath, fileName, probeDepth, varargin );

if strcmp(p.Results.fileName(1,end-3:end),'.mat')
    fileList = dir(p.Results.fileName);
else
    fileList = dir(strcat(p.Results.fileName,'.mat'));
end

S=load(fileList.name);    % load the binSpkCount structure
%S= S.(p.Results.varName); % rename the structure simply as S (get rid of the original varName)

reachOn = find(p.Results.reachWinEdges==p.Results.reachOnTime); % reach-On (start) point within the reachWinEdges
tagStimOn = find(p.Results.tagWinEdges==p.Results.tagOnTime);   % tagStim-On (start) point within the tagWinEdges

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
    %MeanTagStmOn(u)   = mean(S.tagLaser.SpkCountMatZ{u}(uv.tagTimeLogic));    % mean FR during tagging period
    %MeanPreTagStmOn(u)= mean(S.tagLaser.SpkCountMatZ{u}(uv.tagPreTimeLogic)); % mean FR during pre-tagging period
    
    SpkReach = full(cell2mat(S.reach.SpkCountMat{u}));  % spike count mat during the whole reach period
    MeanFR   = nanmean(SpkReach(:))*1000/(length(p.Results.reachWinEdges)/1000);   % mean FR during the entire reach period
    
    stimE.FRidx(u,1) = MeanFR > p.Results.lowFRcut; % in case mean FR > 0.5 Hz
end
clearvars u

%% Save stimE
stimE.S = S;
stimE.p = p; 

saveName = strcat(p.Results.fileName,'_stimE');
save(saveName,'stimE')

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

sites = cell2mat(S.reach.Site); % site IDs

laserVSreach=stimE.meanLaser-stimE.meanReach;
laserVSreach(:,2)=1:length(laserVSreach);
laserVSreach=sortrows(laserVSreach,1);

fig2 = imecOpt3GeomColorMapZ( sites(laserVSreach(stimE.FRidx,2)), laserVSreach(stimE.FRidx,1), 'rb', true, p.Results.colorAxis, p.Results.probeDepth ); % imecOpt3GeomColorMapZ( imecSites, Z, colorScheme, drawAllSites, colorAxis, varargin )
print(fig2, strcat(p.Results.fileName,'_stimE_siteColorMap'), '-dpdf'); % print figure as pdf

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
        
        parse(p, filePath, fileName, probeDepth, vargs{:})
        
    end

end

