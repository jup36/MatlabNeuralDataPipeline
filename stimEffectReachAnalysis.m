function [stimR] = stimEffectReachAnalysis( filePath, varargin )
%'stimEffectReach' examines the effect of laser stim on the reach behavior (kinematics). 
%  It generates/stores a structure 'stimR' which contains relevant quantities such as max velocity, position.  

cd(filePath)

p = parse_input_stimEffectReach( filePath, varargin );
%p = parse_input_stimEffectReach( filePath, {} ); % use this when running line-by-line

fileList = dir('BehVariables*'); % get the behavioral data BehVariables*.mat 

%% Preprocessing of the behavioral data
for f = 1:length(fileList) % increment files
    load(fileList(f).name,'pos1','positionData','ts','reach0'); % load reach info
    
    %% Detect reach Start, Stop and Position for laser and pseudo-laser trials
    if p.Results.reachBeforeLastReward
        positionData = positionData(:,1:ts.reward(end));
    else
    end    
    
    % laser trials
    [ stimR.laserReachStart{f,1}, stimR.laserReachStop{f,1}, ~, stimR.laserReachPos{f,1}, stimR.laserReachValDurIdx{f,1} ] = getPseudoReachTimesForStimTrials( positionData, ts.stmLaser, p.Results.minReachDurCut ); % detect pseudo-reach trajectories around the time of laser stimulation deliveries. 
    % pseudo-laser trials
    [ stimR.pLaserReachStart{f,1}, stimR.pLaserReachStop{f,1}, ~, stimR.pLaserReachPos{f,1}, stimR.pLaserReachValDurIdx{f,1} ] = getPseudoReachTimesForStimTrials( positionData, ts.pseudoLaser, p.Results.minReachDurCut ); %  detect the pseudo-reach trajectories around the time of pseudo Laser TTL pulses (without actual laser delivery)
    
    %% get the StmReachIdx by matching ts.stmLaser and ts.stmLaserReach 
    % index for stmReaches (laserReaches actually satisfying reach criterion)
    stimR.stimReachIdx{f,1} = zeros(length(ts.stmLaser),1);
    for r = 1:length(ts.stmLaserReach) % increment stmReach Trials
        stimR.stimReachIdx{f,1}(find(ts.stmLaser==ts.stmLaserReach(r),1)) = true; % get the tempStmReachIdx
    end
    
    %% Quantify reachSumAmp (AUC), decimated reachPos1, reachMaxPos, reachMaxVel
    % make the select reach and pseudoLaser indices 
    if p.Results.useAllTrials % in case using all trials
        stimR.selReachIdx{f,1} = ones(length(ts.reachStart),1);
        stimR.selLaserIdx{f,1} = ones(length(stimR.laserReachStart{f,1}),1);
        stimR.selPLaserIdx{f,1} = ones(length(stimR.pLaserReachStart{f,1}),1);
    else % in case selecting trials
        stimR.selReachIdx{f,1} = ts.reachStart >= ts.reward(p.Results.selectTrials(1)) & ts.reachStart <= ts.reward(p.Results.selectTrials(end)); 
        stimR.selLaserIdx{f,1} = stimR.laserReachStart{f,1} >= ts.reward(p.Results.selectTrials(1)) & stimR.laserReachStart{f,1} <= ts.reward(p.Results.selectTrials(end)); 
        stimR.selPLaserIdx{f,1} = stimR.pLaserReachStart{f,1} >= ts.reward(p.Results.selectTrials(1)) & stimR.pLaserReachStart{f,1} <= ts.reward(p.Results.selectTrials(end)); 
        %stimR.selReachIdx{f,1} = ts.reachStart >= ts.reward(50) & ts.reachStart <= ts.reward(100); 
        %stimR.selLaserIdx{f,1} = stimR.laserReachStart{f,1} >= ts.reward(50) & stimR.laserReachStart{f,1} <= ts.reward(100); 
        %stimR.selPLaserIdx{f,1} = stimR.pLaserReachStart{f,1} >= ts.reward(50) & stimR.pLaserReachStart{f,1} <= ts.reward(100); 
    end
    
    %% Quantify all reaches
    stimR.reachPos{f,1} = pos1; % store reach position in stimR structure
    stimR.reachDur{f,1} = (ts.reachStop - ts.reachStart)'; % get reach durations
   
    for i = 1:min(size(pos1,1),length(ts.reachStart)) % increment all reaches
        currReach0Idx = ts.reachStart(i):ts.reachStop(i); % current reach0 index (points between reachStart and reachStop)
        if currReach0Idx(end)<length(reach0)
           stimR.reachSumAmp{f,1}(i,1) = sum(reach0(1,currReach0Idx),2); % reach position AUC (integration)
           stimR.reachPos1Dcm{f,1}(i,:) = decimate(pos1(i,:),10); % decimate reach position data 
           stimR.reachMaxPos{f,1}(i,1) = max(reach0(1,currReach0Idx),[],2); % max position each reach
           stimR.reachMaxVel{f,1}(i,1) = max(diff(reach0(1,currReach0Idx)),[],2); % max velocity each reach 
        end
    end
    clearvars i currReach0Idx
    
    % Descriptive stats
    [stimR.reachDur{f,2}(1,1), stimR.reachDur{f,2}(1,2), stimR.reachDur{f,2}(1,3)] = meanstdsem(stimR.reachDur{f,1}); 
    [stimR.reachSumAmp{f,2}(1,1), stimR.reachSumAmp{f,2}(1,2), stimR.reachSumAmp{f,2}(1,3)] = meanstdsem(stimR.reachSumAmp{f,1}); 
    [stimR.reachMaxPos{f,2}(1,1), stimR.reachMaxPos{f,2}(1,2), stimR.reachMaxPos{f,2}(1,3)] = meanstdsem(stimR.reachMaxPos{f,1}); 
    [stimR.reachMaxVel{f,2}(1,1), stimR.reachMaxVel{f,2}(1,2), stimR.reachMaxVel{f,2}(1,3)] = meanstdsem(stimR.reachMaxVel{f,1}); 
    
    % select trials for reaches
    if p.Results.useAllTrials==false
        stimR.selReachDur{f,1} = stimR.reachDur{f,1}(stimR.selReachIdx{f,1},1); % trials x 1
        stimR.selReachSumAmp{f,1} = stimR.reachSumAmp{f,1}(stimR.selReachIdx{f,1},1); % trials x 1
        stimR.selReachPos1Dcm{f,1} = stimR.reachPos1Dcm{f,1}(stimR.selReachIdx{f,1},:); % trials x decimated time points
        stimR.selReachMaxPos{f,1} = stimR.reachMaxPos{f,1}(stimR.selReachIdx{f,1},1); % trials x 1
        stimR.selReachMaxVel{f,1} = stimR.reachMaxVel{f,1}(stimR.selReachIdx{f,1},1); % trials x 1
        
        % Descriptive stats
        [stimR.selReachDur{f,2}(1,1), stimR.selReachDur{f,2}(1,2), stimR.selReachDur{f,2}(1,3)] = meanstdsem(stimR.selReachDur{f,1}); 
        [stimR.selReachSumAmp{f,2}(1,1), stimR.selReachSumAmp{f,2}(1,2), stimR.selReachSumAmp{f,2}(1,3)] = meanstdsem(stimR.selReachSumAmp{f,1}); 
        [stimR.selReachMaxPos{f,2}(1,1), stimR.selReachMaxPos{f,2}(1,2), stimR.selReachMaxPos{f,2}(1,3)] = meanstdsem(stimR.selReachMaxPos{f,1}); 
        [stimR.selReachMaxVel{f,2}(1,1), stimR.selReachMaxVel{f,2}(1,2), stimR.selReachMaxVel{f,2}(1,3)] = meanstdsem(stimR.selReachMaxVel{f,1}); 
    end
    
    %% Quantify laser reaches
    stimR.laserReachDur{f,1} = (stimR.laserReachStop{f,1} - stimR.laserReachStart{f,1})'; % laser reach duration
    for i = 1:length(stimR.laserReachStart{f,1}) % increment laser trials
        currReach0Idx = stimR.laserReachStart{f,1}(i):stimR.laserReachStop{f,1}(i); % current reach0 index (points between laser reachStart and reachStop)
        stimR.laserReachSumAmp{f,1}(i,1) = sum(reach0(1,currReach0Idx),2); % laser reach position AUC (integration)
        stimR.laserReachPos1Dcm{f,1}(i,:) = decimate(stimR.laserReachPos{f,1}(i,:),10); % decimate reach position data of laser trials
        stimR.laserReachMaxPos{f,1}(i,1) = max(reach0(1,currReach0Idx));   % max position each laser reach
        stimR.laserReachMaxVel{f,1}(i,1) = max(diff(reach0(1,currReach0Idx)),[],2); % max velocity each laser reach
    end
    clearvars i currReach0Idx
    
    % Descriptive stats
    [stimR.laserReachDur{f,2}(1,1), stimR.laserReachDur{f,2}(1,2), stimR.laserReachDur{f,2}(1,3)] = meanstdsem(stimR.laserReachDur{f,1}); 
    [stimR.laserReachSumAmp{f,2}(1,1), stimR.laserReachSumAmp{f,2}(1,2), stimR.laserReachSumAmp{f,2}(1,3)] = meanstdsem(stimR.laserReachSumAmp{f,1}); 
    [stimR.laserReachMaxPos{f,2}(1,1), stimR.laserReachMaxPos{f,2}(1,2), stimR.laserReachMaxPos{f,2}(1,3)] = meanstdsem(stimR.laserReachMaxPos{f,1}); 
    [stimR.laserReachMaxVel{f,2}(1,1), stimR.laserReachMaxVel{f,2}(1,2), stimR.laserReachMaxVel{f,2}(1,3)] = meanstdsem(stimR.laserReachMaxVel{f,1});
    
    % select trials for laser reaches
    if p.Results.useAllTrials==false
        stimR.laserSelValTrialIdx{f,1} = stimR.selLaserIdx{f,1} & stimR.laserReachValDurIdx{f,1}';
        stimR.selLaserReachDur{f,1} = stimR.laserReachDur{f,1}(stimR.laserSelValTrialIdx{f,1},1); % trials x 1
        stimR.selLaserReachSumAmp{f,1} = stimR.laserReachSumAmp{f,1}(stimR.laserSelValTrialIdx{f,1},1); % 1 x trials
        stimR.selLaserReachPos1Dcm{f,1} = stimR.laserReachPos1Dcm{f,1}(stimR.laserSelValTrialIdx{f,1},:); % trials x decimated time points
        stimR.selLaserReachMaxPos{f,1} = stimR.laserReachMaxPos{f,1}(stimR.laserSelValTrialIdx{f,1},1); % trials x 1
        stimR.selLaserReachMaxVel{f,1} = stimR.laserReachMaxVel{f,1}(stimR.laserSelValTrialIdx{f,1},1); % trials x 1
    
        % Descriptive stats
        [stimR.selLaserReachDur{f,2}(1,1), stimR.selLaserReachDur{f,2}(1,2), stimR.selLaserReachDur{f,2}(1,3)] = meanstdsem(stimR.selLaserReachDur{f,1});
        [stimR.selLaserReachSumAmp{f,2}(1,1), stimR.selLaserReachSumAmp{f,2}(1,2), stimR.selLaserReachSumAmp{f,2}(1,3)] = meanstdsem(stimR.selLaserReachSumAmp{f,1}); 
        [stimR.selLaserReachMaxPos{f,2}(1,1), stimR.selLaserReachMaxPos{f,2}(1,2), stimR.selLaserReachMaxPos{f,2}(1,3)] = meanstdsem(stimR.selLaserReachMaxPos{f,1}); 
        [stimR.selLaserReachMaxVel{f,2}(1,1), stimR.selLaserReachMaxVel{f,2}(1,2), stimR.selLaserReachMaxVel{f,2}(1,3)] = meanstdsem(stimR.selLaserReachMaxVel{f,1});    
    end
    
    %% Quantify pseudo-laser reaches
    stimR.pLaserReachDur{f,1} = (stimR.pLaserReachStop{f,1} - stimR.pLaserReachStart{f,1})'; % pseudo-laser reach duration
    for i = 1:length(stimR.pLaserReachStart{f,1}) % increment pseudo-laser trials
        currReach0Idx = stimR.pLaserReachStart{f,1}(i):stimR.pLaserReachStop{f,1}(i); % current reach0 index (points between pseudo-laser reachStart and reachStop)
        stimR.pLaserReachSumAmp{f,1}(i,1) = sum(reach0(1,currReach0Idx),2); % pseudo-laser reach position AUC (integration)
        stimR.pLaserReachPos1Dcm{f,1}(i,:) = decimate(stimR.pLaserReachPos{f,1}(i,:),10); % decimate reach position data of pseudo-laser trials
        stimR.pLaserReachMaxPos{f,1}(i,1) = max(reach0(1,currReach0Idx)); % max position each pseudo-laser reach
        stimR.pLaserReachMaxVel{f,1}(i,1) = max(diff(reach0(1,currReach0Idx)),[],2); % max velocity each pseudo-laser reach
    end
    clearvars i currReach0Idx
    
    % Descriptive stats
    [stimR.pLaserReachDur{f,2}(1,1), stimR.pLaserReachDur{f,2}(1,2), stimR.pLaserReachDur{f,2}(1,3)] = meanstdsem(stimR.pLaserReachDur{f,1}); 
    [stimR.pLaserReachSumAmp{f,2}(1,1), stimR.pLaserReachSumAmp{f,2}(1,2), stimR.pLaserReachSumAmp{f,2}(1,3)] = meanstdsem(stimR.pLaserReachSumAmp{f,1}); 
    [stimR.pLaserReachMaxPos{f,2}(1,1), stimR.pLaserReachMaxPos{f,2}(1,2), stimR.pLaserReachMaxPos{f,2}(1,3)] = meanstdsem(stimR.pLaserReachMaxPos{f,1}); 
    [stimR.pLaserReachMaxVel{f,2}(1,1), stimR.pLaserReachMaxVel{f,2}(1,2), stimR.pLaserReachMaxVel{f,2}(1,3)] = meanstdsem(stimR.pLaserReachMaxVel{f,1});
    
    % select trials for pseudo-laser reaches
    if p.Results.useAllTrials==false
        stimR.selpLaserReachDur{f,1} = stimR.pLaserReachDur{f,1}(stimR.selPLaserIdx{f,1},1);      % trials x 1
        stimR.selpLaserReachSumAmp{f,1} = stimR.pLaserReachSumAmp{f,1}(stimR.selPLaserIdx{f,1},1); % 1 x trials
        stimR.selpLaserReachPos1Dcm{f,1} = stimR.pLaserReachPos1Dcm{f,1}(stimR.selPLaserIdx{f,1},:); % trials x decimated time points
        stimR.selpLaserReachMaxPos{f,1} = stimR.pLaserReachMaxPos{f,1}(stimR.selPLaserIdx{f,1},1); % trials x 1
        stimR.selpLaserReachMaxVel{f,1} = stimR.pLaserReachMaxVel{f,1}(stimR.selPLaserIdx{f,1},1); % trials x 1
    
        % Descriptive stats
        [stimR.selpLaserReachDur{f,2}(1,1), stimR.selpLaserReachDur{f,2}(1,2), stimR.selpLaserReachDur{f,2}(1,3)] = meanstdsem(stimR.selpLaserReachDur{f,1}); 
        [stimR.selpLaserReachSumAmp{f,2}(1,1), stimR.selpLaserReachSumAmp{f,2}(1,2), stimR.selpLaserReachSumAmp{f,2}(1,3)] = meanstdsem(stimR.selpLaserReachSumAmp{f,1}); 
        [stimR.selpLaserReachMaxPos{f,2}(1,1), stimR.selpLaserReachMaxPos{f,2}(1,2), stimR.selpLaserReachMaxPos{f,2}(1,3)] = meanstdsem(stimR.selpLaserReachMaxPos{f,1}); 
        [stimR.selpLaserReachMaxVel{f,2}(1,1), stimR.selpLaserReachMaxVel{f,2}(1,2), stimR.selpLaserReachMaxVel{f,2}(1,3)] = meanstdsem(stimR.selpLaserReachMaxVel{f,1});
    end
end
clearvars f

%% save stimR
if ismac
    delimitIdx = strfind(filePath,'/');
elseif ispc
    delimitIdx = strfind(filePath,'\');
else
    disp('Platform not supported')
end

saveName = strcat(filePath(delimitIdx(end)+1:end),'_stimR');
save(saveName,'stimR','p')

%% Plot with Gramm
stimR.expLabel = [zeros(length(stimR.reachMaxPos{1}),1)+1; zeros(length(stimR.pLaserReachMaxPos{1}),1)+2; zeros(length(stimR.laserReachMaxPos{1}),1)+3]; % zeros(length(maxPosConStimReaches),1)+3]; % experiment label
stimR.maxPosArr = [ stimR.reachMaxPos{1}; stimR.pLaserReachMaxPos{1}; stimR.laserReachMaxPos{1} ]; % concatenated max amp (pos1) data
stimR.maxVelArr = [ stimR.reachMaxVel{1}; stimR.pLaserReachMaxVel{1}; stimR.laserReachMaxVel{1} ]; % concatenated max vel (diff(pos1)) data

color = [zeros(length(stimR.reachMaxPos{1}),1); zeros(length(stimR.pLaserReachMaxPos{1}),1)+1; zeros(length(stimR.laserReachMaxPos{1}),1)+2]; % to mark the stim reach trials with a different color

clear g
g(1,1)=gramm('x',stimR.expLabel,'y',stimR.maxPosArr,'color',color);
g(1,2)=gramm('x',stimR.expLabel,'y',stimR.maxVelArr,'color',color);

g(1,1).geom_jitter('width',0.4,'height',0);
g(1,1).set_names('x','label','y','maxPos');
g(1,1).set_title('Max Pos');
 
g(1,2).geom_jitter('width',0.4,'height',0);
g(1,2).set_names('x','label','y','maxVel');
g(1,2).set_title('Max Vel');

fig1 = figure('Position',round([100 100 800 400]./1.5));
g.draw();

print(fig1, saveName, '-dpdf'); % print figure as pdf

%gf.coord_flip();
%figure('Position',[100 100 800 400]);
%gf.draw();

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_stimEffectReach( filePath, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
        % parse input, and extract name-value pairs
        
        default_useAllTrials = true;  % logical to indicate using all trials in analyses
        default_selectTrials = 1:150; % default select trials
        default_minReachDurCut = 50;  % default minimum reach duration cut off (e.g. 50 ms) - to deselect unless the detected reach duration is greater than this value
        default_reachBeforeLastReward = true; % detect only the reaches occur before the last reward delivery
        
        p = inputParser; % create parser object
        
        addRequired(p,'filePath');
        addParameter(p,'useAllTrials', default_useAllTrials)
        addParameter(p,'selectTrials', default_selectTrials)
        addParameter(p,'minReachDurCut', default_minReachDurCut)
        addParameter(p,'reachBeforeLastReward', default_reachBeforeLastReward)
        
        parse(p,filePath,vargs{:})
        
    end


end