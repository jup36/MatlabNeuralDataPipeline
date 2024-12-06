function js2p0_tbytSpkHandJsPreprocess_50ms_stim_parse(filePath)
%This is a preprocessing function to first demarcate trials into blocks of
% different joystick load & position combinations using the function 'jkvtBlockParse'.
% Then it gets trial-by-trial binned (e.g. 20-ms bin) spike count matrices aligned to
% an event ('reachStart' or 'pullStart' or 'trJoystickReady') depending on the trial type.
% For each trial, hand/joystick trajectories (from DLC) and joystick kinetics
% (force applied from joystick encoder) are acquired, if available, with the same alignment
% as the spike data but curtailed as those kinematic informations get cut
% off right after the end of execution. Additionally, some useful
% trial-by-trial behavioral variables are computed for future decoding
% anlaysis such as reach angle, max force, initial hand position etc.
% The output structure 'ss', 'jkvt', 'trI' are saved in the same filePath.
% filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';

%% 1. Load data
%clc; clearvars; close all;
%filePath = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles';
cd(filePath)
% neural and behavioral data
spkDir_CtxStr = dir('binSpkCountSTRCTX*');
spkDir_Cg = dir('binSpkCountCg*');

if ~isempty(spkDir_CtxStr)
    load(fullfile(spkDir_CtxStr(1).folder, spkDir_CtxStr(1).name),'spkTimesCell','jkvt')
    spkTimesCell_CtxStr = spkTimesCell; clearvars spkTimesCell
    % cortex/striatum index
    depth = cell2mat(cellfun(@(a) a(2), spkTimesCell_CtxStr(4,:),'un',0))'; % depth from pial surface
    ctxI = depth<1900; % cortex index
    spkTimesCellCTX = spkTimesCell_CtxStr(:,ctxI);
    depthCtx = depth(ctxI);
    strI = depth>2100; % striatum index
    spkTimesCellSTR = spkTimesCell_CtxStr(:,strI);
    depthStr = depth(strI);
end

if ~isempty(spkDir_Cg)
    load(fullfile(spkDir_Cg(1).folder, spkDir_Cg(1).name),'spkTimesCell')
    if ~exist('jkvt', 'var')
        load(fullfile(spkDir_Cg(1).folder, spkDir_Cg(1).name), 'jkvt')
    end

    spkTimesCellCg = spkTimesCell; clearvars spkTimesCell
    depthCg = cell2mat(cellfun(@(a) a(2), spkTimesCellCg(4,:),'un',0))'; % depth from pial surface
end

%% align hand trajectories to neural data
vfT = {jkvt(:).vFrameTime}'; % video frame time
hTrj = {jkvt(:).hTrjF}'; % hand trajectory

sm_kernel = TNC_CreateGaussian(250,30,500,20); % a kernel for smoothing (mu, sigma, time, dT)
binSize = 50; % 50 ms
spkBin = -1000:2000;

%% parse trial/block information
% pull torque & reach position combinations
pTqs = round(unique([jkvt(:).pull_torque])./10); % pull torque list (divided by 10 for ordering)
rPos = round(unique([jkvt(:).reachP1])*10); % reach position list (multiplied by 10 for ordering)
posTq = unique(repmat(pTqs,[length(rPos) 1])+repmat(rPos',[1 length(pTqs)])); % position-torque combinations
posTqCnt = zeros(length(posTq),1); % count occurrence of each combination

%% parse block transition and block type
jkvt = jkvtBlockParse(jkvt);

%% align behavioral and neural trajectories
% determine where to align based on trial types
% Trial type index
trI.spI = cell2mat(cellfun(@(a) strcmpi(a,'sp'),{jkvt(:).trialType},'un',0)); % successful pull
trI.toI = cell2mat(cellfun(@(a) strcmpi(a,'to'),{jkvt(:).trialType},'un',0)); % time out
trI.psI = cell2mat(cellfun(@(a) strcmpi(a,'ps'),{jkvt(:).trialType},'un',0)); % push
trI.pmppI = cell2mat(cellfun(@(a) strcmpi(a,'pmpp'),{jkvt(:).trialType},'un',0)); % premature pull and push
trI.pmI = cell2mat(cellfun(@(a) strcmpi(a,'pm'),{jkvt(:).trialType},'un',0)); % premature pull and push

stimOnRelToTrStart = nanmean(cell2mat(cellfun(@(a, b) a-b, {jkvt.trStart}, {jkvt.stimLaserOn}, 'UniformOutput', false)));

% get unitTimeBin aligned to stimOn (jkvt(t).stimLaserOn) or pseudoStimOn (jkvt(t).pLaserOn)
if ~isfield(jkvt, 'pLaserOn')
    stimI = ~isnan([jkvt.stimLaserOn]);
    %stimOnRelToTrStart = nanmean(cell2mat(cellfun(@(a, b) a-b, {jkvt.trStart}, {jkvt.stimLaserOn}, 'UniformOutput', false)));
    for jj = find(~stimI)
        jkvt(jj).pLaserOn = jkvt(jj).trStart-stimOnRelToTrStart;
    end
end

for t = 1:size(jkvt,2)
    %% get a time point to align depending on the trial type
    ss(t).trialType = jkvt(t).trialType;
    % get time point to which neural data and hand trajectories are aligned
    if strcmpi(jkvt(t).trialType, 'sp') % a successful pull trial
        if ~isempty(jkvt(t).rStartToPull) % if there's rStartToPull
            ss(t).timeAlign = jkvt(t).rStartToPull; % align to reach start to pull
            ss(t).evtAlign  = 'rStart';
            if ~isempty(jkvt(t).rStopToPull)
                ss(t).rEnd = jkvt(t).rStopToPull; % reachStop
            elseif ~isempty(jkvt(t).pullStops)
                ss(t).rEnd = jkvt(t).pullStops;   % instead use pullStop
            else
                ss(t).rEnd = jkvt(t).trEnd;
            end
        elseif ~isempty(jkvt(t).pullStarts) && ~isempty(jkvt(t).pullStops)
            ss(t).timeAlign = jkvt(t).pullStarts;
            ss(t).rEnd = jkvt(t).trJsReady+jkvt(t).movKins.pullStop; % reachStop
            ss(t).evtAlign  = 'pStart';
        else
            ss(t).timeAlign = jkvt(t).trJsReady; % if no reachStart detected, just align to the joystick ready
            ss(t).evtAlign  = 'trJsReady';
            ss(t).rEnd = jkvt(t).trEnd;
        end
    else % not a success trial
        if ~isempty(jkvt(t).hTrjRstart) && ~isempty(jkvt(t).hTrjRstop) % if there's a detected reachStart align to that
            ss(t).unsuccessReachLogic = true; % to mark unseccessful reach attempts
            ss(t).timeAlign = jkvt(t).vFrameTime(jkvt(t).hTrjRstart(end));
            ss(t).evtAlign  = 'rStart';
            ss(t).rEnd = jkvt(t).vFrameTime(min(jkvt(t).hTrjRstop(end),length(jkvt(t).vFrameTime))); % reachStop
        else
            ss(t).timeAlign = jkvt(t).trJsReady; % if no reachStart detected, just align to the joystick ready
            ss(t).evtAlign  = 'trJsReady';
            ss(t).rEnd = jkvt(t).trEnd;
        end
    end

    % get spike time bins and binned spike count matrices
    ss(t).spkTimeBins = ss(t).timeAlign + (spkBin(1):binSize:spkBin(end-1));
    %spkTime1msBins = ss(t).timeAlign + (spkBin(1):spkBin(end-1));

    if exist('spkTimesCellCTX')==1
        ss(t).unitTimeBCtx = psthBINcellPerTrial(spkTimesCellCTX, ss(t).timeAlign, binSize, [abs(spkBin(1)) abs(spkBin(end))]); % binned spikeCounts aligned to this trial
    end
    if exist('spkTimesCellSTR')==1
        ss(t).unitTimeBStr = psthBINcellPerTrial(spkTimesCellSTR, ss(t).timeAlign, binSize, [abs(spkBin(1)) abs(spkBin(end))]); % binned spikeCounts aligned to this trial
    end
    if exist('spkTimesCellCg')==1
        ss(t).unitTimeBCg = psthBINcellPerTrial(spkTimesCellCg, ss(t).timeAlign, binSize, [abs(spkBin(1)) abs(spkBin(end))]); % binned spikeCounts aligned to this trial
    end

    % align to laserOn, pLaserOn, or reachPrep (2s reach preparatory period)
    if ~isnan(jkvt(t).stimLaserOn) % stim trial
        if trI.spI(t)
            ss(t).rStartRelToStim = jkvt(t).rStartToPull-jkvt(t).stimLaserOn; 
        else
            ss(t).rStartRelToStim = jkvt(t).stimLaserOff-jkvt(t).stimLaserOn; 
        end
        
        if exist('spkTimesCellCTX')==1
            ss(t).utbCtxStimAlign = psthBINcellPerTrial(spkTimesCellCTX, jkvt(t).stimLaserOn, binSize, [1000 4000]); % binned spikeCounts aligned to this trial
            ss(t).utbCtxStimAlign1msBin = psthBINcellPerTrial1msBin(spkTimesCellCTX, jkvt(t).stimLaserOn, 1, [1000 4000]); % binned spikeCounts aligned to this trial
        end
        if exist('spkTimesCellSTR')==1
            ss(t).utbStrStimAlign = psthBINcellPerTrial(spkTimesCellSTR, jkvt(t).stimLaserOn, binSize, [1000 4000]); % binned spikeCounts aligned to this trial
            ss(t).utbStrStimAlign1msBin = psthBINcellPerTrial1msBin(spkTimesCellSTR, jkvt(t).stimLaserOn, 1, [1000 4000]); % binned spikeCounts aligned to this trial
        end
        if exist('spkTimesCellCg')==1
            ss(t).utbCgStimAlign = psthBINcellPerTrial(spkTimesCellCg, jkvt(t).stimLaserOn, binSize, [1000 4000]); % binned spikeCounts aligned to this trial
            ss(t).utbCgStimAlign1msBin = psthBINcellPerTrial1msBin(spkTimesCellCg, jkvt(t).stimLaserOn, 1, [1000 4000]); % binned spikeCounts aligned to this trial
        end
    else % control trial
        if isempty(jkvt(t).rStartToPull)
            takePstimTrI = true;
        elseif jkvt(t).rStartToPull-jkvt(t).pLaserOn>3000 % ensure to exclude trials where reach initiated
            takePstimTrI = true;
        end
        takePstimTrI = takePstimTrI && ~trI.toI(t);

        if takePstimTrI
            if exist('spkTimesCellCTX')==1
                ss(t).utbCtxPstimAlign = psthBINcellPerTrial(spkTimesCellCTX, jkvt(t).pLaserOn, binSize, [1000 4000]); % binned spikeCounts aligned to this trial
            end
            if exist('spkTimesCellSTR')==1
                ss(t).utbStrPstimAlign = psthBINcellPerTrial(spkTimesCellSTR, jkvt(t).pLaserOn, binSize, [1000 4000]); % binned spikeCounts aligned to this trial
            end
            if exist('spkTimesCellCg')==1
                ss(t).utbCgPstimAlign = psthBINcellPerTrial(spkTimesCellCg, jkvt(t).pLaserOn, binSize, [1000 4000]); % binned spikeCounts aligned to this trial
            end
        end
    end


    %% stim trial info
    if isfield(jkvt,'stimLaserOn')
        if ~isnan(jkvt(t).stimLaserOn) && ~isnan(jkvt(t).stimLaserOff)
            ss(t).tLaserStart = jkvt(t).stimLaserOn;
            ss(t).tLaserStop = jkvt(t).stimLaserOff;
            ss(t).spkTimeBlaserI = ss(t).tLaserStart<=ss(t).spkTimeBins & ss(t).spkTimeBins<=ss(t).tLaserStop;
        end
    end
    fprintf('processed trial %d\n', t) % report trial progression
end

%% z-score normalization
ctx_base_mean = nanmean(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBCtx}, 'UniformOutput', false))), 2); 
ctx_base_std = nanstd(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBCtx}, 'UniformOutput', false))), 0, 2); 

str_base_mean = nanmean(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBStr}, 'UniformOutput', false))), 2); 
str_base_std = nanstd(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBStr}, 'UniformOutput', false))), 0, 2); 

%cg_base_mean = nanmean(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBCg}, 'UniformOutput', false))), 2); 
%cg_base_std = nanstd(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBCg}, 'UniformOutput', false))), 0, 2); 

for t = 1:length(ss)
    if ~isempty(ss(t).utbCtxStimAlign)
        ss(t).utbCtxStimAlignZ = (full(ss(t).utbCtxStimAlign)-repmat(ctx_base_mean, 1, size(ss(t).utbCtxStimAlign, 2)))...
            ./repmat(ctx_base_std, 1, size(ss(t).utbCtxStimAlign, 2)); 
    end
    if ~isempty(ss(t).utbStrStimAlign)
        ss(t).utbStrStimAlignZ = (full(ss(t).utbStrStimAlign)-repmat(str_base_mean, 1, size(ss(t).utbStrStimAlign, 2)))...
            ./repmat(str_base_std, 1, size(ss(t).utbStrStimAlign, 2)); 
    end
%     if ~isempty(ss(t).utbCgStimAlign)   
%         ss(t).utbCgStimAlignZ = (full(ss(t).utbCgStimAlign)-repmat(cg_base_mean, 1, size(ss(t).utbCgStimAlign, 2)))...
%             ./repmat(cg_base_std, 1, size(ss(t).utbCgStimAlign, 2)); 
%     end
end


% transfer block information from jkvt to ss
[ss(:).blNumber] = deal(jkvt(:).blNumber);
[ss(:).blType] = deal(jkvt(:).blType);
[ss(:).blShiftLogic] = deal(jkvt(:).blShiftLogic);

% variables to be saved to build trial-type classifiers
ss_trialType = {ss.trialType};
ss_blNumber = {ss.blNumber};
jkvt_rewarded = {jkvt.rewarded};

saveName = filePath(end-19:end-9);

% % Save relevant variable for classification analysis
% if isfield(ss, 'unitTimeBCtx')
%     ss_unitTimeBCtx = {ss.unitTimeBCtx};
%     save(fullfile(filePath, strcat('unitTimeBCtx', '_', saveName)), 'ss_unitTimeBCtx')
% end
% 
% if isfield(ss, 'unitTimeBStr')
%     ss_unitTimeBStr = {ss.unitTimeBStr};
%     save(fullfile(filePath, strcat('unitTimeBStr', '_', saveName)), 'ss_unitTimeBStr')
% end
% 
% if isfield(ss, 'unitTimeBCg')
%     ss_unitTimeBCg = {ss.unitTimeBCg};
%     save(fullfile(filePath, strcat('unitTimeBCg', '_', saveName)), 'ss_unitTimeBCg')
% end
% 
% if isfield(ss, 'trialType')
%     ss_trialType = {ss.trialType};
%     save(fullfile(filePath, strcat('trialType', '_', saveName)), 'ss_trialType')
% end
% 
% if isfield(ss, 'blNumber')
%     ss_blNumbs = {ss.blNumber};
%     save(fullfile(filePath, strcat('blockNums', '_', saveName)), 'ss_blNumbs')
% end

save(fullfile(filePath, strcat('js2p0_tbytSpkHandJsTrjBin_50ms_stimParse_',saveName)),'ss','jkvt','spkTimesCell*','depth*')

% save(fullfile(filePath,strcat('js2p0_tbytSpkHandJsTrjBin_',saveName)),'depthCtx','depthStr','-append')
% save(fullfile('/Users/parkj/Dropbox (HHMI)/j2p0_dataShare/js2p0_tbytSpkHandJsTrjBin_WR40_081919.mat'),'depthCtx','depthStr','-append')
%load(fullfile(filePath,strcat('js2p0_tbytSpkHandJsTrjBin_',saveName)),'ss','jkvt','trI','spkTimesCell')


%% %%%%%%%%%%%%%%%%%%%
%%% Helper function %%
%%%%%%%%%%%%%%%%%%%%%%
% parse trial/block information
    function jkvt = jkvtBlockParse(jkvt)
        tqd = diff([jkvt(1).pull_torque, jkvt(:).pull_torque]'); % torque change
        p1d = diff([jkvt(1).reachP1, jkvt(:).reachP1]'); % position 1 change

        blNumb = 1;  % block number
        blType = []; % block type
        % detect and parse block shifts across trials
        for tt = 1:size(jkvt,2)
            % first trial
            if tt==1
                jkvt(tt).blNumber = blNumb;
                jkvt(tt).blType = 'first';
                jkvt(tt).blShiftLogic = true;
            end
            % second trial and onward
            if tt>=2
                if tqd(tt)==0 && p1d(tt)==0 % if there's no change
                    jkvt(tt).blNumber = blNumb; % remains the same
                    jkvt(tt).blShiftLogic = false; % not shifted
                    jkvt(tt).blType = []; % block shift type
                else % if there's any change, classify the block shift type
                    blNumb = blNumb+1;
                    jkvt(tt).blNumber = blNumb; % remains the same
                    jkvt(tt).blShiftLogic = true; % shifted
                    % parse torque shift
                    if tqd(tt)>0
                        tqi = 'tqUp';
                    elseif tqd(tt)<0
                        tqi = 'tqDown';
                    elseif tqd(tt)==0
                        tqi = 'tqSame';
                    end
                    % parse position shift
                    if p1d(tt)>0
                        p1i = 'rightward';
                    elseif p1d(tt)<0
                        p1i = 'leftward';
                    elseif p1d(tt)==0
                        if jkvt(tt).reachP1==min([jkvt(:).reachP1])
                            p1i = 'sameleft';
                        elseif jkvt(tt).reachP1==max([jkvt(:).reachP1])
                            p1i = 'sameright';
                        end
                    end
                    jkvt(tt).blType = strcat(tqi,p1i); % block shift type
                end
            end
        end
    end

% get interpolated hand trajectory
    function [intV,xq] = interpsm(x,v)
        xq = x(1):x(end); % new timescale
        inthTrjF = @(a) interp1(x,a,xq); % interpolation function
        vC = mat2cell(v,[ 1 1 1 ], size(v,2)); % convert to cell
        intVC = cellfun(@(a) inthTrjF(a), vC, 'un', 0); % interpolated hTrj cell

        lengthIntVC = unique(cellfun(@length, intVC));

        % to use sg filter get sgfiltFramelen
        if lengthIntVC >= 201
            sgfiltFramelen = 101;
        elseif lengthIntVC < 201
            if mod(lengthIntVC,2)==0
                sgfiltFramelen = lengthIntVC-1; % the frame length for sg filter needs to be an odd number
            else
                sgfiltFramelen = lengthIntVC;
            end
        end
        intVCf = cellfun(@(a) sgolayfilt(a,3,sgfiltFramelen), intVC, 'un', 0); % filtered hTjr cell
        intV = cell2mat(intVCf);
    end

% compute reach angle
    function reachAngle = computeReachAngle( hTrjC, jsXYpos )
        %Angle is calculated relative to the straight line projected from the
        % initial hand position
        distToJs = cellfun(@(a) sum((a-repmat(jsXYpos,1,size(a,2))).^2,1), hTrjC, 'un', 0);  % compute the distance from the jsXY
        for j = 1:length(distToJs)
            [~,tmpMinI] = min(distToJs{j});
            u=hTrjC{j}(:,tmpMinI)-hTrjC{j}(:,1); % reach vector
            %v=[hTrjC{j}(1,tmpMinI);hTrjC{j}(2,1)]-hTrjC{j}(:,1); % reference vector
            v=[hTrjC{j}(1,1);hTrjC{j}(2,tmpMinI)]-hTrjC{j}(:,1); % reference vector
            reachAngle(j) = angleTwoVectors(u,v); % get the angle between reach and reference vectors
        end
    end

    function [targinref] = targetintoreftime(A,timeA,timeRef, varargin)
        %This function takes a times series 'A' across a time window ('timeA') and
        % plug it into the new time window 'timeRef'.
        % 1-ms time resolution is assumed.
        p = targetintoreftime( A,timeA,timeRef, varargin );

        Aval = A(timeRef(1)<=timeA&timeA<=timeRef(end)); % valid portion of A to be placed in timeRef

        if p.Results.fillin == 0
            targinref = zeros(size(A,1),size(timeRef,2));
        elseif isnan(p.Results.fillin)
            targinref = nan(size(A,1),size(timeRef,2));
        end

        tRef1 = find(timeRef>=timeA(1),1,'first'); % 1st point
        if ~isempty(tRef1) && ~isempty(Aval)
            targinref(tRef1:tRef1+min(length(targinref),length(Aval))-1) = Aval;
        end

        function p = targetintoreftime( A,timeA,timeRef, vargs )
            % parse input, and extract name-value pairs
            default_fillin = 0; % by default do not re-read the raw bin file, if done already

            p = inputParser; % create parser object
            addRequired(p,'A');
            addRequired(p,'timeA');
            addRequired(p,'timeRef');
            addParameter(p,'fillin',default_fillin);

            parse(p,A,timeA,timeRef,vargs{:})
        end

    end

end
