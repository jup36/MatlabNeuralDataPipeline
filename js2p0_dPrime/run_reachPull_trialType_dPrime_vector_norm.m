%% 'glm_reachPull_dPrime.m' computes sensitivity index (d') that indicates
% the degree to which a neuron discriminates its firing rate as a function
% of a task variable (e.g. left vs. right, low vs. high load, left-low vs. rest).
% d' is computed as as the difference between two averaged responses, divided by their standard deviation.
% Statistical significance can be determined based on the confidence
% interval obtained by shuffling trial Ids.

function run_reachPull_trialType_dPrime_vector_norm(filePath)
%% 0. Load data
cd(filePath)
wr = strfind(filePath, 'WR');
saveName = filePath(wr:wr+10); 

spkDir = dir('**/*binSpkCountSTRCTX*.mat');
spkDir_Cg = dir('**/*binSpkCountCg*.mat');

%load(fullfile('D:\Junchol_Data\JS2p0\collectData','a2dColorMap.mat'), 'colormap2D') % 2d colorMap for the scatter plot

if ~isempty(spkDir) || ~isempty(spkDir_Cg)
    
    if ~isempty(spkDir)
        load(fullfile(spkDir(1).folder, spkDir(1).name),'spkTimesCell','jkvt')
    end
    
    if ~isempty(spkDir_Cg)
        spkTimesCellCg = load(fullfile(spkDir_Cg(1).folder, spkDir_Cg(1).name),'spkTimesCell');
        spkTimesCellCg = spkTimesCellCg.spkTimesCell; 
        if ~exist('jkvt', 'var')
            load(fullfile(spkDir_Cg(1).folder, spkDir_Cg(1).name),'jkvt');
        end
    end
       
    %% 1. Task parameters
    % task parameter
    p.DT = 20; % 20 ms (width of timeBin)
    p.filter_sigma = 100; % ms
    p.window = 1000; %5000; % ms
    p.cut = 4 * p.filter_sigma / p.DT;
    r2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
    
    if ~isempty(spkDir)
        p.isStr = cell2mat(spkTimesCell(5,:));
    end 
    
    % code parameter
    PLOT = false;
    CALC_PRM = true;
    
    %% 2. Behavioral data
    % detect events and continuous joystick pull trajectory
    %rwdTimes = cell2mat(cellfun(@(a) a(~isnan(a)), {jkvt(:).rewardT}, 'un', 0));
    trEnds = cell2mat(cellfun(@(a) a(~isnan(a)), {jkvt(:).trEnd}, 'un', 0));
    sessionDur = trEnds(end)+10*1000; % cut 10 seconds after the last reward
    session_range = [1 sessionDur];
    clip = @(t) t(session_range(1) <= t & t <= session_range(2)) - session_range(1);
    time_bin = (0:p.DT:sessionDur)';
    n_bin = length(time_bin) - 1;
    
    % trial indices
    [posTqC,posTqTypes] = posTrqCombinations( jkvt );
    rStartI = cellfun(@(c) ~isempty(c), {jkvt(:).rStartToPull})';
    pStartI = cellfun(@(c) ~isempty(c), {jkvt(:).pullStarts})';
    leI = cell2mat(cellfun(@(a) contains(a,'p1'), posTqC, 'un', 0 )) & pStartI;
    riI = cell2mat(cellfun(@(a) contains(a,'p2'), posTqC, 'un', 0 )) & pStartI;
    loI = cell2mat(cellfun(@(a) contains(a,'t1'), posTqC, 'un', 0 )) & pStartI;
    hiI = cell2mat(cellfun(@(a) contains(a,'t2'), posTqC, 'un', 0 )) & pStartI;
    leloI = leI & loI;
    lehiI = leI & hiI;
    riloI = riI & loI;
    rihiI = riI & hiI;
    rwI = [jkvt(:).rewarded]';
    
    % if reach Start is missing, put estimated values based on the pullStarts with (-200 ms offset)
    if sum(~rStartI & pStartI)>=1
        tempMissingRstarts = find(~rStartI & pStartI);
        for tt = 1:length(tempMissingRstarts)
            jkvt(tempMissingRstarts(tt)).rStartToPull = jkvt(tempMissingRstarts(tt)).pullStarts - 200;
        end
    end
    rStartI = cellfun(@(c) ~isempty(c), {jkvt(:).rStartToPull})';
    assert(unique(pStartI(rStartI))==true)
    
    % To save task variables
    tV.sessionDur = sessionDur;
    tV.clip = clip;
    tV.time_bin = time_bin;
    tV.evtRstart = [jkvt(:).rStartToPull];
    tV.evtRwd = [jkvt(:).trEnd]+1000;
    tV.posTqC = posTqC;
    tV.posTqTypes = posTqTypes;
    tV.rStartI = rStartI;
    tV.pStartI = pStartI;
    tV.leI = leI;
    tV.riI = riI;
    tV.loI = loI;
    tV.hiI = hiI;
    tV.rwI = rwI;
    tV.leloI = leloI;
    tV.lehiI = lehiI;
    tV.riloI = riloI;
    tV.rihiI = rihiI;
    
    %% compute d'
    % cortico-striatal
    if ~isempty(spkDir)  
        for i_cell = 1:size(spkTimesCell,2)
            spike_time = clip(spkTimesCell{1,i_cell});
            n_spike = length(spike_time);
            [spike_bin,~,spike_binT] = histcounts(spike_time, time_bin);
            spike_rate = sum(spike_bin) / (sessionDur/(1000)); % spikes/Sec
            spike_bin_conv = basis.normal_filter((spike_bin.*(1000/p.DT))', p.filter_sigma, p.DT); % p.filter_sigma = 100 ms
            
            if spike_rate>=0.5
                % align to the task event (reachStart to pull)
                [spike_rStart_count, time_tr_org] = alignToTaskEvent(tV.evtRstart, spike_time, p.DT, p.DT, p.window + 4 * p.filter_sigma);
                spike_rStart = spike_rStart_count'.*(1000./p.DT); % timeBin-by-trial, spike rate (Hz)
                in_t = 1 + p.cut:length(time_tr_org) - p.cut;
                time_tr = time_tr_org(in_t);
                
                spike_rStart_conv0 = basis.normal_filter(spike_rStart, p.filter_sigma, p.DT);
                spike_rStart_conv = spike_rStart_conv0(in_t,:);
                
                % d' trial types
                dPrmC{i_cell,1} = dPrime_trialType_vector_norm_stat(spike_rStart_conv, [tV.leloI(pStartI), tV.lehiI(pStartI), tV.riloI(pStartI), tV.rihiI(pStartI)], [-1 -1 1 1; -1 1 -1 1]);
            else
                dPrmC{i_cell,1} = [];  
            end
            fprintf('processed cell # %d\n', i_cell) % report unit progression
        end
        if ~exist(fullfile(filePath,'Matfiles'),'dir')
            mkdir(fullfile(filePath,'Matfiles'))
        end
        save(fullfile(filePath, 'Matfiles', strcat('glm_dPrime_vec_norm_',saveName)),'dPrmC','tV','p')
    end
    % Cg
    if ~isempty(spkDir_Cg)  
        for i_cell = 1:size(spkTimesCellCg,2)
            spike_time = clip(spkTimesCellCg{1,i_cell});
            n_spike = length(spike_time);
            [spike_bin,~,spike_binT] = histcounts(spike_time, time_bin);
            spike_rate = sum(spike_bin) / (sessionDur/(1000)); % spikes/Sec
            spike_bin_conv = basis.normal_filter((spike_bin.*(1000/p.DT))', p.filter_sigma, p.DT); % p.filter_sigma = 100 ms
            
            if spike_rate>=0.5
                % align to the task event (reachStart to pull)
                [spike_rStart_count, time_tr_org] = alignToTaskEvent(tV.evtRstart, spike_time, p.DT, p.DT, p.window + 4 * p.filter_sigma);
                spike_rStart = spike_rStart_count'.*(1000./p.DT); % timeBin-by-trial, spike rate (Hz)
                in_t = 1 + p.cut:length(time_tr_org) - p.cut;
                time_tr = time_tr_org(in_t);
                
                spike_rStart_conv0 = basis.normal_filter(spike_rStart, p.filter_sigma, p.DT);
                spike_rStart_conv = spike_rStart_conv0(in_t,:);
                
                % d' trial types
                dPrmC_Cg{i_cell,1} = dPrime_trialType_vector_norm_stat(spike_rStart_conv, [tV.leloI(pStartI), tV.lehiI(pStartI), tV.riloI(pStartI), tV.rihiI(pStartI)], [-1 -1 1 1; -1 1 -1 1]);
            else
                dPrmC_Cg{i_cell,1} = []; 
            end
            fprintf('processed cell # %d\n', i_cell) % report unit progression
        end
        if ~exist(fullfile(filePath,'Matfiles'),'dir')
            mkdir(fullfile(filePath,'Matfiles'))
        end
        save(fullfile(filePath, 'Matfiles', strcat('glm_dPrime_Cg_vec_norm_',saveName)),'dPrmC_Cg','tV','p')
    end 
end
end


%% individual unit psth aligned to a task event
% % sort trials
% tq = round([jkvt.pull_torque]./10); % torque
% tqT = tq(S.trI)'; tqT(:,2) = 1:length(S.trI); tqT(:,3) = S.trI;  sortByTq = sortrows(tqT,1);
% ps = [jkvt.reachP1]; % position
% psT = ps(S.trI)'; psT(:,2) = 1:length(S.trI); psT(:,3) = S.trI;  sortByPs = sortrows(psT,1);
% tqpsT(:,1) = tqT(:,1)+psT(:,1); tqpsT(:,2) = 1:length(S.trI); tqpsT(:,3) = S.trI; sortByTqPs = sortrows(tqpsT,1); % torque position combo
%
% %load('binSpkCountSTRCTXWR40_082019.mat', 'rStartToPull')
% %S = rStartToPull
% cellI_stc = 131;
% [cellI] = getUnitIdSfromSpkTimesCell(spkTimesCell, S, cellI_stc);
% %cellI = 11; % 102 79, 80, 35, 77, 85
% thisUnitSpkTimes = S.SpkTimes{cellI};
% individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTq, cellI, [3e3 2e3], [1e3 1e3]);
% individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByPs, cellI, [3e3 2e3], [1e3 1e3]);
% individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTqPs, cellI, [3e3 2e3], [1e3 1e3]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [binSpkMatOut,binEdgesOut] = alignToTaskEvent(evt, spikeTime, binSize, stepSize, oneWindow)
binEdgesOut = -oneWindow:binSize:oneWindow;
bin1msEdge = cellfun(@(a) max(a-oneWindow,0):a+oneWindow+1, num2cell(evt),'un', 0);
bin1msCount = cellfun(@(a) histcounts(spikeTime,a), bin1msEdge, 'un',0);
binSpkMatOut = bin1msSpkCountMat(cell2mat(bin1msCount'),binSize,stepSize);
end


function [unitIdSTC] = getUnitIdSpkTimesCell(spkTimesCellIn, sIn, unitIdS)
%unit Ids might not match between 'spkTimesCell' and the structure 'S' (e.g. S.unitTimeTrial)
% this function identifies the unit from 'S' in 'spkTimesCell' and outputs
% it as 'unitidSTC'

geomI = cell2mat(cellfun(@(a) isequal(a,sIn.geometry{unitIdS}), spkTimesCellIn(4,:), 'un', 0));
wfI = cell2mat(cellfun(@(a) isequal(a,sIn.meanWF{unitIdS}), spkTimesCellIn(6,:), 'un', 0));

unitIdSTC = find(geomI & wfI);

end


function [unitIdS] = getUnitIdSfromSpkTimesCell(spkTimesCellIn, sIn, unitIdC)
%unit Ids might not match between 'spkTimesCell' and the structure 'S' (e.g. S.unitTimeTrial)
% this function identifies the unit from 'S' in 'spkTimesCell' and outputs
% it as 'unitidSTC'

geomI = cell2mat(cellfun(@(a) isequal(a,spkTimesCellIn{4,unitIdC}), sIn.geometry, 'un', 0));
wfI = cell2mat(cellfun(@(a) isequal(a,spkTimesCellIn{6,unitIdC}), sIn.meanWF, 'un', 0));

unitIdS = find(geomI & wfI);
end


function [winModelRate] = alignGlmOutToTaskEvent(evt, time_bin, modelRate, binSize, oneWindow)
windowBins = [floor(-oneWindow/binSize)+1,ceil(oneWindow/binSize)];
binEvt = find(histcounts(evt,time_bin)'); % identify bins per event
windowBound = arrayfun(@(a) a+windowBins,binEvt,'un',0);
for t = 1:length(windowBound)
    if windowBound{t}(1)>=1 && windowBound{t}(2)<=length(modelRate)
        winModelRateC{t} = modelRate(windowBound{t}(1):windowBound{t}(2));
    else
        winModelRateC{t} = NaN(sum(abs(windowBins))+1,1);
    end
end
%winModelRateC = cellfun(@(a) modelRate(a(1):a(2)), windowBound, 'un',0);
winModelRate = cell2mat(winModelRateC);
end


function [trainC,testC,type] = crossValTrainTestSetsEachTrialType(typeC, fold, varargin)
%This function takes trial type information (e.g. four types (two position-by-two torques)
% and generates N-fold cross-validation train and test sets. An additional
% logic can be input for further selection of trials as a varargin.
% use e.g.1: [trainC,testC] = crossValTrainTestSetsEachTrialType(posTqC, 5, pStartI)
% use e.g.2:

useAllTrainTrs = true;

if nargin == 2
    addLogic = true(size(typeC,1),1);
elseif nargin == 3
    addLogic = varargin{1};
end

assert(length(addLogic)==length(typeC))

type = unique(typeC);
for ty = 1:length(type)
    typeS{1,ty} = find(cell2mat(cellfun(@(a) strcmpi(a,type{ty}) , typeC, 'un', 0)) & addLogic);
    typeS{1,ty} = randsample(typeS{1,ty}, length(typeS{1,ty})); % randomize the order
end

trialN = min(cellfun(@length,typeS));
testN = floor(trialN*(1/fold));
trainN = floor(trialN*((fold-1)/fold));

testC  = cell(length(type),fold); % test sets per type and fold
trainC = cell(length(type),fold); % train sets per type and fold

for ff = 1:fold
    tempTest = (ff-1)*testN+1:(ff-1)*testN+testN; % test trials of this fold (common across trial types)
    for tp = 1:length(type)
        if useAllTrainTrs
            tempTrain = find(~ismember(1:length(typeS{tp}),tempTest)); % sample train trials of this trial type of this fold
        else
            tempTrain = randsample(find(~ismember(1:length(typeS{tp}),tempTest)),trainN); % sample train trials of this trial type of this fold
        end
        testC{tp,ff}  = sort(typeS{tp}(tempTest));
        trainC{tp,ff} = sort(typeS{tp}(tempTrain));
    end
end
end


function [testDatC] = getTestDatCell(dataSet,indices)
for ii = 1:length(indices)
    testDatC{ii} = dataSet(:,indices(ii));
end
end


% compute the probability of being a left-target trial
function [prob1, post1, post2] = posteriorProb2(model1_lambda, model2_lambda, testDat)
post1 = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(testDat), model1_lambda, 'un', 0));
post2 = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(testDat), model2_lambda, 'un', 0));
prob1 = post1./(post1+post2);
end


function [excludeLogic] = excludeSegments(dat,pointsToExclude,oneWindow)
% this function takes a whole segment and excludes
%subsegment-windows specified by the points and window.
% dat = X;
% pointsToExclude = testPtsToExclude_bin;
% oneWindow = 2000/p.DT;

if size(dat,1)<size(dat,2)
    dat = dat';
end

windowsToExclude = arrayfun(@(a) a-oneWindow:a+oneWindow, pointsToExclude, 'un', 0);
excludeLogicC = cellfun(@(a) ismember(1:size(dat,1),a), windowsToExclude, 'un', 0)';
excludeLogic = sum(cell2mat(excludeLogicC))>0;

if size(excludeLogic,1)<size(excludeLogic,2)
    excludeLogic = excludeLogic';
end

end


function [r2_psth] = getR2psth(evtInMs, timeBin, rateFit, rateDat, binSize, oneWin)

rr2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
model_psth = func.alignGlmOutToTaskEvent(evtInMs, timeBin, rateFit, binSize, oneWin); % fitted rate aligned to reachStart
spike_psth = func.alignGlmOutToTaskEvent(evtInMs, timeBin, rateDat, binSize, oneWin); % fitted rate aligned to reachStart
r2_psth = max(rr2(reshape(spike_psth,[],1),reshape(model_psth,[],1)),0);

end


function dPrmRez = dPrime_trialType(spike_time_trial,ttI,twoDconvMat)
%spike_time_trial = spike_rStart_conv;
%ttI = [tV.leloI(pStartI), tV.lehiI(pStartI), tV.riloI(pStartI), tV.rihiI(pStartI)];
%twoDconvMat = [-1 -1 1 1; -1 1 -1 1]; % column(1-4): left-low, left-high, right-low, right-high

nTr = size(ttI,1);

%% trial-shuffled distribution of d'
for tt = 1:size(ttI,2)
    for rs = 1:1000 % 1000 resampling
        tmpTr = randsample(nTr,sum(ttI(:,tt))); % random sample trials
        % random sampled spike_time_trial mat
        tmpMat = spike_time_trial(:,tmpTr); % spike (random sample trials)
        tmpNoMat = spike_time_trial(:,~ismember(1:nTr,tmpTr)); % spike (rest)
        rs_dPrm_1d{tt}(:,rs) = (nanmean(tmpMat,2)-nanmean(tmpNoMat,2))./(nanstd(tmpMat,[],2)+nanstd(tmpNoMat,[],2));
        rs_dPrm_2d{tt}(:,:,rs) = twoDconvMat(:,tt)*rs_dPrm_1d{tt}(:,rs)'; % projection onto the 2d space
    end
end

rs_dPrm_2d_X = nansum(cell2mat(reshape(cellfun(@(a) squeeze(a(1,:,:)), rs_dPrm_2d, 'un', 0),1,1,[])),3); % X-dim (left-right)
rs_dPrm_2d_Y = nansum(cell2mat(reshape(cellfun(@(a) squeeze(a(2,:,:)), rs_dPrm_2d, 'un', 0),1,1,[])),3); % Y-dim (low-high)

rs_dPrm_2d_XY = [reshape(rs_dPrm_2d_X,[],1),reshape(rs_dPrm_2d_Y,[],1)];

% resampled/shuffled bivariate distribution mean, cov to get the pdf
dPrmRez.mu_shuff_2d_XY = nanmean(rs_dPrm_2d_XY,1);
dPrmRez.cov_shuff_2d_XY = cov(rs_dPrm_2d_XY);

% cdf
[Xs,Ys] = meshgrid(linspace(-3,3,600)',linspace(-3,3,600)');
XYs = [Xs(:),Ys(:)];
p_cdf0 = mvncdf(XYs, dPrmRez.mu_shuff_2d_XY, dPrmRez.cov_shuff_2d_XY);
%p_cdf = reshape(p_cdf0,600,600);

%% actual d' and p-values based on the cdf of the suffled distribution
for tt = 1:size(ttI,2) % trial Types (e.g. lelo, lehi, rilo, rihi)
    tmpMat = spike_time_trial(:,ttI(:,tt));
    tmpNoMat = spike_time_trial(:,~ttI(:,tt));
    dPrm_1d{tt} = (nanmean(tmpMat,2)-nanmean(tmpNoMat,2))./(nanstd(tmpMat,[],2)+nanstd(tmpNoMat,[],2)+eps);
    dPrm_2d{tt} = twoDconvMat(:,tt)*dPrm_1d{tt}';
end

dPrmRez.trj2d = sum(cell2mat(reshape(dPrm_2d,1,1,[])),3)'; % sum across all trial-type scores
dist2d = @(a) sqrt(a(:,1).^2+a(:,2).^2);
dPrmRez.trjDistZero = dist2d(dPrmRez.trj2d); % distance from zero

for j = 1:size(dPrmRez.trj2d,1)
    [~,minDistI] = min(dist2d(dPrmRez.trj2d(j,:)-XYs)); % locate each point of the 2d trajectory within the grid
    dPrmRez.dist2d_cdf_p(j,1) = min(p_cdf0(minDistI),1-p_cdf0(minDistI));
end

%dPrm_2d_distZero(dist2d_cdf_p>0.01,1)=NaN;

%[~,dPrm_sig_maxDistI] = max(dPrm_2d_distZero);

end


function dPrmRez = dPrime_trialType_vector_norm_stat(spike_time_trial,ttI,twoDconvMat)
%spike_time_trial = spike_rStart_conv; 
%ttI = [tV.leloI(pStartI), tV.lehiI(pStartI), tV.riloI(pStartI), tV.rihiI(pStartI)]; 
%twoDconvMat = [-1 -1 1 1; -1 1 -1 1]; % column(1-4): left-low, left-high, right-low, right-high

nTr = size(ttI,1);  
nTr_sample = round(nTr/4); 

distf = @(a, b) sqrt(a^2 + b^2); 

%% trial-shuffled distribution of d' 
for tt = 1:size(ttI,2)
    for rs = 1:1000 % 1000 resampling
        tmpTr = randsample(nTr, nTr_sample); % randsample(nTr,sum(ttI(:,tt))); % random sample trials 
        % random sampled spike_time_trial mat
        tmpMat = spike_time_trial(:,tmpTr); % spike (random sample trials)
        tmpNoMat = spike_time_trial(:,~ismember(1:nTr,tmpTr)); % spike (rest)
        rs_dPrm_1d{tt}(:,rs) = (nanmean(tmpMat,2)-nanmean(tmpNoMat,2))./(nanstd(tmpMat,[],2)+nanstd(tmpNoMat,[],2)); 
        rs_dPrm_2d{tt}(:,:,rs) = twoDconvMat(:,tt)*rs_dPrm_1d{tt}(:,rs)'; % projection onto the 2d space
    end
end

rs_dPrm_2d_X = nansum(cell2mat(reshape(cellfun(@(a) squeeze(a(1,:,:)), rs_dPrm_2d, 'un', 0),1,1,[])),3); % X-dim (left-right)
rs_dPrm_2d_Y = nansum(cell2mat(reshape(cellfun(@(a) squeeze(a(2,:,:)), rs_dPrm_2d, 'un', 0),1,1,[])),3); % Y-dim (low-high)

max_dist = nan(1000, 1);  
for j = 1:1000 % 1000 resampling
    max_dist(j) = max(cell2mat(arrayfun(@(a, b) distf(a, b), rs_dPrm_2d_X(:,j), rs_dPrm_2d_Y(:,j), 'un', 0))); 
end

%% actual d' and p-values based on the cdf of the suffled distribution
for tt = 1:size(ttI,2) % trial Types (e.g. lelo, lehi, rilo, rihi)
    tmpMat = spike_time_trial(:,ttI(:,tt)); 
    tmpNoMat = spike_time_trial(:,~ttI(:,tt)); 
    dPrm_1d{tt} = (nanmean(tmpMat,2)-nanmean(tmpNoMat,2))./(nanstd(tmpMat,[],2)+nanstd(tmpNoMat,[],2)+eps); 
    dPrm_2d{tt} = twoDconvMat(:,tt)*dPrm_1d{tt}'; 
end

dPrmRez.trj2d = sum(cell2mat(reshape(dPrm_2d,1,1,[])),3)'; % sum across all trial-type scores

dist2d = @(a) sqrt(a(:,1).^2+a(:,2).^2); 
dPrmRez.trjDistZero = dist2d(dPrmRez.trj2d); % distance from zero
dPrmRez.trjDistZero_sigI = dPrmRez.trjDistZero > mean(max_dist) + 2*std(max_dist); 

end
