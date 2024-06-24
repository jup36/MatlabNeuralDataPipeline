
filePath = '/Volumes/Extreme SSD/js2p0/gpfa_model';
load(fullfile(filePath, 'simDataGpfa.mat'), 'trialgroup')

%% preprocessing
trialgroup_spec_C = cell(4, 1);
for i = 1:length(trialgroup.spec)
    trialgroup_spec_C{1}{end+1}= bin1msSpkCountMat(max(trialgroup.spec(i).psth(:, :, 1), 0), 50, 50); % 50-ms bin
    trialgroup_spec_C{2}{end+1}= bin1msSpkCountMat(max(trialgroup.spec(i).psth(:, :, 2), 0), 50, 50);
    trialgroup_spec_C{3}{end+1}= bin1msSpkCountMat(max(trialgroup.spec(i).psth(:, :, 3), 0), 50, 50);
    trialgroup_spec_C{4}{end+1}= bin1msSpkCountMat(max(trialgroup.spec(i).psth(:, :, 4), 0), 50, 50);
end

trialgroup_selt_C = cell(4, 1);
for i = 1:length(trialgroup.selt)
    trialgroup_selt_C{1}{end+1}= bin1msSpkCountMat(max(trialgroup.selt(i).psth(:, :, 1), 0), 50, 50); % 50-ms bin
    trialgroup_selt_C{2}{end+1}= bin1msSpkCountMat(max(trialgroup.selt(i).psth(:, :, 2), 0), 50, 50);
    trialgroup_selt_C{3}{end+1}= bin1msSpkCountMat(max(trialgroup.selt(i).psth(:, :, 3), 0), 50, 50);
    trialgroup_selt_C{4}{end+1}= bin1msSpkCountMat(max(trialgroup.selt(i).psth(:, :, 4), 0), 50, 50);
end

% get datSpec
count = 0;
for j = 1:4
    for jj = 1:10
        count = count + 1;
        datSpec(count).trialId = count;
        datSpec(count).spikes = max(trialgroup_spec_C{j}{jj}, 0);
    end
end

% get datSel
count = 0;
for j = 1:4
    for jj = 1:10
        count = count + 1;
        datSelt(count).trialId = count;
        datSelt(count).spikes = max(trialgroup_selt_C{j}{jj}, 0);
    end
end

%% runGPFA
% specification model
saveNameSpec = fullfile(filePath, 'gpfa_spec');
if exist(saveNameSpec, "dir")~=7
    mkdir(saveNameSpec)
end

gpfaRezSpec = neuralTrajAlreadyBinned(saveNameSpec, datSpec, 'method', 'gpfa', ...
    'xDim', 5, 'kernSDList', 5, 'binWidth', 1); %

% gpfa postprocess
gpfaRezSpec.method = 'gpfa';
% Orthonormalize neural trajectories
[gpfaRezSpec.estParamsOrtho, gpfaRezSpec.seqTrainOrtho] = postprocess(gpfaRezSpec, 'kernSD', 5);

gpfaRezSpec_tt1 = projYtoGetXsm(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSpec(1:10));  % trial type 1
gpfaRezSpec_tt2 = projYtoGetXsm(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSpec(11:20)); % trial type 2
gpfaRezSpec_tt3 = projYtoGetXsm(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSpec(21:30)); % trial type 3
gpfaRezSpec_tt4 = projYtoGetXsm(fullfile(saveNameSpec, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSpec(31:40)); % trial type 4

plot3D_seqCells_color({gpfaRezSpec_tt1, gpfaRezSpec_tt2, gpfaRezSpec_tt3, gpfaRezSpec_tt4}, {[0 0 1], [0 1 1], [1 0 1], [1 0 0]}, 10);

%% selection model
saveNameSelt = fullfile(filePath, 'gpfa_selt');
if exist(saveNameSelt, "dir")~=7
    mkdir(saveNameSelt)
end

gpfaRezSelt = neuralTrajAlreadyBinned(saveNameSelt, datSelt, 'method', 'gpfa', ...
    'xDim', 5, 'kernSDList', 5, 'binWidth', 1); %

% gpfa postprocess
gpfaRezSelt.method = 'gpfa';
% Orthonormalize neural trajectories
[gpfaRezSelt.estParamsOrtho, gpfaRezSelt.seqTrainOrtho] = postprocess(gpfaRezSelt, 'kernSD', 5);

gpfaRezSelt_tt1 = projYtoGetXsm(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSelt(1:10));  % trial type 1
gpfaRezSelt_tt2 = projYtoGetXsm(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSelt(11:20)); % trial type 2
gpfaRezSelt_tt3 = projYtoGetXsm(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSelt(21:30)); % trial type 3
gpfaRezSelt_tt4 = projYtoGetXsm(fullfile(saveNameSelt, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datSelt(31:40)); % trial type 4

plot3D_seqCells_color({gpfaRezSelt_tt1, gpfaRezSelt_tt2, gpfaRezSelt_tt3, gpfaRezSelt_tt4}, {[0 0 1], [0 1 1], [1 0 1], [1 0 0]}, 10);


%% selection model pca
saveNameSeltPca = fullfile(filePath, 'pca_selt');
if exist(saveNameSeltPca, "dir")~=7
    mkdir(saveNameSeltPca)
end

pcaRezSelt = neuralTrajAlreadyBinned(saveNameSeltPca, datSelt, 'method', 'pca', ...
    'xDim', 5, 'kernSDList', 5, 'binWidth', 1); %

pcaXpost1 = zeros(3, 50, 10); 
count = 0; 
for j = 1:10
    count = count+1; 
    pcaXpost1(:, :, count) = pcaRezSelt.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost2 = zeros(3, 50, 10); 
count = 0; 
for j = 11:20
    count = count+1;
    pcaXpost2(:, :, count) = pcaRezSelt.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost3 = zeros(3, 50, 10); 
count = 0; 
for j = 21:30
    count = count+1;
    pcaXpost3(:, :, count) = pcaRezSelt.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpost4 = zeros(3, 50, 10); 
count = 0; 
for j = 31:40
    count = count+1;
    pcaXpost4(:, :, count) = pcaRezSelt.kern.seqTrain(j).xpost(1:3, :); 
end

pcaXpostC = {pcaXpost1, pcaXpost2, pcaXpost3, pcaXpost4}; 

plot3D_pcaSeqCells_color(pcaXpostC, {[0 0 1], [0 1 1], [1 0 1], [1 0 0]}, 10);


{pcaRezSelt.kern.seqTrain.xpost}



plot3D_seqCells_color({gpfaRezSelt_tt1, gpfaRezSelt_tt2, gpfaRezSelt_tt3, gpfaRezSelt_tt4}, {[0 0 1], [0 1 1], [1 0 1], [1 0 0]}, 10);




for f = 1:length(filePath)


    % control trials
    for t = 1:length(ctrlTrId)
        datCtx(t).trialId = ctrlTrId(t);
        datCtx(t).spikes = full(ss(ctrlTrId(t)).unitTimeBCtx);
        datStr(t).trialId = ctrlTrId(t);
        datStr(t).spikes = full(ss(ctrlTrId(t)).unitTimeBStr);
        if isfield(ss, 'unitTimeBCg')
            datCg(t).trialId = ctrlTrId(t);
            datCg(t).spikes = full(ss(ctrlTrId(t)).unitTimeBCg);
        end

        % p-stim-align
        datCtx_pstimAlign(t).trialId = ctrlTrId(t);
        datCtx_pstimAlign(t).spikes = full(ss(ctrlTrId(t)).utbCtxPstimAlign);
        datStr_pstimAlign(t).trialId = ctrlTrId(t);
        datStr_pstimAlign(t).spikes = full(ss(ctrlTrId(t)).utbStrPstimAlign);
        if isfield(ss, 'utbCgPstimAlign')
            datCg_pstimAlign(t).trialId = ctrlTrId(t);
            datCg_pstimAlign(t).spikes = full(ss(ctrlTrId(t)).utbCgPstimAlign);
        end
    end

    % stim trials
    for t = 1:length(stimId)
        datCtx_stim(t).trialId = stimId(t);
        datCtx_stim(t).spikes = full(ss(stimId(t)).unitTimeBCtx);
        datStr_stim(t).trialId = stimId(t);
        datStr_stim(t).spikes = full(ss(stimId(t)).unitTimeBStr);
        if isfield(ss, 'unitTimeBCg')
            datCg_stim(t).trialId = stimId(t);
            datCg_stim(t).spikes = full(ss(stimId(t)).unitTimeBCg);
        end
        % stim trials stim-align
        datCtx_stimAlign(t).trialId = stimId(t);
        datCtx_stimAlign(t).spikes = full(ss(stimId(t)).utbCtxStimAlign);
        datStr_stimAlign(t).trialId = stimId(t);
        datStr_stimAlign(t).spikes = full(ss(stimId(t)).utbStrStimAlign);
        if isfield(ss, 'utbCgStimAlign')
            datCg_stimAlign(t).trialId = stimId(t);
            datCg_stimAlign(t).spikes = full(ss(stimId(t)).utbCgStimAlign);
        end
    end

    %             % successful stim trials
    %             for t = 1:length(stimSpRsId)
    %                 datCtx_stimRs(t).trialId = stimSpRsId(t);
    %                 datCtx_stimRs(t).spikes = full(ss(stimSpRsId(t)).unitTimeBCtx);
    %                 datStr_stimRs(t).trialId = stimSpRsId(t);
    %                 datStr_stimRs(t).spikes = full(ss(stimSpRsId(t)).unitTimeBStr);
    %             end

    % stimFullBt
    for t = 1:length(stimFullBtId)
        datCtx_stimFullBt(t).trialId = stimFullBtId(t);
        datCtx_stimFullBt(t).spikes = full(ss(stimFullBtId(t)).unitTimeBCtx);
        datStr_stimFullBt(t).trialId = stimFullBtId(t);
        datStr_stimFullBt(t).spikes = full(ss(stimFullBtId(t)).unitTimeBStr);
        if isfield(ss, 'unitTimeBCg')
            datCg_stimFullBt(t).trialId = stimId(t);
            datCg_stimFullBt(t).spikes = full(ss(stimFullBtId(t)).unitTimeBCg);
        end
    end

    %% run gpfa on cortex data
    saveNameCtx = fullfile(fileparts(filePath{1}), 'gpfa_ctx');
    if exist(saveNameCtx, "dir")~=7
        mkdir(saveNameCtx)
    end

    % get gpfaRez
    if exist(fullfile(saveNameCtx, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')), 'file')==2
        gpfaRezCtx = load(fullfile(saveNameCtx, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')));
    else
        gpfaRezCtx = neuralTrajAlreadyBinned(saveNameCtx, datCtx, 'method', 'gpfa', ...
            'xDim', 5, 'kernSDList', 5, 'binWidth', 1); % ensure to have the binWidth 1 to prevent further binning
    end

    % postprocess
    gpfaRezCtx.method = 'gpfa';
    % Orthonormalize neural trajectories
    [gpfaRezCtx.estParamsOrtho, gpfaRezCtx.seqTrainOrtho] = postprocess(gpfaRezCtx, 'kernSD', 5);

    %% run gpfa on str data
    saveNameStr = fullfile(fileparts(filePath{1}), 'gpfa_str');
    if exist(saveNameStr, "dir")~=7
        mkdir(saveNameStr)
    end

    % get gpfaRez
    if exist(fullfile(saveNameStr, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')), 'file')==2
        gpfaRezStr = load(fullfile(saveNameStr, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')));
    else
        gpfaRezStr = neuralTrajAlreadyBinned(saveNameStr, datStr, 'method', 'gpfa', ...
            'xDim', 5, 'kernSDList', 5, 'binWidth', 1); % ensure to have the binWidth 1 to prevent further binning
    end

    % postprocess
    gpfaRezStr.method = 'gpfa';
    % Orthonormalize neural trajectories
    [gpfaRezStr.estParamsOrtho, gpfaRezStr.seqTrainOrtho] = postprocess(gpfaRezStr, 'kernSD', 5);

    %% run gpfa on cg data
    if isfield(ss, 'unitTimeBCg')
        saveNameCg = fullfile(fileparts(filePath{1}), 'gpfa_cg');
        if exist(saveNameCg, "dir")~=7
            mkdir(saveNameCg)
        end

        % get gpfaRez
        if exist(fullfile(saveNameCg, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')), 'file')==2
            gpfaRezCg = load(fullfile(saveNameCg, strcat(sprintf('%s_xDim%02d', 'gpfa', 5), '.mat')));
        else
            gpfaRezCg = neuralTrajAlreadyBinned(saveNameCg, datCg, 'method', 'gpfa', ...
                'xDim', 5, 'kernSDList', 5, 'binWidth', 1); % ensure to have the binWidth 1 to prevent further binning
        end

        % postprocess
        gpfaRezCg.method = 'gpfa';
        % Orthonormalize neural trajectories
        [gpfaRezCg.estParamsOrtho, gpfaRezCg.seqTrainOrtho] = postprocess(gpfaRezCg, 'kernSD', 5);
    end

    %% projection of 'test' trials onto the gpfa space (CTX)
    gpfaRezCtx_stim = projYtoGetXsm(fullfile(saveNameCtx, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCtx_stim);
    gpfaRezCtx_stimFullBt = projYtoGetXsm(fullfile(saveNameCtx, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCtx_stimFullBt);
    %gpfaRezCtx_stimRs = projYtoGetXsm(fullfile(saveNameCtx, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCtx_stimRs);

    gpfaRezCtx_stimAlign = projYtoGetXsm(fullfile(saveNameCtx, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCtx_stimAlign);
    gpfaRezCtx_pstimAlign = projYtoGetXsm(fullfile(saveNameCtx, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCtx_pstimAlign);

    % Plot neural trajectories in 3D space
    %plot3D(gpfaRezCtx.seqTrainOrtho, 'xorth', 'dimsToPlot', 1:3);
    %plot3D_seqCells_color({gpfaRezCtx.seqTrainOrtho, gpfaRezCtx_stimFullBt}, {[0 0 1], [1 0 1]}, 4);
    %plot3D_seqCells_color({gpfaRezCtx.seqTrainOrtho, gpfaRezCtx_stimRs}, {[0 0 1], [1 0 1]}, 10);
    %plot3D_seqCells_color({gpfaRezCtx.seqTrainOrtho, gpfaRezCtx_pstimAlign, gpfaRezCtx_stimAlign}, {[0 0 1], [1 0 1], [1 0 0]}, 10);
    %[az, el] = view;
    % print(fullfile(saveNameCtx, 'gpfaRezCtx_repTrj_sp_pstim_stim'), '-dpdf', '-bestfit', '-vector')

    %plot3D_seqCells_color({gpfaRezCtx.seqTrainOrtho, gpfaRezCtx_stimFullBt, gpfaRezCtx_stimAlign}, {[0 0 1], [1 0 1], [1 0 0]}, 10);
    %view(az, el);
    % print(fullfile(saveNameCtx, 'gpfaRezCtx_repTrj_sp_stimFullBt_stim'), '-dpdf', '-bestfit', '-vector')

    [gpfaRezCtx_euclideanDistCum, gpfaRezCtx_euclideanDistMin] = quantify_trjDeviation_seqCells_color(gpfaRezCtx.seqTrainOrtho, ...
        {gpfaRezCtx.seqTrainOrtho, gpfaRezCtx_stimFullBt, gpfaRezCtx_pstimAlign, gpfaRezCtx_stimAlign}, ...
        {[0 0 1], [1 0 1], [0 1 1], [1 0 0]});
    %print(fullfile(saveNameCtx, 'gpfaRezCtx_euclideanDistCum'), '-dpdf', '-bestfit', '-vector')

    save(fullfile(saveNameCtx, 'gpfaRezCtx.mat'), 'gpfaRezCtx*');

    %% projection of 'test' trials onto the gpfa space (STR)
    gpfaRezStr_stim = projYtoGetXsm(fullfile(saveNameStr, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datStr_stim);
    gpfaRezStr_stimFullBt = projYtoGetXsm(fullfile(saveNameStr, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datStr_stimFullBt);
    %gpfaRezStr_stimRs = projYtoGetXsm(fullfile(saveNameStr, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datStr_stimRs);

    gpfaRezStr_stimAlign = projYtoGetXsm(fullfile(saveNameStr, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datStr_stimAlign);
    gpfaRezStr_pstimAlign = projYtoGetXsm(fullfile(saveNameStr, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datStr_pstimAlign);

    % Plot neural trajectories in 3D space
    %plot3D(gpfaRezStr.seqTrainOrtho, 'xorth', 'dimsToPlot', 1:3);
    %plot3D_seqCells_color({gpfaRezStr.seqTrainOrtho, gpfaRezStr_stimFullBt}, {[0 0 1], [1 0 1]}, 4);
    %plot3D_seqCells_color({gpfaRezStr.seqTrainOrtho, gpfaRezStr_stimRs}, {[0 0 1], [1 0 1]}, 10);
    %plot3D_seqCells_color({gpfaRezStr.seqTrainOrtho, gpfaRezStr_pstimAlign, gpfaRezStr_stimAlign}, {[0 0 1], [1 0 1], [1 0 0]}, 10);
    %plot3D_seqCells_color({gpfaRezStr.seqTrainOrtho, gpfaRezStr_stimFullBt, gpfaRezStr_stimAlign}, {[0 0 1], [1 0 1], [1 0 0]}, 10);

    [gpfaRezStr_euclideanDistCum, gpfaRezStr_euclideanDistMin] = quantify_trjDeviation_seqCells_color(gpfaRezStr.seqTrainOrtho, ...
        {gpfaRezStr.seqTrainOrtho, gpfaRezStr_stimFullBt, gpfaRezStr_pstimAlign, gpfaRezStr_stimAlign}, ...
        {[0 0 1], [1 0 1], [0 1 1], [1 0 0]});
    %print(fullfile(saveNameStr, 'gpfaRezStr_euclideanDistCum'), '-dpdf', '-bestfit', '-vector')

    save(fullfile(saveNameStr, 'gpfaRezStr.mat'), 'gpfaRezStr*');

    %% projection of 'test' trials onto the gpfa space (Cg)
    if isfield(ss, 'unitTimeBCg')
        gpfaRezCg_stim = projYtoGetXsm(fullfile(saveNameCg, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCg_stim);
        gpfaRezCg_stimFullBt = projYtoGetXsm(fullfile(saveNameCg, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCg_stimFullBt);
        %gpfaRezCg_stimRs = projYtoGetXsm(fullfile(saveNameCg, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCg_stimRs);

        gpfaRezCg_stimAlign = projYtoGetXsm(fullfile(saveNameCg, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCg_stimAlign);
        gpfaRezCg_pstimAlign = projYtoGetXsm(fullfile(saveNameCg, strcat(sprintf('%s_xDim%02d', 'gpfa', 5))), datCg_pstimAlign);

        % Plot neural trajectories in 3D space
        %plot3D(gpfaRezCg.seqTrainOrtho, 'xorth', 'dimsToPlot', 1:3);
        %plot3D_seqCells_color({gpfaRezCg.seqTrainOrtho, gpfaRezCg_stimFullBt}, {[0 0 1], [1 0 1]}, 4);
        %plot3D_seqCells_color({gpfaRezCg.seqTrainOrtho, gpfaRezCg_stimRs}, {[0 0 1], [1 0 1]}, 10);
        %plot3D_seqCells_color({gpfaRezCg.seqTrainOrtho, gpfaRezCg_pstimAlign, gpfaRezCg_stimAlign}, {[0 0 1], [1 0 1], [1 0 0]}, 10);
        %plot3D_seqCells_color({gpfaRezCg.seqTrainOrtho, gpfaRezCg_stimFullBt, gpfaRezCg_stimAlign}, {[0 0 1], [1 0 1], [1 0 0]}, 10);

        [gpfaRezCg_euclideanDistCum, gpfaRezCg_euclideanDistMin] = quantify_trjDeviation_seqCells_color(gpfaRezCg.seqTrainOrtho, ...
            {gpfaRezCg.seqTrainOrtho, gpfaRezCg_stimFullBt, gpfaRezCg_pstimAlign, gpfaRezCg_stimAlign}, ...
            {[0 0 1], [1 0 1], [0 1 1], [1 0 0]});
        %print(fullfile(saveNameCg, 'gpfaRezCg_euclideanDistCum'), '-dpdf', '-bestfit', '-vector')

        save(fullfile(saveNameCg, 'gpfaRezCg.mat'), 'gpfaRezCg*');
    end
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq = projYtoGetXsm(pathToGpfaRez, dat)
%This utility function handles projection of test trials (each trial dimension: NxT)
% onto the dimensionality-reduction subspace (e.g., gpfa).
% "exactInferenceWithLL.m" is the original code that handles the projection.
%       dif      = bsxfun(@minus, [seq(nList).y], estParams.d); % yDim x sum(T)
%       term1Mat = reshape(CRinv * dif, xDim*T, []); % (xDim*T) x length(nList)

% pathToGpfaRez = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/gpfa_ctx/gpfa_xDim05';
% estParams = gpfaRezCtx.estParams;
% dat: 1xTrials struct array with two fields - trialId, spikes.

% load parameters
load(fullfile(pathToGpfaRez), 'estParams', 'hasSpikesBool');

% get seq and ensure the y-dim matches with C-dim
seq = getSeq(dat, 1);
[yDim, xDim] = size(estParams.C);

if size(estParams.C, 1)==sum(hasSpikesBool)
    if any(hasSpikesBool==false)
        for t = 1:length(seq)
            seq(t).y = seq(t).y(hasSpikesBool, :);
        end
    end
else
    error("Dimensions of C and dat mismatch!")
end

%% Main computations to get the xsm
if estParams.notes.RforceDiagonal
    Rinv     = diag(1./diag(estParams.R));
    logdet_R = sum(log(diag(estParams.R)));
else
    Rinv     = inv(estParams.R);
    Rinv     = (Rinv+Rinv') / 2; % ensure symmetry
    logdet_R = logdet(estParams.R);
end
CRinv  = estParams.C' * Rinv;
CRinvC = CRinv * estParams.C;

Tall = [seq.T];
Tu   = unique(Tall);

T = Tu;

[K_big, K_big_inv, logdet_K_big] = make_K_big(estParams, T);

% There are three sparse matrices here: K_big, K_big_inv, and CRinvC_inv.
% Choosing which one(s) to make sparse is tricky.  If all are sparse,
% code slows down significantly.  Empirically, code runs fastest if
% only K_big is made sparse.
%
% There are two problems with calling both K_big_inv and CRCinvC_big
% sparse:
% 1) their sum is represented by Matlab as a sparse matrix and taking
%    its inverse is more costly than taking the inverse of the
%    corresponding full matrix.
% 2) term2 has very few zero entries, but Matlab will represent it as a
%    sparse matrix.  This makes later computations with term2 ineffficient.

K_big = sparse(K_big);

blah        = cell(1, T);
[blah{:}]   = deal(CRinvC);
%CRinvC_big = blkdiag(blah{:});     % (xDim*T) x (xDim*T)
[invM, logdet_M] = invPerSymm(K_big_inv + blkdiag(blah{:}), xDim,...
    'offDiagSparse', true);

% Note that posterior covariance does not depend on observations,
% so can compute once for all trials with same T.
% xDim x xDim posterior covariance for each timepoint
Vsm = nan(xDim, xDim, T);
idx = 1: xDim : (xDim*T + 1);
for t = 1:T
    cIdx       = idx(t):idx(t+1)-1;
    Vsm(:,:,t) = invM(cIdx, cIdx);
end

% T x T posterior covariance for each GP
VsmGP = nan(T, T, xDim);
idx   = 0 : xDim : (xDim*(T-1));
for i = 1:xDim
    VsmGP(:,:,i) = invM(idx+i,idx+i);
end

% Process all trials with length T
nList    = find(Tall == T);
dif      = bsxfun(@minus, [seq(nList).y], estParams.d); % yDim x sum(T)
term1Mat = reshape(CRinv * dif, xDim*T, []); % (xDim*T) x length(nList)

% Compute blkProd = CRinvC_big * invM efficiently
% blkProd is block persymmetric, so just compute top half
Thalf   = ceil(T/2);
blkProd = zeros(xDim*Thalf, xDim*T);
idx     = 1: xDim : (xDim*Thalf + 1);
for t = 1:Thalf
    bIdx            = idx(t):idx(t+1)-1;
    blkProd(bIdx,:) = CRinvC * invM(bIdx,:);
end
blkProd = K_big(1:(xDim*Thalf), :) *...
    fillPerSymm(speye(xDim*Thalf, xDim*T) - blkProd, xDim, T);
xsmMat  = fillPerSymm(blkProd, xDim, T) * term1Mat; % (xDim*T) x length(nList)

ctr = 1;
for n = nList
    seq(n).xsm   = reshape(xsmMat(:,ctr), xDim, T);
    seq(n).Vsm   = Vsm;
    seq(n).VsmGP = VsmGP;

    ctr = ctr + 1;
end

%% orthonormalize
[Xorth, Corth] = orthogonalize([seq.xsm], estParams.C);
seq = segmentByTrial(seq, Xorth, 'xorth');
estParams.Corth = Corth;

end
