%% load the encoding model data
%This code is modified from /Users/jp3025/Documents/codes/MatlabNeuralDataPipeline/js2p0_collectCode/runGPFA_collect_sim_50msBin_earlyEpoch.m
filePath = '/Users/jp3025/Documents/codes/MatlabNeuralDataPipeline/mesomapping/encodingGngModelSim/simData';
load(fullfile(filePath, 'simData'), 'selector', '-mat');

%% preprocessing
gngGroupC = cell(2, length(selector.psth));

for i = 1:length(selector.psth) % sessions (20)
    goMat = selector.psth{i}(:, :, 1:2:end); % Go trials
    ngMat = selector.psth{i}(:, :, 2:2:end); % Ng trials

    for ii = 1:size(goMat, 3) % iterate trials
        gngGroupC{1, i}{end+1} = binAvg1msSpkCountMat(max(goMat(:, :, ii), 0), 50, 50); % 50-ms bin
    end

    for ii = 1:size(ngMat, 3) % iterate trials
        gngGroupC{2, i}{end+1} = binAvg1msSpkCountMat(max(ngMat(:, :, ii), 0), 50, 50); % 50-ms bin
    end
end

%% run GPFA
% get datSpec
for ss = 1:size(gngGroupC, 2) % iterate session
    count = 0;
    for gg = 1:2
        for tt = 1:10
            count = count + 1;
            sessionDat(count).trialId = count;
            sessionDat(count).spikes = max(gngGroupC{gg, ss}{tt}, 0);
        end
    end

    saveName = fullfile(filePath, 'gng_gpfa');
    %saveName = fullfile(filePath, sprintf('%s_gngSession_%d_xDim%02d', 'gpfa', ss, 5));
    if exist(saveName, "dir")~=7
        mkdir(saveName)
    end

    % this runs GPFA
    gpfaRez{ss} = neuralTrajAlreadyBinned(saveName, sessionDat, 'method', 'gpfa', ...
            'xDim', 5, 'kernSDList', 5, 'binWidth', 1); %

    % gpfa postprocess
    gpfaRez{ss}.method = 'gpfa';
    % Orthonormalize neural trajectories
    [gpfaRez{ss}.estParamsOrtho, gpfaRez{ss}.seqTrainOrtho] = postprocess(gpfaRez{ss}, 'kernSD', 5);

    gpfaRez{ss}.go = projYtoGetXsmC(gpfaRez{ss}, sessionDat(1:10));  % go trials
    gpfaRez{ss}.ng = projYtoGetXsmC(gpfaRez{ss}, sessionDat(11:20)); % ng trials

    %plot3D_seqCells_color({gpfaRez{ss}.go, gpfaRez{ss}.ng}, {[0 1 1], [1 0 1]}, 10);
end

%% plot 
trialsToPlot = [1, 5, 10, 15, 20]; 
% PSTH
% Parameters for min-max normalization and concatenation
goM_all = [];
ngM_all = [];

%% plot PSTH
for ss = 1:length(selector.psth)
    if ismember(ss, trialsToPlot)
        goM = mean(cat(3, gngGroupC{1,ss}{:}), 3); % go trial mean
        ngM = mean(cat(3, gngGroupC{2,ss}{:}), 3); % ng trial mean
        figure; imagesc(goM); clim([1 10]); set(gca, 'TickDir', 'out')
        print(fullfile(filePath, sprintf('gng_simPSTH_goTrial_%d', ss)), '-bestfit', '-vector', '-dpdf'); 
        figure; imagesc(ngM); clim([1 10]); set(gca, 'TickDir', 'out')
        print(fullfile(filePath, sprintf('gng_simPSTH_ngTrial_%d', ss)), '-bestfit', '-vector', '-dpdf'); 
    end
end

% GPFA trajectories
az =  -100.3136;
el = 5.4844;
for ss = 1:size(gngGroupC, 2) % iterate session
    if ismember(ss, trialsToPlot)
        plot3D_seqCells_colorG({gpfaRez{ss}.go, gpfaRez{ss}.ng}, {[0 1 1], [1 0 1]}, 10);
        view([az, el]);
        saveName = sprintf('gpfaFig_session_%d', ss);
        print(fullfile(filePath, saveName), '-bestfit', '-vector', '-dpdf');
    end
end
