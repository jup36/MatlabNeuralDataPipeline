
filePath = {'/Volumes/Extreme SSD/js2p0/WR37_022119', ... % Dual recording without silencing, Trj checked (07/04/22)
    '/Volumes/Extreme SSD/js2p0/WR37_022219', ... % Dual recording without silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR37_022619', ... % Cg recording contra-Cg silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR37_022719', ... % B Only contra-Cg delayed silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR38_052119', ... % Dual recording without silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR38_052219', ... % Dual recording without silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR38_052319', ... % Cg recording contra-Cg silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR38_052419', ... % Corticostriatal recording M1 silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR39_091019', ... % B Only contra-Cg delayed silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR39_091119', ... % B Only contra-Cg delayed silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR39_100219', ... % Dual recording with contra Cg silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR39_100319', ... % B Only contra-Cg delayed silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR40_081919', ... % Dual recording with contra Cg silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR40_082019', ... % Dual recording with contra Cg silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR44_031020'};    % Dual recording with contra Cg delayed silencing, Trj checked

figSaveDir = '/Volumes/Extreme SSD/js2p0/collectData/collectFigure/';

speedfunc = @(a) sqrt(sum(a.^2, 1));

%% prepare color
pastel2 = slanCM('Pastel2', 10);
% plotColorListWithNumbers(pastel2);

pastel1 = slanCM('Pastel1', 10);
% plotColorListWithNumbers(pastel1);

%% Main
for f = 1:length(filePath)
    cd(filePath{f})
    file = dir('**/*js2p0_tbytSpkHandJsTrjBin_WR*.mat');
    load(fullfile(file.folder, file.name),'ss', 'jkvt')
    rez.b_id{f} = [ss.blNumber];
    [~, fileName] = fileparts(filePath{f});

    if isfield(jkvt, 'stimLaserOn')
        %% define inactivation trials
        trI.stimI{f} = ~cellfun(@isnan, {jkvt.stimLaserOn}); % All inactivation trials
        trI.spI{f} = cellfun(@(a) strcmpi(a, 'sp'), {ss.trialType}); % success trials
        trI.rsAlignI{f} = cellfun(@(a) strcmpi(a, 'rStart'), {ss.evtAlign}); % reachStart aligned trials
        trI.spNoStimI{f} = trI.spI{f} & ~trI.stimI{f}; %& trI.rsAlignI{f};
        trI.spNoStimRsAlignI{f} = trI.spI{f} & ~trI.stimI{f} & trI.rsAlignI{f};
        trI.stimSpRsI{f} = trI.stimI{f} & trI.spI{f} & trI.rsAlignI{f}; % successful reach-aligned inactivation trials
        trI.stimBtI{f} = zeros(length(ss), 1);
        trI.stimFullBtI{f} = zeros(length(ss), 1);

        for i = 1:length(ss)
            if ~isempty(ss(i).spkTimeBlaserI) && ~isempty(ss(i).hTrjB)
                ss(i).hTrjBlaserI = ss(i).spkTimeBlaserI(:, 1:min(size(ss(i).hTrjB, 2), length(ss(i).spkTimeBlaserI)));
                if any(ss(i).hTrjBlaserI)
                    if trI.stimSpRsI{f}(i)
                        if find(ss(i).hTrjBlaserI, 1, 'first')<50 && find(ss(i).hTrjBlaserI, 1, 'last')>51
                            trI.stimBtI{f}(i) = true;
                            % hand velocity of trials with full
                            % inactivation
                            if trI.rsAlignI{f}(i) && sum(ss(i).hTrjBlaserI(51:end))==size(ss(i).hVelB(:, 51:end), 2)
                                trI.stimFullBtI{f}(i) = true;
                                ss(i).hVelB_fullBt = ss(i).hVelB(:, 51:end);
                            end
                        end
                    end
                end
            end
        end

        trI.stimBtId{f} = find(trI.stimBtI{f}); % breakthrough trial ID

        %% descriptive stats
        % success rate, breakthrough rate
        firstSp = find(trI.spI{f}, 1, 'first');
        lastSp = find(trI.spI{f}, 1, 'last');
        valI = 1:length(trI.spI{f}) >= firstSp & 1:length(trI.spI{f}) <= lastSp;

        rez.successRate{f} = sum(trI.spI{f}(valI & ~trI.stimI{f}))/sum(valI & ~trI.stimI{f})*100;
        rez.breakThroughRate{f} = length(trI.stimBtId{f})/sum(trI.stimI{f})*100;

        % reaction time
        rez.rtSp{f} = [jkvt(trI.spNoStimI{f}).rt];    % RT: successful no-stim trials
        rez.rtSpRsAlign{f} = [jkvt(trI.spNoStimRsAlignI{f}).rt];    % RT: successful, reachStart-aligned trials
        rez.rtStimBt{f} = [jkvt(trI.stimBtId{f}).rt];  % RT: breakthrough trials
        if isfield(trI, 'stimFullBtI') && any(trI.stimFullBtI{f}==1)
            rez.rtStimFullBt{f} = [jkvt(trI.stimFullBtI{f}==1).rt];  % RT: breakthrough trials
        end

        % response duration
        rez.rdSp{f} = cell2mat(cellfun(@(a, b) a-b, {jkvt.rStopToPull}, {ss.timeAlign}, 'UniformOutput', false));   % RD: successful no-stim trials
        rez.rdSpRsAlign{f} = cell2mat(cellfun(@(a, b) a-b, {jkvt(trI.spNoStimI{f}).rStopToPull}, {ss(trI.spNoStimI{f}).timeAlign}, 'UniformOutput', false));
        rez.rdStimBt{f} = cell2mat(cellfun(@(a, b) a-b, {jkvt(trI.stimBtId{f}).rStopToPull}, {ss(trI.stimBtId{f}).timeAlign}, 'UniformOutput', false));
        if isfield(trI, 'stimFullBtI') && any(trI.stimFullBtI{f})
            rez.rdStimFullBt{f} = cell2mat(cellfun(@(a, b) a-b, {jkvt(trI.stimFullBtI{f}==1).rStopToPull}, {ss(trI.stimFullBtI{f}==1).timeAlign}, 'UniformOutput', false));
        end

        % reach angle
        rez.absRchAng{f} = abs([ss(trI.spNoStimI{f}).rchAngDeg]); % abs RA: successful no-stim trials
        rez.absRchAngId{f} = getTrialId(fileName, trI.spNoStimI{f}, rez.b_id{f}, {ss.rchAngDeg}, rez.absRchAng{f});
        rez.rchAng{f} = [ss(trI.spNoStimI{f}).rchAngDeg]; % RA: all
        rez.absRsAlignRchAng{f} = abs([ss(trI.spNoStimRsAlignI{f}).rchAngDeg]);     % abs RA: successful no-stim trials
        rez.absRsAlignRchAngId{f} = getTrialId(fileName, trI.spNoStimRsAlignI{f}, rez.b_id{f}, {ss.rchAngDeg}, rez.absRsAlignRchAng{f});
        rez.rsAlignRchAng{f} = [ss(trI.spNoStimRsAlignI{f}).rchAngDeg];
        rez.absRchAngStimBt{f} = abs([ss(trI.stimBtId{f}).rchAngDeg]); % abs RA: breakthrough trials
        rez.rchAngStimBt{f} = [ss(trI.stimBtId{f}).rchAngDeg]; % RA: breakthrough trials
        if isfield(trI, 'stimFullBtI') && any(trI.stimFullBtI{f}==1)
            rez.absRchAngStimFullBt{f} = abs([ss(trI.stimFullBtI{f}==1).rchAngDeg]); % abs RA: breakthrough trials
            rez.absRchAngStimFullBtId{f} = getTrialId(fileName, trI.stimFullBtI{f}==1, rez.b_id{f}, {ss.rchAngDeg}, rez.absRchAngStimFullBt{f});
            rez.rchAngStimFullBt{f} = [ss(trI.stimFullBtI{f}==1).rchAngDeg]; % abs RA: breakthrough trials
            rez.rchAngStimFullBtId{f} = getTrialId(fileName, trI.stimFullBtI{f}==1, rez.b_id{f}, {ss.rchAngDeg}, rez.rchAngStimFullBt{f});
        end

        % tortuosity
        rez.tort{f} = [ss(trI.spNoStimI{f}).hXYtort];     % abs RA: successful no-stim trials
        rez.rsAlignTort{f} = [ss(trI.spNoStimRsAlignI{f}).hXYtort];     % abs RA: successful no-stim trials
        rez.tortStimBt{f} = [ss(trI.stimBtId{f}).hXYtort]; % abs RA: breakthrough trials
        if isfield(trI, 'stimFullBtI') && any(trI.stimFullBtI{f}==1)
            rez.tortStimFullBt{f} = [ss(trI.stimFullBtI{f}==1).hXYtort]; % abs RA: breakthrough trials
        end

        % hand trajectories
        rez.trjNoStim{f} = {ss(trI.spNoStimI{f}).hTrjB};
        rez.trjNoStimRsAlign{f} = {ss(trI.spNoStimRsAlignI{f}).hTrjB};
        rez.trjStimBt{f} = {ss(trI.stimBtId{f}).hTrjB};
        if isfield(trI, 'stimFullBtI') && any(trI.stimFullBtI{f}==1)
            rez.trjStimFullBt{f} = {ss(trI.stimFullBtI{f}==1).hTrjB};
        end

        % Assuming ss(1).hVelB is your 3-by-100 matrix
        rez.spdNoStim{f} = cellfun(@(a) speedfunc(a), {ss(trI.spNoStimI{f}).hVelB}, 'UniformOutput', false);
        rez.spdNoStimRsAlign{f} = cellfun(@(a) speedfunc(a), {ss(trI.spNoStimRsAlignI{f}).hVelB}, 'UniformOutput', false);
        rez.spdStimBt{f} = cellfun(@(a) speedfunc(a), {ss(trI.stimBtId{f}).hVelB}, 'UniformOutput', false);

        % speed cut baseline time bins
        rez.spdNoStim{f} = cellfun(@(a) a(51:end), rez.spdNoStim{f}, 'UniformOutput', false);
        rez.spdNoStimRsAlign{f} = cellfun(@(a) a(51:end), rez.spdNoStimRsAlign{f}, 'UniformOutput', false);
        rez.spdStimBt{f} = cellfun(@(a) a(51:end), rez.spdStimBt{f}, 'UniformOutput', false);

        % max speed
        rez.maxSpdNoStim{f} = cell2mat(cellfun(@nanmax, rez.spdNoStim{f}, 'UniformOutput', false));
        rez.maxSpdNoStimRsAlign{f} = cell2mat(cellfun(@nanmax, rez.spdNoStimRsAlign{f}, 'UniformOutput', false));
        rez.maxSpdStimBt{f} = cell2mat(cellfun(@nanmax, rez.spdStimBt{f}, 'UniformOutput', false));

        % mean speed
        rez.meanSpdNoStim{f} = cell2mat(cellfun(@nanmean, rez.spdNoStim{f}, 'UniformOutput', false));
        rez.meanSpdNoStimRsAlign{f} = cell2mat(cellfun(@nanmean, rez.spdNoStimRsAlign{f}, 'UniformOutput', false));
        rez.meanSpdStimBt{f} = cell2mat(cellfun(@nanmean, rez.spdStimBt{f}, 'UniformOutput', false));

        % examine speed in the full breakthrough trials
        if isfield(trI, 'stimFullBtI') && any(trI.stimFullBtI{f}==1) && isfield(ss, 'hVelB_fullBt')
            rez.spdStimFullBt{f} = cellfun(@(a) speedfunc(a), {ss.hVelB_fullBt}, 'UniformOutput', false);
            rez.maxSpdStimFullBt{f} = cell2mat(cellfun(@nanmax, rez.spdStimFullBt{f}, 'UniformOutput', false));  % max speed breakthrough with full inactivation
            rez.meanSpdStimFullBt{f} = cell2mat(cellfun(@nanmean, rez.spdStimFullBt{f}, 'UniformOutput', false));  % mean speed breakthrough with full inactivation
            rez.meanSpdStimFullBt{f} = rez.meanSpdStimFullBt{f}(~isnan(rez.meanSpdStimFullBt{f})); 
        end

        % binned force (forceB)
        rez.forceBNoStim{f} = {ss(trI.spNoStimI{f}).forceB};
        rez.forceBNoStimRsAlign{f} = {ss(trI.spNoStimRsAlignI{f}).forceB};
        rez.forceBNoStimRsAlignId{f} = getTrialId(fileName, trI.spNoStimRsAlignI{f}, rez.b_id{f}, {ss.forceB}, rez.forceBNoStimRsAlign{f});
        rez.forceBstimBt{f} = {ss(trI.stimBtId{f}).forceB};
        if isfield(trI, 'stimFullBtI') && any(trI.stimFullBtI{f}==1)
            rez.forceBstimFullBt{f} = {ss(trI.stimFullBtI{f}==1).forceB};
            rez.forceBstimFullBtId{f} = getTrialId(fileName, trI.stimFullBtI{f}==1, rez.b_id{f}, {ss.forceB}, rez.forceBstimFullBt{f});
        end

        % max force
        rez.maxForceNoStim{f} = cell2mat(cellfun(@(a) nanmax(abs(a)), rez.forceBNoStim{f}, 'UniformOutput', false));
        rez.maxForceNoStimRsAlign{f} = cell2mat(cellfun(@(a) nanmax(abs(a)), rez.forceBNoStimRsAlign{f}, 'UniformOutput', false));
        rez.maxForceStimBt{f} = cell2mat(cellfun(@(a) nanmax(abs(a)), rez.forceBstimBt{f}, 'UniformOutput', false));

        % mean force
        rez.meanForceNoStim{f} = cell2mat(cellfun(@(a) nanmean(abs(a)), rez.forceBNoStim{f}, 'UniformOutput', false));
        rez.meanForceNoStimRsAlign{f} = cell2mat(cellfun(@(a) nanmean(abs(a)), rez.forceBNoStimRsAlign{f}, 'UniformOutput', false));
        rez.meanForceStimBt{f} = cell2mat(cellfun(@(a) nanmean(abs(a)), rez.forceBstimBt{f}, 'UniformOutput', false));

        % examine force in the full breakthrough trials
        if isfield(trI, 'stimFullBtI') && any(trI.stimFullBtI{f}==1) && isfield(ss, 'hVelB_fullBt')
            rez.maxForceStimFullBt{f} = cell2mat(cellfun(@(a) nanmax(abs(a)), rez.forceStimFullBt{f}, 'UniformOutput', false));  % max speed breakthrough with full inactivation
            rez.meanForceStimFullBt{f} = cell2mat(cellfun(@(a) nanmean(abs(a)), rez.forceStimFullBt{f}, 'UniformOutput', false));  % mean speed breakthrough with full inactivation
        end
    end
    fprintf('processed file #%d\n', f)
end

%% Initial hand position
rezCol.trjNoStimRsAlignInitX = cell(1, length(rez.trjNoStimRsAlign));
rezCol.trjStimFullBtInitX = cell(1, length(rez.trjStimFullBt));

for f = 1:length(rez.trjNoStimRsAlign)
    if ~isempty(rez.trjNoStimRsAlign{f})
        % block ID to identify left vs right trials
        bIdNoStimRsAlign = rez.b_id{f}(trI.spNoStimRsAlignI{f});
        bIdNoStimRsAlignLeftLogic = ismember(bIdNoStimRsAlign, [1, 2, 5, 6]);
        % noStimRsAlign Init X
        rezCol.trjNoStimRsAlignInitX{1, f} = cell2mat(cellfun(@(a) nanmedian(a(1, 1:30)), rez.trjNoStimRsAlign{f}, 'UniformOutput', false));
        % noStimRsAlign Init Y
        rezCol.trjNoStimRsAlignInitY{1, f} = cell2mat(cellfun(@(a) nanmedian(a(2, 1:30)), rez.trjNoStimRsAlign{f}, 'UniformOutput', false));
        % median X subtraction
        medX = nanmedian(rezCol.trjNoStimRsAlignInitX{1, f});
        rezCol.trjNoStimRsAlignInitX{1, f} = rezCol.trjNoStimRsAlignInitX{1, f}-medX;
        % median Y subtraction
        medY = nanmedian(rezCol.trjNoStimRsAlignInitY{1, f});
        rezCol.trjNoStimRsAlignInitY{1, f} = rezCol.trjNoStimRsAlignInitY{1, f}-medY;
        % left trials
        rezCol.trjNoStimRsAlignInitXY_LTrial{1, f} = [rezCol.trjNoStimRsAlignInitX{1, f}(bIdNoStimRsAlignLeftLogic); ...
            rezCol.trjNoStimRsAlignInitY{1, f}(bIdNoStimRsAlignLeftLogic)];
        % right trials
        rezCol.trjNoStimRsAlignInitXY_RTrial{1, f} = [rezCol.trjNoStimRsAlignInitX{1, f}(~bIdNoStimRsAlignLeftLogic); ...
            rezCol.trjNoStimRsAlignInitY{1, f}(~bIdNoStimRsAlignLeftLogic)];
    end

    if ~isempty(rez.trjStimFullBt{f})
        % stimFullBt block ID to identify left vs right trials
        bIdStimFullBt = rez.b_id{f}(trI.stimFullBtI{f}==1);
        bIdStimFullBtLeftLogic = ismember(bIdStimFullBt, [1, 2, 5, 6]);
        % stimFullBt Init X
        rezCol.trjStimFullBtInitX{1, f} = cell2mat(cellfun(@(a) nanmedian(a(1, 1:30)), rez.trjStimFullBt{f}, 'UniformOutput', false));
        % stimFullBt Init Y
        rezCol.trjStimFullBtInitY{1, f} = cell2mat(cellfun(@(a) nanmedian(a(2, 1:30)), rez.trjStimFullBt{f}, 'UniformOutput', false));
        % median X subtraction
        rezCol.trjStimFullBtInitX{1, f} = rezCol.trjStimFullBtInitX{1, f}-medX;
        % median Y subtraction
        rezCol.trjStimFullBtInitY{1, f} = rezCol.trjStimFullBtInitY{1, f}-medY;
        % left trials
        rezCol.trjStimFullBtInitXY_LTrial{1, f} = [rezCol.trjNoStimRsAlignInitX{1, f}(bIdStimFullBtLeftLogic); ...
            rezCol.trjNoStimRsAlignInitY{1, f}(bIdStimFullBtLeftLogic)];
        % right trials
        rezCol.trjStimFullBtInitXY_RTrial{1, f} = [rezCol.trjNoStimRsAlignInitX{1, f}(~bIdStimFullBtLeftLogic); ...
            rezCol.trjNoStimRsAlignInitY{1, f}(~bIdStimFullBtLeftLogic)];
    end
end
clearvars f

%% Success rate & breakthrough rate
rezCol.successRate = cell2mat(rez.successRate);
rezCol.breakThroughRate = cell2mat(rez.breakThroughRate);

[~, stat.sRate.p, ~, stat.sRate.stats] = ttest2(rezCol.successRate, rezCol.breakThroughRate);

nameSes = cellfun(@(a) extractNameFromPath(a), filePath, 'UniformOutput', false);
name = unique(cellfun(@(a) a(1:4), nameSes, 'UniformOutput', false));

for f = 1:length(name)
    mLogic = cellfun(@(a) contains(a, name{f}), nameSes);
    rezCol.sRatePerM{f, 1} = name{f};
    rezCol.sRatePerM{f, 2} = nanmean(cell2mat(rez.successRate(mLogic)));
    rezCol.BtRatePerM{f, 1} = name{f};
    rezCol.BtRatePerM{f, 2} = nanmean(cell2mat(rez.breakThroughRate(mLogic)));
end

[~, stat.sRatePerM.p, ~, stat.sRatePerM.stats] = ttest2([rezCol.sRatePerM{:, 2}]', [rezCol.BtRatePerM{:, 2}]');

fig_sRatePerM = scatter_row_by_row([[rezCol.sRatePerM{:, 2}]', [rezCol.BtRatePerM{:, 2}]'], [0 90]);
xlim([0.5 2.5])
set(gca, 'TickDir', 'out')
print(fig_sRatePerM, fullfile(figSaveDir, 'sRate_btRate_control_inactivation_trials'), '-dpdf', '-vector')

%% RT (reaction time) stat
rezCol.controlRt = cell2mat(rez.rtSp);
rezCol.controlRt(rezCol.controlRt<50) = NaN;
rezCol.controlRsAlignRt = cell2mat(rez.rtSpRsAlign);
rezCol.controlRsAlignRt(rezCol.controlRsAlignRt<50) = NaN;
rezCol.stimBtRt = cell2mat(rez.rtStimBt);
rezCol.stimBtRt(rezCol.stimBtRt<50) = NaN;
rezCol.stimFullBtRt = cell2mat(rez.rtStimFullBt);
rezCol.stimFullBtRt(rezCol.stimFullBtRt<50) = NaN;

%[~, stat.rt.p, ~, stat.rt.stats] = ttest2(rezCol.controlRt, rezCol.stimBtRt);
[~, stat.rt.p, ~, stat.rt.stats] = ttest2(rezCol.controlRt(rezCol.controlRt<10000), rezCol.stimBtRt(rezCol.stimBtRt<10000));

histogramProbDensity({rezCol.controlRt, rezCol.stimBtRt}, {pastel1(3, :), pastel2(5, :)}, 0.4, {100, 25})
xlim([0 10000])
%set(gca, 'YScale', 'log');
%set(gca, 'TickDir', 'out')

histogramProbDense({rezCol.controlRt(rezCol.controlRt<10000)./100, rezCol.stimBtRt(rezCol.stimBtRt<10000)./100}, {pastel1(3, :), pastel2(5, :)}, 0.4);
%integralValue = trapz(rtPdfPoints{2}, rtPdfValues{2});

set(gca, 'TickDir', 'out')
%set(gca, 'YScale', 'log');
xlim([0 100])
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'stimNoStimRTpdf'), '-dpdf', '-vector')

%% RD (response duration) stat
rezCol.controlRd = cell2mat(rez.rdSp);
rezCol.stimBtRd = cell2mat(rez.rdStimBt);
rezCol.controlRd = cell2mat(rez.rdSp);
rezCol.controlRsAlignRd = cell2mat(rez.rdSpRsAlign);
rezCol.stimBtRd = cell2mat(rez.rdStimBt);
rezCol.stimFullBtRd = cell2mat(rez.rdStimFullBt);

[~, stat.rd.p, ~, stat.rd.stats] = ttest2(rezCol.controlRd, rezCol.stimBtRd);

[rdPdfValues, rdPdfPoints] = histogramProbDense({rezCol.controlRd./10, rezCol.stimBtRd./10}, {pastel1(3, :), pastel2(5, :)}, 0.4);
integralValue = trapz(rdPdfPoints{1}, rdPdfValues{1});
xlim([0 200])
set(gca, 'TickDir', 'out')
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'stimNoStimRDpdf'), '-dpdf', '-vector')

%% reach angle
rezCol.controlRchAng = cell2mat(rez.absRchAng);
rezCol.controlRchAng(rezCol.controlRchAng==0) = NaN;

rezCol.controlRchAngRaw = cell2mat(rez.rchAng);
rezCol.controlRchAngRaw(rezCol.controlRchAngRaw==0) = NaN;

rezCol.controlRsAlignRchAng = cell2mat(rez.absRsAlignRchAng);
rezCol.controlRsAlignRchAng(rezCol.controlRsAlignRchAng==0) = NaN;
rezCol.controlRsAlignRchAngTrId = combineCellsIntoCell(rez.absRsAlignRchAngId); 

rezCol.controlRsAlignRchAngRaw = cell2mat(rez.rsAlignRchAng);
rezCol.controlRsAlignRchAngRaw(rezCol.controlRsAlignRchAngRaw==0) = NaN;

rezCol.stimBtRchAng = cell2mat(rez.absRchAngStimBt);
rezCol.stimBtRchAng(rezCol.stimBtRchAng==0) = NaN;

rezCol.stimBtRchAngRaw = cell2mat(rez.rchAngStimBt);
rezCol.stimBtRchAngRaw(rezCol.stimBtRchAngRaw==0) = NaN;

rezCol.stimFullBtRchAng = cell2mat(rez.absRchAngStimFullBt);
rezCol.stimFullBtRchAng(rezCol.stimFullBtRchAng==0) = NaN;
rezCol.stimFullBtRchAngTrId = combineCellsIntoCell(rez.absRchAngStimFullBtId); 

rezCol.stimFullBtRchAngRaw = cell2mat(rez.rchAngStimFullBt);
rezCol.stimFullBtRchAngRaw(rezCol.stimFullBtRchAngRaw==0) = NaN;

%[~, stat.rchAng.p, ~, stat.rchAng.stats] = ttest2(rezCol.controlRchAng, rezCol.stimBtRchAng);
[~, stat.rchAng.p, ~, stat.rchAng.stats] = ttest2(rezCol.controlRchAng, rezCol.stimFullBtRchAng);

%% max speed
rezCol.controlMaxSpd = cell2mat(rez.maxSpdNoStim);
rezCol.controlRsAlignMaxSpd = cell2mat(rez.maxSpdNoStimRsAlign);
rezCol.stimBtMaxSpd = cell2mat(rez.maxSpdStimBt);
rezCol.stimFullBtMaxSpd = cell2mat(rez.maxSpdStimFullBt);

[~, stat.maxSpd.p, ~, stat.maxSpd.stats] = ttest2(rezCol.controlRsAlignMaxSpd, rezCol.stimFullBtMaxSpd);

[maxSpdPdfValues, maxSpdPdfPoints] = histogramProbDense({rezCol.controlRsAlignMaxSpd, rezCol.stimFullBtMaxSpd}, {pastel1(3, :), pastel2(5, :)}, 0.4);
xlim([0 120])
set(gca, 'TickDir', 'out')
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'stimNoStimMaxSpdPdf'), '-dpdf', '-vector')

%% mean speed
rezCol.controlMeanSpd = cell2mat(rez.meanSpdNoStim);
rezCol.controlRsAlignMeanSpd = cell2mat(rez.meanSpdNoStimRsAlign);
rezCol.stimBtMeanSpd = cell2mat(rez.meanSpdStimBt);
rezCol.stimFullBtMeanSpd = cell2mat(rez.meanSpdStimFullBt);

[~, stat.meanSpd.p, ~, stat.meanSpd.stats] = ttest2(rezCol.controlRsAlignMeanSpd, rezCol.stimBtMeanSpd);

[meanSpdPdfValues, meanSpdPdfPoints] = histogramProbDense({rezCol.controlRsAlignMeanSpd, rezCol.stimBtMeanSpd}, {pastel1(3, :), pastel2(5, :)}, 0.4);
integralValue = trapz(meanSpdPdfPoints{1}, meanSpdPdfValues{1});
xlim([0 40])
set(gca, 'TickDir', 'out')
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'stimNoStimMeanSpdPdf'), '-dpdf', '-vector')

%% tortuosity
rezCol.controlTort = cell2mat(rez.tort);
rezCol.controlRsAlignTort = cell2mat(rez.rsAlignTort);
rezCol.stimBtTort = cell2mat(rez.tortStimBt);
rezCol.stimFullBtTort = cell2mat(rez.tortStimFullBt);

[~, stat.tort.p, ~, stat.tort.stats] = ttest2(rezCol.controlRsAlignTort, rezCol.stimBtTort);

[tortPdfValues, tortPdfPoints] = histogramProbDense({rezCol.controlRsAlignTort.*10, rezCol.stimBtTort.*10}, {pastel1(3, :), pastel2(5, :)}, 0.4);
integralValue = trapz(tortPdfPoints{1}, tortPdfValues{1});
xlim([0 250])
set(gca, 'TickDir', 'out')
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'stimNoStimTortPdf'), '-dpdf', '-vector')

%% max force
rezCol.controlMaxForce = cell2mat(rez.maxForceNoStim);
rezCol.controlRsAlignMaxForce = cell2mat(rez.maxForceNoStimRsAlign);
rezCol.controlRsAlignForceId = combineCellsIntoCell(rez.forceBNoStimRsAlignId); 

rezCol.stimBtMaxForce = cell2mat(rez.maxForceStimBt);
rezCol.stimFullBtMaxForce = cell2mat(rez.maxForceStimFullBt);
rezCol.stimFullBtForceId = combineCellsIntoCell(rez.forceBstimFullBtId); 

[~, stat.maxForce.p, ~, stat.maxForce.stats] = ttest2(rezCol.controlRsAlignMaxForce, rezCol.stimFullBtMaxForce);
[forcePdfValues, forcePdfPoints] = histogramProbDense({rezCol.controlRsAlignMaxForce, rezCol.stimFullBtMaxForce}, {pastel1(3, :), pastel2(5, :)}, 0.4);
xlim([0 300])
set(gca, 'TickDir', 'out')
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'stimNoStimForcePdf'), '-dpdf', '-vector')

%% mean force
rezCol.controlMeanForce = cell2mat(rez.meanForceNoStim);
rezCol.controlRsAlignMeanForce = cell2mat(rez.meanForceNoStimRsAlign);
rezCol.stimBtMeanForce = cell2mat(rez.meanForceStimBt);
rezCol.stimFullBtMeanForce = cell2mat(rez.meanForceStimFullBt);

%% save
%save(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI')
% load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI')

%% Visualization
% histogram distribution with trial fraction not trial count on as the quantity

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [combinedIds, valTrialId, valBlockId] = getTrialId(fileName, trialI, blockI, datCell, refDat)
% datCell = {ss.rchAngDeg};
% trialI = trI.spNoStimI{f};
% blockI = rez.b_id{f};
% refDat = rez.absRchAng{f};

valDatI = ~cellfun(@isempty, datCell); 

valTrialI = valDatI(:) & trialI(:);
valTrialId = find(valTrialI);
assert(sum(valTrialI)==length(refDat))
valBlockId = blockI(valTrialI);
combinedIds = cell(sum(valTrialI), 1);

for i = 1:length(valBlockId)
    switch valBlockId(i)
        case {1, 5}
            combinedIds{i} = strcat(fileName, '_', sprintf("tr%d_bl%d", valTrialId(i), valBlockId(i)), '_', 'lelo');
        case {2, 6}
            combinedIds{i} = strcat(fileName, '_', sprintf("tr%d_bl%d", valTrialId(i), valBlockId(i)), '_', 'lehi');
        case {3, 7}
            combinedIds{i} = strcat(fileName, '_', sprintf("tr%d_bl%d", valTrialId(i), valBlockId(i)), '_', 'rilo');
        case {4, 8}
            combinedIds{i} = strcat(fileName, '_', sprintf("tr%d_bl%d", valTrialId(i), valBlockId(i)), '_', 'rihi');
    end
end
end
