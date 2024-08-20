
filePaths = {'/Volumes/Extreme SSD/js2p0/WR37_022119', ... % Dual recording without silencing, Trj checked (07/04/22)
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

sessionOfInterest = {'WR38_052119', 'WR38_052419', 'WR39_100219', 'WR40_081919', 'WR40_082019', 'WR44_031020'};
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI')

for f = 1:length(filePath)
    clearvars dat*
    [~, fileName] = fileparts(filePaths{f});
    if ismember(fileName, sessionOfInterest)
        filePath = GrabFiles_sort_trials('js2p0_tbytSpkHandJsTrjBin_50ms_stimPstimPrepExtWoTo_', 0, ...,
            cellfun(@(a) fullfile(a, 'Matfiles'), filePaths(f), 'UniformOutput', false));
        [~, header] = fileparts(filePaths{f}); 
        figSaveDir = fullfile(filePaths{f}, 'Matfiles', 'Figure');

        if ~isempty(filePath)
            load(filePath{1}, 'ss') % load the data structure
            %% get datCtx and datStr
            assert(size(ss, 2)==length(trI.stimI{f}))

            % trial ids
            ctrlTrId = find(~trI.stimI{f} & trI.spI{f});
            stimFullBtId = find(trI.stimFullBtI{f});
            stimId = find(trI.stimI{f});
            stimSpRsId = find(trI.stimSpRsI{f});
            blockIC = {ss.blNumber};
            leLiI = find(cellfun(@(a) ismember(a, [1, 5]), blockIC));
            leHiI = find(cellfun(@(a) ismember(a, [2, 6]), blockIC));
            riLiI = find(cellfun(@(a) ismember(a, [3, 7]), blockIC));
            riHiI = find(cellfun(@(a) ismember(a, [4, 8]), blockIC));

            %% utbCtx
            if isfield(ss, 'unitTimeBCtx')
                utbCtx_ctrlTrsC = cellfun(@full, {ss(ctrlTrId).unitTimeBCtx}, 'UniformOutput', false);
                utbCtx_ctrlTrs = cat(3, utbCtx_ctrlTrsC{:});
                utbCtx_ctrlTrsZ = zscoreNormalizationUnitTimeTrial(utbCtx_ctrlTrs, 15);

                utbCtx_ctrlTrsZ_leLi = utbCtx_ctrlTrsZ(:, :, ismember(ctrlTrId, leLiI));
                utbCtx_ctrlTrsZ_leHi = utbCtx_ctrlTrsZ(:, :, ismember(ctrlTrId, leHiI));
                utbCtx_ctrlTrsZ_riLi = utbCtx_ctrlTrsZ(:, :, ismember(ctrlTrId, riLiI));
                utbCtx_ctrlTrsZ_riHi = utbCtx_ctrlTrsZ(:, :, ismember(ctrlTrId, riHiI));

                figure;
                imagesc(nanmean(utbCtx_ctrlTrsZ_leLi, 3));
                clim([-2 4]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbCtx_ctrlTrsZ_leLi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbCtx_ctrlTrsZ_leHi, 3));
                clim([-2 4]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbCtx_ctrlTrsZ_leHi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbCtx_ctrlTrsZ_riLi, 3));
                clim([-2 4]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbCtx_ctrlTrsZ_riLi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbCtx_ctrlTrsZ_riHi, 3));
                clim([-2 4]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbCtx_ctrlTrsZ_riHi"), '-dpdf', '-vector')
                close;

                [utbCtx_ctrlTrsZ_leLi_Mean, ~, utbCtx_ctrlTrsZ_leLi_Sem] = meanstdsem(smooth2a(nanmean(utbCtx_ctrlTrsZ_leLi, 3), 0, 1));
                [utbCtx_ctrlTrsZ_leHi_Mean, ~, utbCtx_ctrlTrsZ_leHi_Sem] = meanstdsem(smooth2a(nanmean(utbCtx_ctrlTrsZ_leHi, 3), 0, 1));
                [utbCtx_ctrlTrsZ_riLi_Mean, ~, utbCtx_ctrlTrsZ_riLi_Sem] = meanstdsem(smooth2a(nanmean(utbCtx_ctrlTrsZ_riLi, 3), 0, 1));
                [utbCtx_ctrlTrsZ_riHi_Mean, ~, utbCtx_ctrlTrsZ_riHi_Sem] = meanstdsem(smooth2a(nanmean(utbCtx_ctrlTrsZ_riHi, 3), 0, 1));

                figure;
                meanSemErrorBarOverTime([utbCtx_ctrlTrsZ_leLi_Mean; utbCtx_ctrlTrsZ_leHi_Mean; utbCtx_ctrlTrsZ_riLi_Mean; utbCtx_ctrlTrsZ_riHi_Mean], ...
                    [utbCtx_ctrlTrsZ_leLi_Sem; utbCtx_ctrlTrsZ_leHi_Sem; utbCtx_ctrlTrsZ_riLi_Sem; utbCtx_ctrlTrsZ_riHi_Sem], ...
                    linspace(-1, 2, 60), [248,129,88; 232,40,41; 58,84,162; 51,103,154]./255);
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "meanSemCtx_ctrlTrsZ_trialTypes"), '-dpdf', '-vector')
                close
            end

            %% utbStr
            if isfield(ss, 'unitTimeBStr')
                utbStr_ctrlTrsC = cellfun(@full, {ss(ctrlTrId).unitTimeBStr}, 'UniformOutput', false);
                utbStr_ctrlTrs = cat(3, utbStr_ctrlTrsC{:});
                utbStr_ctrlTrsZ = zscoreNormalizationUnitTimeTrial(utbStr_ctrlTrs, 15);

                utbStr_ctrlTrsZ_leLi = utbStr_ctrlTrsZ(:, :, ismember(ctrlTrId, leLiI));
                utbStr_ctrlTrsZ_leHi = utbStr_ctrlTrsZ(:, :, ismember(ctrlTrId, leHiI));
                utbStr_ctrlTrsZ_riLi = utbStr_ctrlTrsZ(:, :, ismember(ctrlTrId, riLiI));
                utbStr_ctrlTrsZ_riHi = utbStr_ctrlTrsZ(:, :, ismember(ctrlTrId, riHiI));

                figure;
                imagesc(nanmean(utbStr_ctrlTrsZ_leLi, 3));
                clim([-2 4]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbStr_ctrlTrsZ_leLi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbStr_ctrlTrsZ_leHi, 3));
                clim([-2 4]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbStr_ctrlTrsZ_leHi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbStr_ctrlTrsZ_riLi, 3));
                clim([-2 4]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbStr_ctrlTrsZ_riRi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbStr_ctrlTrsZ_riHi, 3));
                clim([-2 4]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbStr_ctrlTrsZ_riHi"), '-dpdf', '-vector')
                close;

                [utbStr_ctrlTrsZ_leLi_Mean, ~, utbStr_ctrlTrsZ_leLi_Sem] = meanstdsem(smooth2a(nanmean(utbStr_ctrlTrsZ_leLi, 3), 0, 1));
                [utbStr_ctrlTrsZ_leHi_Mean, ~, utbStr_ctrlTrsZ_leHi_Sem] = meanstdsem(smooth2a(nanmean(utbStr_ctrlTrsZ_leHi, 3), 0, 1));
                [utbStr_ctrlTrsZ_riLi_Mean, ~, utbStr_ctrlTrsZ_riLi_Sem] = meanstdsem(smooth2a(nanmean(utbStr_ctrlTrsZ_riLi, 3), 0, 1));
                [utbStr_ctrlTrsZ_riHi_Mean, ~, utbStr_ctrlTrsZ_riHi_Sem] = meanstdsem(smooth2a(nanmean(utbStr_ctrlTrsZ_riHi, 3), 0, 1));

                figure;
                meanSemErrorBarOverTime([utbStr_ctrlTrsZ_leLi_Mean; utbStr_ctrlTrsZ_leHi_Mean; utbStr_ctrlTrsZ_riLi_Mean; utbStr_ctrlTrsZ_riHi_Mean], ...
                    [utbStr_ctrlTrsZ_leLi_Sem; utbStr_ctrlTrsZ_leHi_Sem; utbStr_ctrlTrsZ_riLi_Sem; utbStr_ctrlTrsZ_riHi_Sem], ...
                    linspace(-1, 2, 60), [248,129,88; 232,40,41; 58,84,162; 51,103,154]./255);
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "meanSemStr_ctrlTrsZ_trialTypes"), '-dpdf', '-vector')
                close
            end

            %% utbCg
            if isfield(ss, 'unitTimeBCg')
                utbCg_ctrlTrsC = cellfun(@full, {ss(ctrlTrId).unitTimeBCg}, 'UniformOutput', false);
                utbCg_ctrlTrs = cat(3, utbCg_ctrlTrsC{:});
                utbCg_ctrlTrsZ = zscoreNormalizationUnitTimeTrial(utbCg_ctrlTrs, 15);

                utbCg_ctrlTrsZ_leLi = utbCg_ctrlTrsZ(:, :, ismember(ctrlTrId, leLiI));
                utbCg_ctrlTrsZ_leHi = utbCg_ctrlTrsZ(:, :, ismember(ctrlTrId, leHiI));
                utbCg_ctrlTrsZ_riLi = utbCg_ctrlTrsZ(:, :, ismember(ctrlTrId, riLiI));
                utbCg_ctrlTrsZ_riHi = utbCg_ctrlTrsZ(:, :, ismember(ctrlTrId, riHiI));

                figure;
                imagesc(nanmean(utbCg_ctrlTrsZ_leLi, 3));
                clim([-2 3]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbCg_ctrlTrsZ_leLi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbCg_ctrlTrsZ_leHi, 3));
                clim([-2 3]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbCg_ctrlTrsZ_leHi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbCg_ctrlTrsZ_riLi, 3));
                clim([-2 3]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbCg_ctrlTrsZ_riRi"), '-dpdf', '-vector')
                close;

                figure;
                imagesc(nanmean(utbCg_ctrlTrsZ_riHi, 3));
                clim([-2 3]); colormap('hot')
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "utbCg_ctrlTrsZ_riHi"), '-dpdf', '-vector')
                close;

                [utbCg_ctrlTrsZ_leLi_Mean, ~, utbCg_ctrlTrsZ_leLi_Sem] = meanstdsem(smooth2a(nanmean(utbCg_ctrlTrsZ_leLi, 3), 0, 1));
                [utbCg_ctrlTrsZ_leHi_Mean, ~, utbCg_ctrlTrsZ_leHi_Sem] = meanstdsem(smooth2a(nanmean(utbCg_ctrlTrsZ_leHi, 3), 0, 1));
                [utbCg_ctrlTrsZ_riLi_Mean, ~, utbCg_ctrlTrsZ_riLi_Sem] = meanstdsem(smooth2a(nanmean(utbCg_ctrlTrsZ_riLi, 3), 0, 1));
                [utbCg_ctrlTrsZ_riHi_Mean, ~, utbCg_ctrlTrsZ_riHi_Sem] = meanstdsem(smooth2a(nanmean(utbCg_ctrlTrsZ_riHi, 3), 0, 1));

                figure;
                meanSemErrorBarOverTime([utbCg_ctrlTrsZ_leLi_Mean; utbCg_ctrlTrsZ_leHi_Mean; utbCg_ctrlTrsZ_riLi_Mean; utbCg_ctrlTrsZ_riHi_Mean], ...
                    [utbCg_ctrlTrsZ_leLi_Sem; utbCg_ctrlTrsZ_leHi_Sem; utbCg_ctrlTrsZ_riLi_Sem; utbCg_ctrlTrsZ_riHi_Sem], ...
                    linspace(-1, 2, 60), [248,129,88; 232,40,41; 58,84,162; 51,103,154]./255);
                set(gca, 'TickDir', 'out')
                print(fullfile(figSaveDir, "meanSemCg_ctrlTrsZ_trialTypes"), '-dpdf', '-vector')
                close
            end

        end
    end
    clearvars dat*
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
