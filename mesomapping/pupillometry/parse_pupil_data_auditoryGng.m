function [ellipse_area, ellipse_center] = parse_pupil_data_auditoryGng(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623';
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';

%% whereabouts
% animal ID
mID = regexp(filePath, 'DA\d{1,3}', 'match', 'once');
[~, header] = fileparts(filePath); 

filePath_pupil = GrabFiles_sort_trials('pupilcropped', 0, {filePath});

% csv files
filePath_csv = dir(fullfile(filePath_pupil{1}, '*.csv'));
filePath_csv_folder = filePath_csv(1).folder;
filePath_csv_nameC = {filePath_csv.name};
filePath_csv_nameC_id = cellfun(@(a) str2double(regexp(a, [mID '_(\d{1,3})'], 'tokens', 'once')), filePath_csv_nameC, 'UniformOutput', true);
[~, sortI] = sort(filePath_csv_nameC_id);
filePath_csv_nameC = filePath_csv_nameC(sortI);

% mp4 files
filePath_mp4 = GrabFiles_sort_trials('cropped.mp4', 0, filePath_pupil);

% get tbytDat
tbytDatPath = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
load(fullfile(tbytDatPath{1}), 'tbytDat')

% get tbytDat_parseGng
tbytDatDatPath = GrabFiles_sort_trials('tbytDat_parseGng', 0, {fullfile(filePath, 'Matfiles')});
tbytDat = load(fullfile(tbytDatDatPath{1}), 'tbytDat');
tbytDat = tbytDat.('tbytDat');

% get the pixel to mm conversion constant
if exist(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj', 'pixelToMm.mat'), 'file')==2
    load(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj', 'pixelToMm.mat'), 'pixelToMm')
else
    pixelToMm = unitConversionPixelToMilli(tbytDat(1).vidFile); % ensure to do calibration on raw video before cropping
end

%fR = 200; % the faceCam frame rate: 200 Hz

%% parse
ellipse_areaC = cell(numel(tbytDat), 1);
vidReader = VideoReader(filePath_mp4{1}); % open the first video to get some info

%
lickPeth = -2:0.02:2; % time relative to the hit lick (first lick that triggers the reward)
cropLick = @(b) b(b >= 2 & b <= 6); % time relative to cue onset (t=0)
cueAlignedTs = -1:0.02:6;


for f = 1:numel(filePath_csv_nameC)
    csvTab = readDlcCsv(fullfile(filePath_csv_folder, filePath_csv_nameC{f}));
    csvTab = removevars(csvTab, 'bodyparts_coords'); % drop redundant info
    newVarNames = {'L_x', 'L_y', 'L_like', ...
        'LD_x', 'LD_y', 'LD_like', ...
        'D_x', 'D_y', 'D_like', ...
        'DR_x', 'DR_y', 'DR_like', ...
        'R_x', 'R_y', 'R_like', ...
        'RV_x', 'RV_y', 'RV_like', ...
        'V_x', 'V_y', 'V_like', ...
        'VL_x', 'VL_y', 'VL_like'};
    csvTab.Properties.VariableNames = newVarNames;

    [ellipse_areaC{f, 1}, ellipse_centerC{f, 1}, eyeOpenPercC{f, 1}] = parse_csv_pupil_coordinates_simple(csvTab, vidReader.height, pixelToMm);

    fprintf('Pupil data in csv file #%d is parsed.\n', f);
end

centerEllipse = nanmedian(cell2mat(ellipse_centerC));
ellipse_centerC = cellfun(@(a) a-centerEllipse, ellipse_centerC, 'UniformOutput', false);

trI.hitTrials = []; % hit trial count
trI.faTrials = [];
trI.missTrials = [];
trI.crTrials = [];

% align pupil area data to task events
for t = 1:numel(tbytDat)
    %fT = selectMiddlePoints(tbytDat(t).faceCamTrel);
    %cueFrameI = cropFrames(fT);

    faceCamI = tbytDat(t).cmosExpTrainI;
    faceCamTs = evtInS.faceCam(evtInS.faceCam(:, 2)==faceCamI, 1);

    correspondingFrameI = interp1(faceCamTs, 1:length(faceCamTs), tbytDat(t).resampledVidFrameTs, 'nearest');

    if correspondingFrameI(end)<length(faceCamTs)
        tbytDat(t).pupil_areaC_cueAlign{1, 1} = tbytDat(t).resampledVidFrameTs-tbytDat(t).evtOn;
        tbytDat(t).pupil_areaC_cueAlign{1, 2} = ellipse_areaC{faceCamI}(correspondingFrameI);

        tbytDat(t).pupil_areaC_cueAlign_itp{1, 1} = cueAlignedTs;
        tbytDat(t).pupil_areaC_cueAlign_itp{1, 2} = interp1(tbytDat(t).pupil_areaC_cueAlign{1, 1}, tbytDat(t).pupil_areaC_cueAlign{1, 2}, cueAlignedTs, 'linear', 'extrap');
        
        % for normalization
        base_pupil_area = tbytDat(t).pupil_areaC_cueAlign_itp{1, 2}(cueAlignedTs<0); 
        mean_base_pupil_area = nanmean(base_pupil_area); 
        std_base_pupil_area = nanstd(base_pupil_area); 

        tbytDat(t).pupil_areaC_cueAlign_itp{1, 3} = (tbytDat(t).pupil_areaC_cueAlign_itp{1, 2}-mean_base_pupil_area)./std_base_pupil_area; % z-scored pupil size

        % hit
        if ~isempty(tbytGng(t).hitLicks) % take hit aligned 
            trI.hitTrials = [trI.hitTrials; t];
            cropLicks.hitLicks{t} = cropLick(tbytDat(t).Lick-tbytDat(t).evtOn);
            hitLickTs = cueAlignedTs-cropLicks.hitLicks{t}(1); % timestamps relative to the hit lick
            
            if max(hitLickTs)>=max(lickPeth) && min(hitLickTs)<=min(lickPeth) 
                hitLickFrameI = interp1(hitLickTs, 1:length(hitLickTs), lickPeth, 'nearest');
                
                pupil_area_lickAlign = tbytDat(t).pupil_areaC_cueAlign_itp{1, 2}(hitLickFrameI); 
                
                if sum(isnan(pupil_area_lickAlign))==0
                    tbytDat(t).pupil_areaC_hitLickAlign_itp{1, 1} = lickPeth; 
                    tbytDat(t).pupil_areaC_hitLickAlign_itp{1, 2} = interp1WithNaNs(hitLickTs(hitLickFrameI), pupil_area_lickAlign, lickPeth);    
                    tbytDat(t).pupil_areaC_hitLickAlign_itp{1, 3} = (tbytDat(t).pupil_areaC_hitLickAlign_itp{1, 2}-mean_base_pupil_area)./std_base_pupil_area; % z-scored pupil size; 
                    %tbytDat(t).pupil_areaC_hitLickAlign_itp{1, 2} = interp1(hitLickTs(hitLickFrameI), pupil_area_lickAlign, lickPeth, 'linear', 'extrap');
                end
            end
        end

        % miss (just align to cue onset)
        if tbytDat(t).rewardTrI==1 && isempty(tbytGng(t).hitLicks)
            trI.missTrials = [trI.missTrials; t];
        end

        % false alarm
        if tbytDat(t).punishTrI==1 && ~isempty(tbytGng(t).faLicks)
            trI.faTrials = [trI.faTrials; t];
            cropLicks.faLicks{t} = cropLick(tbytDat(t).Lick-tbytDat(t).evtOn);
            faLickTs = cueAlignedTs-cropLicks.faLicks{t}(1); % timestamps relative to the fa lick
            
            if max(faLickTs)>=max(lickPeth) && min(faLickTs)<=min(lickPeth)
                faLickFrameI = interp1(faLickTs, 1:length(faLickTs), lickPeth, 'nearest');
                
                pupil_area_lickAlign = tbytDat(t).pupil_areaC_cueAlign_itp{1, 2}(faLickFrameI); 

                tbytDat(t).pupil_areaC_faLickAlign_itp{1, 1} = lickPeth; 
                tbytDat(t).pupil_areaC_faLickAlign_itp{1, 2} = interp1WithNaNs(faLickTs(faLickFrameI), pupil_area_lickAlign, lickPeth);     
                tbytDat(t).pupil_areaC_faLickAlign_itp{1, 3} = (tbytDat(t).pupil_areaC_faLickAlign_itp{1, 2}-mean_base_pupil_area)./std_base_pupil_area; % z-scored pupil size;     
              
            end
        end

        % correct rejection (just align to cue onset)
        if tbytDat(t).punishTrI==1 && isempty(tbytGng(t).faLicks)
            trI.crTrials = [trI.crTrials; t];
        end
    end
    fprintf('Pupil data in trial #%d is parsed.\n', t);
end

%%
save(fullfile(filePath, strcat(header, '_tbytDat_pupillometry')), 'tbytDat')

%% plot ellipse area trajectories
colorList = slanCM('pastel1', 10);
% Hit trials
[pup.meanHitPupilArea, ~, pup.semHitPupilArea] = meanstdsem(cell2mat(cellfun(@(a) a{1, 3}, {tbytDat(trI.hitTrials).pupil_areaC_hitLickAlign_itp}, 'UniformOutput', false)')); 
plotMeanSemColor(pup.meanHitPupilArea, pup.semHitPupilArea, -2:0.02:2, colorList(3, :), {'pupil size aligned to hit licks'});
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'norm_pupil_area_hit_trials.pdf'), '-dpdf', '-vector', '-bestfit')

% Hit trials (cue aligned)
[pup.meanHitPupilAreaCueAlign, ~, pup.semHitPupilAreaCueAlign] = meanstdsem(cell2mat(cellfun(@(a) a{1, 3}, {tbytDat(trI.hitTrials).pupil_areaC_cueAlign_itp}, 'UniformOutput', false)')); 
x_pupil_areaC_cueAlign_itp = tbytDat(trI.hitTrials(1)).pupil_areaC_cueAlign_itp{1}; 
plotMeanSemColor(pup.meanHitPupilAreaCueAlign, pup.semHitPupilAreaCueAlign, x_pupil_areaC_cueAlign_itp, colorList(3, :), {'pupil size cue-aligned in hit trials'});
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:6, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'norm_pupil_area_hit_trials_cueAlign.pdf'), '-dpdf', '-vector', '-bestfit')

% FA trials
[pup.meanFaPupilArea, ~, pup.semFaPupilArea] = meanstdsem(cell2mat(cellfun(@(a) a{1, 3}, {tbytDat(trI.faTrials).pupil_areaC_faLickAlign_itp}, 'UniformOutput', false)')); 
plotMeanSemColor(pup.meanFaPupilArea, pup.semFaPupilArea, -2:0.02:2, colorList(1, :), {'pupil size aligned to fa licks'});
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'norm_pupil_area_fa_trials.pdf'), '-dpdf', '-vector', '-bestfit')

% FA trials (cue aligned)
[pup.meanFaPupilAreaCueAlign, ~, pup.semFaPupilAreaCueAlign] = meanstdsem(cell2mat(cellfun(@(a) a{1, 3}, {tbytDat(trI.faTrials).pupil_areaC_cueAlign_itp}, 'UniformOutput', false)')); 
%x_pupil_areaC_cueAlign_itp = tbytDat(trI.faTrials(1)).pupil_areaC_cueAlign_itp{1}; 
plotMeanSemColor(pup.meanFaPupilAreaCueAlign, pup.semFaPupilAreaCueAlign, x_pupil_areaC_cueAlign_itp, colorList(1, :), {'pupil size cue-aligned in FA trials'});
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:6, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'norm_pupil_area_FA_trials_cueAlign.pdf'), '-dpdf', '-vector', '-bestfit')

% CR trials
[pup.meanCrPupilArea, ~, pup.semCrPupilArea] = meanstdsem(cell2mat(cellfun(@(a) a{1, 3}, {tbytDat(trI.crTrials).pupil_areaC_cueAlign_itp}, 'UniformOutput', false)'));
plotMeanSemColor(pup.meanCrPupilArea, pup.semCrPupilArea, x_pupil_areaC_cueAlign_itp, colorList(2, :), {'pupil size aligned to cr licks'});
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:6, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'norm_pupil_area_cr_trials.pdf'), '-dpdf', '-vector', '-bestfit')

% Miss trials
[pup.meanMissPupilArea, ~, pup.semMissPupilArea] = meanstdsem(cell2mat(cellfun(@(a) a{1, 3}, {tbytDat(trI.missTrials).pupil_areaC_cueAlign_itp}, 'UniformOutput', false)'));
plotMeanSemColor(pup.meanMissPupilArea, pup.semMissPupilArea, x_pupil_areaC_cueAlign_itp, colorList(4, :), {'pupil size aligned to miss licks'});
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:6, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'norm_pupil_area_miss_trials.pdf'), '-dpdf', '-vector', '-bestfit')

% Hit, Miss, FA, CR trials (cue aligned)
plotMeanSemColormap([pup.meanHitPupilAreaCueAlign; pup.meanMissPupilArea; pup.meanFaPupilAreaCueAlign; pup.meanCrPupilArea], ...
    [pup.semHitPupilAreaCueAlign; pup.semMissPupilArea; pup.semFaPupilAreaCueAlign; pup.semCrPupilArea], ...
    x_pupil_areaC_cueAlign_itp, [colorList(3, :); colorList(4, :); colorList(1, :); colorList(2, :)], {'Hit', 'Miss', 'FA', 'CR'});
xlabel('Time (s)'); xlim([-1 2]); ylabel('Pupil area'); set(gca, 'XTick', -10:1:10, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_hit_miss_fa_cr_trials.pdf'), '-dpdf', '-vector', '-bestfit')

% Hit, Miss trials (cue aligned)
plotMeanSemColormap([pup.meanHitPupilAreaCueAlign; pup.meanMissPupilArea], ...
    [pup.semHitPupilAreaCueAlign; pup.semFaPupilAreaCueAlign], ...
    x_pupil_areaC_cueAlign_itp, [colorList(3, :); colorList(4, :)], {'Hit', 'Miss'});
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -10:1:10, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_hit_miss_trials_cueAligned.pdf'), '-dpdf', '-vector', '-bestfit')

% FA, CR trials (cue aligned)
plotMeanSemColormap([pup.meanFaPupilAreaCueAlign; pup.meanCrPupilArea], ...
    [pup.semFaPupilAreaCueAlign; pup.semCrPupilArea], ...
    x_pupil_areaC_cueAlign_itp, [colorList(1, :); colorList(2, :)], {'FA', 'CR'});
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -10:1:10, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_FA_CR_trials_cueAligned.pdf'), '-dpdf', '-vector', '-bestfit')

% %% plot licks
% rasterPlotCell({cropLicks.hitLicks(trI.hitTrials), cropLicks.missLicks(missTrials), cropLicks.faLicks(faTrials), cropLicks.crLicks(crTrials)}, 1)
% print(fullfile(filePath, 'Figure', 'lick_rasters_hit_miss_fa_cr_trials_cueAligned.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% rasterHistogramPlotCell({cropLicks.hitLicks(trI.hitTrials), cropLicks.faLicks(faTrials), cropLicks.crLicks(crTrials)}, -2:0.05:2)
% ylim([0 8])
% axis tight; grid on
% print(fullfile(filePath, 'Figure', 'lick_histogram_hit_fa_cr_trials_cueAligned.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% %% plot ellipse area trajectories of example trials
% figure; plot(-2:0.05:2-0.05, full(cell2mat(ellipse_areaC_stimAlign_binned(22:31))))
% set(gca, 'TickDir', 'out');
% xlabel('Time relative to stim onset');
% ylabel('Pupil Size');
% print(fullfile(filePath, 'Figure', 'pupilSizeTraj_exampleTrials'), '-dpdf', '-vector');
% 
% %% correlation between lick count and pupil size
% ellipse_valI = cell2mat(cellfun(@(a) ~isnan(sum(a)), ellipse_areaC_stimAlign_binned, 'UniformOutput', false));
% lick_nonZeroI = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.stimLick}, 'UniformOutput', false));
% 
% preStimPupSize = cell2mat(cellfun(@(a) mean(a(1:40)), ellipse_areaC_stimAlign_binned, 'UniformOutput', false));
% postStimPupSize = cell2mat(cellfun(@(a) mean(a(41:end)), ellipse_areaC_stimAlign_binned, 'UniformOutput', false));
% periStimPupSize = cell2mat(cellfun(@mean, ellipse_areaC_stimAlign_binned, 'UniformOutput', false));
% stimLickCount = cell2mat(cellfun(@numel, {tbytDat.stimLick}, 'un', 0));
% 
% [corrRez.periStimLickRho, corrRez.periStimLickP] = corr(periStimPupSize', stimLickCount');
% [corrRez.preStimLickRho, corrRez.preStimLickP] = corr(preStimPupSize', stimLickCount');
% [corrRez.postStimLickRho, corrRez.postStimLickP] = corr(postStimPupSize', stimLickCount');
% 
% trialI =  ellipse_valI & lick_nonZeroI;
% [corrRez.periStimLickRho_ValT, corrRez.periStimLickP_ValT] = corr(periStimPupSize(trialI)', stimLickCount(trialI)');
% [corrRez.preStimLickRho_ValT, corrRez.preStimLickP_ValT] = corr(preStimPupSize(trialI)', stimLickCount(trialI)');
% [corrRez.postStimLickRho_ValT, corrRez.postStimLickP_ValT] = corr(postStimPupSize(trialI)', stimLickCount(trialI)');
% 
% [~, periStimPupSize_sortI] = sort(periStimPupSize);
% figure(1); scatter(periStimPupSize', stimLickCount', 50, "blue", "filled");
% set(gca, 'TickDir', 'out');
% xlabel('peri-stim pupil size (mm2)');
% ylabel('lick counts');
% print(fullfile(filePath, 'Figure', 'periStimPupSize_stimLickCount'), '-dpdf', '-vector');
% 
% %% PCA pupil size trajectories
% [pupil_pca.coeff, pupil_pca.score, pupil_pca.latent, pupil_pca.tsquared, pupil_pca.explained, pupil_pca.mu] = pca(full(cell2mat(ellipse_areaC_stimAlign_binned)'));
% 
% % pc1
% [~, pcScore1_sortI] = sort(pupil_pca.score(:, 1), 'descend');
% plot(full(cell2mat(ellipse_areaC_stimAlign_binned(pcScore1_sortI(1)))))
% figure(10)
% plot(pupil_pca.coeff(:, 1))
% print(fullfile(filePath, 'Figure', 'periStimPupSize_pc1'), '-dpdf', '-vector');
% 
% % pc2
% [~, pcScore2_sortI] = sort(pupil_pca.score(:, 2), 'descend');
% plot(full(cell2mat(ellipse_areaC_stimAlign_binned(pcScore2_sortI(1)))))
% figure(11)
% plot(pupil_pca.coeff(:, 2))
% print(fullfile(filePath, 'Figure', 'periStimPupSize_pc2'), '-dpdf', '-vector');
% 
% % pc3
% [~, pcScore3_sortI] = sort(pupil_pca.score(:, 3), 'descend');
% plot(full(cell2mat(ellipse_areaC_stimAlign_binned(pcScore3_sortI(1)))))
% figure(12)
% plot(pupil_pca.coeff(:, 3))
% print(fullfile(filePath, 'Figure', 'periStimPupSize_pc3'), '-dpdf', '-vector');
% 
% % pc4
% [~, pcScore3_sortI] = sort(pupil_pca.score(:, 3), 'descend');
% plot(full(cell2mat(ellipse_areaC_stimAlign_binned(pcScore3_sortI(1)))))
% 
% figure(13)
% plot(pupil_pca.coeff(:, 4))
% print(fullfile(filePath, 'Figure', 'periStimPupSize_pc4'), '-dpdf', '-vector');
% 
% % pc 5
% figure(14)
% plot(pupil_pca.coeff(:, 5))
% print(fullfile(filePath, 'Figure', 'periStimPupSize_pc5'), '-dpdf', '-vector');
% 
% %% eye closure
% eyeClosureC = cellfun(@(a) find(a<0.5), eyeOpenPercC, 'UniformOutput', false);
% figure(20)
% plot(eyeOpenPercC{185}) % an example trial with some eye closure
% 
% set(gca, 'TickDir', 'out');
% xlabel('Frames');
% ylabel('Eye open (% detected dots)');
% print(fullfile(filePath, 'Figure', 'eyeOpenPerc'), '-dpdf', '-vector');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [ellipse_area, faceCamPt, ellipse_center, eyeOpenPer] = parse_csv_pupil_coordinates(csvTable, faceCamPt, height, pixelToMmRatio)
        % match frame numbers
        if length(faceCamPt) > size(csvTable, 1)
            faceCamPt = faceCamPt(1:size(csvTable, 1));
        end

        ellipse_area = nan(size(csvTable, 1), 1);
        eyeOpenPer = nan(size(csvTable, 1), 1);
        ellipse_center = nan(size(csvTable, 1), 2);

        % fit ellipse frame by frame
        for fr = 1:size(csvTable, 1)
            tabFr = csvTable(fr, :);
            X = [];
            Y = [];

            % L coordinates (1)
            if tabFr.L_like >= 0.5
                X(end+1) = tabFr.L_x;
                Y(end+1) = height-tabFr.L_y;
            end
            % LD coordinates (2)
            if tabFr.LD_like >= 0.5
                X(end+1) = tabFr.LD_x;
                Y(end+1) = height-tabFr.LD_y;
            end
            % D coordinates (3)
            if tabFr.D_like >= 0.5
                X(end+1) = tabFr.D_x;
                Y(end+1) = height-tabFr.D_y;
            end
            % DR coordinates (4)
            if tabFr.DR_like >= 0.5
                X(end+1) = tabFr.DR_x;
                Y(end+1) = height-tabFr.DR_y;
            end
            % R coordinates (5)
            if tabFr.R_like >= 0.5
                X(end+1) = tabFr.R_x;
                Y(end+1) = height-tabFr.R_y;
            end
            % RV coordinates (6)
            if tabFr.RV_like >= 0.5
                X(end+1) = tabFr.RV_x;
                Y(end+1) = height-tabFr.RV_y;
            end
            % V coordinates (7)
            if tabFr.V_like >= 0.5
                X(end+1) = tabFr.V_x;
                Y(end+1) = height-tabFr.V_y;
            end
            % VL coordinates (8)
            if tabFr.VL_like >= 0.5
                X(end+1) = tabFr.VL_x;
                Y(end+1) = height-tabFr.VL_y;
            end

            eyeOpenPer(fr) = length(X)/8;

            if length(X) > 5 && length(Y) > 5
                ellipse_fit = fit_ellipse(X, Y);
                ellipse_center(fr, 1) = ellipse_fit.X0_in * pixelToMmRatio;
                ellipse_center(fr, 2) = ellipse_fit.Y0_in * pixelToMmRatio;
                ellipse_area(fr) = ellipse_fit.long_axis * ellipse_fit.short_axis * pi * pixelToMmRatio^2;
            end
        end

    end

    function plotPupilCenter(pupilCenter, sizeRatio, rgbColor)
        if nargin < 2
            sizeRatio = 1.5; % Default size ratio
        end
        if nargin < 3
            rgbColor = [0, 0, 1]; % Default color (blue) if not provided
        end

        % Number of points
        numPoints = size(pupilCenter, 1);

        % Initialize figure
        figure;
        hold on;

        % Size and opacity parameters
        initialSize = 10; % Adjust as needed
        finalSize = initialSize * sizeRatio;
        initialAlpha = 0.01;
        finalAlpha = 0.85;

        % Loop through each point
        for i = 1:numPoints
            % Calculate current size and alpha
            currentSize = initialSize + (finalSize - initialSize) * (i - 1) / (numPoints - 1);
            currentAlpha = initialAlpha + (finalAlpha - initialAlpha) * (i - 1) / (numPoints - 1);

            % Plot circle with specified color and alpha
            scatter(pupilCenter(i, 1), pupilCenter(i, 2), currentSize, 'Filled', ...
                'MarkerFaceColor', rgbColor, 'MarkerFaceAlpha', currentAlpha, ...
                'MarkerEdgeAlpha', currentAlpha);

            % Draw line to next point with the same color and alpha (if not the last point)
            if i < numPoints
                nextAlpha = initialAlpha + (finalAlpha - initialAlpha) * i / (numPoints - 1);
                line([pupilCenter(i, 1), pupilCenter(i+1, 1)], [pupilCenter(i, 2), pupilCenter(i+1, 2)], ...
                    'Color', [rgbColor, nextAlpha], 'LineWidth', 1.5);
            end
        end

        hold off;
        xlabel('X Coordinate');
        ylabel('Y Coordinate');
        title('Pupil Center Movement Over Time');
        set(gca, 'TickDir', 'out')
        % print(fullfile(filePath, 'Figure', 'pupil_center_trial_1'), '-dpdf', '-vector')

    end

end