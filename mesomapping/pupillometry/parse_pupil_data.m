function [ellipse_area, ellipse_center] = parse_pupil_data(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623';
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';

fR = 200; % the faceCam frame rate: 200 Hz

% get tbytDat
tbytDatPath = GrabFiles_sort_trials('tbytDat', 0, {fullfile(filePath, 'Matfiles')});
load(fullfile(tbytDatPath{1}), 'tbytDat')

tbytGngDatPath = GrabFiles_sort_trials('tbytDat_parseGng', 0, {fullfile(filePath, 'Matfiles')});
tbytGng = load(fullfile(tbytGngDatPath{1}), 'tbytDat');
tbytGng = tbytGng.('tbytDat');

% get the pixel to mm conversion constant
if exist(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj', 'pixelToMm.mat'), 'file')==2
    load(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj', 'pixelToMm.mat'), 'pixelToMm')
else
    pixelToMm = unitConversionPixelToMilli(tbytDat(1).vidFile); % ensure to do calibration on raw video before cropping
end

%% parse
ellipse_areaC = cell(numel(tbytDat), 1);
vidReader = VideoReader(fullfile(tbytDat(1).cropVidFile)); % open the first video to get some info

pethTime = [-2 2]; 
pethDur = max(pethTime)-min(pethTime); 
cropFrames = @(a) a >= min(pethTime) & a <= max(pethTime);
cropLick = @(b) b(b >= min(pethTime) & b <= max(pethTime)); 

for t = 1:numel(tbytDat)
    csvTab = readDlcCsv(tbytDat(t).csvFile);
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

    [ellipse_areaC{t, 1}, tbytDat(t).faceCam, ellipse_centerC{t, 1}, eyeOpenPercC{t, 1}] = parse_csv_pupil_coordinates(csvTab, tbytDat(t).faceCam, vidReader.height, pixelToMm);

    fprintf('Pupil data in trial #%d is parsed.\n', t);
end

centerEllipse = nanmedian(cell2mat(ellipse_centerC));
ellipse_centerC = cellfun(@(a) a-centerEllipse, ellipse_centerC, 'UniformOutput', false);

hitTrials = []; % hit trial count
faTrials = []; 
missTrials = []; 
crTrials = []; 

% align pupil area data to task events
for t = 1:numel(tbytGng)
    fT = selectMiddlePoints(tbytGng(t).faceCamTrel);
    cueFrameI = cropFrames(fT);
    if ~isempty(fT)
        % hit
        if ~isempty(tbytGng(t).hitLicks)
            hitTrials = [hitTrials; t]; 
            cropLicks.hitLicks{t} = cropLick(tbytGng(t).Lick-tbytGng(t).stimOn); 
            tempPeth = fT-tbytGng(t).hitLicks(1);
            if min(tempPeth) <= min(pethTime) && max(tempPeth) >= max(pethTime)
                pethFrameI = cropFrames(tempPeth);
                if sum(pethFrameI) < fR*pethDur
                    pethIshort = fR*pethDur-sum(pethFrameI); 
                    pethIsubStart = find(pethFrameI, 1, 'last')+1; 
                    pethFrameI(pethIsubStart:pethIsubStart+pethIshort-1) = true; 
                end
                if sum(isnan(ellipse_areaC{t, 1}(pethFrameI)))==0
                    pup.hitPupilArea{t,1} = ellipse_areaC{t, 1}(pethFrameI)';
                end
                if sum(cueFrameI) < fR*pethDur
                    pethCueIshort = fR*pethDur-sum(cueFrameI); 
                    pethCueIsubStart = find(cueFrameI, 1, 'last')+1; 
                    cueFrameI(pethCueIsubStart:pethCueIsubStart+pethCueIshort-1) = true;
                end
                if sum(isnan(ellipse_areaC{t, 1}(cueFrameI)))==0
                    pup.hitPupilAreaCueAlign{t,1} = ellipse_areaC{t, 1}(cueFrameI)';
                end              
            end
        end
        % miss (just align to cue onset)
        if tbytGng(t).rewardTrI==1 && isempty(tbytGng(t).hitLicks)
            missTrials = [missTrials; t]; 
            cropLicks.missLicks{t} = cropLick(tbytGng(t).Lick-tbytGng(t).stimOn); 
            if min(fT) <= min(pethTime) && max(fT) >= max(pethTime)
                pethFrameI = cropFrames(fT);
                if sum(pethFrameI) < fR*pethDur
                    pethIshort = fR*pethDur-sum(pethFrameI);
                    pethIsubStart = find(pethFrameI, 1, 'last')+1;
                    pethFrameI(pethIsubStart:pethIsubStart+pethIshort-1) = true;
                end
                pup.missPupilArea{t,1} = ellipse_areaC{t, 1}(pethFrameI)';
            end
        end

        % false alarm
        if tbytGng(t).punishTrI==1 && ~isempty(tbytGng(t).faLicks)
            faTrials = [faTrials; t]; 
            cropLicks.faLicks{t} = cropLick(tbytGng(t).Lick-tbytGng(t).stimOn); 
            tempPeth = fT-tbytGng(t).faLicks(1);
            if min(tempPeth) <= min(pethTime) && max(tempPeth) >= max(pethTime)
                pethFrameI = cropFrames(tempPeth);
                if sum(pethFrameI) < fR*pethDur
                    pethIshort = fR*pethDur-sum(pethFrameI); 
                    pethIsubStart = find(pethFrameI, 1, 'last')+1; 
                    pethFrameI(pethIsubStart:pethIsubStart+pethIshort-1) = true; 
                end
                %if sum(isnan(ellipse_areaC{t, 1}(pethFrameI)))==0
                    pup.faPupilArea{t,1} = ellipse_areaC{t, 1}(pethFrameI)';
                %end
                if sum(cueFrameI) < fR*pethDur
                    pethCueIshort = fR*pethDur-sum(cueFrameI); 
                    pethCueIsubStart = find(cueFrameI, 1, 'last')+1; 
                    cueFrameI(pethCueIsubStart:pethCueIsubStart+pethCueIshort-1) = true;
                end
                %if sum(isnan(ellipse_areaC{t, 1}(cueFrameI)))==0
                    pup.faPupilAreaCueAlign{t,1} = ellipse_areaC{t, 1}(cueFrameI)';
                %end  
            end
        end

        % correct rejection (just align to cue onset)
        if tbytGng(t).punishTrI==1 && isempty(tbytGng(t).faLicks)
            crTrials = [crTrials; t]; 
            cropLicks.crLicks{t} = cropLick(tbytGng(t).Lick-tbytGng(t).stimOn); 
            if min(fT) <= min(pethTime) && max(fT) >= max(pethTime)
                pethFrameI = cropFrames(fT);
                if sum(pethFrameI) < fR*pethDur
                    pethIshort = fR*pethDur-sum(pethFrameI); 
                    pethIsubStart = find(pethFrameI, 1, 'last')+1; 
                    pethFrameI(pethIsubStart:pethIsubStart+pethIshort-1) = true; 
                end
                pup.crPupilArea{t, 1} = ellipse_areaC{t, 1}(pethFrameI)';
            end
        end
    end
end

%% plot ellipse area trajectories
% Hit trials
[pup.meanHitPupilArea, ~, pup.semHitPupilArea] = meanstdsem(cell2mat(pup.hitPupilArea)); 
plotMeanSemColor(pup.meanHitPupilArea, pup.semHitPupilArea, -2:0.005:1.995, [23, 245, 30]./255, {'pupil size aligned to hit licks'}); 
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_hit_trials.pdf'), '-dpdf', '-vector', '-bestfit')
% Hit trials (cue aligned)
[pup.meanHitPupilAreaCueAlign, ~, pup.semHitPupilAreaCueAlign] = meanstdsem(cell2mat(pup.hitPupilAreaCueAlign)); 
plotMeanSemColor(pup.meanHitPupilAreaCueAlign, pup.semHitPupilAreaCueAlign, -2:0.005:1.995, [23, 245, 30]./255, {'pupil size cue-aligned in hit trials'}); 
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_hit_trials_cueAlign.pdf'), '-dpdf', '-vector', '-bestfit')

% FA trials
[pup.meanFaPupilArea, ~, pup.semFaPupilArea] = meanstdsem(cell2mat(pup.faPupilArea)); 
plotMeanSemColor(pup.meanFaPupilArea, pup.semFaPupilArea, -2:0.005:1.995, [23, 245, 30]./255, {'pupil size aligned to fa licks'}); 
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_fa_trials.pdf'), '-dpdf', '-vector', '-bestfit')
% FA trials (cue aligned)
[pup.meanFaPupilAreaCueAlign, ~, pup.semFaPupilAreaCueAlign] = meanstdsem(cell2mat(pup.faPupilAreaCueAlign)); 
plotMeanSemColor(pup.meanFaPupilAreaCueAlign, pup.semFaPupilAreaCueAlign, -2:0.005:1.995, [23, 245, 30]./255, {'pupil size cue-aligned in fa trials'}); 
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_fa_trials_cueAlign.pdf'), '-dpdf', '-vector', '-bestfit')

% CR trials 
[pup.meanCrPupilArea, ~, pup.semCrPupilArea] = meanstdsem(cell2mat(pup.crPupilArea)); 
plotMeanSemColor(pup.meanCrPupilArea, pup.semCrPupilArea, -2:0.005:1.995, [23, 245, 30]./255, {'pupil size aligned to cr licks'}); 
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_cr_trials.pdf'), '-dpdf', '-vector', '-bestfit')

% Hit, FA, CR trials
colorList = slanCM('hawaii', 10); 
plotMeanSemColormap([pup.meanHitPupilArea; pup.meanFaPupilArea; pup.meanCrPupilArea], ...
                    [pup.semHitPupilArea; pup.semFaPupilArea; pup.semCrPupilArea], ...
                    -2:0.005:1.995, colorList, {'Hit', 'FA', 'CR'}); 
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_hit_fa_cr_trials.pdf'), '-dpdf', '-vector', '-bestfit')

% Hit, FA, CR trials (cue aligned)
plotMeanSemColormap([pup.meanHitPupilAreaCueAlign; pup.meanFaPupilAreaCueAlign; pup.meanCrPupilArea], ...
                    [pup.semHitPupilAreaCueAlign; pup.semFaPupilAreaCueAlign; pup.semCrPupilArea], ...
                    -2:0.005:1.995, colorList, {'Hit', 'FA', 'CR'}); 
xlabel('Time (s)'); ylabel('Pupil area'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out'); axis tight; grid on
print(fullfile(filePath, 'Figure', 'pupil_area_hit_fa_cr_trials_cueAligned.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot licks
rasterPlotCell({cropLicks.hitLicks(hitTrials), cropLicks.missLicks(missTrials), cropLicks.faLicks(faTrials), cropLicks.crLicks(crTrials)}, 1)
print(fullfile(filePath, 'Figure', 'lick_rasters_hit_miss_fa_cr_trials_cueAligned.pdf'), '-dpdf', '-vector', '-bestfit')

rasterHistogramPlotCell({cropLicks.hitLicks(hitTrials), cropLicks.faLicks(faTrials), cropLicks.crLicks(crTrials)}, -2:0.05:2)
ylim([0 8])
axis tight; grid on
print(fullfile(filePath, 'Figure', 'lick_histogram_hit_fa_cr_trials_cueAligned.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot ellipse area trajectories of example trials
figure; plot(-2:0.05:2-0.05, full(cell2mat(ellipse_areaC_stimAlign_binned(22:31))))
set(gca, 'TickDir', 'out');
xlabel('Time relative to stim onset');
ylabel('Pupil Size');
print(fullfile(filePath, 'Figure', 'pupilSizeTraj_exampleTrials'), '-dpdf', '-vector');

%% correlation between lick count and pupil size
ellipse_valI = cell2mat(cellfun(@(a) ~isnan(sum(a)), ellipse_areaC_stimAlign_binned, 'UniformOutput', false));
lick_nonZeroI = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.stimLick}, 'UniformOutput', false));

preStimPupSize = cell2mat(cellfun(@(a) mean(a(1:40)), ellipse_areaC_stimAlign_binned, 'UniformOutput', false));
postStimPupSize = cell2mat(cellfun(@(a) mean(a(41:end)), ellipse_areaC_stimAlign_binned, 'UniformOutput', false));
periStimPupSize = cell2mat(cellfun(@mean, ellipse_areaC_stimAlign_binned, 'UniformOutput', false));
stimLickCount = cell2mat(cellfun(@numel, {tbytDat.stimLick}, 'un', 0));

[corrRez.periStimLickRho, corrRez.periStimLickP] = corr(periStimPupSize', stimLickCount');
[corrRez.preStimLickRho, corrRez.preStimLickP] = corr(preStimPupSize', stimLickCount');
[corrRez.postStimLickRho, corrRez.postStimLickP] = corr(postStimPupSize', stimLickCount');

trialI =  ellipse_valI & lick_nonZeroI;
[corrRez.periStimLickRho_ValT, corrRez.periStimLickP_ValT] = corr(periStimPupSize(trialI)', stimLickCount(trialI)');
[corrRez.preStimLickRho_ValT, corrRez.preStimLickP_ValT] = corr(preStimPupSize(trialI)', stimLickCount(trialI)');
[corrRez.postStimLickRho_ValT, corrRez.postStimLickP_ValT] = corr(postStimPupSize(trialI)', stimLickCount(trialI)');

[~, periStimPupSize_sortI] = sort(periStimPupSize);
figure(1); scatter(periStimPupSize', stimLickCount', 50, "blue", "filled");
set(gca, 'TickDir', 'out');
xlabel('peri-stim pupil size (mm2)');
ylabel('lick counts');
print(fullfile(filePath, 'Figure', 'periStimPupSize_stimLickCount'), '-dpdf', '-vector');

%% PCA pupil size trajectories
[pupil_pca.coeff, pupil_pca.score, pupil_pca.latent, pupil_pca.tsquared, pupil_pca.explained, pupil_pca.mu] = pca(full(cell2mat(ellipse_areaC_stimAlign_binned)'));

% pc1
[~, pcScore1_sortI] = sort(pupil_pca.score(:, 1), 'descend');
plot(full(cell2mat(ellipse_areaC_stimAlign_binned(pcScore1_sortI(1)))))
figure(10)
plot(pupil_pca.coeff(:, 1))
print(fullfile(filePath, 'Figure', 'periStimPupSize_pc1'), '-dpdf', '-vector');

% pc2
[~, pcScore2_sortI] = sort(pupil_pca.score(:, 2), 'descend');
plot(full(cell2mat(ellipse_areaC_stimAlign_binned(pcScore2_sortI(1)))))
figure(11)
plot(pupil_pca.coeff(:, 2))
print(fullfile(filePath, 'Figure', 'periStimPupSize_pc2'), '-dpdf', '-vector');

% pc3
[~, pcScore3_sortI] = sort(pupil_pca.score(:, 3), 'descend');
plot(full(cell2mat(ellipse_areaC_stimAlign_binned(pcScore3_sortI(1)))))
figure(12)
plot(pupil_pca.coeff(:, 3))
print(fullfile(filePath, 'Figure', 'periStimPupSize_pc3'), '-dpdf', '-vector');

% pc4
[~, pcScore3_sortI] = sort(pupil_pca.score(:, 3), 'descend');
plot(full(cell2mat(ellipse_areaC_stimAlign_binned(pcScore3_sortI(1)))))

figure(13)
plot(pupil_pca.coeff(:, 4))
print(fullfile(filePath, 'Figure', 'periStimPupSize_pc4'), '-dpdf', '-vector');

% pc 5
figure(14)
plot(pupil_pca.coeff(:, 5))
print(fullfile(filePath, 'Figure', 'periStimPupSize_pc5'), '-dpdf', '-vector');

%% eye closure
eyeClosureC = cellfun(@(a) find(a<0.5), eyeOpenPercC, 'UniformOutput', false);
figure(20)
plot(eyeOpenPercC{185}) % an example trial with some eye closure

set(gca, 'TickDir', 'out');
xlabel('Frames');
ylabel('Eye open (% detected dots)');
print(fullfile(filePath, 'Figure', 'eyeOpenPerc'), '-dpdf', '-vector');

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