function [ellipse_area, ellipse_center] = parse_pupil_data(filePath) 
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623/Matfiles'; 

% get tbytDat
tbytDatPath = GrabFiles_sort_trials('tbytDat', 0, {fullfile(filePath)}); 
load(fullfile(tbytDatPath{1}), 'tbytDat')

% get the pixel to mm conversion constant
if exist(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj', 'pixelToMm.mat'), 'file')==2
    load(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj', 'pixelToMm.mat'), 'pixelToMm')
else
    pixelToMm = unitConversionPixelToMilli(tbytDat(1).vidFile); % ensure to do calibration on raw video before cropping
end

%% parse 
ellipse_areaC = cell(numel(tbytDat), 1); 
vidReader = VideoReader(fullfile(tbytDat(1).cropVidFile)); % open the first video to get some info 

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

    [pupilAreaPixel, tbytDat(t).faceCam, pupilEllipseFit, eyeOpenPerc] = parse_csv_pupil_coordinates(csvTab, tbytDat(t).faceCam, vidReader.height); 
    ellipse_areaC{t} = pupilAreaPixel.*(pixelToMm^2); 
    ellipse_center{t} = [pupilEllipseFit.X0_in, pupilEllipseFit.Y0_in]; 
    eyeOpenPercC{t} = eyeOpenPerc; 

    fprintf('Pupil data in trial #%d is parsed.\n', t);
end

% align to stim onset
x = -2:0.005:2-0.005; 
for t = 1:numel(tbytDat)
    [~, stimOnI] = min(abs(tbytDat(t).faceCam-tbytDat(t).stimOn)); 
    ellipse_areaC_stimAlign{t} = ellipse_areaC{t}(stimOnI-2/0.005:stimOnI+2/0.005-1); 
    ellipse_areaC_stimAlign_binned{t} = smooth2a(nanmean(reshape(ellipse_areaC_stimAlign{t}, [], 10), 2), 1, 0); % smooth?
end 

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
function [ellipse_area, faceCamPt, ellipse_fit, eyeOpenPer] = parse_csv_pupil_coordinates(csvTable, faceCamPt, height)
    % match frame numbers
    if length(faceCamPt) > size(csvTable, 1)
        faceCamPt = faceCamPt(1:size(csvTable, 1)); 
    end
    
    ellipse_area = nan(size(csvTable, 1), 1); 
    eyeOpenPer = nan(size(csvTable, 1), 1); 
    
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
            %ellipse_fit_status{fr} = ellipse_fit.status; 
            ellipse_area(fr) = ellipse_fit.long_axis * ellipse_fit.short_axis * pi; 
        end
    end
    
end
