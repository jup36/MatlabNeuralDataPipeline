
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI')

%% prepare color
pastel2 = slanCM('Pastel2', 10); 
%plotColorListWithNumbers(pastel2); 

pastel1 = slanCM('Pastel1', 10); 
%plotColorListWithNumbers(pastel1); 

color_L_seed = pastel1(1,:); 
color_R_seed = pastel2(4,:); 

%% prepare init position data
trjNoStimRsAlignInitXY_LTrial = cell2mat(rezCol.trjNoStimRsAlignInitXY_LTrial); 
trjNoStimRsAlignInitXY_RTrial = cell2mat(rezCol.trjNoStimRsAlignInitXY_RTrial); 

trjNoStimRsAlignInitX_LTrial = trjNoStimRsAlignInitXY_LTrial(1, :); 
trjNoStimRsAlignInitX_RTrial = trjNoStimRsAlignInitXY_RTrial(1, :);

trjNoStimRsAlignInitY_LTrial = trjNoStimRsAlignInitXY_LTrial(2, :); 
trjNoStimRsAlignInitY_RTrial = trjNoStimRsAlignInitXY_RTrial(2, :);

%% prepare reach angle data
controlRiI = cellfun(@(x) contains(x, 'ri'), rezCol.controlRsAlignRchAngTrId);
controlLeI = cellfun(@(x) contains(x, 'le'), rezCol.controlRsAlignRchAngTrId);

medianAngle = nanmean([rezCol.controlRsAlignRchAngRaw, rezCol.stimFullBtRchAngRaw]); % for median subtraction

rezCol.controlRsAlignRchAngRaw = rezCol.controlRsAlignRchAngRaw-medianAngle; % for median subtraction 

controlAngRi = rezCol.controlRsAlignRchAngRaw(controlRiI);

controlAngLe = rezCol.controlRsAlignRchAngRaw(controlLeI); 

%% corr
[corrRez.rho_initX_RA_Ri, corrRez.p_initX_RA_Ri]=corr(trjNoStimRsAlignInitX_RTrial', controlAngRi', 'rows', 'complete'); 
[corrRez.rho_initX_RA_Le, corrRez.p_initX_RA_Le]=corr(trjNoStimRsAlignInitX_LTrial', controlAngLe', 'rows', 'complete'); 
[corrRez.rho_initX_RA_LeRi, corrRez.p_initX_RA_LeRi]=corr([trjNoStimRsAlignInitX_LTrial'; trjNoStimRsAlignInitX_RTrial'], ...
                                                            [controlAngLe'; controlAngRi'], 'rows', 'complete'); 

[corrRez.rho_initY_RA_Ri, corrRez.p_initY_RA_Ri]=corr(trjNoStimRsAlignInitY_RTrial', controlAngRi', 'rows', 'complete'); 
[corrRez.rho_initY_RA_Le, corrRez.p_initY_RA_Le]=corr(trjNoStimRsAlignInitY_LTrial', controlAngLe', 'rows', 'complete'); 
[corrRez.rho_initY_RA_LeRi, corrRez.p_initY_RA_LeRi]=corr([trjNoStimRsAlignInitY_LTrial'; trjNoStimRsAlignInitY_RTrial'], ...
                                                            [controlAngLe'; controlAngRi'], 'rows', 'complete'); 

%% plot
figure; hold on; 
sR = scatter(trjNoStimRsAlignInitX_RTrial', controlAngRi', 30, color_R_seed, 'filled'); 
sR.AlphaData = ones(1, length(controlAngRi)).*0.75;
sR.MarkerFaceAlpha = 'flat';

sL = scatter(trjNoStimRsAlignInitX_LTrial', controlAngLe', 30, color_L_seed, 'filled'); 
sL.AlphaData = ones(1, length(controlAngLe)).*0.75;
sL.MarkerFaceAlpha = 'flat';
set(gca, 'TickDir', 'out')
xlim([-10 10])
ylim([-70 70])
pbaspect([1 1 1])
hold off; 

print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/corr_InitPosX_RA'), '-dpdf', '-vector')

figure; hold on; 
sRY = scatter(trjNoStimRsAlignInitY_RTrial', controlAngRi', 30, color_R_seed, 'filled'); 
sRY.AlphaData = ones(1, length(controlAngRi)).*0.75;
sRY.MarkerFaceAlpha = 'flat';

sLY = scatter(trjNoStimRsAlignInitY_LTrial', controlAngLe', 30, color_L_seed, 'filled'); 
sLY.AlphaData = ones(1, length(controlAngLe)).*0.75;
sLY.MarkerFaceAlpha = 'flat';
set(gca, 'TickDir', 'out')
xlim([-5 5])
ylim([-70 70])
pbaspect([1 1 1])
hold off; 

print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/corr_InitPosY_RA'), '-dpdf', '-vector')


