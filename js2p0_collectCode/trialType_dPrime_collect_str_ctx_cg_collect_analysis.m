%% Load data
load(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData','dPrime_CtxStr_vec_norm_collectRez'),'dPrmRez'); 
dPrmRez_ctx_str = dPrmRez; clearvars dPrmRez

load(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData','dPrime_Cg_vec_norm_collectRez'),'dPrmRez'); 
dPrmRez_cg = dPrmRez; clearvars dPrmRez 

dPrmRez = [dPrmRez_ctx_str, dPrmRez_cg]; clearvars dPrmRez_ctx_str dPrmRez_cg

%% Get colors and projection matrix
load(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData','a2dColorMap.mat'),'colormap2D') % 2d colorMap for the scatter plot
numColors = size(colormap2D,1);
rangeColor = linspace(-2,2,numColors);

% lelo>> le >> lehi >> hi >> rihi >> ri >> rilo >> lo (counter-clockwise from left-low)
q2=sqrt(2);
projMatX = [-1 -q2 -1 0  1 q2  1   0];
projMatY = [-1  0   1 q2 1  0 -1 -q2];

projMat  = [projMatX; projMatY];
timeX = -1000:20:1000;
timeXI = timeX>=-200; % & timeX<1000;

distf = @(a, b) sqrt(a^2 + b^2);

%% sort/organize representative cells of each encoding type, and draw d'-type proportion pie chart
sigI = cell2mat({dPrmRez(:).sigI}');
sigIdx = find(sigI); 
mean_2dPrj = cell2mat({dPrmRez(sigI).meanCoord_2dProj}');
mean_2dPrj(:,3) = sigIdx;

%dPrm_types = unique(max_2dPrj(~isnan(max_2dPrj(:,2)),2));
dPrm_types = mean_2dPrj(:,2);
ttNumbCell = sum(~isnan(dPrm_types));
figure(4); hold on;
xlim([-2 2]);
ylim([-2 2]);
plot([-2 2], [0 0],'k')
plot([0 0], [-2 2],'k')
set(gca,'tickDir','out')
pbaspect([1 1 1])

% sort/organize draw proportion Pie chart of d' types for both regions
for tt = 1:8 % trial types counter-clockwise from left-low
    ttI = dPrm_types==tt; % current trial type logic
    if sum(ttI)>=1
        dPrmTtC{1,tt} = sortrows(mean_2dPrj(dPrm_types==tt,:),-1);
        propol = ceil(sum(ttI)./ttNumbCell*30000);
        currPt = projMat(:,tt);
        [~,trj2dXI] = min(abs(currPt(1)-rangeColor));
        [~,trj2dYI] = min(abs(currPt(2)-rangeColor));
        thisColor = [colormap2D(trj2dXI, trj2dYI, 1), colormap2D(trj2dXI, trj2dYI, 2), colormap2D(trj2dXI, trj2dYI, 3)];
        scatter(currPt(1),currPt(2),propol, 'MarkerEdgeColor', [.5 .5 .5], 'MarkerFaceColor', thisColor, 'LineWidth', .5, 'MarkerFaceAlpha',.75) % plot the proportion of the current trial type
    end
end
hold off
% print(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData/collectFigure','dPrime_vec_norm_Ctx_Str_Cg_proportionPie'),'-dpdf','-painters')
dPrmType_distribution = cell2mat(cellfun(@(a) size(a,1), dPrmTtC, 'un', 0))./ttNumbCell;

