
filePath = {'D:\Junchol_Data\JS2p0\WR37_022119', ... % B Only
    'D:\Junchol_Data\JS2p0\WR37_022219', ... % B Only
    'D:\Junchol_Data\JS2p0\WR37_022619', ... % Cg recording contra-Cg silencing
    'D:\Junchol_Data\JS2p0\WR37_022719', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR38_052119', ... % Dual recording without silencing
    'D:\Junchol_Data\JS2p0\WR38_052219', ... % Dual recording without silencing
    'D:\Junchol_Data\JS2p0\WR38_052319', ... % Cg recording contra-Cg silencing
    'D:\Junchol_Data\JS2p0\WR38_052419', ... % Corticostriatal recording M1 silencing
    'D:\Junchol_Data\JS2p0\WR39_091019', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR39_091119', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR39_100219', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR39_100319', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR40_081919', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR40_082019', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR44_031020'};    % Dual recording with contra Cg delayed silencing

%saveName = {'WR37_022119', 'WR38_052219', 'WR38_052419', 'WR39_100219', 'WR40_081919', 'WR40_082019', 'WR44_031020'};

% completed list:
load(fullfile('D:\Junchol_Data\JS2p0\collectData','a2dColorMap.mat'),'colormap2D') % 2d colorMap for the scatter plot
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

%% plot all cells from both regions
figure(11); hold on;
xlim([-2 2]);
ylim([-2 2]);
plot([-2 2], [0 0],'k')
plot([0 0], [-2 2],'k')
set(gca,'tickDir','out')
pbaspect([1 1 1])

counter = 0;
for f = 1:length(filePath)
    cd(filePath{f})
    spkDir_STRCTX = dir('**/*binSpkCountSTRCTX*.mat');
    dPRM_Dir = dir('**/*glm_dPrime_vec_norm_WR*.mat');
    
    if ~isempty(spkDir_STRCTX) && ~isempty(dPRM_Dir)
        wr = strfind(filePath{f}, 'WR');
        saveName = filePath{f}(wr:wr+10);
        %dPrmDir = dir(fullfile(filePath{f}, strcat('glm_dPrime_*')));
        load(fullfile(dPRM_Dir.folder, dPRM_Dir.name),'dPrmC')
        %spkDir = dir(fullfile(filePath{f}, 'binSpkCountSTRCTX*'));
        load(fullfile(spkDir_STRCTX(1).folder, spkDir_STRCTX(1).name),'spkTimesCell')
        isStr = cell2mat(spkTimesCell(5,:));
        depth = cell2mat(cellfun(@(a) a(2), spkTimesCell(4,:), 'un', 0));
        assert(length(isStr)==length(dPrmC))
        %% get the maximal d' of each cell and draw
        for i_cell = 1:length(dPrmC)
            thisCellId = strcat(saveName,sprintf('_Cell#%d',i_cell));
            if ~isempty(dPrmC{i_cell})
                counter = counter+1;
                s = dPrmC{i_cell}; % structure for this cell
                dPrmRez(counter).cellId = thisCellId; % cellId with fileName and cell# as in the spkTimesCell of binSpkCountSTRCTX*
                dPrmRez(counter).isStr = isStr(i_cell); % isStr
                dPrmRez(counter).depth = depth(i_cell); % recording depth
                dPrmRez(counter).s = s; % original structure
                % find the maximal discriminative encoding (d') point
                sigI = s.trjDistZero_sigI; % index for significance, based on vector norm length
                sigI(~timeXI,1)=false; % to exclude points before reachStart
                %s.trjDistZero(sigI==0) = NaN;
                [~,dPrm_maxDistI] = max(s.trjDistZero);
                % maxCoord
                dPrmRez(counter).maxCoord = s.trj2d(dPrm_maxDistI,:); % the maximal d' point on the 2-d d' space
                [dPrmRez(counter).maxCoord_2dProj(1,1),dPrmRez(counter).maxCoord_2dProj(1,2)] = max(dPrmRez(counter).maxCoord*projMat);
                % meanCoord
                dPrmRez(counter).meanCoord = nanmean(s.trj2d(timeXI,:),1); % the maximal d' point on the 2-d d' space
                [dPrmRez(counter).meanCoord_2dProj(1,1),dPrmRez(counter).meanCoord_2dProj(1,2)] = max(dPrmRez(counter).meanCoord*projMat);
                % sigMeanCoord
                dPrmRez(counter).sigMeanCoord = nanmean(s.trj2d(sigI,:),1); % the maximal d' point on the 2-d d' space
                [dPrmRez(counter).sigMeanCoord_2dProj(1,1), dPrmRez(counter).sigMeanCoord_2dProj(1,2)] = max(dPrmRez(counter).sigMeanCoord*projMat);
                
                % get the point to draw
                %currPt =  dPrmRez(counter).maxCoord;
                currPt =  dPrmRez(counter).meanCoord;
                
                if nansum(sigI) >= 1 % sigI(dPrm_maxDistI)
                    [~,XI] = min(abs(currPt(1)-rangeColor));
                    [~,YI] = min(abs(currPt(2)-rangeColor));
                    thisColor = [colormap2D(XI, YI, 1), colormap2D(XI, YI, 2), colormap2D(XI, YI, 3)];
                    %scatter(currPt(1),currPt(2),70, 'MarkerEdgeColor', [.25 .25 .25], 'MarkerFaceColor', thisColor, 'LineWidth', .5,  'MarkerFaceAlpha',.5)
                    %if sum(~isnan(s.trjDistZero))>=1
                    dPrmRez(counter).sigI = true;
                    scatter(currPt(1),currPt(2),70, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', thisColor, ...
                        'MarkerFaceAlpha',.5)
                else
                    dPrmRez(counter).sigI = false;
                    scatter(currPt(1),currPt(2),10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.25 .25 .25], ...
                        'MarkerFaceAlpha',.75)
                end
            end
            fprintf('processed cell # %d\n', i_cell) % report unit progression
        end
        fprintf('processed file # %d\n', f) % report unit progression
    end
end
%print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure','dPrime_CtxStr_vec_norm'),'-dpdf','-painters')
%save(fullfile('D:\Junchol_Data\JS2p0\collectData','dPrime_CtxStr_vec_norm_collectRez'),'dPrmRez')
%load(fullfile('D:\Junchol_Data\JS2p0\collectData','dPrime_CtxStr_vec_norm_collectRez'),'dPrmRez')

%% draw cortex and striatum d' distribution separately
isStrMat = cell2mat({dPrmRez(:).isStr})';
dPrmRezCtx = dPrmRez(~isStrMat);
dPrmRezStr = dPrmRez(isStrMat);

% draw cortex d' scatter
figure(2); hold on;
xlim([-2 2]);
ylim([-2 2]);
plot([-2 2], [0 0],'k')
plot([0 0], [-2 2],'k')
set(gca,'tickDir','out')
pbaspect([1 1 1])
for cc = 1:length(dPrmRezCtx)
    currPt = dPrmRezCtx(cc).meanCoord;
    if dPrmRezCtx(cc).sigI %
        [~,XI] = min(abs(currPt(1)-rangeColor));
        [~,YI] = min(abs(currPt(2)-rangeColor));
        thisColor = [colormap2D(XI, YI, 1), colormap2D(XI, YI, 2), colormap2D(XI, YI, 3)];
        scatter(currPt(1),currPt(2),70, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', thisColor, ...
            'LineWidth', .5,  'MarkerFaceAlpha',.5)
    else
        scatter(currPt(1),currPt(2),10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.25 .25 .25], ...
            'MarkerFaceAlpha',.75)
    end
end
hold off;
print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure','dPrime_Ctx_scatter_vec_norm'),'-dpdf','-painters')

% draw cortex d' snapshots at different time points
for tp = 1:5:100
    figure(200); hold on;
    xlim([-2 2]);
    ylim([-2 2]);
    plot([-2 2], [0 0],'k')
    plot([0 0], [-2 2],'k')
    set(gca,'tickDir','out')
    pbaspect([1 1 1])
    
    for cc = 1:length(dPrmRezCtx)
        currPt = dPrmRezCtx(cc).s.trj2d(tp,:);
        [~,XI] = min(abs(currPt(1)-rangeColor));
        [~,YI] = min(abs(currPt(2)-rangeColor));
        thisColor = [colormap2D(XI, YI, 1), colormap2D(XI, YI, 2), colormap2D(XI, YI, 3)];        
        if dPrmRezCtx(cc).s.trjDistZero_sigI(tp)==1
            scatter(currPt(1),currPt(2),50, 'MarkerEdgeColor', [.25 .25 .25], 'MarkerFaceColor', thisColor, 'LineWidth', .5,  'MarkerFaceAlpha',.5)
        else
            scatter(currPt(1),currPt(2),25, 'MarkerEdgeColor', [.25 .25 .25], 'MarkerFaceColor', [.25 .25 .25], 'LineWidth', .5,  'MarkerFaceAlpha',.5)
        end
    end
    hText1 = text(0.90,-1.75,strcat(sprintf('%d',timeX(tp)),'ms'),'FontSize',18);
    print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure\dPrm_snapShot',strcat('dPrime_Ctx_scatter_snapShots_',sprintf('%d',timeX(tp)),'ms')),'-dpdf','-painters')
    close all
end
hold off;

% draw striatum d' scatter
figure(3); hold on;
xlim([-2 2]);
ylim([-2 2]);
plot([-2 2], [0 0],'k')
plot([0 0], [-2 2],'k')
set(gca,'tickDir','out')
pbaspect([1 1 1])
for cc = 1:length(dPrmRezStr)
    currPt = dPrmRezStr(cc).meanCoord;
    if dPrmRezStr(cc).sigI %
        [~,XI] = min(abs(currPt(1)-rangeColor));
        [~,YI] = min(abs(currPt(2)-rangeColor));
        thisColor = [colormap2D(XI, YI, 1), colormap2D(XI, YI, 2), colormap2D(XI, YI, 3)];
        scatter(currPt(1),currPt(2),70, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', thisColor, ...
            'LineWidth', .5,  'MarkerFaceAlpha',.5)
    else
        scatter(currPt(1),currPt(2),25, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.25 .25 .25], ...
            'MarkerFaceAlpha',.5)
    end
end
hold off;
print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure','dPrime_Str_scatter_vec_norm'),'-dpdf','-painters')

% draw striatum d' snapshots at different time points
for tp = 1:5:100
    figure(200); hold on;
    xlim([-2 2]);
    ylim([-2 2]);
    plot([-2 2], [0 0],'k')
    plot([0 0], [-2 2],'k')
    set(gca,'tickDir','out')
    pbaspect([1 1 1])
    
    for cc = 1:length(dPrmRezStr)
        currPt = dPrmRezStr(cc).s.trj2d(tp,:);
        [~,XI] = min(abs(currPt(1)-rangeColor));
        [~,YI] = min(abs(currPt(2)-rangeColor));
        thisColor = [colormap2D(XI, YI, 1), colormap2D(XI, YI, 2), colormap2D(XI, YI, 3)];
        
        if dPrmRezStr(cc).s.trjDistZero_sigI(tp)==1
            scatter(currPt(1),currPt(2),50, 'MarkerEdgeColor', [.25 .25 .25], 'MarkerFaceColor', thisColor, 'LineWidth', .5,  'MarkerFaceAlpha',.5)
        else
            scatter(currPt(1),currPt(2),25, 'MarkerEdgeColor', [.25 .25 .25], 'MarkerFaceColor', [.25 .25 .25], 'LineWidth', .5,  'MarkerFaceAlpha',.5)
        end
    end
    hText1 = text(0.90,-1.75,strcat(sprintf('%d',timeX(tp)),'ms'),'FontSize',18);
    print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure\dPrm_snapShot',strcat('dPrime_Str_scatter_snapShots_',sprintf('%d',timeX(tp)),'ms')),'-dpdf','-painters')
    close all
end
hold off;

%% sort/organize representative cells of each encoding type, and draw d'-type proportion pie chart
sigI = cell2mat({dPrmRez(:).sigI}');
%max_2dPrj = cell2mat({dPrmRez(sig_max).maxCoord_2dProj}');
%max_2dPrj(:,3) = 1:size(max_2dPrj,1);
mean_2dPrj = cell2mat({dPrmRez(sigI).meanCoord_2dProj}');
mean_2dPrj(:,3) = 1:size(mean_2dPrj,1);
isStr_collect = cell2mat({dPrmRez(sigI).isStr}');

%dPrm_types = unique(max_2dPrj(~isnan(max_2dPrj(:,2)),2));
dPrm_types = mean_2dPrj(:,2);
ttNumbCell = sum(~isnan(dPrm_types));
ttNumbStr = sum(isStr_collect & ~isnan(dPrm_types));
ttNumbCtx = sum(~isStr_collect & ~isnan(dPrm_types));
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
print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure','dPrime_vec_norm_CtxStr_proportionPie'),'-dpdf','-painters')
dPrmType_CtxStr_distribution = cell2mat(cellfun(@(a) size(a,1), dPrmTtC, 'un', 0))./ttNumbCell;

figure(5); hold on;
xlim([-2 2]);
ylim([-2 2]);
plot([-2 2], [0 0],'k')
plot([0 0], [-2 2],'k')
set(gca,'tickDir','out')
pbaspect([1 1 1])

% sort/organize draw proportion Pie chart of d' types for MCtx
for tt = 1:8 % trial types counter-clockwise from left-low
    ttI = dPrm_types==tt & isStr_collect==0; % current trial type logic
    if sum(ttI)>=1
        dPrmTtC_ctx{1,tt} = sortrows(mean_2dPrj(ttI,:),-1);
        propol = ceil(sum(ttI)./ttNumbCtx*30000);
        currPt = projMat(:,tt);
        [~,trj2dXI] = min(abs(currPt(1)-rangeColor));
        [~,trj2dYI] = min(abs(currPt(2)-rangeColor));
        thisColor = [colormap2D(trj2dXI, trj2dYI, 1), colormap2D(trj2dXI, trj2dYI, 2), colormap2D(trj2dXI, trj2dYI, 3)];
        scatter(currPt(1),currPt(2),propol, 'MarkerEdgeColor', [.5 .5 .5], 'MarkerFaceColor', thisColor, 'LineWidth', .5, 'MarkerFaceAlpha',.75) % plot the proportion of the current trial type
    end
end
hold off
%print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure','dPrime_Ctx_vec_norm_proportionPie'),'-dpdf','-painters')
dPrmType_Ctx_distribution = cell2mat(cellfun(@(a) size(a,1), dPrmTtC_ctx, 'un', 0))./ttNumbCtx;


figure(6); hold on;
xlim([-2 2]);
ylim([-2 2]);
plot([-2 2], [0 0],'k')
plot([0 0], [-2 2],'k')
set(gca,'tickDir','out')
pbaspect([1 1 1])
% sort/organize draw proportion Pie chart of d' types for Str
for tt = 1:8 % trial types counter-clockwise from left-low
    ttI = dPrm_types==tt & isStr_collect==1; % current trial type logic
    if sum(ttI)>=1
        dPrmTtC_str{1,tt} = sortrows(mean_2dPrj(ttI,:),-1);
        propol = ceil(sum(ttI)./ttNumbStr*30000);
        currPt = projMat(:,tt);
        [~,trj2dXI] = min(abs(currPt(1)-rangeColor));
        [~,trj2dYI] = min(abs(currPt(2)-rangeColor));
        thisColor = [colormap2D(trj2dXI, trj2dYI, 1), colormap2D(trj2dXI, trj2dYI, 2), colormap2D(trj2dXI, trj2dYI, 3)];
        scatter(currPt(1),currPt(2),propol, 'MarkerEdgeColor', [.5 .5 .5], 'MarkerFaceColor', thisColor, 'LineWidth', .5, 'MarkerFaceAlpha',.75) % plot the proportion of the current trial type
    end
end
hold off
%print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure','dPrime_Str_vec_norm_proportionPie'),'-dpdf','-painters')
dPrmType_Str_distribution = cell2mat(cellfun(@(a) size(a,1), dPrmTtC_str, 'un', 0))./ttNumbStr;

%save(fullfile('D:\Junchol_Data\JS2p0\collectData','dPrime_CtxStr_collectRez'),'dPrmTtC_str','dPrmTtC_ctx','dPrmTtC','-append')

%%
% figure(7); hold on;
% xlim([-2 2]);
% ylim([-2 2]);
% plot([-2 2], [0 0],'k')
% plot([0 0], [-2 2],'k')
% plot([-2 2], [-2 2],'k')
% plot([-2 2], [2 -2],'k')
% set(gca,'tickDir','out')
% pbaspect([1 1 1])

for tt = 1:8
    figNumb = 10+tt;
    figure(figNumb); hold on;
    xlim([-2 2]);
    ylim([-2 2]);
    plot([-2 2], [0 0],'k')
    plot([0 0], [-2 2],'k')
    plot([-2 2], [-2 2],'k')
    plot([-2 2], [2 -2],'k')
    set(gca,'tickDir','out')
    pbaspect([1 1 1])
    
    cellToDraw = dPrmTtC{1,tt}(1:min(5,size(dPrmTtC{1,tt},1)),3); % draw 10 trials per type
    for cc = 1:length(cellToDraw)
        ttTrj = dPrmRez(cellToDraw(cc)).s.trj2d;
        ttPval = dPrmRez(cellToDraw(cc)).s.dist2d_cdf_p;
        trjDistZero = dPrmRez(cellToDraw(cc)).s.trjDistZero;
        trjDistZero(ttPval>=0.05,1)=NaN; % drop non-significant points
        trjDistZero(~timeXI,1)=NaN; % to exclude points before reachStart
        [~,sig_maxDistI] = max(trjDistZero);
        
        for tm = 2:101 % time bins
            [~,trj2dXI] = min(abs(ttTrj(tm,1)-rangeColor));
            [~,trj2dYI] = min(abs(ttTrj(tm,2)-rangeColor));
            thisColor(tm, :) = [colormap2D(trj2dXI, trj2dYI, 1), colormap2D(trj2dXI, trj2dYI, 2), colormap2D(trj2dXI, trj2dYI, 3)];
            prevPt = ttTrj(tm-1,:);
            currPt = ttTrj(tm,:);
            plot([prevPt(1) currPt(1)],[prevPt(2) currPt(2)],'Color', thisColor(tm, :),'LineWidth',.5)
            if tm == 2
                scatter(currPt(1),currPt(2),70, 'MarkerEdgeColor', thisColor(tm, :), 'MarkerFaceColor', 'none', 'LineWidth', .5);
                %elseif tm == 101
                %    scatter(currPt(1),currPt(2),70, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', thisColor(tm, :), 'LineWidth', .5);
            elseif tm == sig_maxDistI
                scatter(currPt(1),currPt(2),70, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', thisColor(tm, :), 'LineWidth', .5);
            end
            % each point
            %if ttPval(tm,1)<0.01
            %    scatter(currPt(1),currPt(2),10, 'MarkerEdgeColor', [.5 .5 .5], 'MarkerFaceColor', thisColor(tm, :), 'LineWidth', .5);
            %end
        end
    end
    print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure',strcat('dPrime_2dTrj_',sprintf('trialType#%d',tt))),'-dpdf','-painters')
    close all
end

%% create example individual neuronal d' trajectories of each trial type
for tt = 1:8
    if ~isempty(dPrmTtC{tt})
        thisCell = dPrmTtC{tt}(1,3); % representative cell of this dPrm type
        ttTrj1{1,tt}  = dPrmRez(thisCell).s.trj2d;
        ttPval1{1,tt} = dPrmRez(thisCell).s.dist2d_cdf_p;
    end
end

h1=figure; hold on;
xlim([-2 2]);
ylim([-2 2]);
set(gca,'tickDir','out')
pbaspect([1 1 1])
for tm = 2:101 % time bins
    for tt = 1:8
        if ~isempty(ttTrj1{1,tt})
            [~,trj2dXI] = min(abs(ttTrj1{1,tt}(tm,1)-rangeColor));
            [~,trj2dYI] = min(abs(ttTrj1{1,tt}(tm,2)-rangeColor));
            thisColor(tm, :) = [colormap2D(trj2dXI, trj2dYI, 1), colormap2D(trj2dXI, trj2dYI, 2), colormap2D(trj2dXI, trj2dYI, 3)];
            prevPt = ttTrj1{1,tt}(tm-1,:);
            currPt = ttTrj1{1,tt}(tm,:);
            plot([prevPt(1) currPt(1)],[prevPt(2) currPt(2)],'Color', thisColor(tm, :),'LineWidth',.5)
            % each point
            if ttPval1{1,tt}(tm,1)<0.01
                scatter(currPt(1),currPt(2),70, 'MarkerEdgeColor', [.5 .5 .5], 'MarkerFaceColor', thisColor(tm, :), 'LineWidth', .5);
            end
        end
    end
    hText = text(0.9,-1.7,sprintf('time=%d ms',timeX(tm)),'FontSize',14);
    drawnow
    F(tm-1) = getframe(h1);
    delete(hText)
end

video = VideoWriter(fullfile('/Volumes/8TB/Junchol_Data/JS2p0/collectData','exampleTrj_eachTrialType'),'MPEG-4');
video.Quality = 100;
video.FrameRate = 20;
open(video)
writeVideo(video,F)
close(video)
hold off;

%% mean, sem dPrm
for tt = 1:length(projMatX)
    [~,colorXI] = min(abs(projMat(1,tt)-rangeColor));
    [~,colorYI] = min(abs(projMat(2,tt)-rangeColor));
    colorC{tt} = [colormap2D(colorXI, colorYI, 1), colormap2D(colorXI, colorYI, 2), colormap2D(colorXI, colorYI, 3)];
end


for cc = 1:length(dPrmRez)
    tType = dPrmRez(cc).maxCoord_2dProj(2);
    
    if ~isempty(dPrmRez(cc).s) && ~isnan(tType)
        dPrmTrjCell{cc,tType} = (dPrmRez(cc).s.trj2d*projMat(:,tType))';
    end
end

figure; hold on;
for tt = 1:8
    [mdPrm,~,sdPrm] = meanstdsem(cell2mat(dPrmTrjCell(:,tt)));
    plot(mdPrm./max(mdPrm), 'Color', colorC{tt},'LineWidth',2)
end

%%
for tt = 1:8
    [mdPrm_str(tt,:),~,sdPrm_str(tt,:)] = meanstdsem(cell2mat(dPrmTrjCell([dPrmRez(:).isStr]',tt)));
end
% plot pure load
figure; hold on;
%plot(mean(mdPrm_str([4,8],:))./max(mean(mdPrm_str([4,8],:))))
%plot(mean(mdPrm_str([4,8],:)))
plot(mdPrm_str(8,:)./max(mdPrm_str(8,:)))
plot(mdPrm_str(4,:)./max(mdPrm_str(4,:)))

% plot pure direction
%plot(mean(mdPrm([2,6],:))./max(mean(mdPrm([2,6],:))))
plot(mean(mdPrm_str([2,6],:)))

% plot pure mixed
%plot(mean(mdPrm([1,3,5,7],:))./max(mean(mdPrm([1,3,5,7],:))))
plot(mean(mdPrm_str([1,3,5,7],:)))
axis tight
set(gca,'tickDir','out')
print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure','dPrm_neuronType_Str'),'-dpdf','-painters')

%%
for tt = 1:8
    [mdPrm_ctx(tt,:),~,sdPrm_ctx(tt,:)] = meanstdsem(cell2mat(dPrmTrjCell(~[dPrmRez(:).isStr]',tt)));
end
% plot pure load
figure; hold on;
%plot(mean(mdPrm_ctx([4,8],:))./max(mean(mdPrm_ctx([4,8],:))))
plot(mdPrm_ctx(8,:)./max(mdPrm_ctx(8,:)))
plot(mdPrm_ctx(4,:)./max(mdPrm_ctx(4,:)))

% plot pure direction
%plot(mean(mdPrm([2,6],:))./max(mean(mdPrm([2,6],:))))
plot(mean(mdPrm_ctx([2,6],:)))

% plot pure mixed
%plot(mean(mdPrm([1,3,5,7],:))./max(mean(mdPrm([1,3,5,7],:))))
plot(mean(mdPrm_ctx([1,3,5,7],:)))
axis tight
set(gca,'tickDir','out')
print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure','dPrm_neuronType_Ctx'),'-dpdf','-painters')

