function [blockRez] = withinBlockBehAdjust(filePath, figSavePath)
%This function analyzes within-block changes in hand trajectory and force
% trajectories. For each block, the hand trajectories and force trajectories 
% of the first and last 10 trials are plotted, and they are saved in a
% structure named 'blockRez'. To run this, 'js2p0_tbybSpkHandJsPreprocess.m'
% must be run as 'jkvt' and 'ss' in 'js2p0_tbytSpkHandJsTrjBin*.mat' should
% be loaded as inputs. 
% filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';
% figSavePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles/Figure/withinBlockAdjust'; 
% fileLabel = 'WR40_082019'; 

cd(filePath)
fileDir = dir('js2p0_tbytSpkHandJsTrjBin*');
load(fullfile(fileDir(1).folder, fileDir(1).name),'jkvt','ss')

if ~exist(figSavePath,'dir')
    mkdir(fullfile(figSavePath))
end

%% session variables 
rewardI = [jkvt(:).rewarded]; 
blNumb = [ss(:).blNumber]; 
isHTrj = ~cell2mat((cellfun(@isempty,{ss(:).hTrjB},'un',0))); 
%isJsTrj = ~cell2mat((cellfun(@isempty,{ss(:).jTrjB},'un',0))); 
handXY1med = nanmedian(cell2mat(reshape(cellfun(@(a) a(1:2,1),{ss(isHTrj).hTrjB},'un',0),[1,1,sum(isHTrj)])),3);  

%% hand & joystick X-Y trajectories drawn relative to the global reference point (handXY1med)
for b = unique(blNumb)
    trs = find(blNumb==b & isHTrj==1); % all valid trials of the current block
    % draw the first and last 10 trials of the block
    for j = 1:2 % first and last 10 trials
        if j == 1
            trials = trs(1:min(10,length(trs))); % trs(1:floor(length(trs)/2)); % trs(1:10); trials to plot
            figSaveName = sprintf('hand_js_XYtrj_first10ofB%d',b); 
        elseif j == 2
            trials = trs(max(1,end-9):end); %trs(ceil(length(trs)/2):end); % trs(end-9:end); trials to plot
            figSaveName = sprintf('hand_js_XYtrj_last10ofB%d',b);  
        end
        rwdTrialI = rewardI(trials);
        handXyzC = {ss(trials).hTrjB}; %{ss(trials).hTrjB};
        jsXyzC = {ss(trials).jTrjB}; %{ss(trials).jTrjB};
        jsXyzC = jsXyzC(~cellfun(@isempty,jsXyzC));
        
        handXY = cell2mat(reshape(cellfun(@(a) smooth2a(intm(a(1:2,:),100),0,3), handXyzC, 'un', 0),[1,1,length(handXyzC)]));
        jsXY = cell2mat(reshape(cellfun(@(a) smooth2a(intm(a(1:2,:),100),0,3), jsXyzC, 'un', 0),[1,1,length(jsXyzC)]));
        jsXY1 = nanmean([cell2mat({jkvt(trials).jsTreachPosT}),cell2mat({jkvt(trials).jsTreachPosB})],2);
        jsXY1 = jsXY1(1:2,1);
        if ismember(b,[1,2,5,6]); % left blocks
            colorTheme = 'wblue';
        else ismember(b,[3,4,7,8]); % right blocks
            colorTheme = 'wred';
        end
        
        if j == 1 
            [blockRez(b).handXYtrj{1},blockRez(b).jXY1] = plotHandJsTrjXY1(handXY,handXY1med,jsXY,jsXY1,rwdTrialI,colorTheme,figSavePath,figSaveName);
        elseif j == 2
            [blockRez(b).handXYtrj{2}] = plotHandJsTrjXY1(handXY,handXY1med,jsXY,jsXY1,rwdTrialI,colorTheme,figSavePath,figSaveName);
        end
    end
end

%% (reward trials only) hand & joystick X-Y trajectories drawn relative to the global reference point (handXY1med)
for b = unique(blNumb)
    trs = find(blNumb==b & isHTrj==1  & rewardI); % all valid trials of the current block
    % draw the first and last 10 trials of the block
    for j = 1:2 % first and last 10 trials
        if j == 1
            trials = trs(1:min(10,length(trs))); % trs(1:floor(length(trs)/2)); % trs(1:10); trials to plot
            figSaveName = sprintf('rwdHand_js_XYtrj_first10ofB%d',b); 
        elseif j == 2
            trials = trs(max(1,end-9):end); %trs(ceil(length(trs)/2):end); % trs(end-9:end); trials to plot
            figSaveName = sprintf('rwdHand_js_XYtrj_last10ofB%d',b);  
        end
        rwdTrialI = rewardI(trials);
        handXyzC = {ss(trials).hTrjB}; %{ss(trials).hTrjB};
        jsXyzC = {ss(trials).jTrjB}; %{ss(trials).jTrjB};
        jsXyzC = jsXyzC(~cellfun(@isempty,jsXyzC));
        
        handXY = cell2mat(reshape(cellfun(@(a) smooth2a(intm(a(1:2,:),100),0,3), handXyzC, 'un', 0),[1,1,length(handXyzC)]));
        jsXY = cell2mat(reshape(cellfun(@(a) smooth2a(intm(a(1:2,:),100),0,3), jsXyzC, 'un', 0),[1,1,length(jsXyzC)]));
        jsXY1 = nanmean([cell2mat({jkvt(trials).jsTreachPosT}),cell2mat({jkvt(trials).jsTreachPosB})],2);
        jsXY1 = jsXY1(1:2,1);
        if ismember(b,[1,2,5,6]); 
            colorTheme = 'wblue';
        else ismember(b,[3,4,7,8]); 
            colorTheme = 'wred';
        end
        
        if j == 1 
            [blockRez(b).rwdHandXYtrj{1},blockRez(b).rwdJXY1] = plotHandJsTrjXY1(handXY,handXY1med,jsXY,jsXY1,rwdTrialI,colorTheme,figSavePath,figSaveName);
        elseif j == 2
            [blockRez(b).rwdHandXYtrj{2}] = plotHandJsTrjXY1(handXY,handXY1med,jsXY,jsXY1,rwdTrialI,colorTheme,figSavePath,figSaveName);
        end
    end
end

%% force curves aligned to the max (negative) pull force point
for b = unique(blNumb)
    trs = find(blNumb==b & isHTrj==1); % & isJsTrj);
    % draw the first 10 trials of the block
    trials = trs(1:10); %trs(1:floor(length(trs)/2)); % trs(1:10); % the first 10 trials to plot
    forceTrjs = {ss(trials).forceB}; % binned force trajectories
    figSaveName = sprintf('forceTrajectory_b%d_first10',b);
    if ismember(b,[1,2,5,6]);
        colorTheme = 'wblue';
    else ismember(b,[3,4,7,8]);
        colorTheme = 'wred';
    end
    blockRez(b).forceTrj{1} = plotForceAlignToMaxTorque(forceTrjs, figSavePath, figSaveName, 10, colorTheme);
    % draw the last 10 trials of the block
    trials = trs(end-9:end); %trs(ceil(length(trs)/2):end); % trs(end-9:end); % the last 10 trials to plot
    forceTrjs = {ss(trials).forceB}; % binned force trajectories
    figSaveName = sprintf('forceTrajectory_b%d_last10',b);
    blockRez(b).forceTrj{2} = plotForceAlignToMaxTorque(forceTrjs, figSavePath, figSaveName, 10, colorTheme);
end

%% (reward trials only) force curves aligned to the max (negative) pull force point
for b = unique(blNumb)
    trs = find(blNumb==b & isHTrj==1 & rewardI); % & isJsTrj);
    % draw the first 10 trials of the block
    trials = trs(1:min(10,length(trs))); %trs(1:floor(length(trs)/2)); % trs(1:10); % the first 10 trials to plot
    forceTrjs = {ss(trials).forceB}; % binned force trajectories
    figSaveName = sprintf('rwdForceTrajectory_b%d_first10',b);
    if ismember(b,[1,2,5,6]);
        colorTheme = 'wblue';
    else ismember(b,[3,4,7,8]);
        colorTheme = 'wred';
    end
    blockRez(b).rwdForceTrj{1} = plotForceAlignToMaxTorque(forceTrjs, figSavePath, figSaveName, 10, colorTheme);
    % draw the last 10 trials of the block
    trials = trs(max(1,end-9):end); %trs(ceil(length(trs)/2):end); % trs(end-9:end); % the last 10 trials to plot
    forceTrjs = {ss(trials).forceB}; % binned force trajectories
    figSaveName = sprintf('rwdForceTrajectory_b%d_last10',b);
    blockRez(b).rwdForceTrj{2} = plotForceAlignToMaxTorque(forceTrjs, figSavePath, figSaveName, 10, colorTheme);
end

% save 'blockRez'
save(fullfile(fileDir(1).folder, fileDir(1).name),'blockRez','-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force trajectory aligned to the max pull (negative) force point  
function [forceTrjCollect] = plotForceAlignToMaxTorque(forceTrjC, figSavePath, figSaveName, nBackFromMaxForce, colorScheme)
% forceTrjC = {ss(trials1).forceB};
% nBackFromMaxForce = 11; 
% colorScheme = 'wblue'; 
% figSaveName = 'forceTrajectories_b2_first10'; 
nTr = length(forceTrjC); 
forceTrjCollect = nan(nTr,nBackFromMaxForce); 
figure; hold on; 
[cTheme] = TNC_CreateRBColormapJP(nTr*2,colorScheme); % color to assign across trials
c = cTheme(end-nTr+1:end,:); % pick colors from the middle ones
for jj = 1:nTr
    if ~isempty(forceTrjC{jj})
        [minV,minI]=min(forceTrjC{jj}); 
        plot(forceTrjC{jj}(max(1,minI-nBackFromMaxForce+1):minI),'color',c(jj,:))
        scatter(nBackFromMaxForce, minV, 100, 'MarkerEdgeColor', 'none','MarkerFaceColor',c(jj,:),'MarkerFaceAlpha',.8) % draw starting point hTrj
        forceTrjCollect(jj,:)=forceTrjC{jj}(max(1,minI-nBackFromMaxForce+1):minI); 
        %forceTrjCollect(jj,end-length(forceTrjC{jj}(max(1,minI-nBackFromMaxForce+1):minI))+1:end)=forceTrjC{jj}(max(1,minI-nBackFromMaxForce+1):minI); 
    end
end
plot(nanmean(forceTrjCollect),'color',cTheme(end,:),'lineWidth',3)
xlim([1 size(forceTrjCollect,2)+1])
ylim([-120 0])
hold off;
print(fullfile(figSavePath,figSaveName),'-dpdf','-bestfit','-painters')
end 

% plot 2-d (X,Y) hand and joystick trajectories relative to the global reference point (handXY1med)
function [handXYuptoP1,jXY1] = plotHandJsTrjXY1(handXY,handXY1med,jsXY,jsXY1,rwdTrialI,colorTheme,figSavePath,figSaveName)

handXYuptoP1 = cell(size(handXY,3),1); 

ntb = size(handXY,2); % the # of timeBins
ntr = size(handXY,3); % the # of trials 

if ntr>=10
    [cTheme] = TNC_CreateRBColormapJP(ntr*2,colorTheme); % color to assign across trials
    %c = cTheme(max(1,ntr-floor(ntr/2)):max(1,ntr-floor(ntr/2))+ntr-1,:); % pick colors from the middle ones
    c = cTheme(max(end-ntr+1-5,1):end-5,:); % pick colors from the middle ones
else
    [cTheme] = TNC_CreateRBColormapJP(max(ntr*2,10),colorTheme); % color to assign across trials
    %c = cTheme(max(1,ntr-floor(ntr/2)):max(1,ntr-floor(ntr/2))+ntr-1,:); % pick colors from the middle ones
    c = cTheme(max(end-ntr+1-5,1):end-5,:); % pick colors from the middle ones

end

handXY = -handXY; 
handXY1med = -handXY1med;
jsXY = -jsXY;
jsXY1 = -jsXY1; 

%medHandXY1 = nanmedian(handXY(:,1,:),3); % median starting point (reference point)

hXY = handXY-repmat(handXY1med,[1,size(handXY,2),size(handXY,3)]); % normalized by subtracting the median initial hand position 
jXY = jsXY-repmat(handXY1med,[1,size(jsXY,2),size(jsXY,3)]); % normalized by subtracting the median initial hand position ; 
jXY1 = jsXY1-handXY1med; 

% get trajectories by the proximity to the joystick target position
[~,minDjXY1h] = min(sum((hXY-repmat(jXY1,[1,size(hXY,2),size(hXY,3)])).^2),[],2);
p1h = squeeze(minDjXY1h); % nth time bin that is most close to the initial joystick position (putative pullStart point)

[~,minDjXY1j] = min(sum((jXY-repmat(jXY1,[1,size(jXY,2),size(jXY,3)])).^2),[],2);
jp1 = squeeze(minDjXY1j); % nth time bin that is most close to the initial joystick position (putative pullStart point)

figure; hold on;
%scatter(jXY(1), jXY(2), 100, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0 0 0]) % draw joystick 
fstValTr = find(~isempty(hXY));
for jj = 1:size(hXY,3) % trials to draw (just draw them all)
    if rwdTrialI(jj)
        x = hXY(1,1:p1h(jj),jj);
        y = hXY(2,1:p1h(jj),jj);
    else
        x = hXY(1,:,jj);
        y = hXY(2,:,jj);
    end
    handXYuptoP1{jj,1} = [x;y];
    % draw joystick initial position
    if jj == fstValTr
        scatter(jXY1(1), jXY1(2), 200, 'MarkerEdgeColor', 'none','MarkerFaceColor','k','MarkerFaceAlpha',.7) % draw starting point hTrj
    end
    % draw initial and endpoint hand positions
    scatter(x(1), y(1), 50, 'MarkerEdgeColor', 'none','MarkerFaceColor',c(jj,:),'MarkerFaceAlpha',.4) % draw starting point hTrj
    scatter(x(end), y(end), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor',c(jj,:),'MarkerFaceAlpha',.7) % draw last point hTrj
    plot(x,y,'color',c(jj,:),'lineWidth',2)
    %patch([x nan],[y nan],[1:length(y) nan],'FaceColor','none','EdgeColor','interp','lineWidth',2)
    %colormap(c)
    
    % FOR JOYSTICK ALSO JUST DRAW FROM THE CLOSEST POINT FROM THE KNOWN jXY
    if jj<=size(jXY,3)
        jx = jXY(1,jp1(jj):end,jj);
        jy = jXY(2,jp1(jj):end,jj);
    
        plot([jx(1) jx(end)], [jy(1) jy(end)],'k','lineWidth',2)
    end
end
xlim([-12 6])
ylim([-5 15])
pbaspect([1 1 1])
set(gca,'tickDir','out')
colormap(c); colorbar

print(fullfile(figSavePath,figSaveName),'-dpdf','-bestfit','-painters')
hold off;
end

end