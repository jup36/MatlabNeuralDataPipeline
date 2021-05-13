function getHandJsKinematics( filePath ) 
%filePath = 'S:\Junchol_Data\JS2p0\WR40_081419';
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081419';

%This script loads 3-d coordinates of the reaching hand and the joystick saved in 'trj3d.mat', 
% and extracts kinematics of the hand and joystick movements. The process
% comprises 1) defining baseline positions, 2) getting baseline-subtracted
% trajectories and velocities. 'trj3d.mat' has two sets of trajectories
% after stereo-triangulation with reference to the front or side cameras. 
% Joystick trajectories from the front camera tend to be more reliable, thus
% trajectories relative to front camera are used currently. 

%% load 'trj3d.mat'
cd(filePath)
trj3d = dir('**/*trj3d.mat');
T=load(fullfile(trj3d.folder,trj3d.name),'trj3d');
clearvars trj3d
T = T.('trj3d');

hTrjF = {T(:).allPartsMedSgFronXYZ}; % all parts median and Savitzky-Golay filtered hand trajectory front cam
%hTrjS = {T(:).allPartsMedSgSideXYZ}; % all parts median and Savitzky-Golay filtered hand trajectory side cam
jsBot = {T(:).jsBSgFronXYZ}; 
jsTop = {T(:).jsTSgFronXYZ}; 

%% collect jsTime1k_Kinematics data
jsKinFile = dir(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat'));

if length(jsKinFile)==1
    jkvt = load(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat'),'jsTime1k_KV'); % load jsTime1k_KV
else
    disp('Point to "jsTime1k_Kinematics_VideoFiles.mat" with the variable "jsTime1k_KV"!')
    [jsKinFileSelect,jsKinPathSelect] = uigetfile(filePath);
    jkvt = load(fullfile(jsKinPathSelect,jsKinFileSelect),'jsTime1k_KV'); % load jsTime1k_KV
end

jkvt = jkvt.('jsTime1k_KV');
clearvars jsTime1k_KV

% sanity check if the number of trials match
if size(jkvt,2)~=size(T,2)
    error('the trial numbers differ between jsTime1k_KV and trj3d!')
end

%% get the 3d-position/velocity
gmhTrj = nanmedian(cell2mat(hTrjF),2); % global median

% get trajectories relative to the global median first, the preparatory period median can be defined after isolating reach starts
medDst = @(x) sqrt(sum((x-gmhTrj).^2)); % func to get the point-by-point Dstance from the global median position

for t = 1:length(hTrjF)
    if ~isempty(hTrjF{t})
        jkvt(t).hTrjF = hTrjF{t}; 
        jkvt(t).hTrjDstMed = splitapply(medDst, hTrjF{t}, 1:size(hTrjF{t},2)); % point-by-point position relative to the global median
        jkvt(t).hTrjVelMed = diff([jkvt(t).hTrjDstMed(:,1) jkvt(t).hTrjDstMed],1,2)./(4/1000);
        % find the baseline points of low velocity and short-distance from the global median
        [~,sortIhTrjDstMed]=sort(jkvt(t).hTrjDstMed); % sort by distance from the global median
        nearMedLogic = zeros(1,length(jkvt(t).hTrjDstMed));
        nearMedLogic(sortIhTrjDstMed(1:min(250, length(sortIhTrjDstMed))))=1; % taking points not too far from the global median excludes low-velocity points such as time spent with the limb on the JS without moving          
        lowVelLogic = abs(jkvt(t).hTrjVelMed)<5; % low-velocity points 
        
        if sum(nearMedLogic & lowVelLogic)>10
            basePosM = nanmedian(hTrjF{t}(:,nearMedLogic & lowVelLogic), 2); % low velocity (<5 mm/s) points
        else
            basePosM = gmhTrj; % just take the global median points
        end
        baseDst = @(x) sqrt(sum((x-basePosM).^2)); % func to get the point-by-point Dstance from the global median position
        jkvt(t).hTrjDstBase = splitapply(baseDst, hTrjF{t}, 1:size(hTrjF{t},2)); % point-by-point position relative to trial-by-trial baseline
        jkvt(t).hTrjVelBase = diff([jkvt(t).hTrjDstBase(:,1) jkvt(t).hTrjDstBase],1,2)./(4/1000); % point-by-point velocity
        jkvt(t).hTrjDstBaseXyz = hTrjF{t}-basePosM; % normalize by subtracting the median of the baseline position X, Y, Z
        jkvt(t).hTrjVelBaseXyz = diff([jkvt(t).hTrjDstBaseXyz(:,1) jkvt(t).hTrjDstBaseXyz],1,2);
        %elseif isempty(hTrjF{t})
        %    jkvt(t).hTrjDstMed = []; jkvt(t).hTrjVelMed = []; jkvt(t).hTrjDstBase = []; jkvt(t).hTrjVelBase = []; jkvt(t).hTrjDstBaseXyz = []; jkvt(t).hTrjVelBaseXyz = [];
    end
end
clearvars t

%% get the hand velocity and position cut to detect reaches
vel = [jkvt(:).hTrjVelBase]; % baseline normalized 1-d velocity across all time points
sortVel=sort(vel); % sort velocity data
sortVel=sortVel(~isnan(sortVel)); 
velCut = sortVel(round(length(sortVel)*.98)); % set the 98% cutoff for the velocity data (velocity cutoff is required, since reaches are expected to occur at high velocity)

pos = [jkvt(:).hTrjDstBase]; % 1-d position relative to the trial-by-trial baseline
sortPos = sort(pos); % sort position data
sortPos(sortPos<nanmedian(sortPos))=nanmedian(sortPos); % replace less-than median position values with the median position value
posCut = 10; % set the Distance cut off as 10 mm from the baseline sortPos(round(length(sortPos)*.90));

for t = 1:size(jkvt,2)
    [jkvt(t).hTrjRstart, jkvt(t).hTrjRstop] = getReachTimesJs2p0(jkvt(t).hTrjDstBase, jkvt(t).hTrjVelBase, posCut, velCut);
    % posTrace= jkvt(t).hTrjDstBase; velTrace= jkvt(t).hTrjVelBase;
end
clearvars t

%% get the joystick trajectory
% get the valid portions of the joystick trajectories
for tr = 1:length(jsBot)
    if length(jsBot{1,tr})==length(jsTop{1,tr}) % use joystick trajectories from front cameras (more reliable)
        jsValFst = max(T(tr).jsfstPtFB, T(tr).jsfstPtFT);
        frValFst = find(T(tr).vUseFrameIdx,1,'first');
        if jsValFst>frValFst % in most cases, js valid pts are after the first valid frame
            jkvt(tr).jsTrjValB = T(tr).jsBSgFronXYZ(:,jsValFst:end);
            jkvt(tr).jsTrjValT = T(tr).jsTSgFronXYZ(:,jsValFst:end);
        elseif jsValFst<frValFst % in some cases, when corresponding frames started after the js positioining 
            jkvt(tr).jsTrjValB = T(tr).jsBSgFronXYZ;
            jkvt(tr).jsTrjValT = T(tr).jsTSgFronXYZ;
        end  
    end
end
clearvars tr

% the global median joystick coordinates to be used as reference points
gmJsTrjT = nanmedian(cell2mat({jkvt(:).jsTrjValT}),2); % the global median joystick top
gmJsTrjB = nanmedian(cell2mat({jkvt(:).jsTrjValB}),2); % the global median joystick bottom

medDstJsT = @(x) sqrt(sum((x-gmJsTrjT).^2)); % func to get the point-by-point Dstance from the global median position
medDstJsB = @(x) sqrt(sum((x-gmJsTrjB).^2)); % func to get the point-by-point Dstance from the global median position
for tr = 1:size(jkvt,2) 
    % get the distance from the global median 
    if ~isempty(jkvt(tr).jsTrjValT) && ~isempty(jkvt(tr).jsTrjValB)
        % distance
        jkvt(tr).jsDstValGmT = splitapply(medDstJsT, jkvt(tr).jsTrjValT, 1:size(jkvt(tr).jsTrjValT,2)); % point-by-point position relative to the global median; 
        jkvt(tr).jsDstValGmB = splitapply(medDstJsB, jkvt(tr).jsTrjValB, 1:size(jkvt(tr).jsTrjValB,2));
        % velocity
        jkvt(tr).jsVelValGmT = diff([jkvt(tr).jsDstValGmT(:,1) jkvt(tr).jsDstValGmT],1,2)./(4/1000);
        jkvt(tr).jsVelValGmB = diff([jkvt(tr).jsDstValGmB(:,1) jkvt(tr).jsDstValGmB],1,2)./(4/1000);
    end
end
clearvars tr

jsTargetPos = unique([jkvt(:).reachP1]); 
for rp = 1:length(jsTargetPos)
    tmpJsPosIdx = [jkvt(:).reachP1]==jsTargetPos(rp); 
    
    % bind all position/velocity data corresponding to the current joystick position (joystick top) 
    tmpJsTrjValT = [jkvt(tmpJsPosIdx).jsTrjValT]; % Js top Trj
    tmpJsVelValT = [jkvt(tmpJsPosIdx).jsVelValGmT]; % Js top Vel
    tmpJsTrjValTlowVelMedT = nanmedian(tmpJsTrjValT(:,abs(tmpJsVelValT)<1),2); % Js top low-vel points median
    [jkvt(tmpJsPosIdx).jsTreachPosT] = deal(tmpJsTrjValTlowVelMedT); 
    
    % bind all position/velocity data corresponding to the current joystick position (joystick bottom)
    tmpJsTrjValB = [jkvt(tmpJsPosIdx).jsTrjValB]; % Js bottom Trj
    tmpJsVelValB = [jkvt(tmpJsPosIdx).jsVelValGmB]; % Js bottom Vel
    tmpJsTrjValTlowVelMedB = nanmedian(tmpJsTrjValB(:,abs(tmpJsVelValB)<1),2); % Js bottom low-vel points median
    [jkvt(tmpJsPosIdx).jsTreachPosB] = deal(tmpJsTrjValTlowVelMedB); 
    
end
clearvars rp tmp*

save(fullfile(filePath,'jsTime1k_KinematicsTrajectories.mat'), 'jkvt'); % save the raw pixel values
end

function [reachStart, reachStop] = getReachTimesJs2p0( posTrace, velTrace, posCut, velCut )

%posTrace = hTrjSDstBase{7}; 
%velTrace = hTrjSVelBase{7}; 

sampRate = 250; % the sampling rate of videography (default: 250Hz)
timeInt = 1000/sampRate; % time interval in ms 
maxDur = 1500; 
fastPts = find(velTrace>velCut); % points exceeding the velocity cutoff
valFastPts = []; 

%% detect reaches
for i = 1:length(fastPts) % all points exceeding the velocity cutoff
    if isempty(valFastPts) % in case of the 1st high-velocity data point
        if sum(posTrace(fastPts(i):min(fastPts(i)+(200/timeInt), length(posTrace)))>posCut)>(15/timeInt) % this part is just to isolate fast reaches from the short jerks
            valFastPts=[valFastPts, fastPts(i)]; 
        end
    else    % from the 2nd high-velocity data point on
        if (fastPts(i)-fastPts(i-1))>(50/timeInt) % if the current high-velocity point came 80 ms (note that the sampling rate was 250Hz) or more after the previous high-velocity point
            if (fastPts(i)-valFastPts(end)) > (50/timeInt) % ensure the interval between reaches is greater than the set interval threshold (e.g. 2000 ms)
                if sum(posTrace(fastPts(i):min(fastPts(i)+(200/timeInt), length(posTrace)))>posCut)>(15/timeInt) % finally ensure that it is a well-timed reach not a short jerk
                    valFastPts=[valFastPts, fastPts(i)];
                end
            end
        end
    end
end

%% detect reachStart
reachStart=[]; % the goal here is to find the near-zero point that is immediately prior to each high-velocity point  
for i = 1:length(valFastPts)
    currSamp = valFastPts(i); % current valid sample (high-velocity data point) 
    nearZero= find(posTrace(1:currSamp)<posCut*.4,1,'last'); % find the most recent near zero-point
    if nearZero>10
        if currSamp-nearZero>1000  % sometimes the joystick doesn't reset to zero and induces 20 second long reaches. this prevents that by setting the second (additional) posCut
            sortPos2=sort((posTrace(nearZero:currSamp)));
            posCut2 = sortPos2(round(length(sortPos2)*.90)); % reset the posCut2 to be a point immediatly before the reach peak (farthest point)
            nearZero=find(posTrace(1:currSamp)<posCut2,1,'last');
        end
        realStartMinus10=find(velTrace(1:nearZero)<0,10,'last'); % find 10 elements whose velocity is negative prior to the current near zero point  
    else % if there's no near zero point
        realStartMinus10=currSamp; % take the current valid sample (high-velocity data point) as a reach start point
    end
    reachStart = [reachStart, realStartMinus10(1)];
end
clearvars i 

if length(reachStart)>=2 
   reachStart=reachStart([1,find(diff(reachStart)>(100/timeInt))+1]); % ensure that intervals between reaches are greater than 100 ms
end
reachStart=unique(reachStart)+10; 

% detect reachStop
reachStop=[]; 
for i = 1:length(reachStart) 
    currPos = posTrace(reachStart(i):min(reachStart(i)+floor(maxDur/timeInt), length(posTrace))); % take the portion relevant to current reach
    [~,reachPeakPt]=max(currPos(1:floor(length(currPos)/2))); % detect reach peak (the farthest reach point)
    peakToEnd = currPos(reachPeakPt:end); % peak-to-end of reach  
    sortStop=sort(peakToEnd); 
    peakToEndBin = hist(peakToEnd,5); 
    [~,cutoffStop1]=max(peakToEndBin(1:round(length(peakToEndBin)/3)));
    cutoffStop = sortStop(peakToEndBin(cutoffStop1));
    
    tb = find(peakToEnd<cutoffStop); 
    xx = diff(tb)==1; 
    ff = find([false,xx]~=[xx,false]); % spot out the non-consecutive data points 
    gg = find(ff(2:2:end)-ff(1:2:end-1)>=(100/timeInt),1,'first');
    almostEnd=tb(ff(2*gg-1))-1;   % just tranlate g to a point on t  
    if isempty(almostEnd)   % in case there's no such point
        almostEnd=find(peakToEnd<cutoffStop,1);   % just take the point below the reach threshold
        %reachStart(i);
    end
    %thisStop= find(moveavg(diff(peakToEndPos(almostEnd:end)),10)>=0,1);
    reachStop=[reachStop, reachPeakPt+almostEnd+reachStart(i)]; %thisStop+reachPeakTime+almostEnd+reachStart(i)];
end
end
%% examine trajectories by plotting
%rewardTrs = find([jsTime1k_KV(:).rewarded]==1);
% close all;
% trI = 8;
% x = jkvt(trI).hTrjF(1,:); %T(jkvt(trI)).vUseFrameIdx);
% y = jkvt(trI).hTrjF(2,:); %T(jkvt(trI)).vUseFrameIdx);
% z = jkvt(trI).hTrjF(3,:); %T(jkvt(trI)).vUseFrameIdx);
% c = 1:sum(jkvt(trI).vUseFrameIdx); % generate a colormap;
% 
% jsAtTarget = [jkvt(trI).jsTreachPosT jkvt(trI).jsTreachPosB]; % joystick at the reach target position
% jsX = jsAtTarget(1,:); % Js at the reach target position X
% jsY = jsAtTarget(2,:); % Js at the reach target position Y
% jsZ = jsAtTarget(3,:); % Js at the reach target position Z
% 
% figure; hold on; 
% %plot3(gmhTrj(1), gmhTrj(2), abs(gmhTrj(3)), 'ro')
% patch([x nan],[y nan],[z nan],[c nan],'FaceColor','none','EdgeColor','interp')
% %patch([jsX nan],[jsY nan],[jsZ nan],'FaceColor','k')
% plot3(jsX, jsY, jsZ, 'k', 'lineWidth', 3, 'alpha', 0.5)
% 
% %patch([x nan],[z nan],-[y nan],[c nan],'FaceColor','none','EdgeColor','interp')
% xlabel('X(mm)')
% ylabel('Y(mm)')
% zlabel('Z(mm)')
% 
% %patch([x nan],[z nan],abs([y nan]),[c nan],'FaceColor','none','EdgeColor','interp')
% colormap cool
% colorbar
% caxis([1 sum(jkvt(trI).vUseFrameIdx)])
% print(fullfile(filePath,'Figure',sprintf('hTrjDstBaseXyz_tr#%d',trI)),'-dpdf','-painters','-bestfit');
% 

% close all;
% trI = 10;
% x = jkvt(trI).hTrjDstBaseXyz(1,jkvt(10).vUseFrameIdx); %T(jkvt(trI)).vUseFrameIdx);
% y = jkvt(trI).hTrjDstBaseXyz(2,jkvt(10).vUseFrameIdx); %T(jkvt(trI)).vUseFrameIdx);
% z = jkvt(trI).hTrjDstBaseXyz(3,jkvt(10).vUseFrameIdx); %T(jkvt(trI)).vUseFrameIdx);
% c = 1:sum(jkvt(trI).vUseFrameIdx); % generate a colormap;
% 
% figure;
% %plot3(gmhTrj(1), gmhTrj(2), abs(gmhTrj(3)), 'ro')
% patch([x nan],[y nan],[z nan],[c nan],'FaceColor','none','EdgeColor','interp')
% %patch([x nan],[z nan],-[y nan],[c nan],'FaceColor','none','EdgeColor','interp')
% xlabel('X(mm)')
% ylabel('Y(mm)')
% zlabel('Z(mm)')
% 
% %patch([x nan],[z nan],abs([y nan]),[c nan],'FaceColor','none','EdgeColor','interp')
% colormap cool
% colorbar
% caxis([1 sum(jkvt(trI).vUseFrameIdx)])
% print(fullfile(filePath,'Figure',sprintf('hTrjDstBaseXyz_tr#%d',trI)),'-dpdf','-painters','-bestfit');


% % plot/print a trial's X, Y, Z velocity
% xt = (1:size(hTrjVelBSub{7},2))*4;
% hold on; plot(xt,hTrjVelBSub{7}(1,:)); plot(xt,hTrjVelBSub{7}(2,:)); plot(xt,hTrjVelBSub{7}(3,:)); hold off
% xlabel('Time(ms)'); ylabel('Velocity(mm/s)')
% print(fullfile(filePath,'Figure',sprintf('xyzVel_tr#%d',7)),'-dpdf','-painters','-bestfit');
% 
% % plot/print a trial's position/velocity relative to the global median
% figure; plot(xt,hTrjFDstBase{7});
% hold on;
% for ii = 1:length(reachStart)
%     plot(xt(reachStart(ii)), hTrjFDstBase{7}(reachStart(ii)),'or')
% end
% hold off;
% xlabel('Time(ms)'); ylabel('Position relative to baseline (mm)')
% print(fullfile(filePath,'Figure',sprintf('handPosRelToBaseline_tr#%d',7)),'-dpdf','-painters','-bestfit');
% 
% figure; plot(xt,hTrjFVelBase{7});
% hold on;
% for ii = 1:length(reachStart)
%     plot(xt(reachStart(ii)), hTrjFVelBase{7}(reachStart(ii)),'or')
% end
% hold off;
% xlabel('Time(ms)'); ylabel('Velocity (mm/s)')
% print(fullfile(filePath,'Figure',sprintf('handVelRelToBaseline_tr#%d',7)),'-dpdf','-painters','-bestfit');