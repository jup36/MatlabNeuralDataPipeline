
%filePath = 'S:\Junchol_Data\JS2p0\WR40_081419';
filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081419';

%% load 'trj3d.mat'
cd(filePath)
trj3d = dir('**/*trj3d.mat');
T=load(fullfile(trj3d.folder,trj3d.name),'trj3d');
clearvars trj3d
T = T.('trj3d');

hTrjF = {T(:).allPartsMedSgFronXYZ}; % all parts median and Savitzky-Golay filtered hand trajectory front cam
hTrjS = {T(:).allPartsMedSgSideXYZ}; % all parts median and Savitzky-Golay filtered hand trajectory side cam
jsBot = {T(:).jsBs}; 
jsTop = {T(:).jsTs}; 

%% collect jsTime1k_Kinematics data
jsKinFile = dir(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat'));

if length(jsKinFile)==1
    jkv = load(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat'),'jsTime1k_KV'); % load jsTime1k_KV
else
    disp('Point to "jsTime1k_Kinematics_VideoFiles.mat" with the variable "jsTime1k_KV"!')
    [jsKinFileSelect,jsKinPathSelect] = uigetfile(filePath);
    jkv = load(fullfile(jsKinPathSelect,jsKinFileSelect),'jsTime1k_KV'); % load jsTime1k_KV
end

jkv = jkv.('jsTime1k_KV');
clearvars jsTime1k_KV

%% get the 3d-position/velocity
gmhTrj = median(cell2mat(hTrjS),2); % global median

% get trajectories relative to the global median first, the preparatory period median can be defined after isolating reach starts
medDst = @(x) sqrt(sum((x-gmhTrj).^2)); % func to get the point-by-point Dstance from the global median position

for t = 1:length(hTrjS)
    if ~isempty(hTrjS{t})
        jkv(t).hTrjDstMed = splitapply(medDst, hTrjS{t}, 1:size(hTrjS{t},2)); % point-by-point position relative to the global median
        jkv(t).hTrjVelMed = diff([jkv(t).hTrjDstMed(:,1) jkv(t).hTrjDstMed],1,2)./(4/1000);
        % find the baseline points of low velocity and short-distance from the global median
        [~,sortIhTrjDstMed]=sort(jkv(t).hTrjDstMed); % sort by distance from the global median
        nearMedLogic = zeros(1,length(jkv(t).hTrjDstMed));
        nearMedLogic(sortIhTrjDstMed(1:min(250, length(sortIhTrjDstMed))))=1; % taking points not too far from the global median excludes low-velocity points such as time spent with the limb on the JS without moving          
        lowVelLogic = abs(jkv(t).hTrjVelMed)<5; % low-velocity points 
        
        if sum(nearMedLogic & lowVelLogic)>10
            basePosM = median(hTrjS{t}(:,nearMedLogic & lowVelLogic), 2); % low velocity (<5 mm/s) points
        else
            basePosM = gmhTrj; % just take the global median points
        end
        baseDst = @(x) sqrt(sum((x-basePosM).^2)); % func to get the point-by-point Dstance from the global median position
        jkv(t).hTrjDstBase = splitapply(baseDst, hTrjS{t}, 1:size(hTrjS{t},2)); % point-by-point position relative to trial-by-trial baseline
        jkv(t).hTrjVelBase = diff([jkv(t).hTrjDstBase(:,1) jkv(t).hTrjDstBase],1,2)./(4/1000); % point-by-point velocity
        jkv(t).hTrjDstBaseXyz = hTrjS{t}-basePosM; % normalize by subtracting the median of the baseline position X, Y, Z
        jkv(t).hTrjVelBaseXyz = diff([jkv(t).hTrjDstBaseXyz(:,1) jkv(t).hTrjDstBaseXyz],1,2);
        %elseif isempty(hTrjS{t})
        %    jkv(t).hTrjDstMed = []; jkv(t).hTrjVelMed = []; jkv(t).hTrjDstBase = []; jkv(t).hTrjVelBase = []; jkv(t).hTrjDstBaseXyz = []; jkv(t).hTrjVelBaseXyz = [];
    end
end
clearvars t

%% get the velocity and position cut
vel = [jkv(:).hTrjVelBase]; % baseline normalized 1-d velocity across all time points
sortVel=sort(vel); % sort velocity data
velCut = sortVel(round(length(sortVel)*.98)); % set the 98% cutoff for the velocity data (velocity cutoff is required, since reaches are expected to occur at high velocity)

pos = [jkv(:).hTrjDstBase]; % 1-d position relative to the trial-by-trial baseline
sortPos = sort(pos); % sort position data
sortPos(sortPos<median(sortPos))=median(sortPos); % replace less-than median position values with the median position value
posCut = 10; % set the Distance cut off as 10 mm from the baseline sortPos(round(length(sortPos)*.90));

for t = 1:size(jkv,2)
    [jkv(t).hTrjRstart, jkv(t).hTrjRstop] = getReachTimesJs2p0(jkv(t).hTrjDstBase, jkv(t).hTrjVelBase, posCut, velCut);
    % posTrace= jkv(t).hTrjDstBase; velTrace= jkv(t).hTrjVelBase;
end
clearvars t

%% get the joystick trajectory


%% examine trajectories by plotting
%rewardTrs = find([jsTime1k_KV(:).rewarded]==1);
close all;
rwdTrI = 2;
x = T(rewardTrs(rwdTrI)).allPartsMedSgSideXYZ(1,:); %T(rewardTrs(rwdTrI)).vUseFrameIdx);
y = T(rewardTrs(rwdTrI)).allPartsMedSgSideXYZ(2,:); %T(rewardTrs(rwdTrI)).vUseFrameIdx);
z = T(rewardTrs(rwdTrI)).allPartsMedSgSideXYZ(3,:); %T(rewardTrs(rwdTrI)).vUseFrameIdx);
c = 1:sum(T(rewardTrs(rwdTrI)).vUseFrameIdx); % generate a colormap;

figure;
plot3(gmhTrj(1), gmhTrj(2), abs(gmhTrj(3)), 'ro')
patch([x nan],[y nan],[z nan],[c nan],'FaceColor','none','EdgeColor','interp')
%patch([x nan],[z nan],-[y nan],[c nan],'FaceColor','none','EdgeColor','interp')
xlabel('X(mm)')
ylabel('Y(mm)')
zlabel('Z(mm)')

%patch([x nan],[z nan],abs([y nan]),[c nan],'FaceColor','none','EdgeColor','interp')
colormap parula
colorbar
caxis([400 1000])

% plot/print a trial's X, Y, Z velocity
xt = (1:size(hTrjVelBSub{7},2))*4;
hold on; plot(xt,hTrjVelBSub{7}(1,:)); plot(xt,hTrjVelBSub{7}(2,:)); plot(xt,hTrjVelBSub{7}(3,:)); hold off
xlabel('Time(ms)'); ylabel('Velocity(mm/s)')
print(fullfile(filePath,'Figure',sprintf('xyzVel_tr#%d',7)),'-dpdf','-painters','-bestfit');

% plot/print a trial's position/velocity relative to the global median
figure; plot(xt,hTrjSDstBase{7});
hold on;
for ii = 1:length(reachStart)
    plot(xt(reachStart(ii)), hTrjSDstBase{7}(reachStart(ii)),'or')
end
hold off;
xlabel('Time(ms)'); ylabel('Position relative to baseline (mm)')
print(fullfile(filePath,'Figure',sprintf('handPosRelToBaseline_tr#%d',7)),'-dpdf','-painters','-bestfit');

figure; plot(xt,hTrjSVelBase{7});
hold on;
for ii = 1:length(reachStart)
    plot(xt(reachStart(ii)), hTrjSVelBase{7}(reachStart(ii)),'or')
end
hold off;
xlabel('Time(ms)'); ylabel('Velocity (mm/s)')
print(fullfile(filePath,'Figure',sprintf('handVelRelToBaseline_tr#%d',7)),'-dpdf','-painters','-bestfit');


% x = T(rewardTrs(rwdTrI)).allPartsMedSgFron(1,:); %T(rewardTrs(rwdTrI)).vUseFrameIdx);
% y = T(rewardTrs(rwdTrI)).allPartsMedSgFron(2,:); %T(rewardTrs(rwdTrI)).vUseFrameIdx);
% z = T(rewardTrs(rwdTrI)).allPartsMedSgFron(3,:); %T(rewardTrs(rwdTrI)).vUseFrameIdx);
% c = 1:sum(T(rewardTrs(rwdTrI)).vUseFrameIdx); % generate a colormap;
%
% figure;
% patch([x nan],[z nan],abs([y nan]),[c nan],'FaceColor','none','EdgeColor','interp')
% colormap parula
% colorbar
% caxis([400 1000])
%
% print(fullfile(filePath,'Figure',sprintf('tr#%d',rewardTrs(rwdTrI))),'-dpdf','-painters','-bestfit')
% jsTime1k_KV(rewardTrs(rwdTrI)).fVideo
% jsTime1k_KV(rewardTrs(rwdTrI)).sVideo