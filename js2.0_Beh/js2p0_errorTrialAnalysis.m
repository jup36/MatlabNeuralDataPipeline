
filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'; 
%saveName = ''; 

%% 1. Load data
%clc; clearvars; close all;
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';
cd(filePath)
% neural and behavioral data
spkDir = dir('binSpkCountSTRCTX*'); 
load(fullfile(spkDir(1).folder, spkDir(1).name),'spkTimesCell','rStartToPull','jkvt')
S=rStartToPull; clearvars rStartToPull
% behavioral data
%behDir = dir('jsTime1k_KinematicsTrajectories*');
%load(fullfile(behDir(1).folder,fullfile(behDir(1).name)),'jkvt');

%% align hand trajectories to neural data
ts = S.currEvt{1}(:,1); % time stamps (e.g. rStartToPull)
% locate trials map each event (ts) to trials in jkvt
tsMapJkvt1 = find(repmat(ts(1),[size(jkvt,2),1]) <= [jkvt(:).trEnd]',1,'first'); 
tsMapJkvt = arrayfun(@(a) find([jkvt(1:end-1).trEnd]'<= repmat(a,[size(jkvt,2)-1,1]) & repmat(a,[size(jkvt,2)-1,1]) <= [jkvt(2:end).trEnd]')+1, ts(2:end), 'un', 0 ); 
tsJkvtTrs = [{tsMapJkvt1}; tsMapJkvt]; 
tsJkvtTrs = cell2mat(tsJkvtTrs); 

vfT = {jkvt(:).vFrameTime}'; % video frame time
hTrj = {jkvt(:).hTrjF}'; % hand trajectory
%hTrj = {jkvt(:).hTrjDstBaseXyz}'; % baseline-subtracted hand trajectory

spikeB = S.params.binEdges(2:end);
unitTimeTrial = S.unitTimeTrial;

if ~(length(ts)==size(unitTimeTrial,3))
    error('trial numbers do not match!')
end
spk_kernel = TNC_CreateGaussian(250,15,500,1); % mu, sigma, time, dT

binSize = 20; % 20 ms

% cortex/striatum index
depth = cell2mat(cellfun(@(a) a(2), spkTimesCell(4,:),'un',0))'; % depth from pial surface
ctxI = depth<1900; % cortex index
strI = depth>2100; % striatum index

%% parse trial/block information 
% pull torque & reach position combinations
pTqs = round(unique([jkvt(:).pull_torque])./10); % pull torque list (divided by 10 for ordering)
rPos = round(unique([jkvt(:).reachP1])*10); % reach position list (multiplied by 10 for ordering)
posTq = unique(repmat(pTqs,[length(rPos) 1])+repmat(rPos',[1 length(pTqs)])); % position-torque combinations
posTqCnt = zeros(length(posTq),1); % count occurrence of each combination

%% align behavioral and neural trajectories
% determine where to align based on trial types
for t = 1:size(jkvt,2)
    ss(t).trialType = jkvt(t).trialType;
    switch jkvt(t).trialType
        case 'sp' % successful pull 
            % spI = cell2mat(cellfun(@(a) strcmpi(a,'sp'),{jkvt(:).trialType},'un',0)); 
            ss(t).evtAlign = jkvt(t).rStartToPull; 
            ss(t).evtAlignInfo = 'rStartToPull';   
            
            ss(t).unitTimeB = ; % spike structure
        case 'ps' % push error
            if ~isempty(jkvt(t).hTrjRstart)
                ss(t).evtAlign = jkvt(t).hTrjRstart(end);                
            end
            
            
        case 'pmpp' % premature pull and push
            
        case 'pm' % premature pull
            
        case 'to' % time out
    
    end
    
    
end


%% parse trial/block information
function jkvt = jkvtBlockParse(jkvt)
tqd = diff([jkvt(1).pull_torque, jkvt(:).pull_torque]'); % torque change
p1d = diff([jkvt(1).reachP1, jkvt(:).reachP1]'); % position 1 change

blNumb = 1;  % block number
blType = []; % block type
% detect and parse block shifts across trials
for t = 1:size(jkvt,2)
    % first trial
    if t==1
        jkvt(t).blNumber = blNumb;
        jkvt(t).blType = 'first';
        jkvt(t).blShiftLogic = true;
    end
    % second trial and onward 
    if t>=2
        if tqd(t)==0 && p1d(t)==0 % if there's no change
            jkvt(t).blNumber = blNumb; % remains the same
            jkvt(t).blShiftLogic = false; % not shifted
            jkvt(t).blType = []; % block shift type
        else % if there's any change, classify the block shift type
            blNumb = blNumb+1;
            jkvt(t).blNumber = blNumb; % remains the same
            jkvt(t).blShiftLogic = true; % shifted
            % parse torque shift
            if tqd(t)>0
                tqi = 'tqUp';
            elseif tqd(t)<0
                tqi = 'tqDown';
            elseif tqd(t)==0
                tqi = 'tqSame';
            end
            % parse position shift
            if p1d(t)>0
                p1i = 'rightward';
            elseif p1d(t)<0
                p1i = 'leftward';
            elseif p1d(t)==0
                if jkvt(t).reachP1==min([jkvt(:).reachP1])
                    p1i = 'sameleft';
                elseif jkvt(t).reachP1==max([jkvt(:).reachP1])
                    p1i = 'sameright';
                end
            end
            jkvt(t).blType = strcat(tqi,p1i); % block shift type
        end
    end
end
end


