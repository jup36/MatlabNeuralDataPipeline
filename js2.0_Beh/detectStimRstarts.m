
filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'; 

if exist(fullfile(filePath,'jsTime1k_KinematicsTrajectories.mat'),'file')
    S = load(fullfile(filePath,'jsTime1k_KinematicsTrajectories.mat'));  % load jsTime1k_KV.mat
    S = S.('jkvt');
else
    error('No behavioral data found!!!')
end
clearvars 'jkvt'

%% find trials on which laser was turned on during reaching (based on hand trajectories)
sc = 0; 
for t = 1:size(S,2)
    if ~isnan(S(t).stimLaserOn) % stim trial
        
        sc = sc+1;    
        s(sc).tr = t; % stim rStart trial id 
        s(sc).stimLaserOnT = S(t).stimLaserOn; % stim laser on
        s(sc).rewarded = S(t).rewarded; % reward 
        
        s(sc).hTrjDstBase = S(t).hTrjDstBase; % distance relative to the base position 
        s(sc).hTrjDstBaseTrelStim = S(t).vFrameTime-s(sc).stimLaserOnT;  % vTime relative to stim
        
        tmpsI = S(t).vFrameTime(S(t).hTrjRstart) > S(t).stimLaserOn & S(t).vFrameTime(S(t).hTrjRstart) < S(t).stimLaserOff; % check if reach started during laser 
        if sum(tmpsI)>0
           % detect s
           s(sc).sT = S(t).vFrameTime(S(t).hTrjRstart(tmpsI))'; % stim rStart time
           s(sc).sFrameI = S(t).hTrjRstart(tmpsI)'; % stim rStart frame indices
           % detect stimRstop
           s(sc).stimRstopT = S(t).vFrameTime(S(t).hTrjRstop(tmpsI))'; % stim rStop
           s(sc).stimRstopFrameI = S(t).hTrjRstop(tmpsI)'; % stim rStop           
        end
    end
end
clearvars t

figure; hold on; 
tr = s(5).tr; 
%T = S(tr).hTrjDstBaseXyz; % baseline-subtracted hand trajectory
T = S(tr).hTrjF; % baseline-subtracted hand trajectory
plot3( T(1,:), T(2,:), T(3,:) )

tr = s(7).tr; 
%T = S(tr).hTrjDstBaseXyz; % baseline-subtracted hand trajectory
T = S(tr).hTrjF; % baseline-subtracted hand trajectory
plot3( T(1,:), T(2,:), T(3,:) )

tr = s(8).tr; 
%T = S(tr).hTrjDstBaseXyz; % baseline-subtracted hand trajectory
T = S(tr).hTrjF; % baseline-subtracted hand trajectory
plot3( T(1,:), T(2,:), T(3,:) )




