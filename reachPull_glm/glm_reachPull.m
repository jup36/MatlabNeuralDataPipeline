%% 1. Load data
clc; clearvars; close all;

% task parameter
DT = 1; %0.001; % milliseconds

N_pSTART = 10; % reach start, n bumps
RANGE_pSTART = 2000; % ms
BIAS_pSTART = 8;

N_REWARD = 10; % reward
RANGE_REWARD = 3000; % ms
BIAS_REWARD = 8;

N_SPEED = 10; % n bumps
RANGE_SPEED = [0, 40]; % cm/s
%BIAS_SPEED = 10; 

N_H = 10; % n bumps for spike history
RANGE_H = 200; % ms
BIAS_H = 1;

% code parameter
PLOT = true;
CALC_PRM = false;

%% 1. Load data
cd('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles')
% neural data
load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles','binSpkCountCTXWR40_081919.mat'), 'spkTimesCell')

% behavioral data
load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles','jsTime1k_KinematicsTrajectories.mat'),'jkvt','meta')

%% 2. Behavioral data 
% detect events and continuous joystick pull trajectory  
sessionDur = round(str2double(meta.fileTimeSecs)*1000); 
time_bin = (0:DT:sessionDur)';
n_bin = length(time_bin) - 1;

[jkvt,k] = func.evtKinematics( jkvt, sessionDur ); % pullStarts, pullStops, rStartToPull, rStopToPull
k.jsVel = abs(k.jsVel)./10; % convert to positive-valued pull velocities in cm/S

pullStartI = cellfun(@(c) ~isempty(c), {jkvt(:).pullStarts}); 
pullStopI = cellfun(@(c) ~isempty(c), {jkvt(:).pullStops}); 
rewardTimeI = cellfun(@(c) ~isnan(c), {jkvt(:).rewardT}); 

tqList = unique([jkvt(:).pull_torque]);
posList = unique([jkvt(:).reachP1]); 

lowTqI = cellfun(@(c) c==tqList(1), {jkvt(:).pull_torque}); 
highTqI = cellfun(@(c) c==tqList(2), {jkvt(:).pull_torque}); 

leftI = cellfun(@(c) c==posList(1), {jkvt(:).reachP1}); 
rightI = cellfun(@(c) c==posList(2), {jkvt(:).reachP1}); 

leLPull = [jkvt(pullStartI & lowTqI & leftI).pullStarts]; 
leHPull = [jkvt(pullStartI & highTqI & leftI).pullStarts];  
riLPull = [jkvt(pullStartI & lowTqI & rightI).pullStarts]; 
riHPull = [jkvt(pullStartI & highTqI & rightI).pullStarts]; 
reward = [jkvt(rewardTimeI).rewardT]; 

time_leLPull = histcounts(leLPull, time_bin)';
time_leHPull = histcounts(leHPull, time_bin)';
time_riLPull = histcounts(riLPull, time_bin)';
time_riHPull = histcounts(riHPull, time_bin)';
time_reward = histcounts(reward, time_bin)'; 

% bump
[start_base, start_time, start_func] = basis.log_cos(N_pSTART, RANGE_pSTART, DT, BIAS_pSTART, false);
[reward_base, reward_time, reward_func] = basis.log_cos(N_REWARD, RANGE_REWARD, DT, BIAS_REWARD, false);
[speed_basis, speed_time, speed_func] = basis.linear_cos(N_SPEED, RANGE_SPEED, 1, false);

% convolution
X_leLPull = basis.conv(time_leLPull, start_base); % ensure to input histcounted time (not discrete timestamps)
X_leHPull = basis.conv(time_leHPull, start_base);
X_riLPull = basis.conv(time_riLPull, start_base);
X_riHPull = basis.conv(time_riHPull, start_base);
X_reward = basis.conv(time_reward, reward_base); 
X_speed = speed_func(k.jsVel);

%% 3. Spike
% spike
%Direction-tuning: 79, 80, 35, 77, 85
i_cell = 79; % for loop

spike_time = spkTimesCell{1,i_cell};
n_spike = length(spike_time);
spike_bin = histcounts(spike_time, time_bin)';
spike_rate = n_spike / session_duration;

% spike bump
[h_base, h_time] = basis.log_cos(N_H, [DT, RANGE_H], DT, BIAS_H, false);

% convolution
X_h = basis.conv(spike_bin, h_base, h_time(1)/DT);




















