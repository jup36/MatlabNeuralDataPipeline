function runJoystickCSV( filePath )
%runJoystickCSV reads the trial-by-trial summary and the 1-d joystick
%   trajectory (sampled at 100Hz, 10ms resolution) csv files.  
%   Classify trials: 
%       1:  Successful pulls (rewarded)
%       -1: Erroneous pushes 
%       0:  Immature pulls   (before 'go')
%       10: Timeout (20-s timeout)

cd(filePath)

trialSum = dir(fullfile(filePath,'trials.csv')); % get the trial-by-trial summary data 
trialSumTab = readtable(fullfile(filePath,trialSum.name));

trialList = dir(fullfile(filePath,'trial_*')); % get the trial-by-trial data
[~,trialIdx] = sort([trialList.datenum]);      % trial index sorted by time (ascending)

if size(trialSumTab,1)==size(trialList,1)
else
   error('The number of trial-by-trial files do not match!') 
end

% use interpolation to estimate the joystick position at higher temporal resolution
refTime = posixtime(datetime(trialSumTab.trial_start{1},'InputFormat','yyyy-MM-dd-HH-mm-SS')); % use the 1st trialStart time as the reference time of the session
for t = 1:length(trialList) % increment trials 
    currCSV = readtable(fullfile(filePath,trialList(trialIdx(t)).name)); % the current trial's CSV
    bDat.jTraj{t} = currCSV.joystick_position; % joystick trajectory
    bDat.trialType{t} = 
    bDat.
    bDat.
end



end

