function behaviorTimestampsReachDecode(filePath)
%addpath(genpath('/Volumes/RAID2/parkj/MATLAB'))
% get reach properties
load(fullfile(filePath,'BehVariables.mat'),'ts','positionData')
[ reachStartMore, noReachTime, reach0 ] = getReachTimesDecode( positionData );     % all reach traces, aligned to start (pos1), to stop (pos2)

ts.reachStartMore  = reachStartMore;  % reachStart
ts.noReachTime = noReachTime; % noReach time

% Save relevant BehVariables
cd(filePath)
save('BehVariablesReachDecode','positionData','reach0','ts' ) % append the position/velocity data variables

end

