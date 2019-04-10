function [behaviorFile] = TNC_BatchProcessMoverCSVdata(fileDirectory)
% General Information
%
%   Load processing .csv files from the MitoMover Stim task into a large 
%   group variable in the matlab workspace
%       The variable produced will have all the information collected in
%       chronological order- the first file being the oldest file
%       Won't work with earlier MitoMover versions because they don't have
%       stimulation and behavior state data in the files- use 
%
%   Provides some basic information about the reach characteristics. If
%   looking for more task specific analysis, use one of csvLoad files that
%   are task specific 
%
%   Output should be compatible with several other functions
%   
% Inputs
%
%   fileDirectory: path to the files to be run, should be in '' and end 
%                  with a / 
%       Note: file names (animal numbers to be created) should not contain 
%       dashes or underscores 
%
% Outputs 
%
%   behaviorFile: all files in one large workspace variables with the name
%                 specified in the behaviorFile variable
 %%
 % Display what this program will do and the current time

disp('________________Running csvLoad Script_________________________')
disp(' ');
disp(['Timestamp: ' datestr(now)]);
disp(' ');

%allFiles = dir(sprintf('%s*_p.csv',fileDirectory));
allFiles = dir('*_p.csv'); 

% Mouse type must be listed directly after 'Behavior' in the directory

behaviorIndex = strfind(fileDirectory, 'Head Fixed'); % non-applicable to me for now
slashesIndex  = strfind(fileDirectory, '/'); % non-applicable to me for now
slashPostBehav = find(slashesIndex > behaviorIndex, 1, 'first'); % non-applicable to me for now
mousetype  = fileDirectory(slashesIndex(slashPostBehav)+1: slashesIndex(slashPostBehav+1)-1); % non-applicable to me for now    

for i = 1:size(allFiles,1)
    %% File Name and Information
    
    disp(['Directory: ' fileDirectory]);
    fileBase       = allFiles(i).name(1:max(strfind(allFiles(i).name,'_')-1)); % basic part of the filename
    year           = allFiles(i).name(strfind(allFiles(i).name, '201'):strfind(allFiles(i).name, '201')+3); % read out year info from the filename
    
    if strcmp(year, '2012') == 1 || strcmp(year,'2011') == 1
        
        month          = allFiles(i).name(strfind(allFiles(i).name, '201')+4:strfind(allFiles(i).name, '201')+5);
        day            = allFiles(i).name(strfind(allFiles(i).name, '201')+6:strfind(allFiles(i).name, '201')+7);
        session        = allFiles(i).name(4:strfind(allFiles(i).name, '201')-1);
        fileNameBase   = allFiles(i).name(1: 3);

        
    else
        
        month          = allFiles(i).name(strfind(allFiles(i).name, '201')+5:strfind(allFiles(i).name, '201')+6);
        day            = allFiles(i).name(strfind(allFiles(i).name, '201')+8:strfind(allFiles(i).name, '201')+9);
        session        = allFiles(i).name(strfind(allFiles(i).name, '201')+11:max(strfind(allFiles(i).name, '_'))-1);
        fileNameBase   = allFiles(i).name(1: min(strfind(allFiles(i).name, '_'))-1);
    
    end

    mousenumber = fileNameBase(1: nnz(fileNameBase));     

    data.name.year    = year;
    data.name.month   = month;
    data.name.day     = day;
    data.name.session = session; 
    data.name.number  = i;

    data.name.mouseNumber = mousenumber;
    data.name.mouseType   = mousetype; % NA for me now
      
    disp(' ');
    disp(['File ' num2str(i) ' base name: ' fileNameBase]); 
    disp(['Date created: ' year '/' month '/' day]);
    disp(' ');
    
%% Script to load csv behavior data
    clear disp* trial* 

    skipInitTrials = 0;

    fileNameStr = fullfile(fileDirectory,strcat(fileBase,'_p.csv'));
    data.trialParams = dlmread(fileNameStr,',',1+skipInitTrials,0); % 11th column is the column containing information about which trials were stimulated

    fileNameStr = fullfile(fileDirectory,strcat(fileBase, '.csv'));
    data.trialTimes = dlmread(fileNameStr,',',2+skipInitTrials,0);  % time information 
        
    fileNameStr = fullfile(fileDirectory,strcat(fileBase, 'pXY.csv')); 
   
    if exist(fileNameStr,'file')== 1 || exist(fileNameStr,'file')== 2
        
        data.trajectories.contXY = dlmread(fileNameStr,',',2,0);    % extract semicontinuous data into data.trajectories.contXY file
        
    end
    
    % data.trajectories.contXY format
        % data.trajectories.contXY (:,1) = time
        % data.trajectories.contXY (:,2) = x
        % data.trajectories.contXY (:,3) = y
        % data.trajectories.contXY (:,4) = licking data
        % data.trajectories.contXY (:,5) = this must be the behavioral state (paradigmmode)
        % data.trajectories.contXY (:,6) = this must be solenoid (see that it corresponds to the paradigmmode4 in col5)

%% General information about the file

    indArray = 2:size(data.trialTimes,1);

    data.startToEvent  = data.trialTimes(indArray,3) - data.trialTimes(indArray,2);     % time from trial start to threshold cross
    data.endToStart    = data.trialTimes(indArray,2) - data.trialTimes(indArray-1,4);   % inter-trial interval ()
    data.eventToEnd    = data.trialTimes(indArray,4) - data.trialTimes(indArray,3);     % reward delay
    data.threshold.max = max(data.trialParams (:,6));                                   % max threshold
    data.threshold.min = min(data.trialParams(:,6));                                    % min threshold
    numTrials          = size(data.trialTimes,1);                                       % number of total trials
    
  % program starts with ITI- need to remove that with a shift in order to compare continuous data and trial times file
    behavStateDiff      = diff(data.trajectories.contXY(:,5)); % this must be column5 where the behavioral state data are saved
    
    reachTriStartIndex2 = find(behavStateDiff == 2); % this should be able to detect the shift of paradigmmode from 0 to 1 which indicates start of the trial 
    reachTriStartIndex3 = find(behavStateDiff == 3);
    reachTriStartIndex  = [reachTriStartIndex2; reachTriStartIndex3];
    reachTriStartIndex  = sort(reachTriStartIndex);
    contStartTimes      = data.trajectories.contXY(reachTriStartIndex+1,1);  % find the start of trials in the continuous file
    
    % this part simply matches the # of trial and reach
    if nnz(contStartTimes) ~= nnz(data.trialTimes(:,2)) % nnz counts the number of nonzero elements in the matrix 
        
       contLarge = nnz(contStartTimes) > nnz(data.trialTimes(:,2));
       contSmall = nnz(contStartTimes) < nnz(data.trialTimes(:,2));
       if contLarge == 1
           contStartTimes(nnz(data.trialTimes(:,2))+1:end) = []; % curtail the exceeding trial 
       end
       
    end
    
    % why is there shift?
    shift       = contStartTimes - data.trialTimes(:,2);            % data.trialTimes(:,2) contains the trial start time triggered by the controller
    startTimes  = data.trialTimes(:,2) + shift;                     % trial start times
    eventTimes  = data.trialTimes(:,3) + shift;                     % threshold cross times
    rewardTimes = data.trialTimes(:,4) + shift;                     % reward times

       
%     end
%%       
% Do the following for semicontinuous (XY) data files
    
    if exist(fileNameStr,'file')== 1 || exist(fileNameStr,'file')== 2
        
        % Smooth continuous position data
        
        for j=2:3 % two columns containing lever position information
            
            data.trajectories.contXYsmooth(:,j) = sgolayfilt(data.trajectories.contXY(:,j),3,11); 
        
        end
        
        %Convert au's to mm 
        
        beginTrial = find(data.trajectories.contXY(:,6) == 2, 1, 'first');
        
%         data.trajectories.contXYsmooth(beginTrial:end,2) = tan(data.trajectories.contXYsmooth(beginTrial:end,2).*(.1*(pi./180)))*105;
%         data.trajectories.contXYsmooth(beginTrial:end,3) = tan(data.trajectories.contXYsmooth(beginTrial:end,3).*(.1*(pi./180)))*105;
        
        data.trajectories.contXYsmooth(beginTrial:end,2) = tan(data.trajectories.contXYsmooth(beginTrial:end,2).*(.1*(pi./180)))*90; % my joystick is 90 mm long
        data.trajectories.contXYsmooth(beginTrial:end,3) = tan(data.trajectories.contXYsmooth(beginTrial:end,3).*(.1*(pi./180)))*90; % my joystick is 90 mm long
        
        data.trajectories.contXYsmooth(1:beginTrial-1,2) = 0;
        data.trajectories.contXYsmooth(1:beginTrial-1,3) = 0;
        
        % Keep some data the same in columns and don't smooth
        
        data.trajectories.contXYsmooth(:,1)  = data.trajectories.contXY(:,1); % time
        data.trajectories.contXYsmooth(:,7)  = data.trajectories.contXY(:,5); % stimulation time
        data.trajectories.contXYsmooth(:,8)  = data.trajectories.contXY(:,6); % behavioral state
        data.trajectories.contXYsmooth(:,9)  = data.trajectories.contXY(:,4); % licking
        
        
 % Fill in Associated Kinematics
        instTime = diff(data.trajectories.contXYsmooth(:,1));
        
        % find distance of each timepoint
        diffXsq = (diff(data.trajectories.contXYsmooth(:,2))).^2;          % distance in x traveled between time points
        diffYsq = (diff(data.trajectories.contXYsmooth(:,3))).^2;          % distance in y traveled between time points
        
        instDist = sqrt(diffXsq + diffYsq);                                % distance traveled in between time points
        instVel = instDist./instTime;                                      % distance/time between time points
        instAcc = diff(instVel)./instTime(2:numel(instTime));              % velocity/time between time points
        
        data.trajectories.contXYsmooth(:,4) = [0; instDist];               % distance at time points
        data.trajectories.contXYsmooth(:,5) = [0; instVel];                % velocity at time points
        data.trajectories.contXYsmooth(:,6) = [0; 0; instAcc];             % acceleration at time points
        
        data.trajectories.sLeverVm = data.trajectories.contXYsmooth(:,5);  % medfilt1(data.trajectories.contXYsmooth(:,5),5); % smooth velocity
        
%% Isolate reaches from smoothed velocity function
        
        % Threshold velocities; must be >.005mm/ms
%         ABOVE    = find(data.trajectories.sLeverVm > 0.008);
        
        ABOVE    = find(data.trajectories.sLeverVm > 0.005);
        tmp.star = find([10; diff(ABOVE)]>1);
        tmp.enn  = [tmp.star(2:end)-1; size(ABOVE,1)];
        
        tmp.reachIND = [tmp.star tmp.enn];
        
% Delete reaches with approx <100 ms (4 samples = 132ms) above threshold long
remove = [];

for f = 1: size(tmp.reachIND,1)
    
    if tmp.reachIND(f,2) - tmp.reachIND(f,1) + 1 <= 4
        remove = [remove; f];
    end
    
end

tmp.reachIND(remove,:) = []; 

% Merge reaches in the case that next reach start < previous reach stop and
% % Merge reaches where there is approx. <264ms between reaches, given 30Hz
% rate this is about 8 (264ms) samples

remove = [];

for f = 2: size(tmp.reachIND,1)
    
    if ABOVE(tmp.reachIND(f,1)) - ABOVE(tmp.reachIND(f-1,2)) <= 8
        tmp.reachIND(f-1,2) = tmp.reachIND(f,2);
        remove = [remove; f];
    end
    
end

tmp.reachIND(remove,:) = []; 
tmp.star = tmp.reachIND(:,1); 
tmp.enn  = tmp.reachIND(:,2);

% More precisely find start and stops of reaches

        star = [];
        enn = [];
%         for f = 2:size(tmp.star)
%             A = find(data.trajectories.sLeverVm(1: ABOVE(tmp.star(f)),:) <= .0010, 1, 'last');
%             B = find(data.trajectories.sLeverVm(ABOVE(tmp.enn(f)):end,:) <= .0010, 1, 'first') + ABOVE(tmp.enn(f)) - 1;
%             star = [star; A];
%             enn = [enn; B];
%         end
        for f = 1:size(tmp.star)
            
            if f == 1
                clear A B
                A = ABOVE(tmp.star(f));
            elseif f >= 2
                clear startInds stopInds A B
                startInds = beginTrial: ABOVE(tmp.star(f));          
                A = find(data.trajectories.sLeverVm(startInds,:) <= 0.0010, 1, 'last');
                
                if isempty(A) == 1
                    A = beginTrial;
                end
                
                A = startInds(A);
            end
            
            stopInds = ABOVE(tmp.enn(f));
            B = find(data.trajectories.sLeverVm(stopInds:end,:) <= .0010, 1, 'first') + ABOVE(tmp.enn(f)) - 1;
            
            star = [star; A];
            enn = [enn; B];
        end

        if nnz(star) ~= nnz(enn)
            
            if nnz(star) > nnz(enn)
                numDiff = nnz(star)-nnz(enn);
                star((nnz(star)-numDiff+1):end) = [];
            elseif nnz(enn) > nnz(star)
                numDiff = nnz(enn)-nnz(star);
                enn((nnz(enn)-numDiff+1):end) = [];
            end
            
        end
           
reach.timesIND = [star enn];

% Remove duplicates of merged reaches

tmp = diff(reach.timesIND(:,1));
rem = find(tmp == 0);
reach.timesIND(rem,:) = [];
numreach = size(reach.timesIND,1);

% delete reaches that are after the final reward time
reachStarts = data.trajectories.contXYsmooth(reach.timesIND(:,1),1);
delete = find((rewardTimes(end) > reachStarts) == 0);                   
if isempty(delete) == 0
    reach.timesIND(delete, :) = [];
end

%% Information about the extracted reaches 

 for k = 1:size(reach.timesIND,1)
                
    progStart = reach.timesIND(k,1);
    progStop = reach.timesIND(k,2);

        trajAngle   = atan2(data.trajectories.contXYsmooth(progStop,3)-data.trajectories.contXYsmooth(progStart,3),data.trajectories.contXYsmooth(progStop,2)-data.trajectories.contXYsmooth(progStart,2));

        velTraj    = data.trajectories.sLeverVm(progStart : progStop);
        tstmps     = data.trajectories.contXYsmooth(progStart: progStop, 1);
        trajPoints = data.trajectories.contXYsmooth(progStart: progStop, 2:3);

        reach.start(k)   = data.trajectories.contXYsmooth(progStart,1);
        % delete reaches that occur after the last reward is delivered
        
        reach.stop(k)    = data.trajectories.contXYsmooth(progStop,1);
        reach.angle(k,1) = trajAngle;
        
        reach.dist(k,1)  = trapz(velTraj);                                 % average distance
        reach.dist(k,2)  = sum(data.trajectories.contXYsmooth(progStart:progStop,5));  % total distance

        tmp = findpeaks(velTraj);
        reach.numpks(k,1) = numel(tmp);
        reach.dur(k,1)    = reach.stop(k) - reach.start(k);                % reach duration
        reach.vel(k,1)    = max(velTraj);                                  % maximum velocity
        reach.vel(k,2)    = trapz(velTraj) ./ reach.dur(k,1);              % average velocity over reach 
        reach.vel(k,3)    = var(velTraj);                                  % variation in velocity

        reach.acc(k,1) = max(diff(velTraj)./diff(tstmps));                 % maximum acceleration per reach
        reach.acc(k,2) = mean(diff(velTraj)./diff(tstmps));                % mean acceleration
        reach.acc(k,3) = max(diff(velTraj(1:3))./diff(tstmps(1:3)));       % max in first 90 ms of movement

        reach.tort(k,1)  = reach.dist(k,2)./ pdist2([trajPoints(end,1),trajPoints(end,2)],[trajPoints(1,1),trajPoints(1,2)]);

        reach.trial(k,1) = find(rewardTimes >= reach.start(k),1, 'first');
        
 end
 
%  reach.trial = reach.trial + 1;                                            % trial each reach is associated with
 
for k = 1:size(reach.start,2)
    
    trialNum = find(data.trialTimes(:,2) >= reach.start(k) , 1);
    
    if numel(trialNum)==1
        reach.block (k,1) = data.trialParams(trialNum,6);
    else
        reach.block (k,1) = data.trialParams(size(data.trialParams,1),6);
    end
    
end
 
    reach.acc(isinf(reach.acc)) = NaN;
    reach.start = reach.start'; reach.stop = reach.stop'; reach.block = reach.block'; 
    data.reach = reach;
    
    data.reach.times = [reach.start, reach.stop]; 
    
    
%% Generate Reach Analysis

disp(' ');
disp(' ');

maxReach = size(data.reach.vel(:,1),1);

% Per reach analysis
clear forFit;

minThresh = min(data.trialParams(:,6));
maxThresh = max(data.trialParams(:,6));

for reachNum = 1:maxReach;
    
    totalTrials = size(data.trialTimes,1);
    startR  = find(data.trajectories.contXY(:,1) >= data.reach.start(reachNum) , 1);
    stopR   = find(data.trajectories.contXY(:,1) >= data.reach.stop(reachNum) , 1);
    trialNum= find(startTimes >= data.reach.start(reachNum) , 1);
    
        if numel(trialNum)==1
            trialThresh = data.trialParams(trialNum,6);
        else
            trialThresh = data.trialParams(size(data.trialParams,1),6);
        end

        if stopR-startR >= 9
            xVals = sgolayfilt(data.trajectories.contXY(startR:stopR, 2) , 5 , 9);
            yVals = sgolayfilt(data.trajectories.contXY(startR:stopR, 3) , 5 , 9);
        else
            xVals = sgolayfilt(data.trajectories.contXY(startR:stopR, 2) , 5 , 7);
            yVals = sgolayfilt(data.trajectories.contXY(startR:stopR, 3) , 5 , 7);
        end
        
    vVals = data.trajectories.sLeverVm(startR:stopR) .* 100;
    xVals = xVals - xVals(1);
    yVals = yVals - yVals(1);   

%     multiplier = 33;
    multiplier = mean(diff(data.trajectories.contXY(startR:stopR,1)));

    posMult = 100;
    
    try
        k = convhull(xVals,yVals);
    catch ME
        k = 1:numel(xVals);
    end

    hullPntsX = xVals(k);
    hullPntsY = yVals(k);
    hullDists = sqrt( hullPntsX.^2 + hullPntsY.^2 );
    maxDisplace = find(hullDists == max(hullDists),1);                      % index of max displacement
    maxDispTotalTime = find(xVals==hullPntsX(maxDisplace) & yVals==hullPntsY(maxDisplace), 1, 'first');     % index of max displacement in the entire file

    
        changeInX = diff(xVals(1:maxDispTotalTime)).^2;
        changeInY = diff(yVals(1:maxDispTotalTime)).^2;
        changeInPosition = sqrt(changeInX + changeInY);
            
        forFit(reachNum,1)  = sum(changeInPosition).*posMult;                % total displacement            
        forFit(reachNum,2)  = max(hullDists).*posMult;                       % maximum displacement
        forFit(reachNum,3)  = maxDispTotalTime .*multiplier;                 % time from reach start of maximum displacement
        forFit(reachNum,4)  = trialThresh;                                   % trial threshold
        forFit(reachNum,5)  = polyarea( hullPntsX , hullPntsY );             % area around reach trajectory            
        forFit(reachNum,6)  = data.trajectories.contXY(stopR, 1) - data.trajectories.contXY(startR, 1);   % reach duration              
        forFit(reachNum,7)  = sum(changeInPosition)/max(hullDists);          % tortuosity
        forFit(reachNum,8)  = data.trajectories.contXY(startR,1)+ (maxDispTotalTime .* multiplier); % time of maximum displacement
        forFit(reachNum,9)  = data.trajectories.contXY(startR,1);            % time of reach start
        forFit(reachNum,10) = data.trajectories.contXY(stopR,1);             % time of reach stop
        
[rho,h] = corr(forFit(:,4), forFit(:,6));
reachAnalysis.thrDurCorr.rho  = rho; 
reachAnalysis.thrDurCorr.prob = h;

end

reachAnalysis.fitPerReach   = forFit;
reachAnalysis.maxReach      = maxReach;

%% Per reward analysis

clear dispData;

minThresh       = min(data.trialParams(:,6));                              % minimum threshold
maxThresh       = max(data.trialParams(:,6));                              % maximum threshold

for g = 2: nnz(rewardTimes)                                                % use all trials (except 1st)- should remove first 10 trials from other analyses
    
    startR      = find(data.trajectories.contXY(:,1) >= rewardTimes(g-1) , 1);
    stopR       = find(data.trajectories.contXY(:,1) >= rewardTimes(g) , 1);
    trialThresh = data.trialParams(g,6);
    
    yVals       = sgolayfilt(data.trajectories.contXY(startR:stopR, 3) , 5 , 11);
    xVals       = sgolayfilt(data.trajectories.contXY(startR:stopR, 2) , 5 , 11); 
    
    dX = diff(xVals);
    dY = diff(yVals);
    vVals = sqrt( dX.^2 + dY.^2 );
    
    totalDisp                   = trapz(vVals);
    
    dispData.interRewardInt(g-1) = rewardTimes(g) - rewardTimes(g-1);
    dispData.interTrialInt(g-1)  = startTimes(g) - rewardTimes(g-1);
    dispData.disp(g-1)           = totalDisp;
    dispData.sqrDisp(g-1)        = totalDisp^2;
    dispData.vel(g-1)            = max(vVals);   
    dispData.rrate(g-1)          = 1000 ./ ((stopR-startR).*33);
    dispData.thresh(g-1)         = trialThresh;
    
end

[rho,h] = corr(dispData.rrate(10:end)',dispData.disp(10:end)');
reachAnalysis.dispRewCorr.rho  = rho;
reachAnalysis.dispRewCorr.prob = h;

[rho,h] = corr(dispData.rrate(10:end)',dispData.sqrDisp(10:end)');
reachAnalysis.sqrDispRewCorr.rho  = rho;
reachAnalysis.sqrDispRewCorr.prob = h;

reachAnalysis.dispData = dispData;

%% data.timing

% Wait time from reward to next reach
if data.reach.trial(1) == 1     % can't use first reach
    
    firstTriReaches = find (data.reach.trial == 1);
    numFirstReach   = 1;
    prevRewardAll   = data.reach.trial-1;
    previousReward  = prevRewardAll(firstTriReaches(end)+1:end);         % the previous reward
    reachStartTimes = data.reach.start(firstTriReaches(end)+1:end);
    timeFromPrevRew(firstTriReaches,1) = 0;
    
    for w = 1:nnz(previousReward)
        timeFromPrevRew(firstTriReaches(end)+w,1) = reachStartTimes(w) - rewardTimes(previousReward(w));    % of current reach
    end
    
else                            % can use the first reach
   
    numFirstReach   = 0;
    previousReward  = data.reach.trial-1;
    reachStartTimes = data.reach.start;
    timeFromPrevRew = reachStartTimes - rewardTimes(previousReward);
    
end

repTrials = unique(previousReward)+1;                                                                             

for s = 1:nnz(repTrials)
    firstReach (s,2) = find(data.reach.trial == repTrials(s), 1, 'first'); % index of first reach of trial
end

firstReach(:,1) = repTrials-1;                                             % trial number
firstReach(:,3) = timeFromPrevRew(firstReach(:,2));                        % time from previous reward

first10Trials = find(firstReach (:,1) <= 10, 1, 'last');

% first ten trials are not included in the statistics calculations
data.timing.waitITI       = firstReach; 
data.timing.meanWaitITI   = mean(timeFromPrevRew(firstReach(first10Trials+1:end,2)));                      % mean wait to reach from previous reward
data.timing.medianWaitITI = median(timeFromPrevRew(firstReach(first10Trials+1:end,2)));
data.timing.varWaitITI    = var(timeFromPrevRew(firstReach(first10Trials+1:end,2)));
data.timing.rangeWaitITI  = max(timeFromPrevRew(firstReach(first10Trials+1:end,2))) - min(timeFromPrevRew(firstReach(first10Trials+1:end,2)));

 %% Time between unrewarded reaches

 maxReach = nnz(reachAnalysis.fitPerReach(:,1));
 
 for reach = 1:maxReach-1
    timeBtwnReach((reach), 1) = reachAnalysis.fitPerReach(reach+1,9) - reachAnalysis.fitPerReach(reach,8);          %time between max displacement of reach and start of subsequent reach
 end
 
 for trial = 10: nnz(eventTimes)        % only for trials 10 - end
    
     rewardedReaches(trial-9,1) = trial;                                   % trial reach is associated with
     rewardedReaches(trial-9,2) = find(data.reach.start <= eventTimes(trial), 1, 'last');   % reach
         
 end
 
 unrewardedReach = 1:maxReach;
 unrewardedReach(rewardedReaches(2:end,2)) = [];
 unrewardedReach (1:rewardedReaches(1,2)) = [];
 
 if unrewardedReach(end) == maxReach
     unrewardedReach(end) = [];
 end
 
 
 data.timing.waitRD(:, 1) = data.reach.trial(unrewardedReach);             %%%%%%%%% check to make sure accurate
 data.timing.waitRD(:, 2) = unrewardedReach;
 data.timing.waitRD(:, 3) = timeBtwnReach(unrewardedReach);

 data.timing.meanRD = mean(data.timing.waitRD(:,3));
 data.timing.medRD  = median(data.timing.waitRD(:,3));
 data.timing.varRD  = var(data.timing.waitRD(:,3));

%%  Stimulation paradigm and analysis

% stimulation 1: stimulation during ITI
% stimulation 2: stimulation during reach
% stimulation 3: stimulation during reward delay
% stimulation 0: no stimulation

% find stimulation type for each trial and reach that was stimulated
%
% analyze certain aspects of file already collected in the file depending
% on stimulation type

if sum(data.trajectories.contXYsmooth (:,7)) > 1    % if there are ANY trials that are stimulated

  % for trial
    stimulatedTrials         = find(data.trialParams (:, 11) >= 1);         % find stimulated trials
    stimTrialStartTimes      = data.trialTimes(stimulatedTrials,2);% - shift; % trial start of stimulated trials
    stimThreshCrossTimes     = data.trialTimes(stimulatedTrials,3);% - shift; % threshold cross times of stimulated trials 
    stimRewardTimes          = data.trialTimes(stimulatedTrials,4);% - shift; % reward delivery times of stimulated trials
    
    unstimulatedTrials (:,1) =  1:1:size(data.trialParams(:,11),1);
    unstimulatedTrials(stimulatedTrials) = [];                              % find unstimulated trials
    
  % continuous data
    stimIndex = find(diff(data.trajectories.contXYsmooth(:, 7)) >= 1)+1;    % index in continuous data of the stimulation onset 
    stimTimes = data.trajectories.contXYsmooth(stimIndex, 1);               % time of stimulation in the continuous data
    stimType = data.trialParams(stimulatedTrials,11);                       % stimulation types associated with each stimulated trial

  % list all trials as either being stimulated (number 1-3) or unstimulated (0)
    
    stimulationType(unstimulatedTrials,1) = 0;                             
    stimulationType(stimulatedTrials,1)   = stimType;
    
  % check to see if there are any recorded reach starts before first threshold cross- if not, delete that stimulated trial
    
    reachBefFirstStim = find(data.reach.start < stimThreshCrossTimes(1), 1, 'last');
    
    if isempty(reachBefFirstStim) == 1
        stimThreshCrossTimes (1) = []; 
        stimulationType(1)       = [];
        stimTri                  = stimulatedTrials;
        stimulatedTrials(1)      = [];
        stimulatedTrials         = stimulatedTrials - 1;
        disp('Deleted first stimulated trial; No reach prior');
        deleteFirstStim = 1;                                               % if first stimulated trial is deleted, deleteFirstStim is 1
    else
        deleteFirstStim = 0;                                               % if first stimulated trial is NOT deleted, deleteFirstStim is 0
        stimTri                  = stimulatedTrials;
       
    end
     
  % find stimulated reaches- either stim during reach or last reach before stim
  
    if unique(stimulationType(stimulatedTrials)) == 1                      % stimulation type 1
        
        for c = 1: nnz(stimThreshCrossTimes) 
            stimReachIndex(c) = find(data.reach.start < stimThreshCrossTimes(c), 1, 'last'); % last reach before threshold cross 
        end
        
    elseif unique(stimulationType(stimulatedTrials)) == 2                  % stimulation type 2
        
        for c = 1: nnz(stimTimes)
            stimReachIndex(c) = find(data.reach.start < stimTimes(c), 1, 'last'); % stimulated reach
        end

    elseif unique(stimulationType(stimulatedTrials)) == 3                  % stimulation type 3

        for c = 1: nnz(stimThreshCrossTimes) 
            stimReachIndex(c) = find(data.reach.start < stimThreshCrossTimes(c), 1, 'last'); % last reach before threshold cross- rewarded reach
        end
        
    else                                                                   % combination of several different stimulation types 
        stimTypes = unique(stimulationType(stimulatedTrials));
    end 
    
     data.stim.percent    = nnz(stimulatedTrials)/ size(data.trialParams (:, 11),1); % record percentage of trials that are stimulated
     data.stim.plot       = data.trialParams(:,11);                     % binary of trial stim or not to plot
     data.stim.type       = stimulationType;
     data.stim.reachIndex = stimReachIndex; 
     data.stim.trialIndex = stimTri;
     data.stim.time       = data.trajectories.contXY(find(diff(data.trajectories.contXYsmooth(:,7)) == 1)+1,1);
   
else
    data.stim.type = 0;
    data.stim.percent = 0;
    
end


 %% licking analysis
 
%  lickTime = data.trajectories.contXYsmooth((data.trajectories.contXYsmooth(:,9) == 1),1);
%  
%  for k = 10:nnz(rewardTimes)-1                                             % don't use the first 10 trials
%     firstLickPostRew (k-9,1) = find(lickTime < rewardTimes(k), 1, 'last');
%  end
%  
%  % first lick after reward
%   firstLickTimePostRew = lickTime(firstLickPostRew+1);
%   latencyToFirstLick   = firstLickTimePostRew- rewardTimes(10:end-1);
%  
%  % last lick before reward
%  if firstLickPostRew(1) == 1
%      
%     lastLickTimePostRew = lickTime(firstLickPostRew(2:end)-1);
%     latencyBeforeRew    = rewardTimes(2:end-1) - lastLickTimePostRew;
%  
%  else
%      
%      lastLickTimePostRew = lickTime(firstLickPostRew);
%      latencyBeforeRew    = rewardTimes(10:end-1) - lastLickTimePostRew;
%  
%  end
%  
%  data.lick.postRew = latencyToFirstLick;
%  
%  data.lick.befRew  = latencyBeforeRew;
 
%% Write the data structure into large behavior file

    disp(' ');
    disp(' ');
     
    x = struct('data', data, 'reachAnalysis', reachAnalysis);%, 'optimum', optimum);
    behaviorFile(i,1) = x;
    
    %eval([variablename '= x']);  

    disp('Completed');

    clear data tmp reach star enn;
    disp(' Finished with file i' ); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIX THIS DISPLAY THE NUMBER OF FILE COMPLETED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    
disp(' ');

clear data yVals xVals wait2TrialStart wait2RewardDelivery velTraj vVals tstmps trialThresh trialNum trial trajPoints totalTrials trajAngle ...
    totalDisp timeToReward timeToReach stopTime stopR startTimes startTime startRewardR startR rewardTime rewardCloseReach rewardR firstReach firstReachTrial...
    unstimulatedTrials stimulatedTrials stimulationType rewardedReaches previousReward stimulation reachStartWaitIndex timeBtwnReach


end

end

% Stimulation 2 : should be just after the reach start time
        % want to look at WHICH REACHES WERE STIMULATED in the stim2 mode-
        % effect on reaching parameters, of the reward delay post
        % stimulation, which reach within the trial was stimulated, if
        % first and there are subsequent reaches, what are the effects 
        % reward rate
        
        
        % need max number of reaches per trial

% Stimulation 3
        % vigor of next reaches
        % delay time accuracy
        % reward rate
    
    % Stimulation 1 
        % iti reach time accuracy 
        % reward rate
    