fileDirectory = '/Users/parkj/Documents/BehavioralData/373962 #8/032317/';
cd(fileDirectory)

allFiles = dir(sprintf('%s*_p.csv',fileDirectory));

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
    %%data.name.mouseType   = mousetype; % NA for me now
    
    disp(' ');
    disp(['File ' num2str(i) ' base name: ' fileNameBase]);
    disp(['Date created: ' year '/' month '/' day]);
    disp(' ');
    
    %% Script to load csv behavior data
    clear disp* trial*
    
    skipInitTrials = 0;
    
    fileNameStr = [fileDirectory fileBase '_p.csv'];
    data.trialParams = dlmread(fileNameStr,',',1+skipInitTrials,0); % 11th column is the column containing information about which trials were stimulated
    
    fileNameStr = [fileDirectory fileBase '.csv'];
    data.trialTimes = dlmread(fileNameStr,',',2+skipInitTrials,0);  % time information
    
    fileNameStr = [fileDirectory fileBase 'pXY.csv'];
    
    if exist(fileNameStr,'file')== 1 || exist(fileNameStr,'file')== 2 % 'exist' checks existence of variable, script, or file
        
        data.trajectories.contXY = dlmread(fileNameStr,',',2,0);    % extract semicontinuous data into data.trajectories.contXY file
        
    end
    
    % data.trajectories.contXY format
    % data.trajectories.contXY (:,1) = time
    % data.trajectories.contXY (:,2) = x
    % data.trajectories.contXY (:,3) = y
    % data.trajectories.contXY (:,4) = licking data
    % data.trajectories.contXY (:,5) = this must be the behavioral state (paradigmmode)
    % data.trajectories.contXY (:,6) = this must be solenoid (see that it corresponds to the paradigmMode4 in col5)
    
    %% Remove duplicate trials (or I can fix the processing file to not to create duplicates)
    
    
    
    %% General information about the file
    
    indArray = 2:size(data.trialTimes,1); % skipping the 1st row?
    
    % data.trialTimes ( 1: count, 2: trialStart, 3: crossTime, 4: trialEnd, 5: 1stLickTime, 6: valveOpenTime, 7: leftTrial )
    data.startToEvent  = data.trialTimes(indArray,3) - data.trialTimes(indArray,2);     % time from trial start to threshold cross
    data.endToStart    = data.trialTimes(indArray,2) - data.trialTimes(indArray-1,4);   % inter-trial interval ()
    data.eventToEnd    = data.trialTimes(indArray,4) - data.trialTimes(indArray,3);     % reward delay
    data.threshold.max = max(data.trialParams (:,6));                                   % max threshold
    data.threshold.min = min(data.trialParams(:,6));                                    % min threshold
    numTrials          = size(data.trialTimes,1);                                       % number of total trials
    
    % program starts with ITI- need to remove that with a shift in order to compare continuous data and trial times file
    behavStateDiff      = diff(data.trajectories.contXY(:,5)); % this must be column5 where the behavioral state data are saved
    
    reachTriStartIndex2 = find(behavStateDiff == 2); % this should detect the shift of paradigmmode from 0 to 2 which indicates start of the trial
    reachTriStartIndex3 = find(behavStateDiff == 3); % no such case
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
    
    % why is there shift? In each trial, there's time lag between
    % processing and arduino with the time marked in processing has greater
    % values perhaps due to 1) serial communication, and 2) different clocks started at different timing between A and P.
    % Correction for this is made below by adding the lag corresponding to each trial to the time marked in Arduino
    shift       = contStartTimes - data.trialTimes(:,2);            % data.trialTimes(:,2) contains the trial start time triggered by the controller
    startTimes  = data.trialTimes(:,2) + shift;                     % trial start times
    eventTimes  = data.trialTimes(:,3) + shift;                     % threshold cross times
    rewardTimes = data.trialTimes(:,4) + shift;                     % reward times
    
    %% Extract kinematic information from the semicontinuous trajectory
    % Do the following for semicontinuous (XY) data files
    % data.trajectories.contXY (1: time (processing), 2: x traj, 3: y traj, 4: licking, 5: behavioral state, 6: solenoid)
    
    if exist(fileNameStr,'file')== 1 || exist(fileNameStr,'file')== 2
        
        % Smooth continuous position data
        
        for j=2:3 % two columns containing lever position information
            data.trajectories.contXYsmooth(:,j) = sgolayfilt(data.trajectories.contXY(:,j),3,11); % sgolayfilt(x,order,framelen) Savitzky-Golay filtering
        end
        
        % Convert au's to mm
        
        beginTrial = find(data.trajectories.contXY(:,5) == 2, 1, 'first'); % first trial index
        
        data.trajectories.contXYsmooth(beginTrial:end,2) = tan(data.trajectories.contXYsmooth(beginTrial:end,2).*(.1*(pi./180)))*90; % my joystick is 90 mm long, .1 (what is this scalar for?), pi./180 (degree to radian conversion)
        data.trajectories.contXYsmooth(beginTrial:end,3) = tan(data.trajectories.contXYsmooth(beginTrial:end,3).*(.1*(pi./180)))*90; % my joystick is 90 mm long
        
        data.trajectories.contXYsmooth(1:beginTrial-1,2) = 0;
        data.trajectories.contXYsmooth(1:beginTrial-1,3) = 0;
        
        % Keep some data the same in columns and don't smooth
        
        data.trajectories.contXYsmooth(:,1)  = data.trajectories.contXY(:,1); % time
        data.trajectories.contXYsmooth(:,7)  = data.trajectories.contXY(:,6); % solenoid open time
        data.trajectories.contXYsmooth(:,8)  = data.trajectories.contXY(:,5); % behavioral state
        data.trajectories.contXYsmooth(:,9)  = data.trajectories.contXY(:,4); % licking (void in training data)
        
        % Fill in Associated Kinematics
        instTime = diff(data.trajectories.contXYsmooth(:,1)); % instantiation of time
        
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
        
    end
    
    %% Isolate reaches from smoothed velocity function
    
    % Threshold velocities; must be >.005mm/ms
    %         ABOVE    = find(data.trajectories.sLeverVm > 0.008);
    
    ABOVE    = find(data.trajectories.sLeverVm > 0.005); % find the points above the velocity threshold
    tmp.star = find([10; diff(ABOVE)]>1); % > 1 is to find above velocity points belonging to different reaches
    tmp.enn  = [tmp.star(2:end)-1; size(ABOVE,1)]; % tmp.enn corresponds to the next point with the velocity greater than .005
    
    tmp.reachIND = [tmp.star tmp.enn]; % reach indicators
    
    % Delete reaches with approx <100 ms (4 samples = 132ms) above threshold long
    remove = [];
    
    for f = 1: size(tmp.reachIND,1) % increment along all temporary high velocity points
        
        if tmp.reachIND(f,2) - tmp.reachIND(f,1) + 1 <= 4 % two points being within 100ms from each other indicate short reaches (<100ms)
            remove = [remove; f]; %
        end
        
    end
    
    tmp.reachIND(remove,:) = [];
    
    % Merge reaches in the case that next reach start < previous reach stop and
    % Merge reaches where there is approx. <264ms between reaches, given 30Hz
    % rate this is about 8 (264ms) samples
    
    remove = [];
    
    for f = 2: size(tmp.reachIND,1)
        
        if ABOVE(tmp.reachIND(f,1)) - ABOVE(tmp.reachIND(f-1,2)) <= 8 % in case the current reach occurred <264ms after the previous reach end
            tmp.reachIND(f-1,2) = tmp.reachIND(f,2); % merge them together as they're most likely to be the same reach
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
            A = find(data.trajectories.sLeverVm(startInds,:) <= 0.0010, 1, 'last'); % define the start of reach (A) as the last point where the velocity didn't exceed the threshold (0.001)
            
            if isempty(A) == 1
                A = beginTrial;
            end
            
            A = startInds(A);
        end
        
        stopInds = ABOVE(tmp.enn(f));
        B = find(data.trajectories.sLeverVm(stopInds:end,:) <= 0.0010, 1, 'first') + ABOVE(tmp.enn(f)) - 1; % define the end of reach (B) as the first point where the velocity falls back lower than the threshold (0.001)
        
        star = [star; A]; % reach start
        enn = [enn; B];   % reach end
    end
    
    if nnz(star) ~= nnz(enn) % in case there's a mismatch between the numbers of start and end, match the numbers by curtailing the excess
        
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
    
    
end