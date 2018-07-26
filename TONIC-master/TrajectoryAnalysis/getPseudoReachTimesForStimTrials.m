function [ reachStart, reachStop, reachMW, pos1, valReachDurIdx, xpos1, ypos1 ] = getPseudoReachTimesForStimTrials( positionData, stmLaser, minReachDurCut )
%This function detects pseudo-reach trajectories around the time of laser
% stimulation deliveries.
% modified on 6/10/18 to eliminate reaches whose duration < minReachDurCut

%TbetweenR = 2000; % minimum time between reaches, in ms
maxDur=1500;      % max reach duration

xPos = positionData(1,:); % Xposition
yPos = positionData(2,:); % Yposition

reachSG = sgolayfilt(sqrt((yPos-median(yPos)).^2 + (xPos-median(xPos)).^2),3,201); % smoothing with sgolay filter
reachMW = moveavg(sqrt((yPos-median(yPos)).^2 + (xPos-median(xPos)).^2),10);       % smoothing with moving average filter
XposMW  = moveavg(xPos,10); % smoothing with moving average filter
YposMW  = moveavg(yPos,10); % smoothing with moving average filter

vel=[0, diff(reachSG)]; % velocity - differentiation of the position trajectory

sortPos=sort(reachSG);  % sort position data
%sortVel=sort(vel);      % sort velocity data
sortPos(sortPos<median(sortPos))=median(sortPos); % replace less-than median position values with the median position value
cutoff  = sortPos(round(length(sortPos)*.90));    % set the 90% cutoff for the position data
%cutoffd = sortVel(round(length(sortVel)*.98));    % set the 98% cutoff for the velocity data (velocity cutoff is required, since reaches are expected to occur at high velocity)

%% detect reach start
reachStart=[]; % the goal here is to find the near-zero point that is immediately prior to each high-velocity point
for i  = 1:length(stmLaser)
    if stmLaser(i)+1000 < length(reachSG)
        currHighAmpPoints = find(reachSG(1,stmLaser(i):stmLaser(i)+1000)>cutoff) + stmLaser(i);  % current valid sample (high-amplitude data point)
        
        if ~isempty(currHighAmpPoints) % in case there's a point(s) above the position cutoff
            
            currHighAmpPointsIdx = zeros(length(currHighAmpPoints),1);
            for ii = 1:length(currHighAmpPoints)
                if currHighAmpPoints(ii)+200<length(reachSG)
                    currHighAmpPointsIdx(ii,1) = sum(reachSG(1,currHighAmpPoints(ii):currHighAmpPoints(ii)+200)>cutoff)>20; % enforce to cross the threshold for at least 20 ms
                else 
                    currHighAmpPointsIdx(ii,1) = sum(reachSG(1,currHighAmpPoints(ii):end)>cutoff)>20; % enforce to cross the threshold for at least 20 ms
                end
            end
            if ~isempty(find(currHighAmpPointsIdx==1,1))
                currSamp = currHighAmpPoints(find(currHighAmpPointsIdx==1,1));
            else
                [~,currSamp] = max(reachSG(1,stmLaser(i):stmLaser(i)+1000),[],2);
                currSamp = currSamp + stmLaser(i);
            end
            
        else % in case there's no point above the position cutoff
            [~,currSamp] = max(reachSG(1,stmLaser(i):stmLaser(i)+1000),[],2); % just take the max position
            currSamp = currSamp + stmLaser(i);
        end
        
        nearZero= find(reachSG(1:currSamp)<cutoff*.4,1,'last'); % find the most recent near zero-point
        if nearZero>10
            if currSamp-nearZero>1000  % sometimes the joystick doesn't reset to zero and induces 20 second long reaches. this prevents that by setting the second (additional) cutoff
                sortPos2=sort((reachSG(nearZero:currSamp)));
                cutoff2 = sortPos2(round(length(sortPos2)*.90)); % reset the cutoff2 to be a point immediatly before the reach peak (farthest point)
                nearZero=find(reachSG(1:currSamp)<cutoff2,1,'last');
            end
            realStartMinus10=find(vel(1:nearZero)<0,10,'last'); % find 10 elements whose velocity is negative prior to the current near zero point
        else % if there's no near zero point
            realStartMinus10=currSamp; % take the current valid sample (high-velocity data point) as a reach start point
        end
        reachStart = [reachStart, realStartMinus10(1)];
    end
    
    % this will eliminate pseudo starts, you might not want to do eliminate them
    pos1=[];    % reach position data aligned to reachStart
    xpos1=[];   % x position data aligned to reachStart normalized to the x position at the reach start
    ypos1=[];   % y position data aligned to reachStart normalized to the y position at the reach start
    for i = reachStart
        if i > 1000
            if i+1499 > length(reachMW) % in case the reach trajectory end exceeds the length of the time series (this must happen only for the last reach)
                pos1 =[pos1; nan(1,size(pos1,2))]; % position data corresponding to each reach from 200-ms before to 1500-ms after each reach start
                xpos1=[xpos1; nan(1,size(xpos1,2))];
                ypos1=[ypos1; nan(1,size(ypos1,2))];
                
            else % in case the reach trajectory is within the length of the time series
                pos1  = [pos1; reachMW(i-200:i+1499)-reachMW(i)]; % position data corresponding to each reach from 200-ms before to 1500-ms after each reach start
                xpos1 = [xpos1; XposMW(i-200:i+1499)-XposMW(i)];  % x position data aligned to reachStart normalized to the x position at the reach start
                ypos1 = [ypos1; YposMW(i-200:i+1499)-YposMW(i)];  % y position data aligned to reachStart normalized to the y position at the reach start
            end
        end
    end
    % figure; imagesc(pos1); colorbar
    % title('start times, aligned at 200ms')
    
    %% detect reach stop
    reachStop=[]; % the goal is to find the reach below the threshold for a certain duration of time (150 ms)
    for i = 1:length(reachStart)
        
        if reachStart(i)+maxDur < length(reachMW)
            currPos= reachMW(reachStart(i):reachStart(i)+maxDur); % take the smoothed position data corresponding to each reach
            
            [~,reachPeakTime]=max(currPos(1:700));   % detect reach peak (the farthest reach point)
            peakToEndPos=currPos(reachPeakTime:end); % peak-to-end of reach
            sortStop=sort(peakToEndPos);
            binSize=5;  % 10 may be better
            jsBinned=hist(peakToEndPos,binSize); %jsBinned=hist(peakToEndPos,binSize); % jsBinned=hist(peakToEndPos,min(peakToEndPos):binSize:max(peakToEndPos));
            [~,cutoffStop1]=max(jsBinned(1:round(length(jsBinned)/3)));
            cutoffStop = sortStop(jsBinned(cutoffStop1)); % min(peakToEndPos) + sortStop(jsBinned(cutoffStop1)); % the way I think the cutoffStop should be defined!
            %cutoffStop=cutoffStop1*binSize+10+min(peakToEndPos); % cutoffStop1
            
            N=150; % # of ms consecutively below threshold needed to detect a stable reach stop point
            t = find(peakToEndPos<cutoffStop);  % find the data points below the cutoffStop
            x = diff(t)==1;            % spot out the points consecutively below the cutoffStop
            f = find([false,x]~=[x,false]); % spot out the non-consecutive data points
            g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first'); % find the first data point, at which the interval between non-consecutive points is greater than 150 ms, as this means that there are 150 or more points consecutively under the reachstop threshold
            almostEnd=t(f(2*g-1))-1;   % just tranlate g to a point on t
            if isempty(almostEnd)    % in case there's no such point
                almostEnd=find(peakToEndPos<cutoffStop,1);   % just take the point below the reach threshold
                %reachStart(i);
            end
            %thisStop= find(moveavg(diff(peakToEndPos(almostEnd:end)),10)>=0,1);
            reachStop=[reachStop, reachPeakTime+almostEnd+reachStart(i)]; %thisStop+reachPeakTime+almostEnd+reachStart(i)];
            
        else
            reachStart = reachStart(1:end-1); 
        end
    
    end
        
end

    % eliminate reaches shorter than the minReachDurCut
    valReachDurIdx = reachStop-reachStart>minReachDurCut;
    reachStart = reachStart(1,valReachDurIdx);
    reachStop = reachStop(1,valReachDurIdx);
    pos1  = pos1(valReachDurIdx,:);
    xpos1 = xpos1(valReachDurIdx,:);
    ypos1 = ypos1(valReachDurIdx,:);


end

