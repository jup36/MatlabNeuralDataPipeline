function [ reachStart, reachStop, reachMW, pos1, pos2, xpos1, ypos1, xpos2, ypos2 ] = getReachTimesJP( positionData )
%Extracts best reach start times from .ns2 data unless BabData set to 1 for
%processed ContData files. 
% Can be made to extract from other data types. Assumes 1kHz sampling
%   Assumes  X and Y position are channels 2 and 3 (line 13)

%outputs:
%  reachStart and reachStop are vectors of reach times
%  reachMW is the amplitude readout of the whole session
%  pos1 = all reach traces, aligned to start
%  pos2 = all reach traces, aligned to stop
% 
%  sgolay injects a negative dip before a reach as an artifact of the filter.
%  we can use that for reach detection, but makes for poor position traces
%  sgolay filter is a convolution filter which fits successive sub-sets of
%  adjacent data points with a low-degree polynomial by the method of
%  linear least squares. 

TbetweenR = 2000; % minimum time between reaches, in ms
maxDur=1500;      % max reach duration

xPos = positionData(1,:); % Xposition 
yPos = positionData(2,:); % Yposition

reachSG = sgolayfilt(sqrt((yPos-median(yPos)).^2 + (xPos-median(xPos)).^2),3,201); % smoothing with sgolay filter
reachMW = moveavg(sqrt((yPos-median(yPos)).^2 + (xPos-median(xPos)).^2),10);       % smoothing with moving average filter
XposMW  = moveavg(xPos,10); % smoothing with moving average filter
YposMW  = moveavg(yPos,10); % smoothing with moving average filter

vel=[0, diff(reachSG)]; % velocity - differentiation of the position trajectory

sortPos=sort(reachSG);  % sort position data
sortVel=sort(vel);      % sort velocity data 
sortPos(sortPos<median(sortPos))=median(sortPos); % replace less-than median position values with the median position value  
cutoff  = sortPos(round(length(sortPos)*.90));    % set the 90% cutoff for the position data
cutoffd = sortVel(round(length(sortVel)*.98));    % set the 98% cutoff for the velocity data (velocity cutoff is required, since reaches are expected to occur at high velocity)
allValSamps1=find(vel>cutoffd);   % velocities exceeding the cutoff
allValSamps1=allValSamps1(allValSamps1<(length(reachSG)-300)); % this line is just to exclude the data points at the tail
allValidSamps=[];

for jj = 1:length(allValSamps1)   % all points exceeding the position cutoff
    if length(allValidSamps)==0   % in case of the 1st high-velocity data point
        if sum(reachSG(allValSamps1(jj):allValSamps1(jj)+200)>cutoff)>20 % this part is just to sort out high-velocity reaches from the short jerks which is not really a reach
           allValidSamps=[allValidSamps, allValSamps1(jj)];              % if it satisfies the condition for a legitimate reach, then register it as a reach
        end
    else    % from the 2nd high-velocity data point on 
        if (allValSamps1(jj)-allValSamps1(jj-1))>100 % if the current high-velocity point came 100 ms or more after the previous high-velocity point; this will prevent the redundant counting from a continued reach
            if (allValSamps1(jj)-allValidSamps(end)) > TbetweenR % this is to ensure that the interval between reaches is greater than the set interval threshold (e.g. 2000 ms)
                if sum(reachSG(allValSamps1(jj):allValSamps1(jj)+200)>cutoff)>20 % finally ensure that it is a well-timed reach not a short jerk
                    allValidSamps=[allValidSamps, allValSamps1(jj)];             % if it satisfies conditions for a legitimate reach, then register it as a reach
                end
            end
        end
    end
end

%% detect reach start
reachStart=[]; % the goal here is to find the near-zero point that is immediately prior to each high-velocity point  
for i  = 1:length(allValidSamps)
    currSamp = allValidSamps(i);                       % current valid sample (high-velocity data point) 
    nearZero= find(reachSG(1:currSamp)<cutoff*.4,1,'last'); % find the most recent near zero-point
    if nearZero>10
        if currSamp-nearZero>1000  % sometimes the joystick doesn't reset to zero and induces 20 second long reaches. this prevents that by setting the second (additional) cutoff
            sortPos2=sort((reachSG(nearZero:currSamp)));
            cutoff2 = sortPos2(round(length(sortPos2)*.90)); % reset the cutoff2 to be a point immediatly before the reach peak (farthest point)
            nearZero=find(reachSG(1:currSamp)<cutoff2,1,'last');
        end
        realStartMinus10=find(vel(1:nearZero)<0,10,'last'); % find 10 elements whose velocity is negative prior to the current near zero point  
        if isempty(realStartMinus10)
            realStartMinus10=currSamp; % take the current valid sample (high-velocity data point) as a reach start point
        end
    else % if there's no near zero point
        realStartMinus10=currSamp; % take the current valid sample (high-velocity data point) as a reach start point
    end
    reachStart = [reachStart, realStartMinus10(1)];
end

reachStart=reachStart([1,find(diff(reachStart)>TbetweenR)+1]); % ensure that intervals between reaches are greater than 2 sec
reachStart=unique(reachStart)+10; % unique ensures that there's no redundant reaches

% this will eliminate pseudo starts, you might not want to do eliminate them
pos1=[];    % reach position data aligned to reachStart 
xpos1=[];   % x position data aligned to reachStart normalized to the x position at the reach start
ypos1=[];   % y position data aligned to reachStart normalized to the y position at the reach start
for i = reachStart
    if i > 1000
        if i+1499 > length(reachMW) % in case the reach trajectory end exceeds the length of the time series (this must happen only for the last reach)
            pos1=[pos1; nan(1,200+1499+1)]; % just put NaNs of the same length
            xpos1=[xpos1; nan(1,200+1499+1)]; % just put NaNs of the same length
            ypos1=[ypos1; nan(1,200+1499+1)]; % just put NaNs of the same length
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
    end
    
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
    if length(almostEnd)==0    % in case there's no such point
        almostEnd=find(peakToEndPos<cutoffStop,1);   % just take the point below the reach threshold
        %reachStart(i);
    end
    %thisStop= find(moveavg(diff(peakToEndPos(almostEnd:end)),10)>=0,1);
    reachStop=[reachStop, reachPeakTime+almostEnd+reachStart(i)]; %thisStop+reachPeakTime+almostEnd+reachStart(i)];
end

pos2=[];  % reach position data aligned to reachStop 
xpos2=[]; % x position data aligned to reachStart normalized to the x position at the reachStop
ypos2=[]; % y position data aligned to reachStart normalized to the y position at the reachStop
for i =reachStop
    if i>1500 && i < (length(reachMW)-200)
        pos2  = [pos2; reachMW(i-1499:i+200)-reachMW(i)];
        xpos2 = [xpos2; XposMW(i-1499:i+200)-XposMW(i)];
        ypos2 = [ypos2; YposMW(i-1499:i+200)-YposMW(i)];
    else
        pos2  =[pos2; nan(1,200+1499+1)]; % just put NaNs of the same length
        xpos2 =[xpos2; nan(1,200+1499+1)]; % just put NaNs of the same length
        ypos2 =[ypos2; nan(1,200+1499+1)]; % just put NaNs of the same length
        
    end
end

% figure; imagesc(pos2); colorbar
% title('stop times, aligned at 1500ms')

% figure; plot(reachMW);
% hold on; plot(reachStart, reachMW(reachStart),'g*');
% hold on; plot(reachStop, reachMW(reachStop),'r*');
% title('JS position with start and stop times');

% Trial-by-trial inspection 
% reachNumb = 152; % trial
% figure; plot(pos1(reachNumb,:)); hold on; 
% plot(200, 0, 'g*'); % the reach start point is 200th on the pos1, and set to 0 (line 82)
% plot(200 + reachStop(reachNumb) - reachStart(reachNumb), reachMW(reachStop(reachNumb)) - reachMW(reachStart(reachNumb)), 'r*'); % the reach start point is 200th on the pos1, and set to 0 (line 82)

end