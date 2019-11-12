function [reachStart, reachStop] = getReachTimesJs2p0( posTrace, velTrace, posCut, velCut )

%posTrace = hTrjSDstBase{7}; 
%velTrace = hTrjSVelBase{7}; 

sampRate = 250; % the sampling rate of videography (default: 250Hz)
timeInt = 1000/sampRate; % time interval in ms 
maxDur = 1500; 
fastPts = find(velTrace>velCut); % points exceeding the velocity cutoff
valFastPts = []; 

%% detect reaches
for i = 1:length(fastPts) % all points exceeding the velocity cutoff
    if isempty(valFastPts) % in case of the 1st high-velocity data point
        if sum(posTrace(fastPts(i):min(fastPts(i)+(200/timeInt), length(posTrace)))>posCut)>(15/timeInt) % this part is just to isolate fast reaches from the short jerks
            valFastPts=[valFastPts, fastPts(i)]; 
        end
    else    % from the 2nd high-velocity data point on
        if (fastPts(i)-fastPts(i-1))>(50/timeInt) % if the current high-velocity point came 80 ms (note that the sampling rate was 250Hz) or more after the previous high-velocity point
            if (fastPts(i)-valFastPts(end)) > (50/timeInt) % ensure the interval between reaches is greater than the set interval threshold (e.g. 2000 ms)
                if sum(posTrace(fastPts(i):min(fastPts(i)+(200/timeInt), length(posTrace)))>posCut)>(15/timeInt) % finally ensure that it is a well-timed reach not a short jerk
                    valFastPts=[valFastPts, fastPts(i)];
                end
            end
        end
    end
end

%% detect reachStart
reachStart=[]; % the goal here is to find the near-zero point that is immediately prior to each high-velocity point  
for i = 1:length(valFastPts)
    currSamp = valFastPts(i); % current valid sample (high-velocity data point) 
    nearZero= find(posTrace(1:currSamp)<posCut*.4,1,'last'); % find the most recent near zero-point
    if nearZero>10
        if currSamp-nearZero>1000  % sometimes the joystick doesn't reset to zero and induces 20 second long reaches. this prevents that by setting the second (additional) posCut
            sortPos2=sort((posTrace(nearZero:currSamp)));
            posCut2 = sortPos2(round(length(sortPos2)*.90)); % reset the posCut2 to be a point immediatly before the reach peak (farthest point)
            nearZero=find(posTrace(1:currSamp)<posCut2,1,'last');
        end
        realStartMinus10=find(velTrace(1:nearZero)<0,10,'last'); % find 10 elements whose velocity is negative prior to the current near zero point  
    else % if there's no near zero point
        realStartMinus10=currSamp; % take the current valid sample (high-velocity data point) as a reach start point
    end
    reachStart = [reachStart, realStartMinus10(1)];
end
clearvars i 

if length(reachStart)>=2 
   reachStart=reachStart([1,find(diff(reachStart)>(100/timeInt))+1]); % ensure that intervals between reaches are greater than 100 ms
end
reachStart=unique(reachStart)+10; 

% detect reachStop
reachStop=[]; 
for i = 1:length(reachStart) 
    currPos = posTrace(reachStart(i):min(reachStart(i)+floor(maxDur/timeInt), length(posTrace))); % take the portion relevant to current reach
    [~,reachPeakPt]=max(currPos(1:floor(length(currPos)/2))); % detect reach peak (the farthest reach point)
    peakToEnd = currPos(reachPeakPt:end); % peak-to-end of reach  
    sortStop=sort(peakToEnd); 
    peakToEndBin = hist(peakToEnd,5); 
    [~,cutoffStop1]=max(peakToEndBin(1:round(length(peakToEndBin)/3)));
    cutoffStop = sortStop(peakToEndBin(cutoffStop1));
    
    tb = find(peakToEnd<cutoffStop); 
    xx = diff(tb)==1; 
    ff = find([false,xx]~=[xx,false]); % spot out the non-consecutive data points 
    gg = find(ff(2:2:end)-ff(1:2:end-1)>=(100/timeInt),1,'first');
    almostEnd=tb(ff(2*gg-1))-1;   % just tranlate g to a point on t  
    if isempty(almostEnd)   % in case there's no such point
        almostEnd=find(peakToEnd<cutoffStop,1);   % just take the point below the reach threshold
        %reachStart(i);
    end
    %thisStop= find(moveavg(diff(peakToEndPos(almostEnd:end)),10)>=0,1);
    reachStop=[reachStop, reachPeakPt+almostEnd+reachStart(i)]; %thisStop+reachPeakTime+almostEnd+reachStart(i)];
end
