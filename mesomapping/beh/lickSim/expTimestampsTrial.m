function timeStamps = expTimestampsTrial(expd, minMaxTime)
minTime = min(minMaxTime);
maxTime = max(minMaxTime);

curTime = minTime; % starting point
intvs = []; % generated intervals
while curTime <= maxTime
    intv = random(expd, 1, 1);
    curTime = curTime + intv;
    if curTime <= maxTime
        intvs = [intvs; intv];
    end
end
% generate cumulative timestamps
timeStamps = cumsum(intvs) + minTime;
end