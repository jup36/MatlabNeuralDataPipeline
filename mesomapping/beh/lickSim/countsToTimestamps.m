function timeStamps = countsToTimestamps(simulCounts, minMaxTime, expd)
timeStamps = [];
intervals = random(expd, sum(simulCounts), 1);
intSum = cumsum(intervals);

if isempty(intSum)
    timeStamps = [];
else
    if intSum(end) > diff(minMaxTime)
        intervals = intervals(intSum <= diff(minMaxTime));

        if isempty(intervals)
            intervals = 0;
        else
            for i = 1:100
                intervals = [intervals; random(expd, 1)];
                intSum = cumsum(intervals);
                intervals = intervals(intSum <= diff(minMaxTime));
            end
            if length(intervals) > length(simulCounts)
                intervals = intervals(1:length(simulCounts));
            end
        end
    end
    timeStamps = minMaxTime(1) + intSum;
    timeStamps = timeStamps(timeStamps < minMaxTime(2));
end
end