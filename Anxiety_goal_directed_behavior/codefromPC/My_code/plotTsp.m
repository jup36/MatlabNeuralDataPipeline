function fig = plotTsp(tsp, t0, t1, msize)
%% raster plot given the time of spikes
numTrial = length(tsp);
[minT, maxT] = minmaxT(tsp, numTrial);
if ~exist('t0', 'var')
    t0 = floor(minT*10)/10;
end

if ~exist('t1', 'var')
    t1 = ceil(maxT*10)/10;
end

if ~exist('msize', 'var')
    msize = 6;
end

fig = figure;
axis([t0, t1, 0, numTrial]);
hold on;
for m=1:numTrial
    temp = tsp{m};
    plot(temp, (temp>0)*m, '.', 'markersize',msize);
end
xlabel('Time (sec.)');
ylabel('# Neuron');

function [minT, maxT] = minmaxT(tsp, numTrial)
%% find the minimum and maximum t in this time series
minV = zeros(numTrial, 1);
maxV = zeros(numTrial, 1);
for m=1:numTrial
    if isempty(tsp{m})
        minV(m) = nan;
        maxV(m) = nan;
    else
        minV(m) = min(tsp{m});
        maxV(m) = max(tsp{m});
    end
end

minT = min(minV);
maxT = max(maxV);