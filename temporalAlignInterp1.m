function  [interpolatedC, sparseTime] = temporalAlignInterp1(dataC, timeC)
% timeC = rwdDffC_motor(:, 2); 
% dataC = rwdDffC_motor(:, 1); 

% Identify empty cells
%dataI = find(~cellfun(@isempty, dataC));
timeI = find(~cellfun(@isempty, timeC));

% Assuming timestamps for all trials are stored in a cell array: timeC
minTime = min(cell2mat(cellfun(@nanmin, timeC(timeI), 'UniformOutput', false)));
maxTime = max(cell2mat(cellfun(@nanmax, timeC(timeI), 'UniformOutput', false)));

% Define a high-resolution time grid
numHighResPoints = 10000;  % or any other desired resolution
highResTime = linspace(floor(minTime), ceil(maxTime), numHighResPoints);
sparseTime = round(minTime*10)/10:0.05:round(maxTime*10)/10;

% Assuming data for all trials are stored in a cell array: dataC
numTrials = length(dataC);

% Preallocate memory for interpolated data
interpolatedC = cell(numTrials, 1);

for t = 1:numTrials
    if ~isempty(dataC{t})
        tmphigh = interp1(timeC{t}, dataC{t}, highResTime, 'linear', 'extrap');
        interpolatedC{t, 1} = interp1(highResTime, tmphigh, sparseTime, 'nearest');
    end
end
