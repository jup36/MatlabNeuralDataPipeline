function [interpolatedC, sparseTime] = temporalAlignInterp3D(dataC, timeC, timeInt)
% timeC = cell array of time vectors for each trial
% dataC = cell array of 3D matrices for each trial, with the 3rd dimension being time
% timeInt = desired time interval for the output sparse time grid

% Identify non-empty time cells
timeI = find(~cellfun(@isempty, timeC));

% Determine the overall time range across all trials
minTime = min(cell2mat(cellfun(@nanmin, timeC(timeI), 'UniformOutput', false)));
maxTime = max(cell2mat(cellfun(@nanmax, timeC(timeI), 'UniformOutput', false)));

% Define the sparse time grid
sparseTime = round(minTime*10)/10:timeInt:round(maxTime*10)/10; 

% Number of trials
numTrials = length(dataC);

% Preallocate memory for interpolated data
interpolatedC = cell(numTrials, 1);

for t = 1:numTrials
    if ~isempty(dataC{t}) && ~isempty(timeC{t})
        % Get the dimensions of the 3D data matrix
        [rows, cols, ~] = size(dataC{t});
        
        % Preallocate a matrix for the interpolated data
        interpolatedData = zeros(rows, cols, length(sparseTime));
        
        % Interpolate each pixel independently
        for r = 1:rows
            for c = 1:cols
                % Extract the time series for this pixel
                pixelTimeSeries = squeeze(dataC{t}(r, c, :));
                
                % Interpolate the time series to the sparse time grid
                interpolatedData(r, c, :) = interp1(timeC{t}, pixelTimeSeries, sparseTime, 'linear', NaN);
            end
        end
        
        % Assign the interpolated data to the output cell array
        interpolatedC{t} = interpolatedData;
    end
end
end
