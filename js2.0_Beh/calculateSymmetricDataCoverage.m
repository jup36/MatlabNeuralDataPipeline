function dataCoverage = calculateSymmetricDataCoverage(data, step)
    % Calculate symmetric data coverage ranges centered at zero.
    %
    % Parameters:
    %   data - An array containing the data points.
    %   step - The step size to increase the range symmetrically.
    %
    % Returns:
    %   dataCoverage - A matrix where each row represents the range as
    %                  [lower_bound, upper_bound, coverage percentage].

    % Initialize variables
    maxDataPoint = max(abs(data));  % Maximum absolute value in data for range determination
    currentRange = 0;               % Start with zero range
    coverageIndex = 0;              % Index for filling dataCoverage
    coveredData = 0;                % Initialize covered data count
    totalDataPoints = numel(data);  % Total number of data points
    
    % Pre-allocate a large matrix for ranges and coverage - assume we won't have more than 1000 steps
    dataCoverage = zeros(1000, 3);

    % Loop to expand range until all data points are covered
    while currentRange < maxDataPoint
        currentRange = currentRange + step;  % Increment range by step size
        coverageIndex = coverageIndex + 1;   % Increment index
        
        % Calculate the number of data points within the current range
        numDataInRange = sum(data >= -currentRange & data <= currentRange);
        percentageCoverage = 100 * numDataInRange / totalDataPoints;  % Calculate coverage percentage
        
        % Store the range and coverage in the matrix
        dataCoverage(coverageIndex, :) = [-currentRange, currentRange, percentageCoverage];
        
        % Update covered data count and check if all data is covered
        if numDataInRange == totalDataPoints
            break;  % Exit loop if all data points are included in the current range
        end
    end
    
    % Trim unused pre-allocated rows
    dataCoverage = dataCoverage(1:coverageIndex, :);
end
