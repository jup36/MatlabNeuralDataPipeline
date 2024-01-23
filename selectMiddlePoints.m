function middlePoints = selectMiddlePoints(array)
    % Calculate the length of the array
    arrayLength = length(array);

    % Calculate the number of points to exclude from each end
    % This calculation ensures that the total number of excluded points is 6
    numPointsToExclude = ceil((arrayLength - 1600) / 2);

    % Calculate the starting and ending indices
    startIndex = numPointsToExclude + 1;
    endIndex = arrayLength - numPointsToExclude;

    % Select the middle points
    middlePoints = array(startIndex:endIndex);
end
