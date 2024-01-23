function middleArrayI = middleIndex(originalArrayLength, numMiddlePoints)
%This function gets the index to select the middle portion of
%   the originalArrayLength that corresponds to
%   numMiddlePoints.

if originalArrayLength >= numMiddlePoints
    % originalArrayLength = length(targetFaceTs);
    % Calculate start and end indices
    startIndex = ceil((originalArrayLength - numMiddlePoints) / 2) + 1;
    endIndex = startIndex + numMiddlePoints - 1;

    middleArrayI = zeros(originalArrayLength, 1);

    % Select the middle points
    middleArrayI(startIndex:endIndex)=1;
else
    error("The specified middle points exceeds the original array length!")
end
end
