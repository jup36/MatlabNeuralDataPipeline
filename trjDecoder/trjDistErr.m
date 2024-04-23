function error = trjDistErr(trj, trjEst)
%trj = rTrj;
%trjEst = rTrjEst;

% Calculate the squared differences
sqDiffs = (trj - trjEst).^2;

% Take the square root to get the Euclidean distance
euD = sqrt(sqDiffs);

% Compute basic statistics
error.mean = mean(euD, 2);
error.median = median(euD, 2);
error.std = std(euD, 0, 2);
end
