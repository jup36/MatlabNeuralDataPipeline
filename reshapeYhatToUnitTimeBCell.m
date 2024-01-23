function unitTimeBCell = reshapeYhatToUnitTimeBCell(Yhat_concat, numbUnit, numbTime, numbTrial)
% reshape to a 3D array
rsArray = reshape(Yhat_concat, [numbTime, numbUnit, numbTrial]); 

% transpose each matrix
ts_rsArray = permute(rsArray, [2, 1, 3]); 

% convert to cell array 
unitTimeBCell = squeeze(mat2cell(ts_rsArray, numbUnit, numbTime, ones(1, numbTrial))); 


end