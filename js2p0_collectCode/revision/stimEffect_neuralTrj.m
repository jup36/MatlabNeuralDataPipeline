% 

filePath = '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles/js2p0_tbytSpkHandJsTrjBin_50ms_WR40_082019.mat'; 
load(fullfile(filePath), 'ss', 'jkvt')

%% get indices
% valid reach and pull trials (trials with valid hand trajectories)
valRchTrI = cellfun(@(a) ~isempty(a), {ss.spkRchIdx}); 
valPullTrI = cellfun(@(a) ~isempty(a), {ss.spkPullIdx}); 

% stim trials
stimTrI = cellfun(@(a) ~isempty(a), {ss.spkTimeBlaserI}); 
valStimTrI = stimTrI & valRchTrI & valPullTrI; 
valStimTrId = find(valStimTrI);  

% noStim trials
noStimTrI = cellfun(@(a) isempty(a), {ss.spkTimeBlaserI});
valNoStimTrI = noStimTrI & valRchTrI & valPullTrI;
valNoStimTrId = find(valNoStimTrI); 

stbLaserIC = {ss(valStimTrI).spkTimeBlaserI}; % spikeTimeBin laser index 
stbRchIC = {ss(valStimTrI).spkRchIdx}; % spikeTimeBin reach index
stbPullIC = {ss(valStimTrI).spkPullIdx}; % spikeTimeBin pull index

% check overlap between stim and reach and pull
laserRchIC = cellfun(@(a, b) a & b, stbLaserIC, stbRchIC, 'UniformOutput', false); 
cellfun(@sum, laserRchIC)

laserPullIC = cellfun(@(a, b) a & b, stbLaserIC, stbPullIC, 'UniformOutput', false); 
cellfun(@sum, laserPullIC)

%% RRR with trials without silencing
% build X (n x p matrix containing the residual activity of the source)
% population: M1) and Y (n x q matrix containing the residual activity of the target population: STR) matrices. 
[X.concat, X.numbUnit, X.numbTime, X.numbTrial] = concatUnitTimeBCell({ss(valNoStimTrId).unitTimeBCtx}); 
[Y.concat, Y.numbUnit, Y.numbTime, Y.numbTrial] = concatUnitTimeBCell({ss(valNoStimTrId).unitTimeBStr}); 

Xc = {ss(valNoStimTrId).unitTimeBCtx}; 
Yc = {ss(valNoStimTrId).unitTimeBCtx}; 
% run cross-validated (trial-shuffled) rrr with 10 dimensions and 10 folds 
rrrCv = reducedRankRegressCrossVal(Xc, Yc, 10, 10, true); 








[X.concat_stim, X.numbUnit_stim, X.numbTime_stim, X.numbTrial_stim] = concatUnitTimeBCell({ss(valStimTrId).unitTimeBCtx}); 
[Y.concat_stim, Y.numbUnit_stim, Y.numbTime_stim, Y.numbTrial_stim] = concatUnitTimeBCell({ss(valStimTrId).unitTimeBStr}); 

assert(size(X.concat, 1)==size(Y.concat, 1)); % they must have the same number of data points

% run RRR
params.dims = [1:20, size(Y.concat, 2)]; % test dimensions (end: full)
[B, B_, V] = ReducedRankRegress(Y.concat, X.concat, params.dims, 'RIDGEINIT', true, 'SCALE', true); 
% Note: B = Bfull*V(:, 1:dim(1))*V(:, 1:dim(1))'; B_ = Bfull*V, where V = pca(Yhat). 

B_rs = reshape(B, size(B, 1), Y.numbUnit, length(params.dims)); % Stack weight matrices

% get Yhat
yhatC = getYhatStackedB(X.concat, B_rs); % yhat for reduced rank 
yhatC_stim = getYhatStackedB(X.concat_stim, B_rs); % yhat for reduced rank in stim trials


% calculate R2
r2C = cellfun(@(a) calculateR2(Y.concat, a), yhatC, 'UniformOutput', false); 
figure; plot(cell2mat(squeeze(r2C))) % plot r2 versus # of dimensions

r2C_stim = cellfun(@(a) calculateR2(Y.concat_stim, a), yhatC_stim, 'UniformOutput', false); 
figure; plot(cell2mat(squeeze(r2C_stim))) % plot r2 versus # of dimensions

% reshape yhatC back to neuron by timebin dims
yhatC_rs = cellfun(@(a) reshapeYhatToUnitTimeBCell(a, Y.numbUnit, Y.numbTime, Y.numbTrial), yhatC, 'UniformOutput', false); 

% reshape Y and convert to a cell
YC = reshapeYhatToUnitTimeBCell(Y.concat, Y.numbUnit, Y.numbTime, Y.numbTrial); 

% calculate R2 trial by trial
tbytR2 = cell2mat(cellfun(@(a, b) calculateR2(a, b),  YC, yhatC_rs{10}, 'UniformOutput', false)); 



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [concatMat, numbUnit, numbTime, numbTrial] = concatUnitTimeBCell(unitTimeBCell) 
numbTrial = length(unitTimeBCell); 

sizeCell = cellfun(@size, unitTimeBCell, 'UniformOutput', false); 
sizeMat = cell2mat(sizeCell(:)); 
assert(length(unique(sizeMat(:, 1)))==1); % ensure that the number of units match
assert(length(unique(sizeMat(:, 2)))==1); % ensure that the number of time points match

numbUnit = unique(sizeMat(:, 1)); 
numbTime = unique(sizeMat(:, 2)); 

concatMat = full(cell2mat(cellfun(@(a) a', unitTimeBCell, 'UniformOutput', false)'));
end

function unitTimeBCell = reshapeYhatToUnitTimeBCell(Yhat_concat, numbUnit, numbTime, numbTrial)
% reshape to a 3D array
rsArray = reshape(Yhat_concat, [numbTime, numbUnit, numbTrial]); 

% transpose each matrix
ts_rsArray = permute(rsArray, [2, 1, 3]); 

% convert to cell array 
unitTimeBCell = mat2cell(ts_rsArray, numbUnit, numbTime, ones(1, numbTrial)); 

end

function R2 = calculateR2(Y, Y_hat)
    % Y is the matrix of actual values
    % Y_hat is the matrix of predicted values

    % Ensure Y and Y_hat are the same size
    if size(Y) ~= size(Y_hat)
        error('Y and Y_hat must be the same size');
    end

    % Calculate the total sum of squares (SST)
    SST = sum((Y - mean(Y, 'all')).^2, 'all');

    % Calculate the residual sum of squares (SSR)
    SSR = sum((Y - Y_hat).^2, 'all');

    % Calculate R^2
    R2 = 1 - (SSR / SST);
end

function YhatC = getYhatStackedB(X, stackedB)
%this function computes Yhat by multiplying X with weight matrices stacked
% over along the 3rd dimension. The stacked B is assumed to have
% dimensions: # of source units by # of target units by # of stacked weight
% matrices

% YhatM
YhatC = cell(1, 1, size(stackedB, 3)); 

% get intercept
if size(stackedB, 1)-size(X, 2)==0
    intercept = []; 
elseif size(stackedB, 1)-size(X, 2)==1
    intercept = ones(size(X, 1), 1); 
else
    error("Input matrix dimensions do not make sense!")
end
    

for jj = 1:size(stackedB, 3)
    yhat = [intercept X]*stackedB(:, :, jj); 
    YhatC{1, 1, jj} = yhat; 
end

end 



