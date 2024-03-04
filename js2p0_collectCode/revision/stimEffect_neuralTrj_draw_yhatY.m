%The purpose of this code is to generate a plot that compares the predicted
% and actual activity of the target population (e.g., striatum).  

filePath = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/rrrRezCV_stimPstim_multiDim_WR40_081919_Folds10.mat'; 
figPath = fullfile(fileparts(filePath), 'Figure'); 
load(fullfile(filePath), 'rrrCv', 'rrrRezR2')

rrDim = 10; % the number of dimensions to work with

yfC = []; 
yhatfC = []; 

% organize yfC and yhatfC, each cell must have neuron-by-time bin matrix
for ff = 1:length(rrrCv)
    % get y and yhat of this folder
    if size(rrrCv(ff).Yhat_test, 3) >= 10
        yf = rrrCv(ff).Ycc_test; 
        yhatf = rrrCv(ff).Yhat_test(:, :, rrDim); 
    else 
        yf = rrrCv(ff).Ycc_test; 
        yhatf = rrrCv(ff).Yhat_test; 
    end
    assert(isequal(size(yhatf), size(yf))); 
    
    % reshape yf and yhatf should be in data points by neuron dimensions
    yfC = [yfC; reshapeYhatToUnitTimeBCell(yf, rrrCv(ff).numbUnitY, rrrCv(ff).numbTime, rrrCv(ff).numbTrial_test)]; 
    yhatfC = [yhatfC; reshapeYhatToUnitTimeBCell(yhatf, rrrCv(ff).numbUnitY, rrrCv(ff).numbTime, rrrCv(ff).numbTrial_test)]; 
end


sqe = cellfun(@(a, b) sum(sum(sqrt((a-b).^2))), yfC, yhatfC); 
[sqe, sqeI] = sortrows(sqe, 1); 

% smooth Y
yfC_sm = cellfun(@(a) smooth2a(a, 0, 1), yfC, 'UniformOutput', false); 

% get prediction from fitted yhat in discrete numbers 
yhatfC_rs = cellfun(@(a) poissrndSpikeGenerator(a), yhatfC, 'UniformOutput', false); 
yhatfC_rs_sm = cellfun(@(a) smooth2a(a, 0, 1), yhatfC_rs, 'UniformOutput', false); 

% image yhat and y
h_y = plotSpikeCountMatrices(yfC_sm(sqeI(1:20)), [0 3], 0.0015); 
print(h_y, fullfile(figPath, 'y_example_trials'), '-vector', '-dpdf', '-bestfit')

h_yhat = plotSpikeCountMatrices(yhatfC_rs_sm(sqeI(1:20)), [0 3], 0.0015); 
print(h_yhat, fullfile(figPath, 'yhat_example_trials'), '-vector', '-dpdf', '-bestfit')


function yhatUnitByTimeRepPoiss = poissrndSpikeGenerator(yhatUnitByTime)
    % replace negative values in yhatUnitByTime to zeros.  
    yhatUnitByTime(yhatUnitByTime<0)=0; 
    % poiss random sampling
    yhatUnitByTimeRepPoiss = poissrnd(yhatUnitByTime); 
end


function unitTimeBCell = reshapeYhatToUnitTimeBCell(Yhat_concat, numbUnit, numbTime, numbTrial)
% Orient the input matrix as numbUnit x numbDataPoints (time x trial)
% because that orientation is required for proper reshaping.
if size(Yhat_concat, 1) == numbUnit
    Yhat = Yhat_concat;
elseif size(Yhat_concat, 2) == numbUnit
    Yhat = Yhat_concat';
end

% reshape to a 3D array
rsArray = reshape(Yhat, [numbUnit, numbTime, numbTrial]);

% convert to cell array
unitTimeBCell = squeeze(mat2cell(rsArray, numbUnit, numbTime, ones(1, numbTrial)));


end
