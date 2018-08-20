function [ nTrjCorrCellR, nTrjCorrCellP ] = corrBtwNeuralTrjTwoPopsTimeLag( filePath, nTrjFile1, nTrjFile2, varName )
%This function computes trial-to-trial correlations of PC scores of all pairwise combinations of PC dimensions across all time bins 
% comprising the activity covariance of two neural populations. 
% The matlab function 'corr' operates on all columnwise combinations by
% default. 

cd(filePath)

% get the nTrj 
nTrj1 = load(nTrjFile1,varName); % load the nTrj structure
nTrj1 = nTrj1.(varName);         % rename the structure as nTrj1
nTrjMat1 = nTrj1.trjMat; % nTrj mat

nTrj2 = load(nTrjFile2,varName); % load the nTrj structure
nTrj2 = nTrj2.(varName);         % rename the structure as nTrj2
nTrjMat2 = nTrj2.trjMat; % nTrj mat

% Check if the dimensions of the nTrj matrices match! 
if ~isequal(size(nTrjMat1),size(nTrjMat2))
    error('Dimensions of the two nTrj matrices do not match!!!')
end

%% get the corr dim-by-dim bin-by-bin
nTrjCorrCellR = cell(5,5); % cell array to store the timelagged correlation matrices (rho values)
nTrjCorrCellP = cell(5,5); % cell array to store the timelagged correlation matrices (p values)

for d = 1:size(nTrjMat1,1) % increment dimension (e.g. 5-d)
    tempTrj1 = squeeze(nTrjMat1(d,:,:))'; % trial-by-timeBin mat for each dim
    for dd = 1:size(nTrjMat2,1) % increment dimension  
         tempTrj2 = squeeze(nTrjMat2(d,:,:))'; % trial-by-timeBin mat for each dim                        
         [nTrjCorrCellR{d,dd},nTrjCorrCellP{d,dd}] = corr(tempTrj1,tempTrj2); 
         %tempCorrMat2 = corr(reshape(repmat(tempTrj1,size(tempTrj2,2),1),size(tempTrj1,1),[]), repmat(tempTrj2,1,size(tempTrj1,2)));      
    end
end

end

