function [ nTrjCorrCell11, nTrjCorrCell22, nTrjCorrCell12 ] = corrBtwNeuralTrjTwoPopsTimeLag( filePath, nTrjFile1, nTrjFile2, varName )
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

%% get the corr dim-by-dim bin-by-bin bewteen PCs within each population 
% neural population1
nTrjCorrCell11 = cell(5,5,2); % cell array to store the timelagged correlation matrices (rho and p values)

for d = 1:size(nTrjMat1,1) % increment dimension (e.g. 5-d)
    tempTrj1 = squeeze(nTrjMat1(d,:,:))'; % trial-by-timeBin mat for each dim
    for dd = 1:size(nTrjMat1,1) % increment dimension  
         tempTrj2 = squeeze(nTrjMat1(dd,:,:))'; % trial-by-timeBin mat for each dim                        
         [nTrjCorrCell11{d,dd,1},nTrjCorrCell11{d,dd,2}] = corr(tempTrj1,tempTrj2); 
         %tempCorrMat2 = corr(reshape(repmat(tempTrj1,size(tempTrj2,2),1),size(tempTrj1,1),[]), repmat(tempTrj2,1,size(tempTrj1,2)));      
    end
end
clearvars temp*

% neural population2
nTrjCorrCell22 = cell(5,5,2); % cell array to store the timelagged correlation matrices (rho and p values)

for d = 1:size(nTrjMat2,1) % increment dimension (e.g. 5-d)
    tempTrj1 = squeeze(nTrjMat2(d,:,:))'; % trial-by-timeBin mat for each dim
    for dd = 1:size(nTrjMat2,1) % increment dimension  
         tempTrj2 = squeeze(nTrjMat2(dd,:,:))'; % trial-by-timeBin mat for each dim                        
         [nTrjCorrCell22{d,dd,1},nTrjCorrCell22{d,dd,2}] = corr(tempTrj1,tempTrj2); 
         %tempCorrMat2 = corr(reshape(repmat(tempTrj1,size(tempTrj2,2),1),size(tempTrj1,1),[]), repmat(tempTrj2,1,size(tempTrj1,2)));      
    end
end
clearvars temp*

%% get the corr dim-by-dim bin-by-bin between PCs between the two populations
nTrjCorrCell12 = cell(5,5,2); % cell array to store the timelagged correlation matrices (rho and p values)

for d = 1:size(nTrjMat1,1) % increment dimension (e.g. 5-d)
    tempTrj1 = squeeze(nTrjMat1(d,:,:))'; % trial-by-timeBin mat for each dim
    for dd = 1:size(nTrjMat2,1) % increment dimension  
         tempTrj2 = squeeze(nTrjMat2(dd,:,:))'; % trial-by-timeBin mat for each dim                        
         [nTrjCorrCell12{d,dd,1},nTrjCorrCell12{d,dd,2}] = corr(tempTrj1,tempTrj2); 
         %tempCorrMat2 = corr(reshape(repmat(tempTrj1,size(tempTrj2,2),1),size(tempTrj1,1),[]), repmat(tempTrj2,1,size(tempTrj1,2)));      
    end
end

% figure; 
% hold on; plot(1:100, 1:100, ':r', 'LineWidth', 2)
% imagesc(nTrjCorrCell12{1,2,1}); pbaspect([1 1 1 ])
% plot(1:100, 1:100, ':r', 'LineWidth', 2)
% axis tight; 
% set(gca, 'TickDir', 'out')
% hold off; 
% 
% figure; 
% hold on; plot(1:100, 1:100, ':r', 'LineWidth', 2)
% imagesc(nTrjCorrCell12{3,2,1}); pbaspect([1 1 1 ])
% plot(1:100, 1:100, ':r', 'LineWidth', 2)
% axis tight; 
% set(gca, 'TickDir', 'out')
% hold off; 
% 
% figure; 
% hold on; plot(1:100, 1:100, ':r', 'LineWidth', 2)
% imagesc(nTrjCorrCell22{1,2,1}); pbaspect([1 1 1 ])
% plot(1:100, 1:100, ':r', 'LineWidth', 2)
% axis tight; 
% set(gca, 'TickDir', 'out')
% hold off; 

end

