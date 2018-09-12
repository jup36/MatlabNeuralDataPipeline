function [ nTrjCorrCell11, nTrjCorrCell22, nTrjCorrCell12 ] = corrBtwNeuralTrjTwoPopsTimeLag( filePath, nTrjFile1, nTrjFile2, varName, neuralRegion1, neuralRegion2, eventName, cAxisValue )
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
nTrjCorrCell11Name = strcat('nTrjCorr',neuralRegion1,neuralRegion1,eventName,'Rho');
formatSpecTitle = 'PC%d%d'; 
for d = 1:size(nTrjMat1,1) % increment dimension (e.g. 5-d)
    tempTrj1 = squeeze(nTrjMat1(d,:,:))'; % trial-by-timeBin mat for each dim
    for dd = 1:size(nTrjMat1,1) % increment dimension  
         tempTrj2 = squeeze(nTrjMat1(dd,:,:))'; % trial-by-timeBin mat for each dim                        
         [nTrjCorrCell11{d,dd,1},nTrjCorrCell11{d,dd,2}] = corr(tempTrj1,tempTrj2);     
         pcId = sprintf(formatSpecTitle,d,dd); 
         titleImg = [nTrjCorrCell11Name, pcId]; 
         surfMatrix( filePath, nTrjCorrCell11{d,dd,1}, nTrjCorrCell11{d,dd,2}, cAxisValue, 'parula', titleImg )
    end
end
clearvars temp*

% neural population2
nTrjCorrCell22 = cell(5,5,2); % cell array to store the timelagged correlation matrices (rho and p values)
nTrjCorrCell22Name = strcat('nTrjCorr',neuralRegion2,neuralRegion2,eventName,'Rho');

for d = 1:size(nTrjMat2,1) % increment dimension (e.g. 5-d)
    tempTrj1 = squeeze(nTrjMat2(d,:,:))'; % trial-by-timeBin mat for each dim
    for dd = 1:size(nTrjMat2,1) % increment dimension  
         tempTrj2 = squeeze(nTrjMat2(dd,:,:))'; % trial-by-timeBin mat for each dim                        
         [nTrjCorrCell22{d,dd,1},nTrjCorrCell22{d,dd,2}] = corr(tempTrj1,tempTrj2); 
         pcId = sprintf(formatSpecTitle,d,dd);          
         titleImg = [nTrjCorrCell22Name, pcId];    
         surfMatrix( filePath, nTrjCorrCell22{d,dd,1}, nTrjCorrCell22{d,dd,2}, cAxisValue, 'parula', titleImg )   
    end
end
clearvars temp*

%% get the corr dim-by-dim bin-by-bin between PCs between the two populations
nTrjCorrCell12 = cell(5,5,2); % cell array to store the timelagged correlation matrices (rho and p values stacked in the 3rd dimension)
nTrjCorrCell12Name = strcat('nTrjCorr',neuralRegion1,neuralRegion2,eventName,'Rho');
nTrjCorrCell1minus2 = cell(5,5,1); % subtract corr coeffs of the bottom triangle from the upper triangle to infer directionality
nTrjCorrCell2minus1 = cell(5,5,1); % subtract corr coeffs of the upper triangle from the bottom triangle to infer directionality 

for d = 1:size(nTrjMat1,1) % increment dimension (e.g. 5-d)
    tempTrj1 = squeeze(nTrjMat1(d,:,:))'; % trial-by-timeBin mat for each dim
    for dd = 1:size(nTrjMat2,1) % increment dimension  
         tempTrj2 = squeeze(nTrjMat2(dd,:,:))'; % trial-by-timeBin mat for each dim                        
         [nTrjCorrCell12{d,dd,1},nTrjCorrCell12{d,dd,2}] = corr(tempTrj1,tempTrj2); 
         pcId = sprintf(formatSpecTitle,d,dd);          
         titleImg = [nTrjCorrCell12Name, pcId];    
         surfMatrix( filePath, nTrjCorrCell12{d,dd,1}, nTrjCorrCell12{d,dd,2}, cAxisValue, 'parula', titleImg ) 
         
         if d <=3 && dd <=3
            titleImg12 = [nTrjCorrCell12Name, pcId, 'Pop1MinusPop2']; 
            titleImg21 = [nTrjCorrCell12Name, pcId, 'Pop2MinusPop1']; 
            
            triu(nTrjCorrCell12{d,dd,1})-triu(nTrjCorrCell12{d,dd,1}'); % directionality from pop1 to pop2
            triu(nTrjCorrCell12{d,dd,1}')-triu(nTrjCorrCell12{d,dd,1}); % directionality from pop2 to pop1
         end
         
    end
end


end

