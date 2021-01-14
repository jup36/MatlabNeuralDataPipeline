function [trainC,testC,type] = crossValTrainTestSetsEachTrialType(typeC, fold, varargin) 
%This function takes trial type information (e.g. four types (two position-by-two torques)
% and generates N-fold cross-validation train and test sets. An additional
% logic can be input for further selection of trials as a varargin. 
% use e.g.1: [trainC,testC] = crossValTrainTestSetsEachTrialType(posTqC, 5, pStartI) 
% use e.g.2: 

useAllTrainTrs = true; 

if nargin == 2
    addLogic = true(size(typeC,1),1); 
elseif nargin == 3
    addLogic = varargin{1};     
end

assert(length(addLogic)==length(typeC))

type = unique(typeC); 
for ty = 1:length(type)
    typeS{1,ty} = find(cell2mat(cellfun(@(a) strcmpi(a,type{ty}) , typeC, 'un', 0)) & addLogic);  
    typeS{1,ty} = randsample(typeS{1,ty}, length(typeS{1,ty})); % randomize the order
end

trialN = min(cellfun(@length,typeS)); 
testN = floor(trialN*(1/fold)); 
trainN = floor(trialN*((fold-1)/fold)); 

testC  = cell(length(type),fold); % test sets per type and fold
trainC = cell(length(type),fold); % train sets per type and fold

for f = 1:fold
    tempTest = (f-1)*testN+1:(f-1)*testN+testN; % test trials of this fold (common across trial types)     
    for tt = 1:length(type)    
        if useAllTrainTrs
            tempTrain = find(~ismember(1:length(typeS{tt}),tempTest)); % sample train trials of this trial type of this fold
        else
            tempTrain = randsample(find(~ismember(1:length(typeS{tt}),tempTest)),trainN); % sample train trials of this trial type of this fold
        end
        testC{tt,f}  = sort(typeS{tt}(tempTest)); 
        trainC{tt,f} = sort(typeS{tt}(tempTrain)); 
    end
end 
end