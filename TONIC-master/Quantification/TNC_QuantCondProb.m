function [pCondDist] = TNC_QuantCondProb(dataArray,actual,possible,merNum,shuffNum)
% FUNCTION DETAILS: Function assesses the probability of detecting a given word (MERNUM-mer) prior to a given condition (ACTUAL) in the dataset (DATAARRAY). Significance is assessed by comparing this number of observations to SHUFFNUM random shufflings and resamplings of the data.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

%% first the actual data:
indices = find(dataArray==actual);

count = 1;
for k = 1:length(indices)
    if indices(k)>merNum+1
        obsWordsData(count,:) = dataArray(1,indices(k)-1:-1:indices(k)-merNum);
        count = count + 1;
    end    
end

for l = 1:merNum
    distWordsData(l,:) = hist(obsWordsData(:,l)',1:possible);
end

pCondDist.obsWords = obsWordsData;
pCondDist.distWords = distWordsData ./ count;


%% Shuffled data for comparison
% create the large shuffled dataset
for j=1:shuffNum
    shuffData(j,:) = dataArray(randperm(length(dataArray)));
end

% for each dataset (including actual) walk through and find the conditional probability distribution for the given X-mer
for i=1:shuffNum

    indices = find(shuffData(i,:)==actual);

    for k = 1:length(indices)
        if indices(k)>merNum+1
            obsWords(count,:) = shuffData(i,indices(k)-1:-1:indices(k)-merNum);
            count = count + 1;
        end
    end

    
    for l = 1:merNum
        distWords(l,:) = hist(obsWords(:,l)',1:possible);
    end
    
    if i==1
        tmpDistWords = (distWords ./ count);
    else
        tmpDistWords = tmpDistWords + (distWords ./ count); 
    end
    
    clear obsWords distWords
    count = 1;
    
end

pCondDist.distShuff = tmpDistWords ./ shuffNum;

%% Plot the results:
figure(1)
subplot(121)
imagesc(pCondDist.distShuff,[0 1])
title('P({x,x,x,x}|3} [Shuffled]')
subplot(122)
imagesc(pCondDist.distWords,[0 1])
title('P({x,x,x,x}|3} [Data]')
