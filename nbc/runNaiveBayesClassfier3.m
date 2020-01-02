function [decodeRezAll, posteriorsAll] = runNaiveBayesClassfier3(cvItr, numbTr, predictor1, predictor2, predictor3, unitIdx, numbDropCells)

label = [ones(numbTr,1); ones(numbTr,1)*2; ones(numbTr,1)*3];

decodeRezAll  = [];
%cvItr = 100; % cross-validation
count = 0;
numbUnit = size(predictor1,1);
trainLambda = zeros(numbUnit,size(predictor1,2),length(unique(label)));

if sum(unitIdx)>0 % when dropping cells (leave-cells-out)
    assert(sum(unitIdx)>numbDropCells)
    dropCells = find(unitIdx);
    trainLambda = zeros(numbUnit-numbDropCells,size(predictor1,2),length(unique(label))); 
end

for cv = 1:cvItr % run 2-fold cross-validation
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    train = zeros(numbTr,1);
    train(randsample(numbTr, round(numbTr*.9)),1) = 1;
    train = train==1;
    test = ~train;
    unitI = ones(size(unitIdx,1),1);     
    if sum(unitIdx)>0
        unitI(dropCells(randsample(length(dropCells), numbDropCells,false)),1)=0;
    end
    unitI = logical(unitI); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if sum(unitIdx)>0
        predictor11 = predictor1(unitI,:,:);
        predictor22 = predictor2(unitI,:,:);
        predictor33 = predictor3(unitI,:,:);
    else 
        predictor11 = predictor1;
        predictor22 = predictor2;
        predictor33 = predictor3;
        
    end
    
    %% Train
    trainLambda(:,:,1) = squeeze(nanmean(predictor11(:,:,train),3));
    trainLambda(:,:,2) = squeeze(nanmean(predictor22(:,:,train),3));
    trainLambda(:,:,3) = squeeze(nanmean(predictor33(:,:,train),3));
    
    %% Test
    testSet = zeros(size(predictor11,1),size(predictor11,2),3*sum(test));
    testSet(:,:,1:sum(test)) = predictor11(:,:,test);
    testSet(:,:,sum(test)+1:2*sum(test)) = predictor22(:,:,test);
    testSet(:,:,2*sum(test)+1:3*sum(test)) = predictor33(:,:,test);
    
    testC = [ones(sum(test),1);ones(sum(test),1)*2; ones(sum(test),1)*3];
    
    rsTestSet = reshape(testSet, size(testSet,1)*size(testSet,2)*size(testSet,3),1,1);
    
    % get posteriors
    for c = 1:length(unique(label))
        rsTrainLambda = reshape(repmat(trainLambda(:,:,c),[1,1,size(testSet,3)]), size(testSet,1)*size(testSet,2)*size(testSet,3),1,1);
        tmpPosterior = poisspdf(rsTestSet, rsTrainLambda);
        posteriorC{c} = reshape(tmpPosterior,size(testSet,1),size(testSet,2),size(testSet,3));
    end
    
    decodeRez = zeros(size(testC,1), size(predictor11,2));
    for t = 1:length(testC)
        count = count + 1;
        pUnitTimeLable = cell2mat(reshape(cellfun(@(c) squeeze(c(:,:,t)), posteriorC, 'Un', 0)',1,1,[])) + .0001; % get log probability
        [~,decodeLabel] = max(sum(log(pUnitTimeLable),1),[],3); % get sum of the log prob. across units and argmax across posteriors
        decodeRez(t,:) = decodeLabel == repmat(testC(t),1,size(decodeLabel,2)); % compare decoded vs actual labels
        posteriorsAll{count} = posteriorC{testC(t)}(:,:,t); % posterior probability for the actual class
    end
    
    decodeRezAll{cv} = decodeRez;
    fprintf('finished crossVal fold %d\n', cv);
end