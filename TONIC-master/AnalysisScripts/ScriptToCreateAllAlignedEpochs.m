
%% ANALYSIS PARAMETERS
% Define the smoothing to apply to the timeseries data
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 3;
% currParams.smthParams.decay    = 10;
% currParams.smthParams.decay    = 20;
% currParams.smthParams.decay    = 50;
currParams.filter.causal = 0;

if currParams.filter.causal
    [currParams.filter.kernel]  = TNC_CreateCausalKernel(currParams.smthParams.rise,currParams.smthParams.decay,1);
else
    [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
end

currParams.winParams.prior     = 2e3;
currParams.winParams.after     = 3e3; 

currParams.popVec.winSize       = 1;

currParams.stableTrials         = 10; %derived from summary behavior data

PopData.currParams = currParams

%% ALIGNED RASTERS

disp(['___________________________________________________'])
disp(['STARTED aligning all raster plots and updating the PopData structure...']);

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
    
    numUnits = size(PopData.session(i).unit,2);
    countUnits = countUnits+numUnits;

    if numUnits>=2 && strcmp(PopData.session(i).sessClass,'learning')
        numPairwise = numPairwise + sum(1:1:numUnits-1); % or +1 to count 
    end

    for j = 1:numUnits

        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,round(PopData.session(i).unit(j).ts)+1) = 1;

        [respCS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,currParams.stableTrials,PopData.session(i).events.CS.ts,[currParams.winParams.prior,currParams.winParams.after],1,1);
        PopData.session(i).unit(j).respCS.raster            = respCS.raster;
        PopData.session(i).unit(j).respCS.boxcar            = respCS.image.boxcar;

        tmpSmooth = conv(delta,currParams.filter.kernel,'same');
        [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,currParams.stableTrials,PopData.session(i).events.CS.ts,[currParams.winParams.prior,currParams.winParams.after],0,1);
        PopData.session(i).unit(j).respCSS.image.psthAVG    = respCSS.image.psthAVG;
        PopData.session(i).unit(j).respCSS.image.psthZ      = respCSS.image.psthZ;
        PopData.session(i).unit(j).respCSS.image.psthZe     = respCSS.image.psthZe;

        k=k+1;

    end
    
    disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' of ' num2str(NumSessions) ' | ' PopData.session(i).sessId '_' PopData.session(i).sessClass])
    
end
k = k-1;
k = k./2;
disp(['___________________________________________________'])
disp(['COMPLETED ... Total units: ' num2str(k) ' | ' num2str(countUnits) ' | Pairwise: ' num2str(numPairwise)]);
        
%% TIMESTAMP ONLY DATA FOR TAHL AND GERT
countUnits  = 0;
numUnits=0;
countSessions = 0;
NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions

    % check if learning
    if strcmp(PopData.session(i).sessClass,'extinction')

        numUnits = size(PopData.session(i).unit,2);
        countUnits = countUnits+numUnits;
        countSessions = countSessions+1;

        for j = 1:numUnits
            
            TahlGert.unit(k).ts = PopData.session(i).unit(j).ts;
            k=k+1;

        end
        
    end

    disp(['Completed unit: ' num2str(k) ' of ' num2str(numUnits) ' ... sessions: ' num2str(i) ' of ' num2str(NumSessions)])

end

k = k-1;
disp(['___________________________________________________'])
disp(['COMPLETED ... Total units: ' num2str(k) ' | ' num2str(countUnits) '...in ' num2str(countSessions) ' sessions.']);

%% EXTRACT EXAMPLE RASTER PLOTS AND WRITE AS HDF5 DATA
for i = 1:663
    
%     currData = allCSaligned.sorted(:,i);
    currData    = finalAlign.CS.sorted(:,i);
    currData2   = finalAlign.CS.sorted(:,i);

    if i==1
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../testExportAN1.h5',name,currData);

        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../testExportAN2.h5',name,currData2);
    else
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../testExportAN1.h5',name,currData,'WriteMode','append');
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../testExportAN2.h5',name,currData2,'WriteMode','append');
    end 
        
end

    % also write out the cluster indices:
    name = sprintf('/Indices%g',i-1);
    hdf5write('../../testExportAN1.h5',name,allCSaligned.finalInds,'WriteMode','append');
    name = sprintf('/Ids%g',i-1);
    hdf5write('../../testExportAN1.h5',name,allCSaligned.finalIds,'WriteMode','append');
    
%% compile behavior data for figure 1

%% find example cells for clusters for figure 1
disp('Trying to find best example cell indices...');
% 
% allCSaligned.clusters.numClusts = max(allCSaligned.finalIds);
% 
% for m=1:allCSaligned.clusters.numClusts
% 
%     currClust = find(allCSaligned.finalIds==m);
%     disp(['Members of cluster ' num2str(m) ': ' num2str(length(currClust))]);
% 
%     currClustInds = allCSaligned.finalInds(currClust);
%     
%     tmpMean = allCSaligned.psthZ(:,currClustInds);
%     clustAvg(:,m) = mean(tmpMean,2);
% 
%     for n=1:length(currClustInds)
% %         score(n) = dot(allCSaligned.psthZ(:,currClust(n)),clustAvg(:,m));
%         score(n) = corr(allCSaligned.psthZ(2000:6000,currClustInds(n)),clustAvg(2000:6000,m));
%     end
%     
%     bestFitInd = find(score==max(score),1);
%     
%     % store the indices to the best examples
%     allCSaligned.exemplars(m) = currClustInds(bestFitInd);
%     
%     disp(['For cluster ' num2str(m) ' the best match is in clusterPosition: ' num2str(currClust(bestFitInd)) ' totalIndex: ' num2str(currClustInds(bestFitInd))]);
%     clear score;
% end

% new version uses cascading sorting, two stage sorted indices are:

numClasses = size(trackSort.resort,2);

for m=1:numClasses
    
    clear clustAvg score
    
    numInds = length(trackSort.resort(m).Inds);
    subInds = round(numInds./5);
    lowClustInds    = trackSort.resort(m).Inds(1:subInds);
    midClustInds    = trackSort.resort(m).Inds(3.*subInds:4.*subInds);
    hiClustInds     = trackSort.resort(m).Inds(numInds-subInds:numInds);
    
    tmpMean = allCSaligned.psthZ(:,lowClustInds);
    clustAvg(:,m) = mean(tmpMean,2);
    figure(9); plot(clustAvg(:,m)); drawnow; title(num2str(m));
    for n=1:length(lowClustInds)
        score(n,1) = dot(allCSaligned.psthZ(:,lowClustInds(n)),clustAvg(:,m));
    end
    bestFitInd = find(score(:,1)==max(score(:,1)),1);
    allCSaligned.exemplars(m,1) = lowClustInds(bestFitInd);
    clear clustAvg
    disp(['For cluster ' num2str(m) ' the best match is in clusterPosition: ' num2str(lowClustInds(bestFitInd))]);
    
    tmpMean = allCSaligned.psthZ(:,midClustInds);
    clustAvg(:,m) = mean(tmpMean,2);
    figure(9); plot(clustAvg(:,m)); drawnow;
    for n=1:length(midClustInds)
        score(n,2) = dot(allCSaligned.psthZ(:,midClustInds(n)),clustAvg(:,m));
    end
    bestFitInd = find(score(:,2)==max(score(:,2)),1);
    allCSaligned.exemplars(m,2) = midClustInds(bestFitInd);
    clear clustAvg
    disp(['For cluster ' num2str(m) ' the best match is in clusterPosition: ' num2str(midClustInds(bestFitInd))]);

    tmpMean = allCSaligned.psthZ(:,hiClustInds);
    clustAvg(:,m) = mean(tmpMean,2);
    figure(9); plot(clustAvg(:,m)); drawnow;
    for n=1:length(hiClustInds)
        score(n,3) = dot(allCSaligned.psthZ(:,hiClustInds(n)),clustAvg(:,m));
    end
    bestFitInd = find(score(:,3)==max(score(:,3)),1);
    allCSaligned.exemplars(m,3) = hiClustInds(bestFitInd);
    clear clustAvg    
    disp(['For cluster ' num2str(m) ' the best match is in clusterPosition: ' num2str(hiClustInds(bestFitInd))]);
    

end

%% prep example raster plots and psths

% cs aligned timestamps for spikes are stored in:
% PopData.session(1).unit(1).respCS.raster.trial(1).ts

% iterate through each exemplar
numExamples = size(allCSaligned.exemplars,1);

% m = 1; % for the low class
% m = 2; % for the mid class
m = 3; % for the hi class

for j=1:numExamples

    currIndex = allCSaligned.exemplars(j,m);

    disp([num2str(j) '...' num2str(m)]);
    
    % grab the session and unit information
    currSess = allCSaligned.sess(currIndex);
    currUnit = allCSaligned.unit(currIndex);

    % index into the PopData and concatenate the timestamps
    numTrials = size(PopData.session(currSess).unit(currUnit).respCS.raster.trial,2);
    
    for k=1:numTrials
    
        % create a equal sized trial id vector
        thisTrialTS     = PopData.session(currSess).unit(currUnit).respCS.raster.trial(k).ts;
        numStamps       = size(thisTrialTS,1);
        thisTrialIDS    = ones(numStamps,1).*k;
        lastTrialArray  = [thisTrialIDS,thisTrialTS];
        
        if k==1
            finalArrayTS = thisTrialTS;            
            finalArrayID = thisTrialIDS;            
        else
            finalArrayTS = [finalArrayTS;thisTrialTS];
            finalArrayID = [finalArrayID;thisTrialIDS];
        end
        
    end
    
    figure(1); subplot(4,1,1:3);
    plot(finalArrayTS,finalArrayID,'k.');
    subplot(4,1,4);
    plot(-1000:3000, allCSaligned.psthZ(:,currIndex),'k');
    drawnow; pause(2);
    
    % save to an hdf5 structure for reading/plotting in igor
    nameA = sprintf('/spikeTimesX%g',j-1);
    nameB = sprintf('/spikeTimesY%g',j-1);
    nameC = sprintf('/psth%g',j-1);

    if j==1
        hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameA,finalArrayTS);
        hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameB,finalArrayID,'WriteMode','append');
        hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameC,allCSaligned.psthZ(:,currIndex),'WriteMode','append');
    else
        hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameA,finalArrayTS,'WriteMode','append');
        hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameB,finalArrayID,'WriteMode','append');
        hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameC,allCSaligned.psthZ(:,currIndex),'WriteMode','append');
    end
    
end

%% extract the population mean licking behavior and an example behavior
disp(['___________________________________________________'])
disp(['STARTED aligning behavior and updating the PopData structure...']);

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
validBehavSess = [];
k=1;

for i=1:NumSessions
    
        disp(['Session ' num2str(i)]);

        numStamps = numel(PopData.session(i).events.EL.ts);
        
        if numStamps>1
            
            PopData.session(i).behavior.found  = 1;

            delta = zeros(1,ceil(PopData.session(i).events.EL.ts(numStamps)));
            delta(1,round(PopData.session(i).events.EL.ts)+1) = 1;

            % use a longer window for licking behvior to tolerate long delays to retrieval:
            [respCS] = TNC_AlignRasters(delta,PopData.session(i).events.EL.ts,currParams.stableTrials,PopData.session(i).events.CS.ts,[2000,8000],1,0);
            PopData.session(i).behavior.raster  = respCS.raster;
            PopData.session(i).behavior.psthAVG = respCS.image.psthAVG;
            PopData.session(i).behavior.psthSEM = respCS.image.psthSEM;
            
            validBehavSess = [validBehavSess,i];
            
        else
            
            PopData.session(i).behavior.raster  = 0;
            PopData.session(i).behavior.psthAVG = 0;
            PopData.session(i).behavior.psthSEM = 0;
            
        end
    
end

PopData.validBehavSess = validBehavSess;

disp(['___________________________________________________'])

%% compile and smooth the mean licking behavior and an example behavior
disp(['___________________________________________________'])
disp(['STARTED aligning behavior and updating the PopData structure...']);
clear lickPSTHavg lickPSTHsm
NumSessions = size(PopData.session,2);
numTrials = 0;
k=1;

for i=1:NumSessions
    
    if strcmp(PopData.session(i).sessClass,'learning')

        if PopData.session(i).behavior.found
            disp(['Session ' num2str(i)]);
            numTrials           = numTrials + size(PopData.session(i).behavior.raster.trial,2);
            tAVG                = PopData.session(i).behavior.psthAVG';
            lickPSTHavg(k,:)    = tAVG ./ trapz(tAVG);
            sAVG                = conv(PopData.session(i).behavior.psthAVG',currParams.filter.kernel,'same');
            lickPSTHsm(k,:)     = sAVG ./ trapz(sAVG);
            k = k+1;
        end
    end
    
end
k=k-1;

allLicks.raw    = lickPSTHavg.*1000;
allLicks.smooth = lickPSTHsm.*1000;
allLicks.ravg    = mean(lickPSTHavg.*1000,1);
allLicks.savg    = mean(lickPSTHsm.*1000,1);
allLicks.serr    = std(lickPSTHsm.*1000,0,1)./sqrt(k-1);

disp(['Aligned ' num2str(numTrials) ' trials from ' num2str(k) ' sessions of licking behavior.']);

figure(2);
subplot(121); imagesc(allLicks.smooth,[0 1]);
subplot(122); plot(-2000:3000,allLicks.savg,'k');

%% Create a US aligned sorted psth plot
clear finalAlign;
countUnits  = 0;
NumSessions = size(PopData.session,2);
k=1;

countUnits = 0;
NumSessions = size(PopData.session,2)
k=1;

for i=1:NumSessions
    % check if learning
    if strcmp(PopData.session(i).sessClass,'learning')
        offset = PopData.session(i).USfiat;

        USpos = (offset.*1000) + 2000;
        disp(['Session ' num2str(i) '... US offset ' num2str(offset) ]);    
        numUnits = size(PopData.session(i).unit,2);
        countUnits = countUnits+numUnits;
        for j = 1:numUnits
            if isnan(PopData.session(i).unit(j).isi.hist.logCount(1,1))
                disp('Invalid ISI distribution found');
            else
                finalAlign.CS.psthZ(:,k) = PopData.session(i).unit(j).respCSS.image.psthZ(1000:4000)';
                finalAlign.US.psthZ(:,k) = PopData.session(i).unit(j).respCSS.image.psthZ(USpos-1500:USpos+3500)';
                finalAlign.sess(1,k)  = i;
                finalAlign.unit(1,k)  = j;
                k=k+1;
            end
        end
    end
    
end

k = k-1;

disp(['Total units: ' num2str(k) ' | ' num2str(countUnits)])

