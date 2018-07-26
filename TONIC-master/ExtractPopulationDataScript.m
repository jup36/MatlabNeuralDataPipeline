% Extract data from the population structure 
% 
% analyses to perform:
% 
%     quantifying the population data for significant modulation to:
%         cs learning - cs extinction
%         reward delivery [retrieval]
%         delay period [projection?]
%         
%     classifying cells/responses
%         isi distributions
%         joe style analysis of different psths for both cs and reward delivery alignment
%         waveform classifier projection onto three frequency band sinc functions
% 

%% Quantify [spike count] CS responses
countUnits  = 0;
NumSessions = size(PopData.session,2);
k           = 1;

timeWindow(1,:) = 2000:2060; % early response 1
timeWindow(2,:) = 2000:2060; % early response 2
timeWindow(3,:) = 2060:2150; % dopamine response
timeWindow(4,:) = 2500:4500; % trace period
timeWindow(5,:) = 5500:7500; % consumption

for i=1:NumSessions
    % check if learning
    if strcmp(PopData.session(i).sessClass,'learning')
        numUnits = size(PopData.session(i).unit,2);
        countUnits = countUnits+numUnits;

        for j = 1:numUnits
            % count the number of spikes in a specific window for every trial
            for m = 1:numWindows
                
            % store the histogram of spike counts for learning
            end
            
            k = k+1;
        end        
    end    
end

k = k-1;
disp(['Total learning units: ' num2str(k) ' | ' num2str(countUnits)])

%%%%%%%%% Repeat for extinction
countUnits  = 0;
k = 1;

% From the ordered collection of histograms calculate the ROC

% Run ROC and other significance tests

k = k-1;
disp(['Total extinction units: ' num2str(k) ' | ' num2str(countUnits)])

%% Plot evolution of population vector over time

% calculate the pca of the response space to improve clustering performance
justPhasicL = allCSaligned.sort.susPeak.sortedZ(1950:2550,:);
justPhasicE = allCSEXaligned.sort.susPeak.sortedZ(1950:2550,:);
covMat2  = cov(justPhasic);
[v,d]   = eig(covMat2);
tmp = size(v);

figure(1);
subplot(5,1,1:4);
imagesc(justPhasicL',[-10 10]);


for m = 1:size(justPhasicL,1)
    popVec(m,1) = dot(justPhasicL(m,:),v(:,tmp(1,2)-0));
    popVec(m,2) = dot(justPhasicL(m,:),v(:,tmp(1,2)-1));
    popVec(m,3) = dot(justPhasicL(m,:),v(:,tmp(1,2)-2));
    popVecE(m,1) = dot(justPhasicE(m,:),v(:,tmp(1,2)-0));
    popVecE(m,2) = dot(justPhasicE(m,:),v(:,tmp(1,2)-1));
    popVecE(m,3) = dot(justPhasicE(m,:),v(:,tmp(1,2)-2));
end
m
figure(2);clf;
% subplot(3,3,1:6);
plot3(popVec(:,1),popVec(:,2),popVec(:,3),'k','LineWidth',1); grid on; hold on;
plot3(popVecE(:,1),popVecE(:,2),popVecE(:,3),'r','LineWidth',1);
plot3([popVec(62,1),popVecE(62,1)],[popVec(62,2),popVecE(62,2)],[popVec(62,3),popVecE(62,3)],'bo-','LineWidth',2);
plot3([popVec(95,1),popVecE(95,1)],[popVec(95,2),popVecE(95,2)],[popVec(95,3),popVecE(95,3)],'co-','LineWidth',2);
plot3([popVec(122,1),popVecE(122,1)],[popVec(122,2),popVecE(122,2)],[popVec(122,3),popVecE(122,3)],'go-','LineWidth',2);

popVecDist = sqrt( (popVec(:,1)-popVecE(:,1)).^2 + (popVec(:,2)-popVecE(:,2)).^2 + (popVec(:,3)-popVecE(:,3)).^2);
figure(1);
subplot(5,1,5);
semilogy(popVecDist,'k-','LineWidth',2); axis([0 600 1 1000]);

% subplot(3,3,7:9);
% plot(1:663,v(:,tmp(1,2)-0),'r',1:663,v(:,tmp(1,2)-1),'k',1:663,v(:,tmp(1,2)-2),'b');
% subplot(3,3,8);
% plot(v(:,tmp(1,2)-1));
% subplot(3,3,9);
% plot(v(:,tmp(1,2)-2));
% hold off;
% plot(popVec(:,3),popVec(:,2),'k','LineWidth',2);

%% CREATE ALLCSALIGNED STRUCTURE
clear allCSaligned
countUnits = 0;
NumSessions = size(PopData.session,2)
k=1;

for i=1:NumSessions
    % check if learning
    if strcmp(PopData.session(i).sessClass,'learning')

        if i<10
            disp(['Session 00' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
        elseif i<100
            disp(['Session 0' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
        else
            disp(['Session ' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
        end
            
        numUnits = size(PopData.session(i).unit,2);
        countUnits = countUnits+numUnits;

        for j = 1:numUnits
            if isnan(PopData.session(i).unit(j).isi.hist.logCount(1,1))
                disp('Invalid ISI distribution found');
            else
                allCSaligned.psthS(:,k)     = PopData.session(i).unit(j).respCSS.image.psthAVG'.*1000;
                allCSaligned.psthZ(:,k)     = PopData.session(i).unit(j).respCSS.image.psthZ';
                allCSaligned.psthZe(:,k)    = PopData.session(i).unit(j).respCSS.image.psthZe';
                allCSaligned.boxcar(:,k)    = PopData.session(i).unit(j).respCS.boxcar'.*100;
                allCSaligned.sess(1,k)      = i;
                allCSaligned.unit(1,k)      = j;
                allCSaligned.id(1,k).sessID = PopData.session(i).sessId;
                allCSaligned.id(1,k).name   = PopData.session(i).unit(j).name;
                k=k+1;
            end
        end
    end
    
end

k = k-1;

disp(['Total units: ' num2str(k) ' | ' num2str(countUnits)])

numUnits = countUnits;

figure(6);
subplot(121);
[mapName] = TNC_CreateRBColormap(1024,'rb');
imagesc(allCSaligned.psthZ',[-15 15]);
colormap(mapName);
xlabel('Time (ms)');
ylabel('Cell index');
title('CS Aligned ZDF');

subplot(122);
imagesc(allCSaligned.boxcar',[-100 100]);
colormap(mapName);
xlabel('Time (ms)');
ylabel('Cell index');
title('CS Aligned Boxcar');

%% PSTH classifier CS
% New approach: try a cascade of clustering:
%     Step 1: cluster based upon responses during the CS (Dopamine vs. Non-Dopamine)
%     Step 2: separate out the groups into 3 super respond, sus

% clear pcas;
% 
% numClusts = 3;
% allCSaligned.clusters.numClusts = numClusts;
% pcaClustering = 1;
% window = 2000:2500;
% 
% if pcaClustering
%     % calculate the pca of the response space to improve clustering performance
%     covMat  = cov(allCSaligned.psthZ(window,:)');
%     [v,d]   = eig(covMat);
%     PopData.clustering.varExp = cumsum(diag(d))./sum(diag(d));
%     tmp = size(v);
%     pcas(1,1:size(v,1)) = v(:,tmp(1,2)-0)';
%     pcas(2,1:size(v,1)) = v(:,tmp(1,2)-1)';
%     pcas(3,1:size(v,1)) = v(:,tmp(1,2)-2)';
%     pcas(4,1:size(v,1)) = v(:,tmp(1,2)-3)';
%     pcas(5,1:size(v,1)) = v(:,tmp(1,2)-4)';
%     pcas(6,1:size(v,1)) = v(:,tmp(1,2)-5)';
% 
%     Clustering.pcas = pcas;
%     
%     % project data onto first 3 pcs
%     for m = 1:size(allCSaligned.psthZ,2)
%         allCSaligned.clustering.proj(m,1) = dot(allCSaligned.psthZ(window,m),v(:,tmp(1,2)-0));
%         allCSaligned.clustering.proj(m,2) = dot(allCSaligned.psthZ(window,m),v(:,tmp(1,2)-1));
%         allCSaligned.clustering.proj(m,3) = dot(allCSaligned.psthZ(window,m),v(:,tmp(1,2)-2));
%         allCSaligned.clustering.proj(m,4) = dot(allCSaligned.psthZ(window,m),v(:,tmp(1,2)-3));
%         allCSaligned.clustering.proj(m,5) = dot(allCSaligned.psthZ(window,m),v(:,tmp(1,2)-4));
%         allCSaligned.clustering.proj(m,6) = dot(allCSaligned.psthZ(window,m),v(:,tmp(1,2)-5));
%     end
% 
%     allCSaligned.clusters.ids1 = kmeans(allCSaligned.clustering.proj,numClusts);
%     
% else
%     
%     allCSaligned.clusters.ids1 = kmeans(allCSaligned.psthZ',numClusts);
% 
% end
% 
% [y,allCSaligned.sort.indices1]  = sort(allCSaligned.clusters.ids1,1,'ascend')
% allCSaligned.sorted1            = allCSaligned.psthZ(:,allCSaligned.sort.indices1);
% 
% figure(6);
% [mapName] = TNC_CreateRBColormap(1024,'rb');
% colormap(mapName);
% subplot(1,5,1:2);
% tmp = max(max(allCSaligned.sorted1));
% imagesc(allCSaligned.sorted1',[-10 10]);
% xlabel('Time (ms)');
% ylabel('Cell index');
% title('CS Aligned PSTH');

% subplot(1,5,3);
% plot(allCSaligned.clusters.ids1(allCSaligned.sort.indices1),701:-1:1,'k.'); axis([0 10 1 701]); axis off;
% 
% colormap(mapName);
% subplot(1,5,4:5);
% imagesc(corr(allCSaligned.sorted),[-1,1]);
% ylabel('Cell index');
% xlabel('Cell index');
% title('CS Aligned Signal Correlations');


% for clusters 2 repeat the pca
tmpDaInds       = find(allCSaligned.clusters.ids1==1);
tmpNonDaInds2    = find(allCSaligned.clusters.ids1==2);
lowClusterPSTH  = allCSaligned.psthZ(:,tmpDaInds);
midCluster      = allCSaligned.psthZ(:,tmpNonDaInds2);
size(midCluster)
window          = 4000:6000;
numClusts       = 12;

allCSaligned.sort.indices1 = tmpDaInds;
disp('Made it through 1 of 3 clustering');

if pcaClustering
    % calculate the pca of the response space to improve clustering performance
    covMat  = cov(midCluster(window,:)');
    [v,d]   = eig(covMat);
    PopData.clustering.varExp2 = cumsum(diag(d))./sum(diag(d));
    tmp = size(v);
    pca2(1,1:size(v,1)) = v(:,tmp(1,2)-0)';
    pca2(2,1:size(v,1)) = v(:,tmp(1,2)-1)';
    pca2(3,1:size(v,1)) = v(:,tmp(1,2)-2)';

    Clustering.pca2 = pca2;
    
    % project data onto first 3 pcs
    for m = 1:size(midCluster,2)
        allCSaligned.clustering.proj2(m,1) = dot(midCluster(window,m),v(:,tmp(1,2)-0));
        allCSaligned.clustering.proj2(m,2) = dot(midCluster(window,m),v(:,tmp(1,2)-1));
        allCSaligned.clustering.proj2(m,3) = dot(midCluster(window,m),v(:,tmp(1,2)-2));
    end

    allCSaligned.clusters.ids2 = kmeans(allCSaligned.clustering.proj2,numClusts);
    
else
    
    allCSaligned.clusters.ids2 = kmeans(midCluster',numClusts);

end

% sort away
[y2,allCSaligned.sort.indices2]  = sort(allCSaligned.clusters.ids2,1,'ascend')
midClusterPSTH                  = midCluster(:,allCSaligned.sort.indices2);

% for clusters 2 repeat the pca
tmpNonDaInds3    = find(allCSaligned.clusters.ids1==3);
lastCluster     = allCSaligned.psthZ(:,tmpNonDaInds3);
window          = 4000:6000;
numClusts       = 5;
size(lastCluster)

disp('Made it through 2 of 3 clustering');

if pcaClustering
    % calculate the pca of the response space to improve clustering performance
    covMat  = cov(lastCluster(window,:)');
    [v,d]   = eig(covMat);
    PopData.clustering.varExp3 = cumsum(diag(d))./sum(diag(d));
    tmp = size(v);
    pca3(1,1:size(v,1)) = v(:,tmp(1,2)-0)';
    pca3(2,1:size(v,1)) = v(:,tmp(1,2)-1)';
    pca3(3,1:size(v,1)) = v(:,tmp(1,2)-2)';

    Clustering.pca2 = pca2;
    
    % project data onto first 3 pcs
    for m = 1:size(lastCluster,2)
        allCSaligned.clustering.proj3(m,1) = dot(lastCluster(window,m),v(:,tmp(1,2)-0));
        allCSaligned.clustering.proj3(m,2) = dot(lastCluster(window,m),v(:,tmp(1,2)-1));
        allCSaligned.clustering.proj3(m,3) = dot(lastCluster(window,m),v(:,tmp(1,2)-2));
    end

    allCSaligned.clusters.ids3 = kmeans(allCSaligned.clustering.proj3,numClusts);
    
else
    
    allCSaligned.clusters.ids3 = kmeans(lastCluster',numClusts);

end

% sort away
[y3,allCSaligned.sort.indices3]  = sort(allCSaligned.clusters.ids3,1,'ascend');
lastClusterPSTH                 = lastCluster(:,allCSaligned.sort.indices3);

disp('Made it through all clustering');
allCSaligned.sorted           = [lowClusterPSTH,midClusterPSTH,lastClusterPSTH];
allCSaligned.finalInds        = [tmpDaInds;tmpNonDaInds2(allCSaligned.sort.indices2);tmpNonDaInds3(allCSaligned.sort.indices3)]
allCSaligned.finalIds         = [ones(size(tmpDaInds));y2+1;y3+13];

figure(6);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
subplot(1,5,1:2);
tmp = max(max(allCSaligned.sorted));
imagesc(allCSaligned.sorted',[-10 10]);
xlabel('Time (ms)');
ylabel('Cell index');
title('CS Aligned PSTH');

subplot(1,5,3);
plot(allCSaligned.finalIds,701:-1:1,'k.'); axis([0 20 1 701]); axis off;

colormap(mapName);
subplot(1,5,4:5);
imagesc(corr(allCSaligned.sorted),[-1,1]);
ylabel('Cell index');
xlabel('Cell index');
title('CS Aligned Signal Correlations');

%% SIMPLE APPROACH - SORT BASED UPON RESPONSE MAGNITUDE
% New approach: try a cascade of clustering:
%     Step 1: cluster based upon responses during the CS (Dopamine vs. Non-Dopamine)
%     Step 2: separate out the groups into 3 super respond, sus

numClusts = 6;
allCSaligned.clusters.numClusts = numClusts;
window = 2000:2200;
window2 = 2500:4000;

responseLatency = zeros(size(allCSaligned.psthZ,2),1);
responsePeak    = zeros(size(allCSaligned.psthZ,2),1);
responseWidth   = zeros(size(allCSaligned.psthZ,2),1);
responseAmp     = zeros(size(allCSaligned.psthZ,2),1);
responseMid     = zeros(size(allCSaligned.psthZ,2),1);

for i = 1:size(allCSaligned.psthZ,2)
    responseSUS(i,1) = trapz(allCSaligned.psthZ(window2,i))./150;
    tmpDATA = allCSaligned.psthZ(window,i);
    peakLoc = find(tmpDATA>3);
    % how do i make sure all of the points are in order? as opposed to skipping points

    if length(peakLoc)>5

        for j=2:length(peakLoc)
            if peakLoc(j) ~= peakLoc(j-1)+1
                peakLoc = peakLoc(1:j-1);
                break;
            end
        end
            
        responseAmp(i)      = trapz(tmpDATA(peakLoc))./length(peakLoc);

        if responseAmp(i)>3
            responseLatency(i)  = find(tmpDATA>2,1);
            responsePeak(i)     = find(tmpDATA==max(tmpDATA),1);
            responseMid(i,1)      = mean(peakLoc);
        else
            responseLatency(i) = 0;
            responsePeak(i) = 0;
            responseAmp(i)  = 0;
        end

    else
        
        peakLoc = find(tmpDATA<-3);
        
        if length(peakLoc)>5
            
            for j=2:length(peakLoc)
                if peakLoc(j) ~= peakLoc(j-1)+1
                    peakLoc = peakLoc(1:j-1);
                    break;
                end
            end

            responseAmp(i)      = trapz(tmpDATA(peakLoc))./length(peakLoc);

            if responseAmp(i)<-2
                responseLatency(i)  = find(tmpDATA<-2,1);
                responsePeak(i)     = find(tmpDATA==min(tmpDATA),1);
                responseMid(i,1)    = mean(peakLoc);
            else
                responseLatency(i) = 0;
                responsePeak(i) = 0;
                responseAmp(i)  = 0;
            end

        end
        
    end
    
    responseWidth(i) = length(peakLoc);
    
end


[y,allCSaligned.sort.amp.inds]  = sort(responseAmp,1,'descend');
allCSaligned.sort.amp.sortedZ   = allCSaligned.psthZ(:,allCSaligned.sort.amp.inds);
allCSaligned.sort.amp.sortedS   = allCSaligned.psthS(:,allCSaligned.sort.amp.inds);

[y,allCSaligned.sort.amp.inds2] = sort(responsePeak,1,'descend');
allCSaligned.sort.amp.sortedZ2  = allCSaligned.psthZ(:,allCSaligned.sort.amp.inds2);
allCSaligned.sort.amp.sortedS2  = allCSaligned.psthS(:,allCSaligned.sort.amp.inds2);
% 
% allCSaligned.sort.km.ids        = kmeans([responseAmp,responseLatency],numClusts);
% [y,allCSaligned.sort.amp.inds3] = sort(allCSaligned.sort.km.ids,1,'ascend');
% allCSaligned.sort.amp.sortedZ3  = allCSaligned.psthZ(:,allCSaligned.sort.amp.inds3);
% allCSaligned.sort.amp.sortedS3  = allCSaligned.psthS(:,allCSaligned.sort.amp.inds3);


figure(6);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
subplot(1,5,1:2);
imagesc(allCSaligned.sort.amp.sortedZ(1500:2650,:)',[-10 10]);
xlabel('Time (ms)');
ylabel('Cell index');
title('CS Aligned PSTH');

% subplot(1,5,3);
% plot(allCSaligned.sort.km.ids(allCSaligned.sort.amp.inds3),701:-1:1,'k.'); axis([0 numClusts+1 1 701]); axis off;

subplot(1,5,4:5);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
% imagesc(corr(allCSaligned.sort.amp.sorted),[-1,1]);
imagesc(allCSaligned.sort.amp.sortedZ2(:,:)',[-10 10]);
ylabel('Cell index');
xlabel('Time (ms)');
title('CS Aligned PSTH');

%% Combine all CS aligned, extinction PSTHs into a single matrix for classification
clear allCSEXaligned
NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
    % check if learning
    if strcmp(PopData.session(i).sessClass,'extinction')
        
        if i<10
            disp(['Session 00' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
        elseif i<100
            disp(['Session 0' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
        else
            disp(['Session ' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
        end
        
        numUnits = size(PopData.session(i).unit,2);

        for j = 1:numUnits
            if isnan(PopData.session(i).unit(j).isi.hist.logCount(1,1))
                disp('Invalid ISI distribution found');
            else
                allCSEXaligned.psthS(:,k)   = PopData.session(i).unit(j).respCSS.image.psthAVG'.*1000;
                allCSEXaligned.psthZ(:,k)   = PopData.session(i).unit(j).respCSS.image.psthZ';
                allCSEXaligned.psthZe(:,k)  = PopData.session(i).unit(j).respCSS.image.psthZe';
                allCSEXaligned.boxcar(:,k)  = PopData.session(i).unit(j).respCS.boxcar'.*100;
                allCSEXaligned.sess(1,k)  = i;
                allCSEXaligned.unit(1,k)  = j;
                allCSEXaligned.id(1,k).sessID = PopData.session(i).sessId;
                allCSEXaligned.id(1,k).name   = PopData.session(i).unit(j).name;

                k=k+1;
            end
        end
    end
    
end

k=k-1;

disp(['Total units: ' num2str(k)])

figure(7);
subplot(121);
[mapName] = TNC_CreateRBColormap(1024,'rb');
imagesc(allCSEXaligned.psthZ',[-15 15]);
colormap(mapName);
xlabel('Time (ms)');
ylabel('Cell index');
title('CSe Aligned ZDF');

subplot(122);
imagesc(allCSEXaligned.boxcar',[-100 100]);
colormap(mapName);
xlabel('Time (ms)');
ylabel('Cell index');
title('CSe Aligned Boxcar');

%% PSTH classifier CS Extinction

% option 1: kmeans clustering of psth
numClusts = 7;
clusterEX = 0;

if clusterEX
    allCSEXaligned.clusters.numClusts = numClusts;
    allCSEXaligned.clusters.ids = kmeans(allCSEXaligned.psthZ',numClusts);
    [y,allCSEXaligned.sort.indices] = sort(allCSEXaligned.clusters.ids,1,'ascend');
    allCSEXaligned.sorted = allCSEXaligned.psthZ(:,allCSEXaligned.sort.indices);
else
    allCSEXaligned.sorted = allCSEXaligned.psthZ(:,allCSaligned.finalInds);
end

figure(7); clf;
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
subplot(121);
tmp = max(max(allCSEXaligned.sorted));
imagesc(allCSEXaligned.sorted',[-10 10]);
xlabel('Time (ms)');
ylabel('Cell index');
title('CS Aligned PSTH');

colormap(mapName);
subplot(122);
imagesc(corr(allCSEXaligned.sorted),[-1,1]);
ylabel('Cell index');
xlabel('Cell index');
title('CS Aligned Pairwise Correlations');

%% Average Kmeans clusters to compare response types for learning and extinction
disp('...CS...');
for m=1:max(allCSaligned.finalIds)
    currClust = find(allCSaligned.finalIds==m);
    disp(['Members of cluster ' num2str(m) ': ' num2str(length(currClust))]);

    tmpMean = allCSaligned.sorted(:,currClust(1):currClust(length(currClust)));
    respGroupCS.avg(:,m) = mean(tmpMean,2);
    respGroupCS.err(:,m) = std(tmpMean,0,2) ./ sqrt(length(currClust));
    tmpMeanEX = allCSEXaligned.psthZ(:,allCSaligned.finalInds(currClust));
    respGroupCSex.avg(:,m) = mean(tmpMeanEX,2);
    respGroupCSex.err(:,m) = std(tmpMeanEX,0,2) ./ sqrt(length(currClust));
end

% disp('...US...');
% for m=1:allUSaligned.clusters.numClusts
%     currClust = find(allUSaligned.clusters.ids==m);
%     disp(['Members of cluster ' num2str(m) ': ' num2str(length(currClust))]);
%     tmpMean = allUSaligned.psthZ(:,currClust);
%     respGroupUS.avg(:,m) = mean(tmpMean,2);
%     respGroupUS.err(:,m) = std(tmpMean,0,2) ./ sqrt(length(currClust));
% end

%% For each response type in learning and extinction take comparison
disp('...CS...');
for m=1:allCSaligned.clusters.numClusts
    currClust = find(allCSaligned.clusters.ids==m);
    disp(['Members of cluster ' num2str(m) ': ' num2str(length(currClust))]);
    tmpMean = allCSaligned.psthZ(:,currClust);
    respGroupCS.avg(:,m) = mean(tmpMean,2);
    respGroupCS.err(:,m) = std(tmpMean,0,2) ./ sqrt(length(currClust));
    clustSizesCS(m) = length(currClust);
end

disp('...CSex...');
for m=1:allCSaligned.clusters.numClusts
    currClust = find(allCSaligned.clusters.ids==m);
    disp(['Members of cluster ' num2str(m) ': ' num2str(length(currClust))]);
    tmpMean = allCSEXaligned.psthZ(:,currClust);
    respGroupCSex.avg(:,m) = mean(tmpMean,2);
    respGroupCSex.err(:,m) = std(tmpMean,0,2) ./ sqrt(length(currClust));
    clustSizesEX(m) = length(currClust);
end

%% PSTH classifier US

% option 1: kmeans clustering of psth
numClusts = 6;
allUSaligned.clusters.numClusts = numClusts;
allUSaligned.clusters.ids = kmeans(allUSaligned.psthZ',numClusts);
[y,allUSaligned.sort.indices] = sort(allUSaligned.clusters.ids,1,'ascend');

allUSaligned.sorted = allUSaligned.psthZ(:,allUSaligned.sort.indices);

figure(7);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
subplot(121);
tmp = max(max(allUSaligned.sorted));
imagesc(allUSaligned.sorted',[-tmp tmp]);
xlabel('Time (ms)');
ylabel('Cell index');
title('US Aligned PSTH');

colormap(mapName);
subplot(122);
imagesc(corr(allUSaligned.sorted),[-1,1]);
ylabel('Cell index');
xlabel('Cell index');
title('US Aligned Pairwise Correlations');

%% Combine all ISI distributions into a single matrix
clear allISI
NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
    % check if learning
    if strcmp(PopData.session(i).sessClass,'extinction')
        numUnits = size(PopData.session(i).unit,2);
        for j = 1:numUnits
            if isnan(PopData.session(i).unit(j).isi.hist.logCount(1,1))
                disp('Invalid ISI distribution found');
            else
                allISI.raw(:,k) = PopData.session(i).unit(j).isi.hist.logCount';
                k=k+1;
            end
        end
    end
    
end
k=k-1;
disp(['Total units: ' num2str(k)])
        
%% ISI classifier

% option 1: kmeans clustering of log-scaled isi pdf
numClusts = 4;
allISI.clusters.numClusts = numClusts;
allISI.clusters.ids = kmeans(allISI.raw',numClusts);
[y,allISI.sort.indices] = sort(allISI.clusters.ids,1,'ascend');

allISI.sorted = allISI.raw(:,allISI.sort.indices);
% allISI.sorted = allISI.raw(:,allCSaligned.sort.indices);
% allISI.sorted.corr = corr(allISI.sorted);

logSpace          = 0:0.1:5;
allISI.logTimes = 10.^logSpace;

for i=1:size(allISI.raw,2)
    % calculate the entropy
    tmp = -1.*(allISI.raw(:,i).*log2(allISI.raw(:,i)));
    tmpInds = find(tmp>0);
    allISI.entropy(i) = sum(tmp(tmpInds));
    tmpMaxLoc = find(allISI.raw(:,i)==max(allISI.raw(:,i)),1);
    allISI.avLogRate(i) = allISI.logTimes(tmpMaxLoc);
end

figure(5);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
subplot(121);
tmp = max(max(allISI.sorted));
imagesc(allISI.sorted',[-tmp tmp]);
xlabel('ISI index');
ylabel('Cell index');

colormap(mapName);
subplot(122);
imagesc(corr(allISI.sorted),[-1,1]);
ylabel('Cell index');
xlabel('Cell index');

% figure(4);
% semilogy(allISI.entropy,allISI.avLogRate,'k.');

% option 2: projection onto first 3 principal components of isi covariance matrix
% allISI.covMat = cov(allISI.sorted);

%% WAVEFORM COLLECTOR AND DISPLAY OF PROPERTIES

NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions

    if strcmp(PopData.session(i).sessClass,'learning')
        numUnits = size(PopData.session(i).unit,2);
        
        for j = 1:numUnits
            wfSize = size(PopData.session(i).unit(j).WfMean);
            if wfSize(1,1)>1
                allWfs.wf(1:wfSize(1,1),k) = PopData.session(i).unit(j).WfMean;
                allWfs.area(k) = trapz(abs(PopData.session(i).unit(j).WfMean));
                allWfs.min(k) = min(PopData.session(i).unit(j).WfMean);
                allWfs.max(k) = max(PopData.session(i).unit(j).WfMean);
                k=k+1;
            end
        end        
    end
    
end

k

% disp(['Total units: ' num2str(k)])
% figure(6);
% plot3(allWfs.area,allWfs.min,allWfs.max,'ko');

%% QUANTIFY contrasts between learning and extinction cell by cell 
figure(5); clf;

thresholdSD = 3;

timeWindow(1).indices = 2000:2030; % early response 1
timeWindow(2).indices = 2030:2060; % early response 2
timeWindow(3).indices = 2060:2150; % dopamine response
timeWindow(4).indices = 2500:3500; % trace period
timeWindow(5).indices = 5000:6000; % consumption

for m = 1:size(timeWindow,2)

    for k = 1:size(allCSaligned.psthZ,2)
        
        quantRespMag.contrast(m).unit(k,1) = (trapz(allCSaligned.psth(timeWindow(m).indices,k))./size(timeWindow(m).indices,2));
        quantRespMag.contrast(m).unit(k,2) = (trapz(allCSEXaligned.psth(timeWindow(m).indices,k))./size(timeWindow(m).indices,2));
        if abs(quantRespMag.contrast(m).unit(k,2)) > thresholdSD || abs(quantRespMag.contrast(m).unit(k,1)) > thresholdSD
            quantRespMag.contrast(m).unit(k,3) = quantRespMag.contrast(m).unit(k,2)-quantRespMag.contrast(m).unit(k,1);
        else
            quantRespMag.contrast(m).unit(k,3) = NaN;
        end
        quantRespMag.contrast(m).unit(k,4) = allCSaligned.clusters.ids(k);
    end
    
    tmpIns = find(1-isnan(quantRespMag.contrast(m).unit(:,3)))
    quantRespMag.valid(m).inds = tmpIns;
    tmpVals     = quantRespMag.contrast(m).unit(quantRespMag.valid(m).inds,3);
    tmpClusts   = quantRespMag.contrast(m).unit(quantRespMag.valid(m).inds,4);
    quantRespMag.valid(m).values = tmpVals;
    quantRespMag.valid(m).clusters = tmpClusts;
    
    histArray = -100:0.5:100;
    figure(5);
    subplot(size(timeWindow,2)+1,1,m);
    hist(quantRespMag.contrast(m).unit(:,3),histArray);
    peakVal = max(hist(quantRespMag.contrast(m).unit(:,3),histArray))+5;
    peakVal = 30;
    hold on; plot([0,0],[0,100],'k--'); plot(median(quantRespMag.valid(m).values),peakVal,'kv','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);plot(mean(quantRespMag.valid(m).values),peakVal,'rv','MarkerSize',10,'LineWidth',2); 
    name = ['Epoch: +' num2str(min(timeWindow(m).indices)-2000) '-' num2str(max(timeWindow(m).indices)-2000) 'ms (CS=0ms)']
    title(name); ylabel('Count'); 
    if m==size(timeWindow,2)
        xlabel('Ext - Acq (Z-score)');
    end
    axis([-20 20 -1 35]);
%     
%     subplot(size(timeWindow,2)+1,1,size(timeWindow,2)+1);
%     byClusterHist = (hist(quantRespMag.valid(m).clusters,1:8));
%     if m==1 
%         plot(byClusterHist./clustSizesCS,'ro-');
%         xlabel('Cluster ID'); ylabel('Fraction of cluster');
%     end
%     if m==2 
%         plot(byClusterHist./clustSizesCS,'mo-');
%     end
%     if m==3 
%         plot(byClusterHist./clustSizesCS,'co-');
%     end
%     if m==4 
%         plot(byClusterHist./clustSizesCS,'bo-');
%     end
%     if m==5 
%         plot(byClusterHist./clustSizesCS,'ko-');
%     end
%     axis([0 9 -0.2 1.2]); grid on;
%     hold on;
end

%% QUANTIFY response to CS in learning and extinction 

thresholdSD = 3;

% timeWindow(1).indices = 2000:2030; % early response 1
% timeWindow(2).indices = 2030:2060; % early response 2
% timeWindow(3).indices = 2060:2150; % dopamine response
% timeWindow(4).indices = 2500:3500; % trace period
% timeWindow(5).indices = 5000:6000; % consumption
timeWindow(1).indices = 1950:2150;
% timeWindow(2).indices = 2060:2160;

    for k = 1:size(allCSaligned.psthZ,2)
        
        tmpCS   = allCSaligned.psthZ(timeWindow(1).indices,k);
        tmpCSe  = allCSEXaligned.psthZ(timeWindow(1).indices,k);
        tmpCSr  = allCSaligned.psth(timeWindow(1).indices,k);
        tmpCSre = allCSEXaligned.psth(timeWindow(1).indices,k);
        
        if max(abs(tmpCS)) > max(abs(tmpCSe))
            if max(tmpCS) > min(abs(tmpCS))
                allSupra = find(tmpCS>3);
                numSupra = length(allSupra);
                if numSupra>5 && length(allSupra)<140
                    responseMag.unit(k,1) = trapz(tmpCS(allSupra(1):allSupra(numSupra)))./numSupra;
                    responseMag.unit(k,2) = trapz(tmpCSe(allSupra(1):allSupra(numSupra)))./numSupra;
                    responseMag.raw(k,1)  = trapz(tmpCSr(allSupra(1):allSupra(numSupra)));
                    responseMag.raw(k,2)  = trapz(tmpCSre(allSupra(1):allSupra(numSupra)));
                    responseMag.peak(k,1) = allSupra(1);
                    responseMag.peak(k,2) = find(tmpCS==max(tmpCS),1);
                else
                    responseMag.unit(k,1) = 0;
                    responseMag.unit(k,2) = 0;
                    responseMag.raw(k,1)  = 0;
                    responseMag.raw(k,2)  = 0;
                    responseMag.peak(k,1) = 0;
                    responseMag.peak(k,2) = 0;
                end                    
            else
%                 allSupra = find(tmpCS<-4);
%                 numSupra = length(allSupra);
%                 if numSupra>1 && allSupra(1)<110
%                     responseMag.unit(k,1) = trapz(tmpCS(allSupra(1):allSupra(numSupra)))./numSupra;
%                     responseMag.unit(k,2) = trapz(tmpCSe(allSupra(1):allSupra(numSupra)))./numSupra;
%                     responseMag.peak(k,1) = allSupra(1);
%                     responseMag.peak(k,2) = find(tmpCS==min(tmpCS),1);
%                 else
                    responseMag.unit(k,1) = 0;
                    responseMag.unit(k,2) = 0;
                    responseMag.raw(k,1)  = 0;
                    responseMag.raw(k,2)  = 0;
                    responseMag.peak(k,1) = 0;
                    responseMag.peak(k,2) = 0;
%                 end                                    
            end
        else
            if max(tmpCSe) > min(abs(tmpCSe))
                allSupra = find(tmpCSe>3);
                numSupra = length(allSupra);
                if numSupra>5 && length(allSupra)<140
                    responseMag.unit(k,1) = trapz(tmpCS(allSupra(1):allSupra(numSupra)))./numSupra;
                    responseMag.unit(k,2) = trapz(tmpCSe(allSupra(1):allSupra(numSupra)))./numSupra;
                    responseMag.raw(k,1)  = trapz(tmpCSr(allSupra(1):allSupra(numSupra)));
                    responseMag.raw(k,2)  = trapz(tmpCSre(allSupra(1):allSupra(numSupra)));
                    responseMag.peak(k,1) = allSupra(1)+10;
                    responseMag.peak(k,2) = find(tmpCSe==max(tmpCSe),1);
                else
                    responseMag.unit(k,1) = 0;
                    responseMag.unit(k,2) = 0;
                    responseMag.raw(k,1)  = 0;
                    responseMag.raw(k,2)  = 0;
                    responseMag.peak(k,1) = 0;
                    responseMag.peak(k,2) = 0;
                end                                    
            else
%                 allSupra = find(tmpCSe<-4);
%                 numSupra = length(allSupra);
%                 if numSupra>1 && allSupra(1)<110
%                     responseMag.unit(k,1) = trapz(tmpCS(allSupra(1):allSupra(numSupra)))./numSupra;
%                     responseMag.unit(k,2) = trapz(tmpCSe(allSupra(1):allSupra(numSupra)))./numSupra;
%                     responseMag.peak(k,1) = allSupra(1);
%                     responseMag.peak(k,2) = find(tmpCSe==min(tmpCSe),1);
%                 else
                    responseMag.unit(k,1) = 0;
                    responseMag.unit(k,2) = 0;
                    responseMag.raw(k,1)  = 0;
                    responseMag.raw(k,2)  = 0;
                    responseMag.peak(k,1) = 0;
                    responseMag.peak(k,2) = 0;
%                 end                                    
            end            
        end
        
%         responseMag.unit(k,1) = (trapz(allCSaligned.psthZ(timeWindow(1).indices,k))./size(timeWindow(1).indices,2));
%         responseMag.unit(k,2) = (trapz(allCSEXaligned.psthZ(timeWindow(1).indices,k))./size(timeWindow(1).indices,2));
% 
%         responseMag.unit(k,3) = (trapz(allCSaligned.psthZ(timeWindow(2).indices,k))./size(timeWindow(2).indices,2));
%         responseMag.unit(k,4) = (trapz(allCSEXaligned.psthZ(timeWindow(2).indices,k))./size(timeWindow(2).indices,2));
% 
%         responseMag.unit(k,1) = (max(allCSaligned.psthZ(timeWindow(1).indices,k)));
%         responseMag.unit(k,2) = (max(allCSEXaligned.psthZ(timeWindow(1).indices,k)));
% 
%         responseMag.unit(k,3) = (max(allCSaligned.psthZ(timeWindow(2).indices,k)));
%         responseMag.unit(k,4) = (max(allCSEXaligned.psthZ(timeWindow(2).indices,k)));
% 
%         responseMag.unit(k,1) = (sum(allCSaligned.psth(timeWindow(1).indices,k)));
%         responseMag.unit(k,2) = (sum(allCSEXaligned.psth(timeWindow(1).indices,k)));

    end

% plot(responseMag.unit(:,3),responseMag.unit(:,4),'ro',responseMag.unit(:,1),responseMag.unit(:,2),'k.',-10:40,-10:40,'r-',-10:40,(-10:40)+3,'r--',-10:40,(-10:40)-3,'r--'); ylabel('Extinction'); xlabel('Learning');
figure(5); clf;
subplot(121)
plot(-10:50,-10:50,'r-',-10:50,(-10:50)+3,'r--',-10:50,(-10:50)-3,'r--',[-10 45],[3 3],'k--',[3 3],[-10 45],'k--',responseMag.unit(:,1),responseMag.unit(:,2),'k.'); ylabel('Extinction'); xlabel('Learning');
box off;axis([-10 45 -10 45]);
subplot(122)
plot(responseMag.peak(:,2)-50,responseMag.unit(:,2)-responseMag.unit(:,1),'k.',[0 200], [3 3], 'r--',[0 200], [-3 -3], 'r--',[0 200], [0 0], 'r-');
xlabel('Latency (ms)'); ylabel('Extinction - Learning');
box off; %axis([0 200 -40 40]);

figure(2); clf;
% plot(-10:50,-10:50,'r-',[-10 45],[3 3],'k--',[3 3],[-10 45],'k--',responseMag.raw(:,1),responseMag.raw(:,2),'k.'); ylabel('Extinction'); xlabel('Learning');
plot(responseMag.peak(:,2)-50,responseMag.raw(:,2)-responseMag.raw(:,1),'k.');
% hist(responseMag.peak(:,2)-50,-50:1:200);
% box off; axis([-20 200 -1 35]);

disp(['Number of phasic responders: ' num2str(length(find(responseMag.peak(:,2)>0)))]);

% plot(responseMag.unit(:,1),responseMag.unit(:,2),'k.',0:20,0:20,'r--'); ylabel('Extinction'); xlabel('Learning');

%% get all US times to create separately aligned PSTH for CS and reward
NumSessions = size(PopData.session,2);

for i=1:NumSessions
    % check if learning
    if strcmp(PopData.session(i).sessClass,'learning')
        disp(['Session ' num2str(i)]);
        if length(PopData.session(i).events.US.ts) > 1
        	rewardOffset(i) = mean(PopData.session(i).events.US.ts(11:21)) - mean(PopData.session(i).events.CS.ts(10:20));
        end
    end
end

%% Multipeak fitting of peak time measuremnt

% Peaks:
% 0.56
% 12.6
% 45.5
% 78.1
% 166.0
% 
% with FWHMs:
% 3.3
% 4.6
% 12.5
% 50.5
% 21.1

amps        = ones(1,5);
peaks       = [0.6,12.6,45.5,78.1,166.0];
fwhmsL      = [3.3, 9.2,20.0,17.5, 41.0];
fwhmsH      = [3.3, 12.9,15.1,46.9,50.0];

figure(6); plot(responsePeak,responseAmp,'k.',peaks,amps,'r.',peaks+fwhmsH,amps,'ro',peaks-fwhmsL,amps,'ro'); drawnow;

% Classify the neurons based upon responseLatency
numCells    = size(responsePeak,1);
[y,inds]    = sort(responsePeak,1,'descend');

numPeaks    = length(peaks);
thresholds  = peaks+fwhmsH;

for i=1:numPeaks
    
    if i>1
        tmpStruct.class(i).inds = find(y>thresholds(i-1) & y<=thresholds(i));
        trackSort.real(i).inds  = inds(tmpStruct.class(i).inds);
    else
        tmpStruct.class(i).inds = find(y<=thresholds(i));
        trackSort.real(i).inds  = inds(tmpStruct.class(i).inds);
    end
    
    latencyGroups(tmpStruct.class(i).inds) = i;
    
end

figure(6); 
subplot(1,4,1:3);
colormap(mapName);
imagesc(allCSaligned.sort.amp.sortedZ2(1500:2650,:)',[-10 10]);
subplot(1,4,4);
plot(latencyGroups,663:-1:1,'k.'); axis([0 numPeaks+1 0 663]); axis off;

%% for each latency peak sort the responses by amplitude

totalGroups = max(latencyGroups);

for i=1:totalGroups

    % find the cells matching the current latency peak
    indsToGroup     = trackSort.real(i).inds;
    tmpAmpData      = responseAmp(indsToGroup);
    tmpSubData      = allCSaligned.psthZ(:,indsToGroup);
    tmpSubDataEX    = allCSEXaligned.psthZ(:,indsToGroup);
    [y,inds]        = sort(tmpAmpData,1,'descend');

    tmpSubDataSort      = tmpSubData(:,inds);
    tmpSubDataEXSort    = tmpSubDataEX(:,inds);
    
    allCSaligned.sort.latency.group(i).indAMP = inds;
    
    % within those cells sort by response amplitude
    if i>1
        allCSaligned.sort.ampPeak.sortedZ = [tmpSubDataSort,allCSaligned.sort.ampPeak.sortedZ];
        allCSEXaligned.sort.ampPeak.sortedZ = [tmpSubDataEXSort,allCSEXaligned.sort.ampPeak.sortedZ];
    else
        allCSaligned.sort.ampPeak.sortedZ = tmpSubDataSort;
        allCSEXaligned.sort.ampPeak.sortedZ = tmpSubDataEXSort;
    end

    figure(3); subplot(1,7,1:3); colormap(mapName);
    imagesc(allCSaligned.sort.ampPeak.sortedZ(1500:2500,:)',[-15 15]); drawnow; pause(1);

    figure(7); subplot(1,7,1:3); colormap(mapName);
    imagesc(allCSEXaligned.sort.ampPeak.sortedZ(1500:2500,:)',[-15 15]); drawnow; pause(1);

    figure(8); 
    imagesc(allCSEXaligned.sort.ampPeak.sortedZ(1500:2500,:)'-allCSaligned.sort.ampPeak.sortedZ(1500:2500,:)',[-15 15]); drawnow; pause(1);
    colormap(mapName);
end

%% for each latency peak sort the responses by the sustained response

totalGroups = max(latencyGroups);

for i=1:totalGroups

    % find the cells matching the current latency peak
    indsToGroup         = trackSort.real(i).inds;
    tmpSusData          = responseSUS(indsToGroup);
    tmpSubData          = allCSaligned.psthZ(:,indsToGroup);
    tmpSubDataS         = allCSaligned.psthS(:,indsToGroup);
    tmpSubDataEX        = allCSEXaligned.psthZ(:,indsToGroup);
    tmpSubDataSEX       = allCSEXaligned.psthS(:,indsToGroup);
    [y,inds]            = sort(tmpSusData,1,'descend');

    tmpSubDataSort      =   tmpSubData(:,inds);
    tmpSubDataEXSort    = tmpSubDataEX(:,inds);
    
    allCSaligned.sort.latency.group(i).indSUS   = inds;

    allCSaligned.sort.latency.group(i).psthZ    = tmpSubData;
    allCSaligned.sort.latency.group(i).psthS    = tmpSubDataS;

    allCSEXaligned.sort.latency.group(i).psthZ  = tmpSubDataEX;
    allCSEXaligned.sort.latency.group(i).psthS  = tmpSubDataSEX;

    figure(4); subplot(5,1,i); imagesc(cov(tmpSubDataSort'),[-1 1]); drawnow;
    allCSaligned.sort.latency.group(i).avg      = mean(tmpSubDataS,2);
    allCSaligned.sort.latency.group(i).err      = std(tmpSubDataS,0,2) ./ sqrt(length(inds)-1);
    allCSEXaligned.sort.latency.group(i).avg    = mean(tmpSubDataSEX,2);
    allCSEXaligned.sort.latency.group(i).err    = std(tmpSubDataSEX,0,2) ./ sqrt(length(inds)-1);
    figure(5); subplot(1,5,i); plot(1:5001,allCSaligned.sort.latency.group(i).avg,'k',1:5001,allCSEXaligned.sort.latency.group(i).avg,'r'); drawnow;

    % within those cells sort by response amplitude
    if i>1
        allCSaligned.sort.susPeak.sortedZ = [tmpSubDataSort,allCSaligned.sort.susPeak.sortedZ];
        allIds = [ones(1,length(indsToGroup)).*i,allIds];
        allCSEXaligned.sort.susPeak.sortedZ = [tmpSubDataEXSort,allCSEXaligned.sort.susPeak.sortedZ];
    else
        allCSaligned.sort.susPeak.sortedZ = tmpSubDataSort;
        allIds = ones(1,length(indsToGroup)).*i;
        allCSEXaligned.sort.susPeak.sortedZ = tmpSubDataEXSort;
    end

    figure(3); subplot(1,7,5:7); colormap(mapName);
    imagesc(allCSaligned.sort.susPeak.sortedZ(1500:4000,:)',[-20 20]); drawnow; pause(1);

    figure(7); subplot(1,7,5:7); colormap(mapName);
    imagesc(allCSEXaligned.sort.susPeak.sortedZ(1500:4000,:)',[-20 20]); drawnow; pause(1);

    disp(num2str(corr(allCSaligned.sort.latency.group(i).indAMP,allCSaligned.sort.latency.group(i).indSUS)));
end

figure(3); subplot(1,7,4);
plot(allIds,663:-1:1,'k.'); axis([0 totalGroups+1 1 663]); axis off;

%% To facilitate checking data print a list of my indices and the names of sessions units:
countUnits  = 0;
NumSessions = size(PopData.session,2);
k           = 1;

fid = fopen('ListAllData.csv','w');
fprintf(fid,['\n________________Listing of All Data_________________________\n']);
fprintf(fid,['Timestamp: ' datestr(now) '\n']);

for i=1:663

    sessInd = allCSaligned.sess(i);
    unitInd = allCSaligned.unit(i);

    fprintf(fid,[num2str(i) ',' num2str(sessInd) ',' num2str(unitInd) ',' PopData.session(sessInd).sessId ',' PopData.session(sessInd).unit(unitInd).name '\n']);

end

fclose(fid);
