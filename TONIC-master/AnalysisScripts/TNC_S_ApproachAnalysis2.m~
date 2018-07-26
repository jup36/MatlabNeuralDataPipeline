% Approach behavior analysis script
%% STANDARD ANALYSIS PARAMETERS
% Define the smoothing to apply to the timeseries data
currParams.smthParams.rise     = 1;
% currParams.smthParams.decay    = 100;
currParams.smthParams.decay    = 1;
% currParams.smthParams.decay    = 10;
% currParams.smthParams.decay    = 25;
% currParams.smthParams.decay    = 60;
currParams.filter.style = 0; % 0: Gaussian; 1: Causal Gaussian; 2: Boxcar

switch currParams.filter.style
    
    case 0
        [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
        
    case 1
        [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
        currParams.filter.kernel(1:currParams.smthParams.decay.*15) = 0;
                
    case 2
        currParams.filter.kernel = zeros(1,currParams.smthParams.decay.*3);
        currParams.filter.kernel(currParams.smthParams.decay:currParams.smthParams.decay.*2) = 1./currParams.smthParams.decay;

end

currParams.winParams.prior     = 1e3;
currParams.winParams.after     = 3e3; 

currParams.popVec.winSize       = 1;

currParams.stableTrials         = 10; %derived from summary behavior data

% hand curated dopamine neuron list:
PopData.daList = [42,67,161,177,187,329,347,503,504,514,522,523,528,529,530,532,539,540,541,544,545,546];

PopData.currParams = currParams;

halfVal = max(currParams.filter.kernel) ./ 2;

figure(1); clf; plot(PopData.currParams.filter.kernel,'ko-','LineWidth',2,'Color',[0 0 0]);

disp(' ');
disp(' ');
disp('________________________________');
disp(' ');
disp('Initialized analysis parameters.');
disp(['Kernel FWHM: ' num2str( find(currParams.filter.kernel>halfVal,1,'last') - find(currParams.filter.kernel>halfVal,1) ) ' ms']);
disp('________________________________');
disp(' ');
disp(' ');

%% ADD NEW DATA TO THE STRUCTURE

% fileName = 'DA-f03-140501trace-002-09.nex';
fileName = 'F06-120501trace-Tstrong-001-02.nex';
clear newSession
[data] = TNC_LoadData(0, 0, fileName);
a = dir([fileName(1:numel(fileName)-4) '*nev'])
[dataB] = TNC_LoadData(0, 0, a.name)

addLogic = 1;
numUnits = size(data.neurons,1);
numSessions = numel(PopData.session)
% numSessions = 180

for pp=1:numUnits
    
    newSession.unit(pp).name    = data.neurons{pp}.name;

    currName                    = data.neurons{pp}.name;    
    newSession.unit(pp).el      = str2num(currName(5:7));
    newSession.unit(pp).un      = strfind('abcdefghijklm',currName(8));
    newSession.unit(pp).uID     = [num2str(newSession.unit(pp).el) '_' num2str(newSession.unit(pp).un)];
    
    newSession.unit(pp).ts      = double(round(data.neurons{pp}.timestamps.*1000));
    newSession.unit(pp).WfMean  = double(mean(data.waves{pp}.waveforms,2));
    newSession.unit(pp).WfStd   = double(std(data.waves{pp}.waveforms,[],2));
    
    newSession.unit(pp).isi     = TNC_QuantISI(newSession.unit(pp).ts);
    newSession.sessClass        = 'learning';
    newSession.sessId           = fileName(1:10);
    
    CSinds                      = find(dataB.Data.Spikes.Electrode==137 & dataB.Data.Spikes.Unit==1);
    newSession.events.CS.ts     = double(round(dataB.Data.Spikes.Timestamps(CSinds)./30));
    newSession.trials           = numel(CSinds);

    ELinds                      = find(dataB.Data.Spikes.Electrode==138);
    newSession.events.EL.ts     = double(round(dataB.Data.Spikes.Timestamps(ELinds)./30));
    
    USinds                      = find(dataB.Data.SerialDigitalIO.UnparsedData==65529);
    newSession.events.US.ts     = double(round(dataB.Data.SerialDigitalIO.TimeStamp(USinds)./30));
    
end

if addLogic
    disp(['Adding the current session to the PopData structure as session ' num2str(numSessions+1)]);
    PopData.session(numSessions+1).unit         = newSession.unit;
    PopData.session(numSessions+1).sessClass    = newSession.sessClass;
    PopData.session(numSessions+1).sessId       = newSession.sessId;
    PopData.session(numSessions+1).events       = newSession.events;
    PopData.session(numSessions+1).trials       = newSession.trials;
    
    PopData.session(numSessions+1)
end

%% EXTRACT LICKING BEHAVIOR DATA MEAN & TRIAL-BY-TRIAL
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
            [respCS] = TNC_AlignRasters(delta,PopData.session(i).events.EL.ts,currParams.stableTrials,PopData.session(i).events.CS.ts,[2e3,1e4],1,0);
            PopData.session(i).behavior.raster  = respCS.raster;
            PopData.session(i).behavior.psthAVG = respCS.image.psthAVG;
            PopData.session(i).behavior.psthSEM = respCS.image.psthSEM;
            
            validBehavSess = [validBehavSess,i];
            
        else
            
            PopData.session(i).behavior.raster  = 0;
            PopData.session(i).behavior.psthAVG = 0;
            PopData.session(i).behavior.psthSEM = 0;
            PopData.session(i).behavior.found   = 0;
            
        end
    
end

PopData.validBehavSess = validBehavSess;

disp(['___________________________________________________'])

%% CREATE LONG, TRIAL by TRIAL PSTHS AND RASTERS FOR ALL LEARNING CELLS

disp(['___________________________________________________'])
disp(['STARTED aligning all raster plots and updating the PopData structure...']);

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
% for i=181;        
    if strcmp(PopData.session(i).sessClass,'learning')

        numUnits = size(PopData.session(i).unit,2);
        countUnits = countUnits+numUnits;
        PopData.session(i).trials = size(PopData.session(i).events.CS.ts,1);
    
        if numUnits>=2 && strcmp(PopData.session(i).sessClass,'learning')
            numPairwise = numPairwise + sum(1:1:numUnits-1); % or +1 to count 
        end

        for j = 1:numUnits

            numStamps = length(PopData.session(i).unit(j).ts);
            delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
            delta(1,round(PopData.session(i).unit(j).ts)+1) = 1;

            [respCS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[1e3,1.1e4],1,1);
            PopData.session(i).unit(j).respCS.raster      = respCS.raster;
            PopData.session(i).unit(j).respCS.psthAVG     = respCS.image.psthAVG;

            tmpSmooth = conv(delta,currParams.filter.kernel,'same');
            [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[1e3,1.1e4],0,1);
            PopData.session(i).unit(j).respCSS.aligned    = respCSS.image.aligned;
            PopData.session(i).unit(j).respCSS.psthZ      = respCSS.image.psthZ;
            PopData.session(i).unit(j).respCSS.noZ        = respCSS.noZ;

            k=k+1;

        end

        disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' of ' num2str(NumSessions) ' | ' PopData.session(i).sessId '_' PopData.session(i).sessClass])
    end    
end
k = k-1;
k = k./2;
disp(['___________________________________________________'])
disp(['COMPLETED ... Total units: ' num2str(k) ' | ' num2str(countUnits) ' | Pairwise: ' num2str(numPairwise)]);

%% CHECK FOR SESSIONS WITH "POPULATION" RECORDINGS
disp(' ');
disp(' ');
disp(' ');

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
k=1; clear sessList;

for i=1:NumSessions
% for i=83;        
    if strcmp(PopData.session(i).sessClass,'learning')

        numUnits = size(PopData.session(i).unit,2);
        countUnits = countUnits+numUnits;
        PopData.session(i).trials = size(PopData.session(i).events.CS.ts,1);
            
        if numUnits>2 && PopData.session(i).behavior.found==1
            if numel(find(PopData.session(i).behavior.flCSon<5e3)) > 40
                sessList(k) = i;
                PopData.sessUsed.sessId(k) = i;
                PopData.sessUsed.numUnits(k) = numel(PopData.session(i).unit);
                k=k+1;
                disp(['Session: ' num2str(i) ' is a candidate session | ' PopData.session(i).sessId])
            end
        end
        

    end    
end
k = k-1;
disp(['___________________________________________________'])
disp([' '])
disp(['COMPLETED ... Valid Sessions: ' num2str(k)]);
disp(['          ... Total units: ' num2str(sum(PopData.sessUsed.numUnits))]);
disp(['          ... Median units/session: ' num2str(median(PopData.sessUsed.numUnits))]);
disp(['          ... Mean units/session: ' num2str(mean(PopData.sessUsed.numUnits)) ' +/- ' num2str(std(PopData.sessUsed.numUnits))]);
disp(['___________________________________________________'])

%% CLEAN UP PROBLEMATIC UNITS IN THE DATA STRUCTURE

% Need to walk through and check for strangely strong zero offset cross correlations

% Store the cleaned, finalized population data structure as a separate
% file

%% START HERE FOR PRE-SAVED DATASTRUCTURE...

%% CREATE A MATRIX OF ALL PETHs

count = 0;
clear allCSresp

for aa = 1:numel(sessList)
% for i = 83;
    i=sessList(aa);

    numUnits    = size(PopData.session(i).unit,2);
    
    for j=1:numUnits
        
        count = count + 1;
        allCSresp.psthZ(count,:) = PopData.session(i).unit(j).respCSS.psthZ;        
        allSpikeRate.avg(count,1) = mean(PopData.session(i).unit(j).respCS.psthAVG(1:1000));
        
    end
                
end


numClasses = 10;


[clustIds , c , sumdist] = kmeans(allCSresp.psthZ(:,1000:4000),numClasses,'Start','cluster');
[vals,inds] = sort(clustIds);


[mapName] = TNC_CreateRBColormap(1024,'cb');

figure(203); 
    subplot(131);
imagesc(allCSresp.psthZ(inds,:),[-2 2]); 
colormap(mapName);

    subplot(132);
imagesc(corr(allCSresp.psthZ(inds,500:2800)'),[-1 1]);
colormap(mapName);

    subplot(133);
imagesc(corr(allCSresp.psthZ),[-1 1]);
colormap(mapName);

%% PCA on time course
 
classOrder = [1 7 10 3 9 4 6 2 5 8];

colorMapNew = [ [numClasses - [1:numClasses]]' zeros(numClasses,1) [1:numClasses]'] ./ numClasses;

[vectors, values] = eig(cov(allCSresp.psthZ(:,1000:4000)));

varCaptured = diag(values)./sum(diag(values));
figure(210); subplot(211); semilogy(50:-1:1,varCaptured(2952:3001),'k.'); ylabel('Variance captured'); xlabel('Dimension');

subplot(212); colormap(colorMapNew);
plot(1:3001,vectors(:,2999:3001));  axis tight; 
sum(varCaptured(2999:3001));

clear projectio*

for i=1:numel(inds)
   
    for j=1:4
        projection(i,j) = dot(allCSresp.psthZ(i,1000:4000)',vectors(:,3001-(j-1)));
    end
    
end

for k=1:numClasses
    theseInds = find(clustIds==k);
    projectionMeans(k,:) = mean(projection(theseInds,:),1);
    projectionErrs(k,:) = std(projection(theseInds,:),0,1) ./ sqrt(numel(theseInds)-1);
end


figure(213); clf;
%     subplot(121); %scatter(projection(:,1),projection(:,2),2,clustIds); colormap(colorMapNew); hold on;
    scatter3(projectionMeans(classOrder,1),projectionMeans(classOrder,2),projectionMeans(classOrder,3),40,[1:numClasses],'filled'); colormap(colorMapNew);
hold on;
    for i = 1:numClasses
        plot3([0 projectionMeans(classOrder(i),1)] , [ 0 projectionMeans(classOrder(i),2)] , [ 0 projectionMeans(classOrder(i),3)],'-','Color',colorMapNew(i,:),'LineWidth',1);
        plot3(projectionMeans(classOrder(i),1)+[projectionErrs(classOrder(i),1) -projectionErrs(classOrder(i),1)] , [ projectionMeans(classOrder(i),2) projectionMeans(classOrder(i),2)] , [ projectionMeans(classOrder(i),3) projectionMeans(classOrder(i),3)],'-','Color',colorMapNew(i,:),'LineWidth',1);
        plot3([ projectionMeans(classOrder(i),1) projectionMeans(classOrder(i),1)] , projectionMeans(classOrder(i),2) + [projectionErrs(classOrder(i),2) -projectionErrs(classOrder(i),2)] , [ projectionMeans(classOrder(i),3) projectionMeans(classOrder(i),3)],'-','Color',colorMapNew(i,:),'LineWidth',1);
        plot3([ projectionMeans(classOrder(i),1) projectionMeans(classOrder(i),1)] , [ projectionMeans(classOrder(i),2) projectionMeans(classOrder(i),2)] , projectionMeans(classOrder(i),3) + [projectionErrs(classOrder(i),3) -projectionErrs(classOrder(i),3)],'-','Color',colorMapNew(i,:),'LineWidth',1);
    end
    plot3(xlim,[0 0],[0 0],'k--');
    plot3([0 0],ylim,[0 0],'k--');
    plot3([0 0],[0 0],zlim,'k--');
    for k=1:numClasses
        if projectionMeans(k,1)>0
            text(projectionMeans(k,1)+3,projectionMeans(k,2),projectionMeans(k,3),num2str(k));
        else
            text(projectionMeans(k,1)-3,projectionMeans(k,2),projectionMeans(k,3),num2str(k));
        end
    end
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
axis off; view([26 26]);

%     subplot(122); %scatter(projection(:,2),projection(:,3),2,clustIds); colormap(colorMapNew);
%     scatter(projectionMeans(:,2),projectionMeans(:,3),50,[1:numClasses],'filled'); colormap(colorMapNew);
% xlabel('PC3');
% ylabel('PC4');

%% Get mean response classes

    figure(205); clf;
    
    numClasses = 10;
    
for j=1:numClasses

    classOrder = [4 1 3 2 5 9 6 8 7 10 ];
    
    i = classOrder(j);
    
    theseInds = find(clustIds==i);
    subsetMatrix = allCSresp.psthZ(theseInds,:);
    allCSresp.classMeans(i,:) = mean(subsetMatrix);
    allCSresp.classErrs(i,:) = std(subsetMatrix,0,1) ./ sqrt(numel(theseInds)-1);
    
%     plot(allCSresp.classMeans(i,500:2800)+(j.*2),'Color',[1-(j./numClasses) 0 (j./numClasses)],'LineWidth',2); hold on;
%     plot(allCSresp.classMeans(i,500:2800)+allCSresp.classErrs(i,500:2800)+(j.*2),'-','Color',[1-(j./numClasses) 0 (j./numClasses)]);
%     plot(allCSresp.classMeans(i,500:2800)-allCSresp.classErrs(i,500:2800)+(j.*2),'-','Color',[1-(j./numClasses) 0 (j./numClasses)]);
%     
    shadedErrorBar(-1000:1.1e4,allCSresp.classMeans(i,:)+(j.*1.5),allCSresp.classErrs(i,:),{'color',[1-(j./numClasses) (j./numClasses).*0.67 (j./numClasses)]}); hold on;
    
    plot([-1e3 1e4],[j.*1.5 j.*1.5],'k--'); axis tight; axis off;
    text(-1500,(j.*1.5)+0.5,[num2str(i) '..' sprintf('%g',numel(theseInds))]);
    
end

%% tmp
disp(' '); clear forDisplay;

% figure(202); clf;
% 
% for i=2:15
%     [clustIds , c , sumdist] = kmeans(allCSresp.psthZ(:,2000:4500),i,'Start','cluster');
%     disp(['Assumed number: ' num2str(i) '... average distance: ' num2str(mean(sumdist).*i)]);
%     forDisplay.x(i-1) = i;
%     forDisplay.y(i-1) = mean(sumdist).*i;
% end

% figure(202); plot(forDisplay.x,forDisplay.y,'ko-'); hold on;

Y = 3:0.25:5;

for i = 1:numel(Y);

   for k = 2:20; %define possible sizes of k

       Color = {'m','b','g','r','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.


       %Kmeans requires rows to contain observations (so cells), columns to contain variables

       FullChar_2 = allCSresp.psthZ(:,2000:4500);

       [idx,C,sumd,D] = kmeans(FullChar_2,k,'replicates',3,'distance','cityblock');


       %find distance from closest cluster centroid

       for j = 1:size(FullChar_2);

           ClostestToCentroid(j,1) = D(j,idx(j,:));

       end

       %compute the distortion - i.e. the normalised square distance between each
       %observation and its clostest cluster center.

       normD = ClostestToCentroid/max(ClostestToCentroid);

       meanDistortion = mean(normD.^2);

       valk(k,1) = meanDistortion;

       YmeanDistortion(i,k) = meanDistortion.^Y(1,i);

   end

end



figure1 = figure('Color',[1 1 1]);
set(gca, 'box','off');


for i = 1:numel(Y);

   [peak,idx] = max(YmeanDistortion(i,:));

   subplot(3,3,i);
   plot(YmeanDistortion(i,:));
   hold on
   plot(idx,peak,'ro');

   xlabel ('Number of Clusters')
   ylabel ('Jump');

   str = sprintf('Y = %f',Y(1,i));

   title(str);
   hold on

   peakKVal(i,1) = idx;

end

%% CHECK FIRST AND LAST LICK DETECTION

disp(['___________________________________________________'])
disp(['STARTED aligning behavior and updating the PopData structure...']);

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
validBehavSess = [];
k=1;
hahnloser = 0; 
numTrialsSpk = 40;
sortedByLat = 1;

% for i=1:NumSessions
% for aa = 1:numel(sessList)
%     i = sessList(aa);
for i=182 % specify the list of sessions to analyze  
    
    if strcmp(PopData.session(i).sessClass,'learning')
        disp(['Session ' num2str(i)]);

        if isfield(PopData.session(i),'behavior')
            if isfield(PopData.session(i).behavior,'found')
                % plot raster of licking data with fl overlaid
                numTrials = size(PopData.session(i).behavior.raster.trial,2);
                
                h=figure(40); clf; 
                set(h,'Color',[1 1 1]);
                disp('Finding Lick Window...');
                for k = 1:numTrials
                    
                    theseTimeStamps = PopData.session(i).behavior.raster.trial(k).ts;
                    firstLickPostCS = find(theseTimeStamps>100,1);
                    slowLicks       = find(diff(theseTimeStamps)>1500);
                    timeOfUS        = find(PopData.session(i).events.US.ts>PopData.session(i).events.CS.ts(k),1);
                        
                    if numel(timeOfUS)>0
                        if (PopData.session(i).events.US.ts(timeOfUS)-PopData.session(i).events.CS.ts(k))<9000
                            PopData.session(i).behavior.USon(1,k)   = PopData.session(i).events.US.ts(timeOfUS)-PopData.session(i).events.CS.ts(k);
                        else
                            disp('No correct US found.');
                            if numel(PopData.session(i).events.US.ts)==1
                                disp('Asserting 2.5 sec');
                                PopData.session(i).behavior.USon(1,k)   = 2500;
                            else
                                disp('Asserting 2 sec');
                                PopData.session(i).behavior.USon(1,k)   = 2000;
                            end
                        end
                    else
                        disp('No US found at all.');
                        if numel(PopData.session(i).events.US.ts)==1
                            disp('Asserting 2.5 sec');
                            PopData.session(i).behavior.USon(1,k)   = 2500;
                        else
                            disp('Asserting 2 sec');
                            PopData.session(i).behavior.USon(1,k)   = 2000;
                        end
                    end
                    
                    if numel(firstLickPostCS)>0
                        disp('Found a postCS lick.');
                        PopData.session(i).behavior.flCSon(1,k) = theseTimeStamps(firstLickPostCS);
                    else
                        disp('Found no postCS lick.');
                        PopData.session(i).behavior.flCSon(1,k) =  9999;
                        PopData.session(i).behavior.flCSoff(1,k) = 10000;
                    end

                    if numel(theseTimeStamps)==0
                        % consider only slow licks that happen after the start of licking
                        disp('No licking. Use default licking estimates.');
                        PopData.session(i).behavior.flCSon(1,k) =  9999;
                        PopData.session(i).behavior.flCSoff(1,k) = 10000;
                    else
                        if numel(firstLickPostCS)>0
                            if numel(slowLicks)>0
                                tmp = find(slowLicks>firstLickPostCS & theseTimeStamps(slowLicks)>PopData.session(i).behavior.USon(1,k),1);
                                if numel(tmp)>0
                                    disp('A qualifying slow lick was found.')
                                    PopData.session(i).behavior.flCSoff(1,k)    = theseTimeStamps(slowLicks(tmp));
                                else
                                    disp('All slow licks came early.')
                                    PopData.session(i).behavior.flCSoff(1,k)    = theseTimeStamps(numel(theseTimeStamps));
                                end
                            else
                                disp('Found no slow lick before end of trial epoch. Using last lick.')
                                if numel(theseTimeStamps)>0
                                    PopData.session(i).behavior.flCSoff(1,k)= theseTimeStamps(numel(theseTimeStamps));
                                else
                                    PopData.session(i).behavior.flCSoff(1,k) = 10000;
                                end
                            end
                        end
                    end

                    if PopData.session(i).behavior.flCSoff(1,k)<PopData.session(i).behavior.USon(1,k)
                        disp('QUALITY CHECK: Last lick looks to be too early. Is there a US?');                        
                        PopData.session(i).behavior.USon(1,k)
                    end

                end

                
                disp('Sorting...');
                if sortedByLat==1
                    [vals,trialInds] = sort(PopData.session(i).behavior.flCSon);
                else
                    trialInds = 1:numTrials;
                end
                
                disp('Plotting routine for behavior...');
                subplot(5,1,1); hold on;
                for k = 1:numTrials
                    thisTrialRow    = trialInds(k);
                    theseTimeStamps = PopData.session(i).behavior.raster.trial(thisTrialRow).ts;
                    validStamps = 1:numel(PopData.session(i).behavior.raster.trial(thisTrialRow).ts);
                    plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*k,'k.','MarkerSize',1);%,'.','Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
%                     plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*k,'b.');%,'.','Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
                end
                
                subplot(5,1,1); 
                title(PopData.session(i).sessId);
                set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
                
                plot(PopData.session(i).behavior.flCSon(trialInds),1:numTrials,'ko',PopData.session(i).behavior.USon(trialInds),1:numTrials,'r-',PopData.session(i).behavior.flCSoff(trialInds),1:numTrials,'bo',[0 0], [-1 numTrialsSpk+1],'b-','LineWidth',1);
                axis([-1000 10000 -1 numTrialsSpk+1]);
                ylabel('Sorted Trials');
                title(i);
                
                if hahnloser
                    disp('Plotting routine for physiology...');
                    disp(' ');
                    numUnits = size(PopData.session(i).unit,2);
                    for m=1:numUnits
                        %thisUnit = unitsByLat(m);
                        thisUnit = m;

                        for k = 1:numTrialsSpk
                            thisTrialRow    = ((m-1).*numTrialsSpk) + (5.*m) + k;
                            theseTimeStamps = PopData.session(i).unit(thisUnit).respCS.raster.trial(trialInds(k)).ts;
                            validStamps     = find(theseTimeStamps>-1000);
                                                        
                            if rem(m,2)
                                subplot(5,1,2:5); hold on;
                                plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
                            else
                                subplot(5,1,2:5);  hold on;
                                plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);            
                            end
                            
                            plot(PopData.session(i).behavior.flCSon(trialInds(1:numTrialsSpk)),((m-1).*numTrialsSpk) + (5.*m) + [1:numTrialsSpk],'k');
                            plot(PopData.session(i).behavior.USon(trialInds(1:numTrialsSpk)),((m-1).*numTrialsSpk) + (5.*m) + [1:numTrialsSpk],'r');
                            plot(PopData.session(i).behavior.flCSoff(trialInds(1:numTrialsSpk)),((m-1).*numTrialsSpk) + (5.*m) + [1:numTrialsSpk],'k.','MarkerSize',3);

                        end
                    end
                end

                subplot(5,1,2:5);
                set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
                plot([0 0], [-1 ((numUnits-1).*numTrialsSpk) + (10.*(numUnits+1)) + k],'b-');
                axis([-1000 10000 -1 ((numUnits).*numTrialsSpk) + (5.*(numUnits+1))]);
                xlabel('Time (ms)');
                ylabel('Sorted Trials per Unit');
                drawnow; %pause();
            end           
        else
            disp('No good behavioral data... skipping');
        end
    end
    
%     pause();
end

%% GET THE DISTRIBUTION OF FIRST LICK LATENCIES
allFL = [];
allLL = [];
clear allFLhist allLLhist
figure(206); clf;

% for i=1:NumSessions
for aa = 1:numel(sessList)
% for i=83 % specify the list of sessions to analyze
  
        i = sessList(aa);
        currFL = PopData.session(i).behavior.flCSon;
        PopDataSub.session(aa).firstLick = currFL;
        currLL = PopData.session(i).behavior.flCSoff;
        PopDataSub.session(aa).lastLick = currLL;
        
        allFL = [allFL currFL];
        allLL = [allLL currLL - currFL];
        
end

bins = 10.^[2:0.025:4];
allFLhist = hist(allFL,bins);
allLLhist = hist(allLL,bins);

allLLhist(1) = 0;

semilogx(bins(1:80),allFLhist(1:80),'k');hold on;
semilogx(bins(1:80),allLLhist(1:80),'b');

%% FIND ALL LARGE SUSTAINED MODULATIONS
winDEL = 2151:2501;
winSUS = 2201:5001;

threshold = 3

respQuant.del.sustain       = zeros(size(allCSaligned.psthZ,2),1);
respQuant.postcs.sustain    = zeros(size(allCSaligned.psthZ,2),1);
% respQuant.perilick.phasic   = zeros(size(allCSaligned.psthZ,2),1);


for i = 1:size(allCSaligned.psthZ,2)
    
    delData = allCSaligned.psthZ(winDEL,i);
    susData = allCSaligned.psthZ(winSUS,i);

    csCheck = find(abs(susData)>threshold);
    
    if numel(csCheck)>50 % ~the s.d. of the smoothing kernel 
        respQuant.del.sustain(i)       = (mean(delData));
        respQuant.postcs.sustain(i)    = (mean(susData));
    else
        respQuant.del.sustain(i)       = 0;
        respQuant.postcs.sustain(i)    = 0;
    end
    
end
figure(41)
subplot(221);
plot(1:size(allCSaligned.psthZ,2),respQuant.del.sustain,'k.',1:size(allCSaligned.psthZ,2),respQuant.postcs.sustain,'ro');
subplot(223);
plot(1:size(allCSaligned.psthZ,2),respQuant.del.sustain-respQuant.postcs.sustain,'ro');
subplot(222);
plot(respQuant.del.sustain,respQuant.postcs.sustain,'k.',[-20 20],[-20 20],'k--',[0 0],[-20 20],'k-',[-20 20],[0 0],'k-');
ylabel('Post CS (500-2500 ms)');xlabel('Delay Period (150-500 ms)');
size(find(abs(respQuant.postcs.sustain)>0))

%% Export Example sustained cells
[values,indices] = sort(respQuant.postcs.sustain,'descend');

figure(42);
subplot(2,2,[1,3])

numInc = 3;
numDec = 662;

exampleSusInc = indices(numInc);
exampleSusDec = indices(numDec);

plot(1:663,values,'k.',numInc,values(numInc),'ro',numDec,values(numDec),'ro');

subplot(2,2,2)
plot(allCSaligned.boxcar(:,exampleSusInc),'k')
% plot(allCSEXaligned.boxcar(:,exampleSusInc),'k')

subplot(2,2,4)
plot(allCSaligned.boxcar(:,exampleSusDec),'k')
% plot(allCSEXaligned.boxcar(:,exampleSusDec),'k')

% export the rasters, psths, and first lick indicators
% iterate through each exemplar
numExamples = 2; m = 2; saveFlag = 0;

for j=1:numExamples

    if j==1
        currIndex = exampleSusInc;
    else
        currIndex = exampleSusDec;
    end

    disp([num2str(j) '...' num2str(m)]);
    
    % grab the session and unit information
    currSess = allCSaligned.sess(currIndex);
    currUnit = allCSaligned.unit(currIndex);
%     currSess = allCSEXaligned.sess(currIndex);
%     currUnit = allCSEXaligned.unit(currIndex);

    % index into the PopData and concatenate the timestamps
    numTrials = size(PopData.session(currSess).unit(currUnit).respCS.raster.trial,2);
    
    % get first lick latency data
    flRaw = PopData.session(currSess).behavior.flCSon;
    % sorted firstLickTimes
    [flSort,flSorti] = sort(flRaw,'ascend');
    offset = find(flSort>501,1);

%         [flSort,flSorti] = sort(1:50,'ascend');
%         offset = 0;
    
    numTrials = 50;
    
    for k=1:numTrials
    
        % sort by the latency to first lick
        currTrial = flSorti(k+offset);
        
        % create a equal sized trial id vector
        thisTrialTS     = PopData.session(currSess).unit(currUnit).respCS.raster.trial(currTrial).ts;
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
    
    figure(1); subplot(4,1,[1:3]);
    plot(flSort(1+offset:numTrials+offset),1:numTrials,'ro',finalArrayTS,finalArrayID,'k.');
    if j==1
        tmpExport(:,1) = flSort(1+offset:numTrials+offset);
        tmpExport(:,2) = 1:numTrials;
    else
        tmpExport(:,3) = flSort(1+offset:numTrials+offset);
        tmpExport(:,4) = 1:numTrials;        
    end
    subplot(4,1,4);
    plot(-2000:3000, allCSaligned.psthZ(:,currIndex),'k');
%     plot(-2000:3000, allCSEXaligned.psthZ(:,currIndex),'k');
    drawnow; pause(1);
    
    if saveFlag
        % save to an hdf5 structure for reading/plotting in igor
        nameA = sprintf('/spikeTimesX%g',j-1);
        nameB = sprintf('/spikeTimesY%g',j-1);
        nameC = sprintf('/psth%g',j-1);
        nameD = sprintf('/psthERR%g',j-1);

        if j==1
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameA,finalArrayTS);
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameB,finalArrayID,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameC,allCSaligned.psthZ(:,currIndex),'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameD,allCSaligned.psthZe(:,currIndex),'WriteMode','append');
        else
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameA,finalArrayTS,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameB,finalArrayID,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameC,allCSaligned.psthZ(:,currIndex),'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameD,allCSaligned.psthZe(:,currIndex),'WriteMode','append');
        end
    end
    
end

%% WARP DATA INTO A CS -> FLICK PHASE

% simply go through and assign a phase to every spike between tone onset and 150 ms after first lick
% for i=1:NumSessions

debugLogic = 0;

for aa = 1:numel(sessList)
    i=sessList(aa);
% for i = 83;
    
    if strcmp(PopData.session(i).sessClass,'learning')
        disp(['Session ' num2str(i)]);

        if isfield(PopData.session(i),'behavior')
            if isfield(PopData.session(i).behavior,'found')

                numUnits    = size(PopData.session(i).unit,2);
                    
                if debugLogic
                    figure(43); clf;
                    figure(44); clf;
                    figure(45); clf;
                    figure(46); clf;
                end
                
                validTrials = find(PopData.session(i).behavior.flCSon>100 & PopData.session(i).behavior.flCSon<2000 & PopData.session(i).behavior.flCSoff>(PopData.session(i).behavior.USon+500) & PopData.session(i).behavior.flCSoff<9500 );
%                 validTrials = find(PopData.session(i).behavior.flCSon>100 & PopData.session(i).behavior.flCSon<1600 & PopData.session(i).behavior.flCSoff>(PopData.session(i).behavior.USon+500));
                numTrials   = numel(validTrials);
                PopData.session(i).behavior.validPhaseTrials = validTrials;

                for m=1:numUnits
                    thisUnit = m

                    for q=1:5
                        allPhase.seg(q).phases = [];
                        allPhase.seg(q).times = [];
                    end
                    allPhase.seg(1).pOff = -0.5;
                    allPhase.seg(2).pOff = 0;
                    allPhase.seg(3).pOff = 0.718;
                    allPhase.seg(4).pOff = 0.718 + 1.077;
                    allPhase.seg(5).pOff = 0.718 + 1.077 + 4.488;
                    
                    allPhase.seg(1).eachDur = [];
                    allPhase.seg(2).eachDur = [];
                    allPhase.seg(3).eachDur = [];
                    allPhase.seg(4).eachDur = [];
                    allPhase.seg(5).eachDur = [];
                    
                    % find appropriate trials where the ts are in sequence
                    % CS on, CS off, FL, US, LL                    

                    
                    % cycle through all valid trials
                    for k = 1:numTrials
                        
                        thisTrial = validTrials(k);
                        theseTimeStamps = PopData.session(i).unit(thisUnit).respCS.raster.trial(thisTrial).ts;
                        
                        size(theseTimeStamps)
                        
                        % warp the baseline data (0.5 radian)
                        validStamps     = find( theseTimeStamps<0 );
                        thisSegDuration = 1000;
                        if numel(validStamps)>0
                            currPhases  = theseTimeStamps(validStamps) .* (0.5 ./ thisSegDuration); % multiply by rads / sec
                            currTimes   = theseTimeStamps(validStamps);
                            if abs(numel(currPhases)-numel(currTimes)) > 0
                                disp('An error has occurred.');
                            end
                            if k==1
                                allPhase.seg(1).phases = currPhases'; % multiply by rads / sec
                                allPhase.seg(1).times = currTimes'; % multiply by rads / sec
                                allPhase.seg(1).pOff = -0.5;
                            else
                                allPhase.seg(1).phases = [allPhase.seg(1).phases currPhases'];
                                allPhase.seg(1).times = [allPhase.seg(1).times currTimes'];
                            end
                            allPhase.seg(1).trial(k).ph = currPhases;
                            allPhase.seg(1).trial(k).ts = currTimes;
                            allPhase.seg(1).trial(k).du = thisSegDuration;
                        else
                            allPhase.seg(1).trial(k).ph = [];
                            allPhase.seg(1).trial(k).ts = [];                            
                            allPhase.seg(1).trial(k).du = thisSegDuration;
                        end                            
                        allPhase.seg(1).eachDur = [allPhase.seg(1).eachDur thisSegDuration];

                        
                        % warp the CS onset -> FL interval (0.718 rads)
                        validStamps     = find( theseTimeStamps>0 & theseTimeStamps<PopData.session(i).behavior.flCSon(1,thisTrial) );
                        thisSegDuration = PopData.session(i).behavior.flCSon(1,thisTrial);                        
                        if numel(validStamps)>0
                            currPhases = theseTimeStamps(validStamps) .* (0.718 ./ thisSegDuration); % multiply by rads / sec
                            currTimes = theseTimeStamps(validStamps); % multiply by rads / sec
                            if k==1
                                allPhase.seg(2).phases = currPhases'; % multiply by rads / sec
                                allPhase.seg(2).times = currTimes'; % multiply by rads / sec
                                allPhase.seg(2).pOff = 0;
                            else
                                allPhase.seg(2).phases = [allPhase.seg(2).phases currPhases'];
                                allPhase.seg(2).times = [allPhase.seg(2).times currTimes'];
                            end
                            allPhase.seg(2).trial(k).ph = currPhases;
                            allPhase.seg(2).trial(k).ts = currTimes;
                            allPhase.seg(2).trial(k).du = thisSegDuration;
                        else
                            allPhase.seg(2).trial(k).ph = [];
                            allPhase.seg(2).trial(k).ts = [];                            
                            allPhase.seg(2).trial(k).du = thisSegDuration;
                        end
                        allPhase.seg(2).eachDur = [allPhase.seg(1).eachDur thisSegDuration];
                        
                        % warp the FL -> US (1.077 rads)
                        validStamps     = find( theseTimeStamps>PopData.session(i).behavior.flCSon(1,thisTrial) & theseTimeStamps<PopData.session(i).behavior.USon(1,thisTrial) );
                        thisSegDuration = PopData.session(i).behavior.USon(1,thisTrial) - PopData.session(i).behavior.flCSon(1,thisTrial);                        
                        if numel(validStamps)>0
                            currPhases = (theseTimeStamps(validStamps)- PopData.session(i).behavior.flCSon(1,thisTrial)) .* (1.077 ./ thisSegDuration); % multiply by rads / sec
                            currTimes = (theseTimeStamps(validStamps)- PopData.session(i).behavior.flCSon(1,thisTrial)); % multiply by rads / sec
                            if k==1
                                allPhase.seg(3).phases = currPhases'; % multiply by rads / sec
                                allPhase.seg(3).times = currTimes'; % multiply by rads / sec
                                allPhase.seg(3).pOff = 0.718;
                            else
                                allPhase.seg(3).phases = [allPhase.seg(3).phases currPhases'];
                                allPhase.seg(3).times = [allPhase.seg(3).times currTimes'];
                            end
                            allPhase.seg(3).trial(k).ph = currPhases+allPhase.seg(3).pOff;
                            allPhase.seg(3).trial(k).ts = currTimes + PopData.session(i).behavior.flCSon(1,thisTrial);
                            allPhase.seg(3).trial(k).du = thisSegDuration;
                        else
                            allPhase.seg(3).trial(k).ph = [];
                            allPhase.seg(3).trial(k).ts = [];                            
                            allPhase.seg(3).trial(k).du = thisSegDuration;
                        end
                        allPhase.seg(3).eachDur = [allPhase.seg(1).eachDur thisSegDuration];

                        % warp the US -> LL interval (4.488 rads)
                        validStamps     = find( theseTimeStamps>PopData.session(i).behavior.USon(1,thisTrial) & theseTimeStamps<PopData.session(i).behavior.flCSoff(1,thisTrial) );
                        thisSegDuration = PopData.session(i).behavior.flCSoff(1,thisTrial)-PopData.session(i).behavior.USon(1,thisTrial);                        
                        if numel(validStamps)>0
                            currPhases = (theseTimeStamps(validStamps)-PopData.session(i).behavior.USon(1,thisTrial)) .* (4.488 ./ thisSegDuration); % multiply by rads / sec
                            currTimes = (theseTimeStamps(validStamps)-PopData.session(i).behavior.USon(1,thisTrial)); % multiply by rads / sec
                            if k==1
                                allPhase.seg(4).phases = currPhases'; % multiply by rads / sec
                                allPhase.seg(4).times = currTimes'; % multiply by rads / sec
                                allPhase.seg(4).pOff = 0.718 + 1.077;
                            else
                                allPhase.seg(4).phases = [allPhase.seg(4).phases currPhases'];
                                allPhase.seg(4).times = [allPhase.seg(4).times currTimes'];
                            end
                            allPhase.seg(4).trial(k).ph = currPhases+allPhase.seg(4).pOff;
                            allPhase.seg(4).trial(k).ts = currTimes + PopData.session(i).behavior.USon(1,thisTrial);
                            allPhase.seg(4).trial(k).du = thisSegDuration;
                        else
                            allPhase.seg(4).trial(k).ph = [];
                            allPhase.seg(4).trial(k).ts = [];                            
                            allPhase.seg(4).trial(k).du = thisSegDuration;
                        end
                        allPhase.seg(4).eachDur = [allPhase.seg(1).eachDur thisSegDuration];

                        % warp the LL -> POST interval (1 radian)
                        validStamps     = find( theseTimeStamps>PopData.session(i).behavior.flCSoff(1,thisTrial) & theseTimeStamps<PopData.session(i).behavior.flCSoff(1,thisTrial)+1000 );
                        thisSegDuration = 1000;
                        
                        if numel(validStamps)>0
                            currPhases = (theseTimeStamps(validStamps)-PopData.session(i).behavior.flCSoff(1,thisTrial)) .* (1 ./ thisSegDuration); % multiply by rads / sec
                            currTimes = (theseTimeStamps(validStamps)-PopData.session(i).behavior.flCSoff(1,thisTrial)); % multiply by rads / sec
                            if k==1
                                allPhase.seg(5).phases = currPhases'; % multiply by rads / sec
                                allPhase.seg(5).times = currTimes'; % multiply by rads / sec
                                allPhase.seg(5).pOff = 0.718 + 1.077 + 4.488;
                            else
                                allPhase.seg(5).phases = [allPhase.seg(5).phases currPhases'];
                                allPhase.seg(5).times = [allPhase.seg(5).times currTimes'];
                            end
                            allPhase.seg(5).trial(k).ph = currPhases+allPhase.seg(5).pOff;
                            allPhase.seg(5).trial(k).ts = currTimes + PopData.session(i).behavior.flCSoff(1,thisTrial);
                            allPhase.seg(5).trial(k).du = thisSegDuration;
                        else
                            allPhase.seg(5).trial(k).ph = [];
                            allPhase.seg(5).trial(k).ts = [];                            
                            allPhase.seg(5).trial(k).du = thisSegDuration;
                        end
                        allPhase.seg(5).eachDur = [allPhase.seg(1).eachDur thisSegDuration];
                        
                    end
                    
                    for ppp=1:5
                        allPhase.seg(ppp).meanDur = mean(allPhase.seg(ppp).eachDur);
                    end
                    
                    PopData.session(i).unit(thisUnit).allPhase = allPhase;
                    
                    if debugLogic

                        figure(43);
                        subplot(2,numUnits,numUnits-m+1);
                        hist(allPhase.seg(1).phases,25);
                        axis tight;
                        subplot(2,numUnits,(numUnits-m+1+numUnits));
                        hist(allPhase.seg(1).times,25);
                        axis tight;
                        disp(['Num stamps in phases: ' num2str(numel(allPhase.seg(1).phases)) ' and in times: ' num2str(numel(allPhase.seg(1).times)) ' ...']);

                        figure(44);
                        subplot(2,numUnits,numUnits-m+1);
                        hist(allPhase.seg(2).phases,25);
                        axis tight;
                        subplot(2,numUnits,(numUnits-m+1+numUnits));
                        hist(allPhase.seg(2).times,25);
                        axis tight;
                        disp(['Num stamps in phases: ' num2str(numel(allPhase.seg(2).phases)) ' and in times: ' num2str(numel(allPhase.seg(2).times)) ' ...']);

                        figure(45);
                        subplot(2,numUnits,numUnits-m+1);
                        hist(allPhase.seg(3).phases,25);
                        axis tight;
                        subplot(2,numUnits,(numUnits-m+1+numUnits));
                        hist(allPhase.seg(3).times,25);
                        axis tight;

                        figure(46);
                        subplot(2,numUnits,numUnits-m+1);
                        hist(allPhase.seg(4).phases,25);
                        axis tight;
                        subplot(2,numUnits,(numUnits-m+1+numUnits));
                        hist(allPhase.seg(4).times,25);
                        axis tight;

                        figure(47);
                        subplot(2,numUnits,numUnits-m+1);
                        hist(allPhase.seg(5).phases,25);
                        axis tight;
                        subplot(2,numUnits,(numUnits-m+1+numUnits));
                        hist(allPhase.seg(5).times,25);
                        axis tight;

                        drawnow;
                    
                    end
                    
                end 
            end
        else
            disp('No good behavioral data... skipping');
        end
    
        drawnow;
    end
end

%% WARP IN TO PRIOR REWARD -> CS SPACE

figure(410); clf;
figure(411); clf;

% for aa = 1:numel(sessList)
for i=181 % specify the list of sessions to analyze

%     i = sessList(aa);
    disp(['Session ' num2str(i)]);
            
    numUnits = size(PopData.session(i).unit,2);
    [colors] = TNC_CreateRBColormap(numUnits,'bo');

    PopData.session(i).trials = size(PopData.session(i).events.CS.ts,1);
          
    for j = 1:numUnits

        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,round(PopData.session(i).unit(j).ts)+1) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');
        
        [respPCS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[8e4,2e3],1,1);
        PopData.session(i).unit(j).respPCS.raster      = respPCS.raster;
        PopData.session(i).unit(j).respPCS.psthAVG     = respPCS.image.psthAVG;
        PopData.session(i).unit(j).respPCS.psthSEM     = respPCS.image.psthSEM;

        k=k+1;

    end

    disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' of ' num2str(NumSessions) ' | ' PopData.session(i).sessId '_' PopData.session(i).sessClass])  

    for pp=1:PopData.session(i).trials
        if pp==1
            PopData.session(i).timeSinceCS(pp) = 8e4;
        else
            PopData.session(i).timeSinceCS(pp) = PopData.session(i).events.CS.ts(pp) - PopData.session(i).events.CS.ts(pp-1);
        end
    end
    
    [tN,tI] = sort(PopData.session(i).timeSinceCS,'descend');
    avgTsCS = mean(PopData.session(i).timeSinceCS)
    
    for j = 1:numUnits
        
        allStamps   = [];
        allRows     = [];
        allCSpX     = [];
        allCSpY     = [];
        deltaP      = zeros(1,6280+ceil(2e6./avgTsCS));
        image       = zeros(PopData.session(i).trials , 6280+ceil(2e6./avgTsCS));
        
        for kk=1:PopData.session(i).trials

            deltaP      = zeros(1,6280+ceil(2e6./avgTsCS));
            currTrial   = tI(kk);        
            theseStamps = PopData.session(i).unit(j).respPCS.raster.trial(currTrial).ts';
            theseRows   = ( kk + ((j-1).*(PopData.session(i).trials + 10))) .* ones(1,numel(theseStamps));

            allStamps   = [allStamps theseStamps];
            allRows     = [allRows theseRows];
            allCSpX     = [allCSpX -(PopData.session(i).timeSinceCS(currTrial))];
            allCSpY     = [allCSpY , kk + ((j-1).*(PopData.session(i).trials + 10))];
            
            % calculate warped stamps
            validStamps = find(theseStamps>-PopData.session(i).timeSinceCS(currTrial) & theseStamps<1);
            validStamps2= find(theseStamps>0);
            PopData.session(i).unit(j).respPCS.raster.trial(currTrial).ph = [(theseStamps(validStamps) .* (6280 ./ PopData.session(i).timeSinceCS(currTrial))) theseStamps(validStamps2).*(1./2e3).*ceil(2e6./avgTsCS)];
            
            deltaP(ceil(theseStamps(validStamps) .* (6280 ./ PopData.session(i).timeSinceCS(currTrial))) + 6280) = (avgTsCS / PopData.session(i).timeSinceCS(currTrial));
            deltaP( ceil(theseStamps(validStamps2).*(1./2e3).*ceil(2e6./avgTsCS)) + 6280 ) = 1;
            
            image(kk,:) =  conv(deltaP,currParams.filter.kernel,'same');
            
        end
        
        PopData.session(i).unit(j).respPCS.WpsthAVG = mean(image,1);
        PopData.session(i).unit(j).respPCS.WpsthSEM = std(image,[],1) ./ sqrt(PopData.session(i).trials-1);
        
        figure(410);
        plot(allStamps,allRows,'.','Color',colors(j,:),'MarkerSize',1); hold on; 
        plot(allCSpX,allCSpY,'k-');
        plot([0 0] , [((j-1).*(PopData.session(i).trials + 10)) PopData.session(i).trials + ((j-1).*(PopData.session(i).trials + 10))] , 'k');
        axis([-8e4 2e3 0 (numUnits-1).*(PopData.session(i).trials + 10)]);

        postCSpnts  = ceil(2e6./avgTsCS);
        
        figure(411); subplot(121);
%         plot(-6279:61,PopData.session(i).unit(j).respPCS.WpsthAVG + (j*0.02),'Color',colors(j,:)); hold on;
        shadedErrorBar(-6279:postCSpnts,PopData.session(i).unit(j).respPCS.WpsthAVG - mean(PopData.session(i).unit(j).respPCS.WpsthAVG) + (j*0.06),PopData.session(i).unit(j).respPCS.WpsthSEM,{'Color',colors(j,:)}); hold on;
        axis([-6280 62 0 (numUnits+1)*0.06]);
         drawnow;  
         if j== numUnits
             plot([0 0], [0 (j+1)*0.06],'k--');
             plot([-5975 -5975], [0 (j+1)*0.06],'k--','Color',[0.25 0.25 0.25]);
         end
         
         subplot(122);
         plot(-8e4:2e3,PopData.session(i).unit(j).respPCS.psthAVG - mean(PopData.session(i).unit(j).respPCS.psthAVG) + (j*0.025),'Color',colors(j,:)); hold on;
        axis([-8e4 2e3 0 (numUnits+1)*0.025]);
         drawnow;  
         if j== numUnits
             plot([0 0], [0 (j+1)*0.025],'k--');
         end
    end
end

%% WARP INTO A SEQUENCE OF PRIOR LL >> CS >> FL >> LL

%% WARP INTO CS >> FL for analysis with Parvez

figure(610); clf;
figure(611); clf;
count = 1;

PopData.ampOfPeak = [];
PopData.timeOfPeak = [];
NumSessions = numel(sessList);
for aa = 1:numel(sessList)
    i = sessList(aa);

%     for i = [103] % specify the list of sessions to analyze

    disp(['Session ' num2str(i)]);
    figure(610); clf;
    figure(611); clf;

    numUnits = size(PopData.session(i).unit,2);
    [colors] = TNC_CreateRBColormap(numUnits,'bo');

    PopData.session(i).trials = size(PopData.session(i).events.CS.ts,1);
          
    for j = 1:numUnits

        numStamps = length(PopData.session(i).unit(j).ts);
        delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,round(PopData.session(i).unit(j).ts)+1) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');
        
        [respJCS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[1e3,6e3],1,1);
        PopData.session(i).unit(j).respJCS.raster      = respJCS.raster;
        PopData.session(i).unit(j).respJCS.psthAVG     = respJCS.image.psthAVG;
        PopData.session(i).unit(j).respJCS.psthSEM     = respJCS.image.psthSEM;
        PopData.session(i).unit(j).respJCS.image       = respJCS.image.aligned;

        k=k+1;

    end

    disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' of ' num2str(NumSessions) ' | ' PopData.session(i).sessId '_' PopData.session(i).sessClass])  

    for pp=1:PopData.session(i).trials
        if pp==1
            PopData.session(i).timeSinceCS(pp) = 8e4;
        else
            PopData.session(i).timeSinceCS(pp) = PopData.session(i).events.CS.ts(pp) - PopData.session(i).events.CS.ts(pp-1);
        end
    end
    
    [tN,tI] = sort(PopData.session(i).behavior.flCSon,'ascend');
    tmp = find(tN<5e3);
    avgCS2FL = mean(tN(tmp));
    validTrials = tI(tmp);
    
    phaseDecay = floor( currParams.smthParams.decay .* (6280/avgCS2FL) ./ 3 );
    [phasekernel]  = TNC_CreateGaussian(phaseDecay.*15,phaseDecay,phaseDecay.*30,1);
    
    clear lat
    
    for j = 1:numUnits
        
        allStamps   = [];
        allRows     = [];
        allCSpX     = [];
        allCSpY     = [];
        lat.flOn    = [];
        lat.flOff   = [];
        lat.preCS   = [];
        deltaP      = zeros(1,6280+2e3);
        imagePh     = zeros(numel(validTrials) , 6280+2e3);
        imageTm     = zeros(numel(validTrials) , numel(PopData.session(i).unit(j).respJCS.image(1,:)));
        
        for kk=1:numel(validTrials)

            deltaP      = zeros(1,6280+2e3);
            currTrial   = validTrials(kk);        
            theseStamps = PopData.session(i).unit(j).respJCS.raster.trial(currTrial).ts';
            theseRows   = ( kk + ((j-1).*(PopData.session(i).trials + 10))) .* ones(1,numel(theseStamps));

            allStamps   = [ allStamps theseStamps ];
            allRows     = [ allRows theseRows ];
            allCSpX     = [ allCSpX PopData.session(i).behavior.flCSon(currTrial) ];
            allCSpY     = [ allCSpY , kk + ((j-1).*(PopData.session(i).trials + 10)) ];
            
            % calculate warped stamps
            validStamps = find(theseStamps>0 & theseStamps<PopData.session(i).behavior.flCSon(currTrial));
            validStamps2= find(theseStamps>PopData.session(i).behavior.flCSon(currTrial) & theseStamps<PopData.session(i).behavior.flCSon(currTrial)+1000);
            validStamps3= find(theseStamps<0);
            
            PopData.session(i).unit(j).respJCS.raster.trial(currTrial).ph = [   theseStamps(validStamps3) , ...                                                                                
                                                                                (theseStamps(validStamps) .* (6280 ./ PopData.session(i).behavior.flCSon(currTrial))) , ...
                                                                                (theseStamps(validStamps2)-PopData.session(i).behavior.flCSon(currTrial)) + 6280];
            
            deltaP( ceil(theseStamps(validStamps3)) + 1e3 ) = 1e3 / 6280;
            deltaP( ceil((theseStamps(validStamps) .* (6280 ./ PopData.session(i).behavior.flCSon(currTrial)))) + 1e3 ) = (avgCS2FL / PopData.session(i).behavior.flCSon(currTrial));
            deltaP( ceil((theseStamps(validStamps2)-PopData.session(i).behavior.flCSon(currTrial)) + 7280) ) = 1e3 / 6280;
            
            imagePh(kk,:)   = conv(deltaP,phasekernel,'same');
            imageTm(kk,:)   = PopData.session(i).unit(j).respJCS.image(currTrial,:);
            lat.flOn(kk)    = PopData.session(i).behavior.flCSon(currTrial);
            lat.flOff(kk)   = PopData.session(i).behavior.flCSoff(currTrial);
            lat.preCS(kk)   = PopData.session(i).timeSinceCS(currTrial);
            
        end
        
        PopData.session(i).unit(j).respJCS.imagePh  = imagePh;
        PopData.session(i).unit(j).respJCS.imageTm  = imageTm;
        PopData.session(i).unit(j).respJCS.lat      = lat;
        
        PopData.session(i).unit(j).respJCS.WpsthAVG = mean(imagePh,1) - mean(mean(imagePh(:,1:1e3),1));
        PopData.session(i).unit(j).respJCS.WpsthSEM = std(imagePh,[],1) ./ sqrt(numel(validTrials)-1);

        PopData.session(i).unit(j).respJCS.TpsthAVG = mean(imageTm,1) - mean(mean(imageTm(:,1:1e3),1));
        PopData.session(i).unit(j).respJCS.TpsthSEM = std(imageTm,[],1) ./ sqrt(numel(validTrials)-1);
        
        if abs(max(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280))) > abs(min(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)))
            PopData.WpsthAVG(count,:) = PopData.session(i).unit(j).respJCS.WpsthAVG ./ abs(max(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)));
            PopData.WpsthSEM(count,:) = PopData.session(i).unit(j).respJCS.WpsthSEM ./ abs(max(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)));

            PopData.TpsthAVG(count,:) = PopData.session(i).unit(j).respJCS.TpsthAVG ./ abs(max(PopData.session(i).unit(j).respJCS.TpsthAVG(1e3:1e3+avgCS2FL)));
            PopData.TpsthSEM(count,:) = PopData.session(i).unit(j).respJCS.TpsthSEM ./ abs(max(PopData.session(i).unit(j).respJCS.TpsthAVG(1e3:1e3+avgCS2FL)));
            
            PopData.timeOfPeak(count,:) = find( PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)==max(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)) , 1 ) .* 1;

            PopData.ampOfPeak(count,1) = (max(mean(imagePh(:,1e3:7280),1)) - mean(mean(imagePh(:,1:1e3)))) ./ (max(mean(imagePh(:,1e3:7280),1)) + mean(mean(imagePh(:,1:1e3))));
            PopData.ampOfPeak(count,2) = (max(mean(imageTm(:,1e3:3e3),1)) - mean(mean(imageTm(:,1:1e3)))) ./ (max(mean(imageTm(:,1e3:3e3),1)) + mean(mean(imageTm(:,1:1e3))));

        else
            PopData.WpsthAVG(count,:) = PopData.session(i).unit(j).respJCS.WpsthAVG ./ abs(min(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)));
            PopData.WpsthSEM(count,:) = PopData.session(i).unit(j).respJCS.WpsthSEM ./ abs(min(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)));

            PopData.TpsthAVG(count,:) = PopData.session(i).unit(j).respJCS.TpsthAVG ./ abs(min(PopData.session(i).unit(j).respJCS.TpsthAVG(1e3:1e3+avgCS2FL)));
            PopData.TpsthSEM(count,:) = PopData.session(i).unit(j).respJCS.TpsthSEM ./ abs(min(PopData.session(i).unit(j).respJCS.TpsthAVG(1e3:1e3+avgCS2FL)));

            PopData.timeOfPeak(count,:) = find( PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)==min(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)) , 1 ) .* -1;
            
            PopData.ampOfPeak(count,1) = (min(mean(imagePh(:,1e3:7280),1)) - mean(mean(imagePh(:,1:1e3)))) ./ (abs(min(mean(imagePh(:,1e3:7280),1))) + mean(mean(imagePh(:,1:1e3))));
            PopData.ampOfPeak(count,2) = (min(mean(imageTm(:,1e3:3e3),1)) - mean(mean(imageTm(:,1:1e3)))) ./ (abs(min(mean(imageTm(:,1e3:3e3),1))) + mean(mean(imageTm(:,1:1e3))));
        end
        
        count = count+1;
        
        figure(610);
        plot(allStamps,allRows,'.','Color',colors(j,:),'MarkerSize',1); hold on; 
        plot(allCSpX,allCSpY,'k-');
        plot([0 0],[min(allCSpY) max(allCSpY)],'k-');
        plot([0 0] , [((j-1).*(PopData.session(i).trials + 10)) PopData.session(i).trials + ((j-1).*(PopData.session(i).trials + 10))] , 'k');
        axis([-1e3 6e3 0 (numUnits+1).*(PopData.session(i).trials + 10)]);
        
        scaler = 0.008;
        figure(611); subplot(121);
        shadedErrorBar(-999:7280,PopData.session(i).unit(j).respJCS.WpsthAVG + (j*scaler),PopData.session(i).unit(j).respJCS.WpsthSEM,{'Color',colors(j,:)}); hold on;
        plot([-1e3 7280] , [(j*scaler) (j*scaler)] , 'k','Color',[0.75 0.75 0.75]);
        axis([-1e3 7280 0 (numUnits+1)*scaler]);
         drawnow;  
         if j== numUnits
             plot([0 0], [0 (j+1)*scaler],'k--');
             plot([6280 6280], [0 (j+1)*scaler],'k--','Color',[0.25 0.25 0.25]);
         end
         
         subplot(122);
        shadedErrorBar(-1e3:6e3,PopData.session(i).unit(j).respJCS.TpsthAVG + (j*0.025),PopData.session(i).unit(j).respJCS.TpsthSEM,{'Color',colors(j,:)}); hold on;
        axis([-1e3 7280 0 (numUnits+1)*scaler]);
        axis([-1e3 6e3 0 (numUnits+1)*0.025]);
         drawnow;  
         if j== numUnits
             plot([0 0], [0 (j+1)*0.025],'k--');
             plot([avgCS2FL avgCS2FL], [0 (j+1)*0.025],'k--','Color',[0.25 0.25 0.25]);
         end
    end
end

[pmMap] = TNC_CreateRBColormap(1024,'cb');
[tpV,tpI] = sort(PopData.timeOfPeak);
tpNeg = find(PopData.timeOfPeak<=0);
tpPos = sort(PopData.timeOfPeak>0);

figure(612); clf;
figure(613); clf;

figure(612);
subplot(121);
imagesc(PopData.WpsthAVG(tpI,:),[-1.5 1.5]); colormap(pmMap);
hold on; plot([1e3 1e3], [0 count],'k--'); plot([7280 7280], [0 count],'k--');% plot(abs(PopData.timeOfPeak(tpI))+1e3, 1:count-1 , 'k-');
set(gca,'TickDir','out');

subplot(122);
imagesc(PopData.TpsthAVG(tpI,:),[-1.5 1.5]); colormap(pmMap);
axis([250 3.5e3 0 count]);
set(gca,'TickDir','out');
    
figure(613); subplot(211);
    shadedErrorBar(-999:7280,mean(PopData.WpsthAVG(tpPos,:),1),std(PopData.WpsthAVG(tpPos,:),[],1)./sqrt(count-1),{'Color',pmMap(1024,:)}); hold on;
    shadedErrorBar(-999:7280,mean(PopData.WpsthAVG(tpNeg,:),1),std(PopData.WpsthAVG(tpNeg,:),[],1)./sqrt(count-1),{'Color',pmMap(1,:)});
    plot([-900 7200] , [0 0], 'k--');
    axis([-750 7030 -0.75 0.75]);
subplot(212);
    shadedErrorBar(-1e3:6e3,mean(PopData.TpsthAVG(tpPos,:),1),std(PopData.TpsthAVG(tpPos,:),[],1)./sqrt(count-1),{'Color',pmMap(1024,:)}); hold on;
    shadedErrorBar(-1e3:6e3,mean(PopData.TpsthAVG(tpNeg,:),1),std(PopData.TpsthAVG(tpNeg,:),[],1)./sqrt(count-1),{'Color',pmMap(1,:)});
    plot([-900 3e3] , [0 0], 'k--');
    axis([-750 2.5e3 -1.5 0.75]);
    
    figure(614);
    PopData.timeOfPeakH = hist(abs(PopData.timeOfPeak),0:250:6280);
    PopData.timeOfPeakHx = [0:250:6280]./1000;
    bar(PopData.timeOfPeakHx,PopData.timeOfPeakH,'k');
    
%     [modMap] = TNC_CreateRBColormap(1024,'bo');    
clear modMap
    modMap(:,1) = 0:0.001:1;
    modMap(:,2) = 0.67:-0.00067:0;
    modMap(:,3) = 1:-0.001:0;
    figure(616); clf;
    scatter(PopData.ampOfPeak(:,1),PopData.ampOfPeak(:,2),10,abs(PopData.timeOfPeak)); colormap(modMap); hold on;
    plot([-1 1] , [-1 1] , 'k--'); plot([0 0] , [-1 1] , 'k-'); plot([-1 1] , [0 0] , 'k-');
    xlabel('Modulation Index (phase)');
    ylabel('Modulation Index (time)');
    axis([-1.1 1.1 -1.1 1.1]);
    set(gca,'TickDir','out');

%% CREATE A HAHNLOSER TYPE RASTER PLOT FOR TIMES AND FOR PHASE
% Toggles between two styles. One is grouped by trial and one is grouped by unit

export =0;

disp(['___________________________________________________'])
disp(['PLOTTING ROUTINE FOR PHASE AND TIME COMPARISONS...']);

countUnits  = 0;
NumSessions = size(PopData.session,2);
validBehavSess = [];
k=1;

% grouping = 1; % groups by unit
grouping = 1; % groups by trial

sortedByLat = 1; % sort the valid trials in order or response latency
% sortedByLat = 0;

% for i=1:NumSessions
for aa = 1:numel(sessList)
% for i=83 % specify the list of sessions to analyze
% for i=109 % specify the list of sessions to analyze
        i = sessList(aa);
        disp(['Session ' num2str(i)]);

        validTrials = PopData.session(i).behavior.validPhaseTrials;
        numTrials = numel(validTrials);
        
        if numTrials==0
            
            disp('Did not find any valid trials');

        else
            
            disp(['Found ' num2str(numTrials) ' valid trials']);
        
            h1=figure(40); clf; 
            set(h1,'Color',[1 1 1]);
            h2=figure(41); clf; 
            set(h2,'Color',[1 1 1]);

            if sortedByLat==1
                disp('Sorting by latency...');
                yLabelStr = '(sorted)';
                latencies = PopData.session(i).behavior.flCSon(validTrials);
                [vals,trialIndsTmp] = sort(latencies);
                trialInds = validTrials(trialIndsTmp);
            else
                disp('No sorting was applied...');
                trialInds = validTrials;
                yLabelStr = '(not sorted)';
            end

            disp('Plotting routine for behavior in time...');
                figure(40);
                subplot(5,1,1); hold on;
                for k = 1:numTrials
                    thisTrialRow    = trialInds(k);
                    theseTimeStamps = PopData.session(i).behavior.raster.trial(thisTrialRow).ts;
                    validStamps = 1:numel(PopData.session(i).behavior.raster.trial(thisTrialRow).ts);
                    plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*k,'k.','MarkerSize',1);%,'.','Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
                end

                subplot(5,1,1); title(['Session ' num2str(i) ' plotted in real time']);
                set(gca,'TickDir','out','TickLength',[0.005 0]); box off;

                plot(PopData.session(i).behavior.flCSon(trialInds),1:numTrials,'ko',PopData.session(i).behavior.USon(trialInds),1:numTrials,'r-',PopData.session(i).behavior.flCSoff(trialInds),1:numTrials,'ko',[0 0], [-1 numTrials+1],'b-','LineWidth',1);
                axis([-1000 10000 -1 numTrials+1]);
                ylabel(['Trials' yLabelStr]);

            disp('Plotting routine for behavior in phase...');
                figure(41);
                subplot(5,1,1); hold on;
                for k = 1:numTrials
                    thisTrialRow    = trialInds(k);
                    theseTimeStamps = PopData.session(i).behavior.raster.trial(thisTrialRow).ts;
                    validStamps = 1:numel(PopData.session(i).behavior.raster.trial(thisTrialRow).ts);
                    plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*k,'k.','MarkerSize',1);%,'.','Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
                end

                subplot(5,1,1); title(['Session ' num2str(i) ' plotted in a warped phase space']);
                set(gca,'TickDir','out','TickLength',[0.005 0]); box off;

                plot(PopData.session(i).behavior.flCSon(trialInds),1:numTrials,'ko',PopData.session(i).behavior.USon(trialInds),1:numTrials,'r-',PopData.session(i).behavior.flCSoff(trialInds),1:numTrials,'ko',[0 0], [-1 numTrials+1],'b-','LineWidth',1);
                axis([-1000 10000 -1 numTrials+1]);
                ylabel(['Trials' yLabelStr]);

            switch grouping
                
                case 1
                    disp('Grouping by units...');
                    disp(' ');
                    numUnits = size(PopData.session(i).unit,2);

                    % __________________________
                    % __________________________
                    % FOR TIMES
                    % __________________________
                    % __________________________
                    
                    figure(40); hold on;
                    
                    allStamps = []; allRows = []; allIDs = [];
                    
                    for m=1:numUnits
                        
                        delta = zeros(1,1.2e4);

                        thisUnit = m;

                        for k = 1:numTrials
                            
                            thisTrialRow    = ((m-1).*numTrials) + (5.*m) + k;
                            theseTimeStamps = [PopData.session(i).unit(thisUnit).allPhase.seg(1).trial(trialIndsTmp(k)).ts,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(2).trial(trialIndsTmp(k)).ts,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(3).trial(trialIndsTmp(k)).ts,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(4).trial(trialIndsTmp(k)).ts,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(5).trial(trialIndsTmp(k)).ts];

                           if numel(theseTimeStamps)>0
                               tmp = find(theseTimeStamps<1.1e4 & theseTimeStamps>-0.5e3);
                           
                               delta(ceil(theseTimeStamps(tmp)+500)) = delta(ceil(theseTimeStamps(tmp)+500)) + 1;
                           
                               if rem(m,2)
                                    subplot(5,1,2:5); hold on;
                                    plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
                               else
                                    subplot(5,1,2:5);  hold on;
                                    plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);            
                               end                            
                            
                               plot(zeros(1,numTrials),((m-1).*numTrials) + (5.*m) + [1:numTrials],'b');
                               plot(PopData.session(i).behavior.flCSon(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');
%                                plot(PopData.session(i).behavior.USon(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'r');
%                                plot(PopData.session(i).behavior.flCSoff(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');
                                
                               if export
                                allStamps   = [allStamps theseTimeStamps'];
                                allRows     = [allRows ones(1,numel(theseTimeStamps)).*thisTrialRow]; 
                                allIDs      = [allIDs ones(1,numel(theseTimeStamps)).*m./numUnits];
                               end
                                
                           end
                           
                        end

                        PopData.session(i).unit(thisUnit).allPhase.timePSTH.delta = delta ./ numTrials;
                        
                    end

                    
                    if export
                        nameA = sprintf('/spikeTimesX%g',j-1);
                        nameB = sprintf('/spikeTimesY%g',j-1);
                        nameC = sprintf('/cellIDs%g',j-1);
                        hdf5write(['MatlabExport/egHahnStyleTime_' num2str(aa) '.h5'],nameA,allStamps);
                        hdf5write(['MatlabExport/egHahnStyleTime_' num2str(aa) '.h5'],nameB,allRows,'WriteMode','append');
                        hdf5write(['MatlabExport/egHahnStyleTime_' num2str(aa) '.h5'],nameC,allIDs,'WriteMode','append');
                    end
                    
                    subplot(5,1,2:5);
                    set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
                    axis([-500 7500 -1 ((numUnits).*numTrials) + (5.*(numUnits+1))]);
                    xlabel('Time (ms)');
                    ylabel('Sorted Trials per Unit');
                    drawnow; %pause();
                    
                    
                    % __________________________
                    % __________________________
                    % FOR PHASES
                    % __________________________
                    % __________________________
                    
                    figure(41); hold on;
                    
                    allStamps = []; allRows = []; allIDs = [];
                    
                    for m=1:numUnits
                        
                        delta = zeros(1,1.2e4);
                        deltaN = zeros(1,1.2e4);
                        
                        thisUnit = m;

                        for k = 1:numTrials
                            
                            thisTrialRow    = ((m-1).*numTrials) + (5.*m) + k;
                            theseTimeStamps = [PopData.session(i).unit(thisUnit).allPhase.seg(1).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(2).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(3).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(4).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(5).trial(trialIndsTmp(k)).ph];

                            tmp = find(theseTimeStamps>-0.5);          
                            delta(ceil((theseTimeStamps(tmp)+0.5).*1000)) = delta(ceil((theseTimeStamps(tmp)+0.5).*1000)) + 1;

                            % normalized delta for rate matching
                            segDurations = [1000 718 1077 4488 1000];
                            for pp=1:5
                                someStamps  = PopData.session(i).unit(thisUnit).allPhase.seg(pp).trial(trialIndsTmp(k)).ph;
                                tmpN        = find(someStamps>-0.5);
                                deltaN(ceil((someStamps(tmpN)+0.5).*1000)) = deltaN(ceil((someStamps(tmpN)+0.5).*1000)) + (PopData.session(i).unit(thisUnit).allPhase.seg(pp).meanDur ./ PopData.session(i).unit(thisUnit).allPhase.seg(pp).trial(trialIndsTmp(k)).du);
                            end
                            
                            
%             validStamps = find(theseStamps>-PopData.session(i).timeSinceCS(currTrial) & theseStamps<1);
%             validStamps2= find(theseStamps>0);
%             PopData.session(i).unit(j).respPCS.raster.trial(currTrial).ph = [(theseStamps(validStamps) .* (6280 ./ PopData.session(i).timeSinceCS(currTrial))) theseStamps(validStamps2).*(1./2e3).*ceil(2e6./avgTsCS)];
%             
%             deltaP(ceil(theseStamps(validStamps) .* (6280 ./ PopData.session(i).timeSinceCS(currTrial))) + 6280) = (avgTsCS / PopData.session(i).timeSinceCS(currTrial));
%             deltaP( ceil(theseStamps(validStamps2).*(1./2e3).*ceil(2e6./avgTsCS)) + 6280 ) = 1;
%             
%             image(kk,:) =  conv(deltaP,currParams.filter.kernel,'same');

                            
                            % plotting
                            if rem(m,2)
                                subplot(5,1,2:5); hold on;
                                plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
                            else
                                subplot(5,1,2:5);  hold on;
                                plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);            
                            end
                            
                            if export
                                allStamps   = [allStamps theseTimeStamps'];
                                allRows     = [allRows ones(1,numel(theseTimeStamps)).*thisTrialRow]; 
                                allIDs      = [allIDs ones(1,numel(theseTimeStamps)).*m./numUnits];
                            end
                            
                            plot(zeros(1,numTrials),((m-1).*numTrials) + (5.*m) + [1:numTrials],'b');
                            plot(ones(1,numTrials).*0.718,((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');
                            plot(ones(1,numTrials).*1.795,((m-1).*numTrials) + (5.*m) + [1:numTrials],'r');
                            plot(ones(1,numTrials).*6.283,((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');

                        end
                        
                        PopData.session(i).unit(thisUnit).allPhase.phasePSTH.delta = delta ./ numTrials;
                        PopData.session(i).unit(thisUnit).allPhase.phasePSTH.deltaN = deltaN ./ numTrials;

                    end

                    
                    if export
                        nameA = sprintf('/spikeTimesX%g',j-1);
                        nameB = sprintf('/spikeTimesY%g',j-1);
                        nameC = sprintf('/cellIDs%g',j-1);
                        hdf5write(['MatlabExport/egHahnStylePhase_' num2str(aa) '.h5'],nameA,allStamps);
                        hdf5write(['MatlabExport/egHahnStylePhase_' num2str(aa) '.h5'],nameB,allRows,'WriteMode','append');
                        hdf5write(['MatlabExport/egHahnStylePhase_' num2str(aa) '.h5'],nameC,allIDs,'WriteMode','append');
                    end

                    
                    subplot(5,1,2:5);
                    set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
                    axis([-0.5 6.785 -1 ((numUnits).*numTrials) + (5.*(numUnits+1))]);
                    xlabel('Phase (rad)');
                    ylabel('Sorted Trials per Unit');
                    drawnow; pause(1);


                    
                case 0
                    
                    % __________________________
                    % __________________________
                    % FOR PHASES
                    % __________________________
                    % __________________________
                    
                    figure(38); clf; hold on;
                    
                    % sort by init response latency to make it easier to see
%                     [vals,sortedLats] = sort(InitRespLatency);                    
                    [vals,sortedLats] = sort(1:numUnits ); 
                    
                    for m=1:numUnits
                        
                        delta = zeros(1,1.2e4);
                        
                        thisUnit = sortedLats(m);

                        for k = 1:numTrials
                            
                            thisTrialRow    = ((numUnits+3).*k) + m;
                            theseTimeStamps = [PopData.session(i).unit(thisUnit).allPhase.seg(1).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(2).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(3).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(4).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(5).trial(trialIndsTmp(k)).ph];
                                           
                            delta(ceil((theseTimeStamps+0.5).*1000)) = delta(ceil((theseTimeStamps+0.5).*1000)) + 1;

% %                             if rem(m,2)
%                                 subplot(5,1,2:5); hold on;
%                                 plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',6,'Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
% %                             else
% %                                 subplot(5,1,2:5);  hold on;
% %                                 plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',6,'Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);            
% %                             end

                            plot(zeros(1,numUnits),((numUnits+3).*k) + [1:numUnits],'b');
                            plot(ones(1,numUnits).*0.718,((numUnits+3).*k)+ [1:numUnits],'k');
                            plot(ones(1,numUnits).*1.795,((numUnits+3).*k) + [1:numUnits],'r');
                            plot(ones(1,numUnits).*6.283,((numUnits+3).*k) + [1:numUnits],'k');

                        end
                        
                        PopData.session(i).unit(thisUnit).allPhase.phasePSTH.delta = delta ./ numTrials;

                    end


                    
                    subplot(5,1,2:5);
                    set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
                    axis([-1 6.285 -1 ((numUnits).*numTrials) + (5.*(numUnits+1))]);
                    xlabel('Phase (rad)');
                    ylabel('Sorted Trials per Unit');
                    drawnow; pause(1);


            end

        end            
end

%% COMPARE PHASE AND TIME SPACE PSTHS
figure(39); clf;
% for i=1:NumSessions
count=1;
InitRespLatency

for aa=1:numel(sessList)
    i=sessList(aa);
% for i=83;    
    figure(39); clf;

        disp(['Session ' num2str(i)]);
        numUnits = size(PopData.session(i).unit,2);

        if size(PopData.session(i).PhaseMat,1) > 1
        
        % sort by init response latency to make it easier to see
    %     [vals,sortedLats] = sort(InitRespLatency); 
        [vals,sortedLats] = sort( 1:numUnits ); 

        PhaseMat = zeros(numUnits,7783);

        for m=1:numUnits

            thisUnit = sortedLats(m);

%             tmpSmoothP = conv(PopData.session(i).unit(thisUnit).allPhase.phasePSTH.delta,currParams.filter.kernel,'same');
            tmpSmoothP = conv(PopData.session(i).unit(thisUnit).allPhase.phasePSTH.deltaN,currParams.filter.kernel,'same');
            PopData.session(i).unit(thisUnit).allPhase.phasePSTH.psth = tmpSmoothP(1:7783);
%             disp('Finished phase psth.');

            tmpSmoothT = conv(PopData.session(i).unit(thisUnit).allPhase.timePSTH.delta,currParams.filter.kernel,'same');
            PopData.session(i).unit(thisUnit).allPhase.timePSTH.psth = tmpSmoothT(1:7783);
%             disp('Finished time psth.');

            % compare to psth with no selection and no warping
            numTrials = size(PopData.session(i).unit(thisUnit).respCS.raster.trial,2);
            delta = zeros(1,11000);

            for k = 1:numTrials
               theseTimeStamps = PopData.session(i).unit(thisUnit).respCS.raster.trial(k).ts;
               tmp = find(theseTimeStamps<9000 & theseTimeStamps>-500);
               delta(ceil(theseTimeStamps(tmp)+500)) = delta(ceil(theseTimeStamps(tmp)+500)) + 1;
            end

            delta = delta ./ numTrials;
            tmpSmoothU = conv(delta,currParams.filter.kernel,'same');
%             disp('Finished U psth.');

%             figure(39); hold on;
%             plot(-1000:6299,(tmpSmoothP./max(tmpSmoothP))+m-1,'Color',[1 0 0]);
%             plot(-1000:6299,(tmpSmoothT./max(tmpSmoothP))+m-1,'Color',[0 0.67 1]);            
%             plot(-1000:6299,(tmpSmoothU./max(tmpSmoothP))+m-1,'Color',[0.5 0.5 0.5]);

%             if mean( tmpSmoothP(4000:6000) ) - mean( tmpSmoothP(1:500) ) < 0
                InitRespLatency(count) = find(tmpSmoothP(500:1218)==max(tmpSmoothP(500:1218)),1) .* 1;
%             else
%                 InitRespLatency(count) = find(tmpSmoothP==max(tmpSmoothP(500:1218)),1) .* 1;                
%             end
            PhaseMat(thisUnit,:) = tmpSmoothP(1:7783) - mean(tmpSmoothP(1:500));
            
            phasePSTHs(count,:)     = ( tmpSmoothP(1:7783)-mean(tmpSmoothP(1:500)) ) ./ std(tmpSmoothP(1:7783));
            phasePSTHsN(count,:)    = phasePSTHs(count,:) ./ max( phasePSTHs(count,:) );
            timePSTHs(count,:)      = ( tmpSmoothT(1:7783)-mean(tmpSmoothT(1:500)) ) ./ std(tmpSmoothT(1:7783));
            timePSTHsN(count,:)     = timePSTHs(count,:) ./ max( phasePSTHs(count,:) );
            rawPSTHs(count,:)       = ( tmpSmoothU(1:7783)-mean(tmpSmoothU(1:500)) ) ./ std(tmpSmoothU(1:7783));
            
            warpEffect(count,1)     = max(tmpSmoothP(500:1218)).*1000;
            warpEffect(count,2)     = max(tmpSmoothU(500:1218)).*1000;
            warpEffect(count,3)     = i;
            warpEffect(count,4)     = thisUnit;
            
            count=count+1;
            

        end

%         plot([0,0],[0,numUnits],'k--');
%         plot([718,718],[0,numUnits],'k--');
%         plot([1795,1795],[0,numUnits],'k--');
%         plot([6283,6283],[0,numUnits],'k--');
%         axis([-1100 6400 -1 numUnits+1]);

        PopData.session(i).PhaseMat = PhaseMat;

%         figure(36);
%         imagesc(corr(PhaseMat),[-1 1]);
% 
%         figure(37);
%         imagesc(corr(PhaseMat'),[-1 1]);
        else
                
            disp('No good behavioral data... skipping');

        end
end

% [mapName] = TNC_CreateRBColormap([min(min(phasePSTHs)).*100 , max(max(phasePSTHs)).*100],'bo');
[mapName] = TNC_CreateRBColormap(4096,'bo');
[amps,inds] = sort( InitRespLatency ); 

figure(3); clf;
subplot(121);
imagesc(timePSTHsN(inds,1:1718),[-1 1]);
    colormap(mapName);
    hold on;
    plot(InitRespLatency(inds)+500,1:numel(InitRespLatency),'k-');
    title('Normalized response in real time');
    ylabel('Unit Index'); xlabel('Time from CS onset (ms)');
subplot(122);
imagesc(phasePSTHsN(inds,1:1718),[-1 1]);
    hold on;
    plot(InitRespLatency(inds)+500,1:numel(InitRespLatency),'k-');
    colormap(mapName);
    title('Normalized response in warped time');
    ylabel('Unit Index'); xlabel('Time from CS onset (ms)');
     
%% For phase PSTH show the clustered plot
    
    numClasses = 10;


[clustIds , c , sumdist] = kmeans(timePSTHs(:,500:2200),numClasses,'Start','cluster');
[vals,inds] = sort(clustIds);


[mapName] = TNC_CreateRBColormap(1024,'mbr');

figure(204); 
    subplot(131);
imagesc(phasePSTHs(inds,:),[-5 5]); 
colormap(mapName);

    subplot(132);
imagesc(timePSTHs(inds,:),[-5 5]);
colormap(mapName);

    subplot(133);
imagesc(corr(phasePSTHs),[-1 1]);
colormap(mapName);

figure(206); clf;
    
classOrder = [1:10];

for j=1:numClasses
    
    i = classOrder(j);
    
    theseInds = find(clustIds==i);
    subsetMatrix = phasePSTHs(theseInds,:);
    allCSresp.PclassMeans(i,:) = mean(subsetMatrix);
    allCSresp.PclassErrs(i,:) = std(subsetMatrix,0,1) ./ sqrt(numel(theseInds)-1);

    subsetMatrix = timePSTHs(theseInds,:);
    allCSresp.TclassMeans(i,:) = mean(subsetMatrix);
    allCSresp.TclassErrs(i,:) = std(subsetMatrix,0,1) ./ sqrt(numel(theseInds)-1);

    subplot(121);
    shadedErrorBar(-499:7283,allCSresp.PclassMeans(i,:)-(j.*3),allCSresp.PclassErrs(i,:),{'color',[1-(j./numClasses) 0 (j./numClasses)]}); hold on;
    plot([-500 7283],[-j.*3 -j.*3],'k--'); axis off;
    text(-1500,-(j.*3)+0.5,[num2str(i) '..' sprintf('%g',numel(theseInds))]);
    
    subplot(122);
    shadedErrorBar(-499:7283,allCSresp.TclassMeans(i,:)-(j.*3),allCSresp.TclassErrs(i,:),{'color',[1-(j./numClasses) 0 (j./numClasses)]}); hold on;    
    plot([-500 7283],[-j.*3 -j.*3],'k--'); axis off;
    text(-1500,-(j.*3)+0.5,[num2str(i) '..' sprintf('%g',numel(theseInds))]);

    
end

%% DIMENSIONALITY REDUCTION OF POPULATION ACTIVITY - SINGLE SESSION
% This segment of the code is meant to calculate two kinds of reduced
% representations:
%   1) What are the minimal characteristic time courses of activation over the approach window?
%   2) What is the trajectory of the population vector over the approach window?

% Define the time window over which to reduce dimensions (window of
% approach seems most applicable)
window = 1000:3500;
widthSig = currParams.smthParams.decay;
thisSession = 83;
sortedByLat = 1;


for aa=1:numel(sessList)

    [PopVec.kernel]  = TNC_CreateGaussian(widthSig.*15,widthSig,widthSig.*30,1);
    
    thisSession = sessList(aa);
    
    PhaseMat = PopData.session(thisSession).PhaseMat;
    numUnits = size(PhaseMat,1);

disp(' ');
disp(' ');
disp('________________________________');
disp('Computing principal components...');
disp(' ');

    % First we consider the basic features of temporal evolution of activity (an analogy to the continuous behavioral data we intend to extract)
    covMatForTimeCourse = cov(PhaseMat(:,window));
    [v,d] = eig(covMatForTimeCourse);
    TimeCourse.pca.vecs = zeros(numel(window),3);
    TimeCourse.pca.vals = diag(d)./sum(diag(d));
    TimeCourse.pca.vecs(:,1) = v(:,numel(TimeCourse.pca.vals));
    TimeCourse.pca.vecs(:,2) = v(:,numel(TimeCourse.pca.vals)-1);
    TimeCourse.pca.vecs(:,3) = v(:,numel(TimeCourse.pca.vals)-2);
    
    TimeCourse.pca.varExp = sum(TimeCourse.pca.vals(numel(TimeCourse.pca.vals)-2:numel(TimeCourse.pca.vals)));
    disp(['Variance explained (time): ' num2str(sum(TimeCourse.pca.vals(numel(TimeCourse.pca.vals)-2:numel(TimeCourse.pca.vals)))) ' | ' num2str(TimeCourse.pca.vals(numel(TimeCourse.pca.vals))) ' ' num2str(TimeCourse.pca.vals(numel(TimeCourse.pca.vals)-1)) ' ' num2str(TimeCourse.pca.vals(numel(TimeCourse.pca.vals)-2)) ]);
    disp(' ');

    % Second we consider the evolution of the population vector using the warped response patterns
    covMatForPopVec = cov(PhaseMat(:,window)');
    [v,d] = eig(covMatForPopVec);
    PopVec.pca.vecs = zeros(size(PhaseMat,1),3);
    PopVec.pca.vals = diag(d)./sum(diag(d));
    PopVec.pca.vecs(:,1) = v(:,numel(PopVec.pca.vals));
    PopVec.pca.vecs(:,2) = v(:,numel(PopVec.pca.vals)-1);
    PopVec.pca.vecs(:,3) = v(:,numel(PopVec.pca.vals)-2);
    PopVec.pca.varExp = sum(PopVec.pca.vals(numel(PopVec.pca.vals)-2:numel(PopVec.pca.vals)));
    disp(['Variance explained (population): ' num2str(sum(PopVec.pca.vals(numel(PopVec.pca.vals)-2:numel(PopVec.pca.vals)))) ' | ' num2str(PopVec.pca.vals(numel(PopVec.pca.vals))) ' ' num2str(PopVec.pca.vals(numel(PopVec.pca.vals)-1)) ' ' num2str(PopVec.pca.vals(numel(PopVec.pca.vals)-2))]);

disp(' ');
disp('Visualizing structure of covariance matrices...');
disp(' ');

    figure(21);
    subplot(2,3,1);
    imagesc(covMatForTimeCourse);
    title('Covariance of PSTH time courses');

    subplot(2,3,2); hold off;
    plot(TimeCourse.pca.vecs(:,1),'Color',[0 0 0],'LineWidth',2); hold on;
    plot(TimeCourse.pca.vecs(:,2),'Color',[1 0 0],'LineWidth',2);
    plot(TimeCourse.pca.vecs(:,3),'Color',[0 0.67 1],'LineWidth',2);
    title('Leading eigenvectors for TimeCourse');

    subplot(2,3,3); hold off; step = 10;
    for j=step+1:step:numel(window)
        plot3(TimeCourse.pca.vecs(j-(step-1):j,1),TimeCourse.pca.vecs(j-(step-1):j,2),TimeCourse.pca.vecs(j-(step-1):j,3),'Color',[1-(j./numel(window)) (0.67.*(j./numel(window))) j./numel(window)],'LineWidth',2); hold on;
    end
    title('Evolution of TimeCourse dimensions');
    grid on;

    subplot(2,3,4);
    imagesc(covMatForPopVec);
    title('Covariance of individual PSTHs');

    subplot(2,3,5); hold off;
    plot(PopVec.pca.vecs(:,1),'Color',[0 0 0],'LineWidth',2); hold on;
    plot(PopVec.pca.vecs(:,2),'Color',[1 0 0],'LineWidth',2);
    plot(PopVec.pca.vecs(:,3),'Color',[0 0.67 1],'LineWidth',2);
    title('Leading eigenvectors for PopVec');

    subplot(2,3,6); hold off;
    % project all units for display
    for i=1:size(PhaseMat,1)
        TimeCourse.proj(1,i) = dot(PhaseMat(i,window),TimeCourse.pca.vecs(:,1));
        TimeCourse.proj(2,i) = dot(PhaseMat(i,window),TimeCourse.pca.vecs(:,2));
        TimeCourse.proj(3,i) = dot(PhaseMat(i,window),TimeCourse.pca.vecs(:,3));

        subplot(2,3,6);

        plot3(TimeCourse.proj(1,i),TimeCourse.proj(2,i),TimeCourse.proj(3,i),'o','Color',[1-(i./numUnits) (0.67.*(i./numUnits)) i./numUnits],'LineWidth',2,'MarkerSize',8); hold on;
    end
    title('Projection onto TimeCourse dimensions');
    grid on;

disp(' ');
disp('Projecting individual trial trajectories in low dimensions...');
disp(' ');
    
    validTrials = PopData.session(thisSession).behavior.validPhaseTrials;
    numTrials = numel(validTrials);

    if sortedByLat==1
        disp('Sorting by latency...');
        yLabelStr = '(sorted)';
        latencies = PopData.session(thisSession).behavior.flCSon(validTrials);
        [vals,trialIndsTmp] = sort(latencies);
        trialInds = validTrials(trialIndsTmp);
    else
        disp('No sorting was applied...');
        trialInds = validTrials;
        yLabelStr = '(not sorted)';
    end

    for k=1:numTrials
                          
        VectorMat = zeros(numUnits,7300); 
        PopVec.trial(k).proj = zeros(3,7300);
        trialIndsTmp(k);
        
        for l=1:numUnits
            
            thisUnit = l;
            
            % collect "phaseStamps"
            theseTimeStamps = [PopData.session(thisSession).unit(thisUnit).allPhase.seg(1).trial(trialIndsTmp(k)).ph,
                               PopData.session(thisSession).unit(thisUnit).allPhase.seg(2).trial(trialIndsTmp(k)).ph,
                               PopData.session(thisSession).unit(thisUnit).allPhase.seg(3).trial(trialIndsTmp(k)).ph,
                               PopData.session(thisSession).unit(thisUnit).allPhase.seg(4).trial(trialIndsTmp(k)).ph];
                           
            ts1msRes = ceil(theseTimeStamps.*1000) + 1000;
            % smooth to create instantaneous rate vectors
            delta = zeros(1,7300);
            delta(1,ts1msRes) = 1;
            tmpSmoothP = conv(delta,PopVec.kernel,'same');

            % compile into a matrix for easy projections
            VectorMat(l,:) = tmpSmoothP;
        end
        
        for m = 1:7300
            PopVec.trial(k).proj(1,m) = dot(VectorMat(:,m),PopVec.pca.vecs(:,1));
            PopVec.trial(k).proj(2,m) = dot(VectorMat(:,m),PopVec.pca.vecs(:,2));
            PopVec.trial(k).proj(3,m) = dot(VectorMat(:,m),PopVec.pca.vecs(:,3));
        end

    end
    
disp(' ');
disp('Visualizing individual trial trajectories...');
disp(' ');

    figure(22); clf; hold on;
    bins = 5;
    binUnits = floor(numTrials./bins);
    PopVec.trial(k).proj = zeros(3,7300); 
    
    for k=1:bins

        tempMat1 = zeros(binUnits,7300);
        tempMat2 = zeros(binUnits,7300);
        tempMat3 = zeros(binUnits,7300);

        for m=1:binUnits
            thisTrial = m + ((k-1)*binUnits);
            tempMat1(m,:) = PopVec.trial(thisTrial).proj(1,:);
            tempMat2(m,:) = PopVec.trial(thisTrial).proj(2,:);
            tempMat3(m,:) = PopVec.trial(thisTrial).proj(3,:);
        end
        
        PopVec.trialBin(k).proj         = zeros(3,7300);
        PopVec.trialBin(k).proj(1,:)    = mean(tempMat1);
        PopVec.trialBin(k).proj(2,:)    = mean(tempMat2);
        PopVec.trialBin(k).proj(3,:)    = mean(tempMat3);
        
        % colorize by fl latency
        figure(22);        
        thisColor = [1-(k./bins) (0.67.*(k./bins)) k./bins];
        window = 1:7000;
        plot3(PopVec.trialBin(k).proj(1,window),PopVec.trialBin(k).proj(2,window),PopVec.trialBin(k).proj(3,window),'-','Color',thisColor,'LineWidth',4); 
        plot3(PopVec.trialBin(k).proj(1,1000),PopVec.trialBin(k).proj(2,1000),PopVec.trialBin(k).proj(3,1000),'o','Color',thisColor,'LineWidth',2,'MarkerSize',10);
        plot3(PopVec.trialBin(k).proj(1,1718),PopVec.trialBin(k).proj(2,1718),PopVec.trialBin(k).proj(3,1718),'s','Color',thisColor,'LineWidth',2,'MarkerSize',10);
        plot3(PopVec.trialBin(k).proj(1,2795),PopVec.trialBin(k).proj(2,2795),PopVec.trialBin(k).proj(3,2795),'d','Color',thisColor,'LineWidth',2,'MarkerSize',10);
        plot3(PopVec.trialBin(k).proj(1,1055:1150),PopVec.trialBin(k).proj(2,1055:1150),PopVec.trialBin(k).proj(3,1055:1150),'k--','LineWidth',4);
        grid on; view([144 62]);
        xlabel('PCA 1');
        ylabel('PCA 2');
        zlabel('PCA 3');
    end
            
    PopData.session(thisSession).PopVec = PopVec;
    PopData.session(thisSession).TimeCourse = TimeCourse;

    clear PopVec TimeCourse;

end    
    
    
    
    
disp(' ');
disp('________________________________');
disp(' ');

%% DIMENSIONALITY REDUCTION OF POPULATION MEAN

figOffset = 0 

% for aa=1:numel(sessList)
for i=83;    
%     i=sessList(aa);

    if size(PopData.session(i).PhaseMat,1)>1
            figure(39+figOffset); clf; 
        [mappedA, mapping] = compute_mapping(phasePSTHs', 'PCA', 3);
%         [mappedA, mapping] = compute_mapping(PopData.session(i).PhaseMat', 'PCA', 3);

        figure(39+figOffset); subplot(4,1,1:3);
        plot3(mappedA(25:500,1),mappedA(25:500,2),mappedA(25:500,3),'color',[0.5 0.5 0.5]);hold on;
        plot3(mappedA(500:1218,1),mappedA(500:1218,2),mappedA(500:1218,3),'color',[0 0.67 1]);
        plot3(mappedA(1218:2295,1),mappedA(1218:2295,2),mappedA(1218:2295,3),'color',[1 0 0]);
        plot3(mappedA(2295:6783,1),mappedA(2295:6783,2),mappedA(2295:6783,3),'color',[0 0 0]);
        plot3(mappedA(6783:7753,1),mappedA(6783:7753,2),mappedA(6783:7753,3),'color',[0.5 0.5 0.5]);

%         title(['Session ' num2str(i) ' | Number of units: ' num2str(numel(PopData.session(i).unit))]);
        title(['All session data']);
        ylabel('PC2'); xlabel('PC1'); zlabel('PC3');

        %     plot(mappedA(25:500,1),mappedA(25:500,2),'.-','color',[0.5 0.5 0.5]);hold on;
        %     plot(mappedA(500:1218,1),mappedA(500:1218,2),'.-','color',[0 0.67 1]);
        %     plot(mappedA(1218:2295,1),mappedA(1218:2295,2),'.-','color',[1 0 0]);
        %     plot(mappedA(2295:6783,1),mappedA(2295:6783,2),'.-','color',[0 0 0]);
        %     plot(mappedA(6783:7753,1),mappedA(6783:7753,2),'.-','color',[0.5 0.5 0.5]);

        figure(39+figOffset); subplot(4,1,4);

        velocity = sqrt( conv(diff(mappedA(:,1)),[0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0],'same').^2 + conv(diff(mappedA(:,2)),[0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0],'same').^2 + conv(diff(mappedA(:,3)),[0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0],'same').^2  );
        % plot(-499:7282,velocity);
        plot(50:450,velocity(50:450,1) - mean(velocity(50:450,1)),'color',[0.5 0.5 0.5]);hold on;
        plot(600:1168,velocity(600:1168,1) - mean(velocity(50:450,1)),'color',[0 0.67 1]);hold on;
        plot(1268:2260,velocity(1268:2260,1) - mean(velocity(50:450,1)),'color',[1 0 0]);hold on;
        plot(2345:6733,velocity(2345:6733,1) - mean(velocity(50:450,1)),'color',[0 0 0]);hold on;
        plot(6848:7733,velocity(6848:7733,1) - mean(velocity(50:450,1)),'color',[0.5 0.5 0.5]);hold on;
        plot([0 7783],[0 0],'k--');
        set(gca,'TickDir','out'); box off; ylabel('Velocity (a.u.)'); xlabel('Phase');
        % plot(mappedA(1218:2295,1),mappedA(1218:2295,2),'.-','color',[1 0 0]);
        % plot(mappedA(2295:6783,1),mappedA(2295:6783,2),'.-','color',[0 0 0]);
        % plot(mappedA(6783:7753,1),mappedA(6783:7753,2),'.-','color',[0.5 0.5 0.5]);

%         pause();
    
     else

        disp('No good phase data... skipping');   
        
    end
end

%% DIMENSIONALITY REDUCTION OF SESSIONS DATA - FAST METHOD
% This segment of the code is meant to calculate two kinds of reduced
% representations:
%   1) What are the minimal characteristic time courses of activation over the approach window?
%   2) What is the trajectory of the population vector over the approach window?

% Define the time window over which to reduce dimensions (window of
% approach seems most applicable)
window = 1000:4000;
widthSig = currParams.smthParams.decay;
modWidth = 4;
sortedByLat = 1;
visualizeTrialwiseProjections = 0;

for aa=1:numel(sessList)
    
    thisSession = sessList(aa);
    
% for aa=1
%     thisSession = 83
    
    [PopVec.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15.*modWidth,currParams.smthParams.decay.*modWidth,currParams.smthParams.decay.*modWidth.*30,1);
    figure(21); clf;
    PhaseMat = PopData.session(thisSession).PhaseMat;
    numUnits = size(PhaseMat,1);
    
    if numUnits>3
        
        validSess(aa) = 1;

disp(' ');
disp(' ');
disp('________________________________');
disp('Computing principal components...');
disp(' ');

    % First we consider the basic features of temporal evolution of activity (an analogy to the continuous behavioral data we intend to extract)
    [mappedA_t, mapping_t] = compute_mapping(PhaseMat', 'PCA', 3);
    TimeCourse.pca.vecs = mapping_t.M;
    TimeCourse.pca.load = mappedA_t;
    
    % Second we consider the evolution of the population vector using the warped response patterns
    covMatForTimeCourse = cov(PhaseMat(:,window));
    [mappedA_u, mapping_u] = compute_mapping(PhaseMat, 'PCA', 3);
    PopVec.pca.vecs = mapping_u.M;
    PopVec.pca.load = mappedA_u;
    
disp(' ');
disp('Visualizing structure of covariance matrices...');
disp(' ');

    subplot(2,1,1);
        plot3(TimeCourse.pca.load(25:500,1),TimeCourse.pca.load(25:500,2),TimeCourse.pca.load(25:500,3),'color',[0.5 0.5 0.5]);hold on;
        plot3(TimeCourse.pca.load(500:1218,1),TimeCourse.pca.load(500:1218,2),TimeCourse.pca.load(500:1218,3),'color',[0 0.67 1]);
        plot3(TimeCourse.pca.load(1218:2295,1),TimeCourse.pca.load(1218:2295,2),TimeCourse.pca.load(1218:2295,3),'color',[1 0 0]);
        plot3(TimeCourse.pca.load(2295:6783,1),TimeCourse.pca.load(2295:6783,2),TimeCourse.pca.load(2295:6783,3),'color',[0 0 0]);
        plot3(TimeCourse.pca.load(6783:7753,1),TimeCourse.pca.load(6783:7753,2),TimeCourse.pca.load(6783:7753,3),'color',[0.5 0.5 0.5]);
        axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]); view([22 62]);
    subplot(2,1,2);
        scatter3(PopVec.pca.load(:,1),PopVec.pca.load(:,2),PopVec.pca.load(:,3),20,1:numUnits); hold on;
        plot3([0 0],[-1 1],[0 0],'k-',[0 0],[0 0],[-1 1],'k-',[-1 1],[0 0],[0 0],'k-'); grid off;

disp(' ');
disp('Projecting individual trial trajectories in low dimensions...');
disp(' ');
    
    validTrials = PopData.session(thisSession).behavior.validPhaseTrials;
    numTrials = numel(validTrials);

    if sortedByLat==1
        disp('Sorting by latency...');
        yLabelStr = '(sorted)';
        latencies = PopData.session(thisSession).behavior.flCSon(validTrials);
        [vals,trialIndsTmp] = sort(latencies);
        trialInds = validTrials(trialIndsTmp);
    else
        disp('No sorting was applied...');
        trialInds = validTrials;
        yLabelStr = '(not sorted)';
    end

    for k=1:numTrials
                          
        VectorMat = zeros(size(PhaseMat)); 
        PopVec.trial(k).proj = zeros(3,size(PhaseMat,2));
        trialIndsTmp(k);
        
        for l=1:numUnits
            
            thisUnit = l;
            
            % collect "phaseStamps"
            theseTimeStamps = [PopData.session(thisSession).unit(thisUnit).allPhase.seg(1).trial(trialIndsTmp(k)).ph,
                               PopData.session(thisSession).unit(thisUnit).allPhase.seg(2).trial(trialIndsTmp(k)).ph,
                               PopData.session(thisSession).unit(thisUnit).allPhase.seg(3).trial(trialIndsTmp(k)).ph,
                               PopData.session(thisSession).unit(thisUnit).allPhase.seg(4).trial(trialIndsTmp(k)).ph,
                               PopData.session(thisSession).unit(thisUnit).allPhase.seg(5).trial(trialIndsTmp(k)).ph];
                           
            ts1msRes = ceil(theseTimeStamps.*1000) + 500;
            % smooth to create instantaneous rate vectors0
            delta = zeros(1,size(PhaseMat,2));
            delta(1,ts1msRes) = 1;
            tmpSmoothP = conv(delta,PopVec.kernel,'same');

            % compile into a matrix for easy projections
            VectorMat(l,:) = (tmpSmoothP - mean(tmpSmoothP));
            
%             if l==8
%                 figure(78); clf;
%                 plot(VectorMat(l,:),'k'); hold on;
%                 plot(PhaseMat(l,:));
%                 pause(0.1);
%             end
            
        end
        
        PopVec.trial(k).proj = VectorMat'*TimeCourse.pca.vecs;
        
    end
    
    if visualizeTrialwiseProjections
disp(' ');
disp('Visualizing individual trial trajectories...');
disp(' ');

            figure(22); clf; hold on;
            bins = 5;
            binUnits = floor(numTrials./bins);
            PopVec.trial(k).proj = zeros(3,size(PhaseMat,2)); 

            for k=1:bins

                tempMat1 = zeros(binUnits,size(PhaseMat,2));
                tempMat2 = zeros(binUnits,size(PhaseMat,2));
                tempMat3 = zeros(binUnits,size(PhaseMat,2));

                for m=1:binUnits
                    thisTrial = m + ((k-1)*binUnits);
                    tempMat1(m,:) = PopVec.trial(thisTrial).proj(1,:);
                    tempMat2(m,:) = PopVec.trial(thisTrial).proj(2,:);
                    tempMat3(m,:) = PopVec.trial(thisTrial).proj(3,:);
                end

                PopVec.trialBin(k).proj         = zeros(3,size(PhaseMat,2));
                PopVec.trialBin(k).proj(1,:)    = mean(tempMat1);
                PopVec.trialBin(k).proj(2,:)    = mean(tempMat2);
                PopVec.trialBin(k).proj(3,:)    = mean(tempMat3);

                % colorize by fl latency
                figure(22);        
                thisColor = [1-(k./bins) (0.67.*(k./bins)) k./bins];
                window = 1:7000;
                plot3(PopVec.trialBin(k).proj(1,window),PopVec.trialBin(k).proj(2,window),PopVec.trialBin(k).proj(3,window),'-','Color',thisColor,'LineWidth',4); 
                plot3(PopVec.trialBin(k).proj(1,1000),PopVec.trialBin(k).proj(2,1000),PopVec.trialBin(k).proj(3,1000),'o','Color',thisColor,'LineWidth',2,'MarkerSize',10);
                plot3(PopVec.trialBin(k).proj(1,1718),PopVec.trialBin(k).proj(2,1718),PopVec.trialBin(k).proj(3,1718),'s','Color',thisColor,'LineWidth',2,'MarkerSize',10);
                plot3(PopVec.trialBin(k).proj(1,2795),PopVec.trialBin(k).proj(2,2795),PopVec.trialBin(k).proj(3,2795),'d','Color',thisColor,'LineWidth',2,'MarkerSize',10);
                plot3(PopVec.trialBin(k).proj(1,1055:1150),PopVec.trialBin(k).proj(2,1055:1150),PopVec.trialBin(k).proj(3,1055:1150),'k--','LineWidth',4);
                grid on; view([144 62]);
                xlabel('PCA 1');
                ylabel('PCA 2');
                zlabel('PCA 3');
            end

            PopData.session(thisSession).PopVec = PopVec;
            PopData.session(thisSession).TimeCourse = TimeCourse;

            
    end
   
    figure(77); clf;
    numTrials = numel(PopVec.trial);
    [mapName] = TNC_CreateRBColormap(numTrials,'bo');
    for jj = 2:10:numTrials
        subplot(ceil(numTrials/20),2,ceil(jj/10));
    %     plot3(PopVec.trial(jj).proj(:,1),PopVec.trial(jj).proj(:,2),PopVec.trial(jj).proj(:,3),'color',mapName(jj,:)); hold on;
            plot3(PopVec.trial(jj).proj(25:500,1),PopVec.trial(jj).proj(25:500,2),PopVec.trial(jj).proj(25:500,3),'color',[0.5 0.5 0.5]);hold on;
            plot3(PopVec.trial(jj).proj(500:1218,1),PopVec.trial(jj).proj(500:1218,2),PopVec.trial(jj).proj(500:1218,3),'color',[0 0.67 1]);
            plot3(PopVec.trial(jj).proj(1218:2295,1),PopVec.trial(jj).proj(1218:2295,2),PopVec.trial(jj).proj(1218:2295,3),'color',[1 0 0]);
            plot3(PopVec.trial(jj).proj(2295:6783,1),PopVec.trial(jj).proj(2295:6783,2),PopVec.trial(jj).proj(2295:6783,3),'color',[0 0 0]);
            plot3(PopVec.trial(jj).proj(6783:7753,1),PopVec.trial(jj).proj(6783:7753,2),PopVec.trial(jj).proj(6783:7753,3),'color',[0.5 0.5 0.5]);
            axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]); view([22 62]);
    end
        
    clear PopVec TimeCourse;

    else
        
        validSess(aa) = 0;
    
    end
end

disp(' ');
disp('________________________________');
disp(' ');

sessList
validSess

%% DIMENSIONALITY REDUCTION OF POPULATION ACTIVITY - PSEUDO-SIMULTANEOUS

%% LOADING OF CONTINUOUS AND BEHAVIORAL DATA
% hard coded at the moment
fileNameStr         = 'DA-Ma-05 100231 trace-d-001.ns4'
Ns4DATA             = openNSx('read',fileNameStr);
session.creation    = Ns4DATA.MetaTags.CreateDateTime;

%% WARPING OF BEHAVIORAL DATA
% uses interp and decimate to resample behavioral data into the warped coordinates

%% CALCULATE POP VECTOR
% as a function of time bin (variable width) or phase bin (variable width)
% if pop activity is necessary then position of popvector in pca space trial by trial should predict differences in response latency
% if pop activity is a good way to pursue then ask the limits of approach behavior (i.e. when does structure appear in the pop vector)
    
binWidth = 1;
numBins = 10000 ./ binWidth;
binMin = 0;
binMax = binMin+binWidth;
sortedByLat = 1;
threeD=0;
clear PopVec
disp('Sorting...');

for aa=1:numel(sessList)

    i = sessList(aa);

    validTrials = PopData.session(i).behavior.validPhaseTrials;
    numTrials = numel(validTrials);
        
    if sortedByLat==1
        [vals,trialInds] = sort(PopData.session(i).behavior.flCSon);
    else
        trialInds = 1:numTrials;
    end

    numUnits = size(PopData.session(i).unit,2);
    disp(['Collecting vectors for ' num2str(numUnits) ' units.']);
    
    for m=1:numUnits

        thisUnit = m;

        for k = 1:numel(trialInds)

            % create smoothed inst rate
            delta = zeros(1,10000);
            tmpTimeStamps = round(PopData.session(i).unit(thisUnit).respCS.raster.trial(trialInds(k)).ts);
%             posStamps = find(tmpTimeStamps>-1000 & tmpTimeStamps<9000);
            theseTimeStamps = tmpTimeStamps(posStamps)+1000;
            delta(ceil(theseTimeStamps)) = 1;
            smthData = conv(delta,currParams.filter.kernel,'same');

    %         binMin = 0;
    %         binMax = binWidth;
    % 
    %         for l = 1:numBins
    % %             PopVec.trials(k).vecData(m,l) = numel(find(theseTimeStamps>binMin & theseTimeStamps<=binMax));
    %             PopVec.trials(k).vecData(m,l) = trapz(smthData(binMin:binMax));
    %             binMin = binMin+binWidth;
    %             binMax = binMax+binWidth;
    %         end

            % create smoothed rate in a trialwise matrix
            PopVec.trials(k).vecData(m,:) = smthData - mean(smthData);

        end

        figure(41); subplot(121);
        imagesc(PopVec.trials(k).vecData,[-0.1 0.1]);
        drawnow;

    end




    % calculate popVector pca for this session
        pcaWin = 1000:6000;
        for k = 1:numTrialsSpk
            covMat = cov(PopVec.trials(k).vecData(:,pcaWin)');
            if k==1
                tmpCov = covMat;
            else
                tmpCov = tmpCov+covMat;
            end
        end
        meanCovMat = tmpCov./numTrialsSpk;
        figure(41); subplot(122);
        imagesc(meanCovMat);
        drawnow;
        
    % slow method
        [v,d] = eig(meanCovMat);
        PopVec.pca.vals = diag(d)./sum(diag(d));
        PopVec.pca.vecs(:,1) = v(:,numel(PopVec.pca.vals));
        PopVec.pca.vecs(:,2) = v(:,numel(PopVec.pca.vals)-1);
        PopVec.pca.vecs(:,3) = v(:,numel(PopVec.pca.vals)-2);
        PopVec.pca.varExp = sum(PopVec.pca.vals(numel(PopVec.pca.vals)-2:numel(PopVec.pca.vals)));
        disp(['Variance explained by first 3 components: ' num2str(sum(PopVec.pca.vals(numel(PopVec.pca.vals)-2:numel(PopVec.pca.vals))))]);
    
    % fast method

    
    % plot the individual trial popVectors
    figure(42); clf; hold on;
    numTrials = 20;
    startTrial = 1;
    pcaTrajWin = 1000:2000;
    for k = startTrial:startTrial+numTrials
        for l = 1:numBins
            PopVec.trials(k).pcaTraj(1,l) = dot(PopVec.trials(k).vecData(:,l),PopVec.pca.vecs(:,1));
            PopVec.trials(k).pcaTraj(2,l) = dot(PopVec.trials(k).vecData(:,l),PopVec.pca.vecs(:,2));
            PopVec.trials(k).pcaTraj(3,l) = dot(PopVec.trials(k).vecData(:,l),PopVec.pca.vecs(:,3));
        end
        index = k-startTrial;
        thisColor = [1-(index./numTrials) 0.67*(index./numTrials) index./numTrials];

        tOfFL = ceil(PopData.session(i).behavior.flCSon(trialInds(k))) + 500;
        pcaTrajWin = 500:tOfFL+1000;

        if k<10
            subplot(211); hold on;
        else
            subplot(212); hold on;
        end

        grid on;

        if threeD==1
            grid on; view([-20 20]);
            xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3');
            plot3(PopVec.trials(k).pcaTraj(1,pcaTrajWin),PopVec.trials(k).pcaTraj(2,pcaTrajWin),PopVec.trials(k).pcaTraj(3,pcaTrajWin),'-','Color',thisColor,'LineWidth',1);
        %     plot3(PopVec.trials(k).pcaTraj(1,pcaTrajWin(1)),PopVec.trials(k).pcaTraj(2,pcaTrajWin(1)),PopVec.trials(k).pcaTraj(3,pcaTrajWin(1)),'^','MarkerSize',10,'Color',thisColor,'LineWidth',2);
            plot3(PopVec.trials(k).pcaTraj(1,pcaTrajWin(500)),PopVec.trials(k).pcaTraj(2,pcaTrajWin(500)),PopVec.trials(k).pcaTraj(3,pcaTrajWin(500)),'d','MarkerSize',10,'Color',thisColor,'LineWidth',2);
            plot3(PopVec.trials(k).pcaTraj(1,pcaTrajWin(numel(pcaTrajWin)-1000)),PopVec.trials(k).pcaTraj(2,pcaTrajWin(numel(pcaTrajWin)-1000)),PopVec.trials(k).pcaTraj(3,pcaTrajWin(numel(pcaTrajWin)-1000)),'s','MarkerSize',10,'Color',thisColor,'LineWidth',2);
        else
            xlabel('PCA 1'); ylabel('PCA 2');
            plot(PopVec.trials(k).pcaTraj(1,pcaTrajWin),PopVec.trials(k).pcaTraj(2,pcaTrajWin),'-','Color',thisColor,'LineWidth',1);
            plot3(PopVec.trials(k).pcaTraj(1,pcaTrajWin(numel(pcaTrajWin))),PopVec.trials(k).pcaTraj(2,pcaTrajWin(numel(pcaTrajWin))),PopVec.trials(k).pcaTraj(3,pcaTrajWin(numel(pcaTrajWin))),'^','MarkerSize',10,'Color',thisColor,'LineWidth',2);
            plot(PopVec.trials(k).pcaTraj(1,pcaTrajWin(500)),PopVec.trials(k).pcaTraj(2,pcaTrajWin(500)),'d','MarkerSize',10,'Color',thisColor,'LineWidth',2);
            plot(PopVec.trials(k).pcaTraj(1,pcaTrajWin(numel(pcaTrajWin)-1000)),PopVec.trials(k).pcaTraj(2,pcaTrajWin(numel(pcaTrajWin)-1000)),'s','MarkerSize',10,'Color',thisColor,'LineWidth',2);
        end
    end

end

%% CALCULATE CROSS CORRELATION HISTOGRAMS (-50 to +50 ms) FOR JENNY

    % Calc the cross correlation for baseline period and CS-FL
    
    % All pairwise correlations
    
%% COMPUTE FEATURE MATRIX FOR PARVEZ
% Proposed features for correlation with reaction time:

for aa=1:numel(sessList)

    thisSession = sessList(aa);

    clear popTraj

    disp('No sorting was applied...');
    validTrials = PopData.session(thisSession).behavior.validPhaseTrials;
    trialInds   = validTrials;
    numTrials   = numel(validTrials);
    numUnits    = size(PopData.session(thisSession).unit,2);
    
    range(1,:) = [1:240] + 1000;
    range(2,:) = [241:480] + 1000;
    range(3,:) = [481:720] + 1000;

    latencies   = PopData.session(thisSession).behavior.flCSon(trialInds);
    PopVec      = PopData.session(thisSession).PopVec;

    % For population metrics cycle through valid trials
    for k=1:numTrials
    % for k=1:1

        q = trialInds(k);

        % pca{1,2,3} loading in evenly divided segments (quintiles?)    
        popTraj(k,1) = mean(PopVec.trial(k).proj(1,range(1,:)));
        popTraj(k,2) = mean(PopVec.trial(k).proj(2,range(1,:)));
        popTraj(k,3) = mean(PopVec.trial(k).proj(3,range(1,:)));

        popTraj(k,4) = mean(PopVec.trial(k).proj(1,range(2,:)));
        popTraj(k,5) = mean(PopVec.trial(k).proj(2,range(2,:)));
        popTraj(k,6) = mean(PopVec.trial(k).proj(3,range(2,:)));

        popTraj(k,7) = mean(PopVec.trial(k).proj(1,range(3,:)));
        popTraj(k,8) = mean(PopVec.trial(k).proj(2,range(3,:)));
        popTraj(k,9) = mean(PopVec.trial(k).proj(3,range(3,:)));

        latency      = round(PopData.session(thisSession).behavior.flCSon(q));
                
        for p=1:numUnits

            % disp('start.')
            
            theseTimeStamps     = PopData.session(thisSession).unit(p).respCS.raster.trial(k).ts;
            validStampsI        = find(theseTimeStamps > -599 & theseTimeStamps < 600);
            validStamps         = round(theseTimeStamps(validStampsI) + 600);
            delta               = zeros(1,1200);
            delta(validStamps)  = 1;
            currPSTH            = conv(delta,currParams.filter.kernel,'same');

            numPnts             = numel([-599:latency]);
            valPStampsI         = find(theseTimeStamps > -599 & theseTimeStamps < latency);
            valPStamps          = round(theseTimeStamps(valPStampsI) + 600);
            deltaP              = zeros(1,numPnts);
            deltaP(valPStamps)  = 1;
            curPPSTH            = conv(deltaP,currParams.filter.kernel,'same');

            thesePhaseStamps    = ceil(PopData.session(thisSession).unit(p).allPhase.seg(2).trial(k).ph.*1000);        
            deltaW              = zeros(1,750);
            
            if numel(thesePhaseStamps) > 0
                deltaW(thesePhaseStamps)= 1;
            end
            
            curWPSTH                = conv(deltaW,currParams.filter.kernel,'same');

            
            % build matrix for corr
            if p==1
                currTrialCorrMat = zeros(numUnits,1200);
                currPhaseCorrMat = zeros(numUnits,750);
            end

            currTrialCorrMat(p,:) = currPSTH;
            currPhaseCorrMat(p,:) = curWPSTH;

            if p==numUnits
                corrStruc = corr(currTrialCorrMat');
                lowInds = find( tril( corrStruc , -1 ) );
                if k==1
                    corrPWise = zeros(numTrials,numel(lowInds));
                end
                corrPWise(k,:) = corrStruc(lowInds);
            end

            if p==numUnits
                corrStrucP = corr(currPhaseCorrMat');
                lowIndsP = find( tril( corrStrucP , -1 ) );
                if k==1
                    corPPWise = zeros(numTrials,numel(lowIndsP));
                end
                corPPWise(k,:) = corrStrucP(lowIndsP);
            end

            baseline    = trapz(currPSTH(1:450))./450;

            response    = trapz(currPSTH(450:1200))./750;
            respRaw     = currPSTH(450:1200);

            respP       = trapz(curPPSTH(450:latency))./(latency+450);
            respPR      = curPPSTH(450:latency);

            respW       = curWPSTH;

            % peak firing rate - raw PSTH and phase PSTH
            peaks(k,p) = max(respRaw)-mean(currPSTH(1:450));
            peaksP(k,p) = max(respPR)-mean(currPSTH(1:450));
            peaksW(k,p) = max(respW);

            % time of peak - raw PSTH and phase PSTH
            peakT(k,p)      = find( respRaw==max(respRaw) , 1);
            peakTP(k,p)     = find( respPR==max(respPR) , 1);
            peakTW(k,p)     = find( respW==max(respW) , 1);

            % time of trough - raw PSTH and phase PSTH
            troughT(k,p)    = find( respRaw==min(respRaw) , 1);
            troughTP(k,p)   = find( respPR==min(respPR) , 1);
            troughTW(k,p)   = find( respW==min(respW) , 1);

            % integrated firing
            resps(k,p) = response-baseline;
            respsP(k,p) = respP-baseline;

            % response duration
            respD(k,p) = numel(find(respRaw > (std(baseline).*3)+baseline));
            respDP(k,p) = numel(find(respPR > (std(baseline).*3)+baseline));
            
%             % peak spikes relative to mean
%             varSpks(k,p) = response - trapz(PopData.session(thisSession).unit(p).respCS.psthAVG());
%             varSpksP(k,p) = response -


        end

    end



    forParvezC = latencies;

    forParvezA = [ popTraj , peaks , troughT, peakT , resps , respD , peaksP , troughTP , peakTP , respsP , respDP , peaksW , peakTW , troughTW , corrPWise , corPPWise ];
    disp('here.');
    forParvezAl= [ 1.*ones(1,size(popTraj,2)) , 2.*ones(1,size(peaks,2)) , 3.*ones(1,size(troughT,2)), 4.*ones(1,size(peakT,2)) , 5.*ones(1,size(resps,2)) , 6.*ones(1,size(respD,2)) , 7.*ones(1,size(peaksP,2)) , 8.*ones(1,size(troughTP,2)) , 9.*ones(1,size(peakTP,2)) , 10.*ones(1,size(respsP,2)) , 11.*ones(1,size(respDP,2)) , 12.*ones(1,size(peaksW,2)) , 13.*ones(1,size(peakTW,2)) , 14.*ones(1,size(troughTW,2)) , 15.*ones(1,size(corrPWise,2)) , 16.*ones(1,size(corPPWise,2)) ];


    forParvezDisp = [ popTraj , peaks , troughT, peakT , resps , respD , peaksP , troughTP , peakTP , respsP , respDP , peaksW , peakTW , troughTW];

    figure(100);
    imagesc(corr(forParvezDisp),[0 1]);

    PopAppData.session(aa).sessNum = thisSession;
    PopAppData.session(aa).numUnits = numUnits;
    PopAppData.session(aa).forParvezC = forParvezC;
    PopAppData.session(aa).forParvezA = forParvezA;
    PopAppData.session(aa).forParvezAl = forParvezAl;
    
    clear popTraj peaks troughT peakT resps respD peaksP troughTP peakTP respsP respDP peaksW peakTW troughTW corrPWise corPPWise

end
% for reference:
%     popTraj , peaks , troughT, peakT , resps , respD , peaksP , troughTP , peakTP , respsP , respDP , peaksW , peakTW , troughTW , corrPWise , corPPWise ];
% 
% 		9           30	     51      72		93     114 	     135         156      177      198      219 	 240     261      282      492         702

%% NEW COMPUTE FEATURE MATRIX FOR PARVEZ
 
% really all we need is the trialwise psths for phase, time, and the behavioral measure of latency
%         PopData.session(i).unit(j).respJCS.imagePh  = imagePh;
%         PopData.session(i).unit(j).respJCS.imageTm  = imageTm;
%         PopData.session(i).unit(j).respJCS.lat      = lat;
 
%=================================
% feature list
%=================================
% mean firing rate (raw)
% max firing rate (raw)
% min firing rate (raw)
% slope of firing rate (raw)
% time of max
% time of min
% time duration of elevated firing
% quartile loading on to PC1-3 (time)
%  
% mean firing rate (warped)
% max firing rate (warped)
% min firing rate (warped)
% slope of firing rate (warped)
% phase of max
% phase of min
% phase duration of elevated firing
% quartile loading on to PC1-3 (warped)
%  
% pairwise correlation time
% pairwise correlation phase
%=================================
 
baseWin     = [1:1e3];
timeWin     = [1e3:2.9e3]; % [1e3:1e3+latencies(k)];
pcWin1      = [1e3:1.5e3];
pcWin2      = [1.5e3:2.2e3];
pcWin3      = [2.2e3:2.9e3];
pcWin1p     = [1e3:3.1e3];
pcWin2p     = [3.1e3:5.2e3];
pcWin3p     = [5.2e3:7.25e3];
phaseWin    = [1e3:7250];
            
for aa = 1:numel(sessList)
    
    i = sessList(aa);
    disp(' ');
    disp('%=====================================');
    disp('%=====================================');
    disp(['Session ' num2str(i)]);
    disp(' ');
% for i=83;
    
    clear *meanRespMat pairwise_phase pairwise_time latencies *_firing_rate pairwise* corrInds
        
    numUnits = size(PopData.session(i).unit,2);
    numTrials = size(PopData.session(i).unit(1).respJCS.imagePh,1);  
    
    for j = 1:numUnits
        TmeanRespMat(j,:) = PopData.session(i).unit(j).respJCS.TpsthAVG(timeWin);
        WmeanRespMat(j,:) = PopData.session(i).unit(j).respJCS.WpsthAVG(phaseWin);
    end
    
    [mappedA_t, mapping_t] = compute_mapping(TmeanRespMat', 'PCA', 3);
    [mappedA_w, mapping_w] = compute_mapping(WmeanRespMat', 'PCA', 3);
    
    figure(1); clf;
    subplot(121); plot(timeWin-2e3,mappedA_t(:,1),'r',timeWin-2e3,mappedA_t(:,2),'k',timeWin-2e3,mappedA_t(:,3),'b');
    subplot(122); plot(phaseWin-2e3,mappedA_w(:,1),'r',phaseWin-2e3,mappedA_w(:,2),'k',phaseWin-2e3,mappedA_w(:,3),'b');
    
    pcompsT = mapping_t.M;
    pcompsP = mapping_w.M;  
          
    for k=1:numTrials
    
        latencies(k)    = PopData.session(i).unit(j).respJCS.lat.flOn(k);
        pc1snippet      = zeros(numUnits,numel(pcWin1));
        pc2snippet      = zeros(numUnits,numel(pcWin2));
        pc3snippet      = zeros(numUnits,numel(pcWin3));
        pc1snippetP      = zeros(numUnits,numel(pcWin1p));
        pc2snippetP      = zeros(numUnits,numel(pcWin2p));
        pc3snippetP      = zeros(numUnits,numel(pcWin3p));
        count           = 1;
            
%         disp(['Trial... ' num2str(k)]);
        
        for j = 1:numUnits
            
            baselineT   = mean(PopData.session(i).unit(j).respJCS.imageTm(k,baseWin));
            baselineP   = mean(PopData.session(i).unit(j).respJCS.imagePh(k,baseWin));
            timeWin     = [1e3:2.9e3]; % [1e3:1e3+latencies(k)];
            phaseWin    = [1e3:7250];            
 
            %=====================================
            % basic stats in time
            %=====================================
            t_base_firing_rate(k,j)     = baselineT;            
            t_mean_firing_rate(k,j)     = mean(PopData.session(i).unit(j).respJCS.imageTm(k,timeWin)) - baselineT;            
            t_max_firing_rate(k,j)      = max(PopData.session(i).unit(j).respJCS.imageTm(k,timeWin)) - baselineT;
            t_min_firing_rate(k,j)      = min(PopData.session(i).unit(j).respJCS.imageTm(k,timeWin)) - baselineT;
                                     pp = polyfit(timeWin,PopData.session(i).unit(j).respJCS.imageTm(k,timeWin),1);
            t_slope_firing_rate(k,j)    = pp(1);
            t_maxt_firing_rate(k,j)     = find(PopData.session(i).unit(j).respJCS.imageTm(k,timeWin) == max(PopData.session(i).unit(j).respJCS.imageTm(k,timeWin)) ,1);
            t_mint_firing_rate(k,j)     = find(PopData.session(i).unit(j).respJCS.imageTm(k,timeWin) == min(PopData.session(i).unit(j).respJCS.imageTm(k,timeWin)) ,1);

            %=====================================
            % basic stats in phase
            %=====================================
            p_base_firing_rate(k,j)     = baselineP;            
            p_mean_firing_rate(k,j)     = mean(PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin)) - baselineP;            
            p_max_firing_rate(k,j)      = max(PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin)) - baselineP;
            p_min_firing_rate(k,j)      = min(PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin)) - baselineP;
                                     pp = polyfit(phaseWin,PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin),1);
            p_slope_firing_rate(k,j)    = pp(1);
            p_maxt_firing_rate(k,j)     = find(PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin) == max(PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin)) ,1);
            p_mint_firing_rate(k,j)     = find(PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin) == min(PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin)) ,1);

            %=====================================
            % duration of elevated firing in time
            %=====================================
            if abs(t_min_firing_rate(k,j)) > abs(t_max_firing_rate(k,j))
                t_dur_firing_rate(k,j) = numel( find( PopData.session(i).unit(j).respJCS.imageTm(k,timeWin)<(baselineT-(0.5.*abs(t_min_firing_rate(k,j)))) ) );
            else
                t_dur_firing_rate(k,j) = numel( find( PopData.session(i).unit(j).respJCS.imageTm(k,timeWin)>(baselineT+(0.5.*abs(t_max_firing_rate(k,j)))) ) );                
            end

            %=====================================
            % duration of elevated firing in phase
            %=====================================
            if abs(p_min_firing_rate(k,j)) > abs(p_max_firing_rate(k,j))
                p_dur_firing_rate(k,j) = numel( find( PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin)<(baselineT-(0.5.*abs(p_min_firing_rate(k,j)))) ) );
            else
                p_dur_firing_rate(k,j) = numel( find( PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin)>(baselineT+(0.5.*abs(p_max_firing_rate(k,j)))) ) );                
            end

            %=====================================
            % pairwise correlations for both time and phase
            %=====================================
            for q = 1:j-1
                
                tmp = corr(PopData.session(i).unit(j).respJCS.imageTm(k,timeWin)',PopData.session(i).unit(q).respJCS.imageTm(k,timeWin)');
                if isnan(tmp)
                    pairwise_time(k,count) = 0;
                else
                    pairwise_time(k,count) = tmp;
                end
                
                tmp = corr(PopData.session(i).unit(j).respJCS.imagePh(k,phaseWin)',PopData.session(i).unit(q).respJCS.imagePh(k,phaseWin)');
                if isnan(tmp)
                    pairwise_phase(k,count) = 0;
                else
                    pairwise_phase(k,count) = tmp;
                end
                
                corrInds(count,1) = j;
                corrInds(count,2) = q;
                count = count+1;
                
            end
            
            %=====================================
            % population vector projected into low dimensions for time
            %=====================================
            pc1snippet(j,:) = PopData.session(i).unit(j).respJCS.imageTm(k,pcWin1);
            pc2snippet(j,:) = PopData.session(i).unit(j).respJCS.imageTm(k,pcWin2);
            pc3snippet(j,:) = PopData.session(i).unit(j).respJCS.imageTm(k,pcWin3);
            
            if j==numUnits
                t_pc1a_firing_rate(k,1)     = mean(pc1snippet'*pcompsT(:,1));
                t_pc1b_firing_rate(k,1)     = mean(pc1snippet'*pcompsT(:,2));
                t_pc1c_firing_rate(k,1)     = mean(pc1snippet'*pcompsT(:,3));

                t_pc2a_firing_rate(k,1)     = mean(pc2snippet'*pcompsT(:,1));
                t_pc2b_firing_rate(k,1)     = mean(pc2snippet'*pcompsT(:,2));
                t_pc2c_firing_rate(k,1)     = mean(pc2snippet'*pcompsT(:,3));

                t_pc3a_firing_rate(k,1)     = mean(pc3snippet'*pcompsT(:,1));
                t_pc3b_firing_rate(k,1)     = mean(pc3snippet'*pcompsT(:,2));
                t_pc3c_firing_rate(k,1)     = mean(pc3snippet'*pcompsT(:,3));             
            end
            
            %=====================================
            % population vector projected into low dimensions for time
            %=====================================
            pc1snippetP(j,:) = PopData.session(i).unit(j).respJCS.imagePh(k,pcWin1p);
            pc2snippetP(j,:) = PopData.session(i).unit(j).respJCS.imagePh(k,pcWin2p);
            pc3snippetP(j,:) = PopData.session(i).unit(j).respJCS.imagePh(k,pcWin3p);
            
            if j==numUnits
                p_pc1a_firing_rate(k,1)     = mean(pc1snippet'*pcompsP(:,1));
                p_pc1b_firing_rate(k,1)     = mean(pc1snippet'*pcompsP(:,2));
                p_pc1c_firing_rate(k,1)     = mean(pc1snippet'*pcompsP(:,3));

                p_pc2a_firing_rate(k,1)     = mean(pc2snippet'*pcompsP(:,1));
                p_pc2b_firing_rate(k,1)     = mean(pc2snippet'*pcompsP(:,2));
                p_pc2c_firing_rate(k,1)     = mean(pc2snippet'*pcompsP(:,3));

                p_pc3a_firing_rate(k,1)     = mean(pc3snippet'*pcompsP(:,1));
                p_pc3b_firing_rate(k,1)     = mean(pc3snippet'*pcompsP(:,2));
                p_pc3c_firing_rate(k,1)     = mean(pc3snippet'*pcompsP(:,3));             
            end
            
        end
        
        
    end
 
    params = [ t_base_firing_rate, t_mean_firing_rate, t_max_firing_rate, t_min_firing_rate, t_slope_firing_rate, t_maxt_firing_rate, t_mint_firing_rate, t_dur_firing_rate, ...
               t_pc1a_firing_rate, t_pc1b_firing_rate, t_pc1c_firing_rate, t_pc2a_firing_rate, t_pc2b_firing_rate, t_pc2c_firing_rate, t_pc3a_firing_rate, t_pc3b_firing_rate, t_pc3c_firing_rate, ...                
               p_base_firing_rate, p_mean_firing_rate, p_max_firing_rate, p_min_firing_rate, p_slope_firing_rate, p_maxt_firing_rate, p_mint_firing_rate, p_dur_firing_rate, ...
               p_pc1a_firing_rate, p_pc1b_firing_rate, p_pc1c_firing_rate, p_pc2a_firing_rate, p_pc2b_firing_rate, p_pc2c_firing_rate, p_pc3a_firing_rate, p_pc3b_firing_rate, p_pc3c_firing_rate, ...                
               pairwise_time, pairwise_phase];

    labels1 = [ 1.*ones(1,size(t_base_firing_rate,2)) ,  2.*ones(1,size(t_mean_firing_rate,2)) ,  3.*ones(1,size(t_max_firing_rate,2)),    4.*ones(1,size(t_min_firing_rate,2)) ,   5.*ones(1,size(t_slope_firing_rate,2)) , 6.*ones(1,size(t_maxt_firing_rate,2)) ,  7.*ones(1,size(t_mint_firing_rate,2)) ,  8.*ones(1,size(t_dur_firing_rate,2)), ... 
               9.*ones(1,size(t_pc1a_firing_rate,2)) , 10.*ones(1,size(t_pc1b_firing_rate,2)) , 11.*ones(1,size(t_pc1c_firing_rate,2)) , 12.*ones(1,size(t_pc2a_firing_rate,2)) , 13.*ones(1,size(t_pc2b_firing_rate,2)) , 14.*ones(1,size(t_pc2c_firing_rate,2)) , 15.*ones(1,size(t_pc3a_firing_rate,2)) , 16.*ones(1,size(t_pc3b_firing_rate,2)), 17.*ones(1,size(t_pc3c_firing_rate,2))];
    labels2 = labels1 + (max(labels1));
    labels3 = [ (max(labels2)+1) .* ones(1,size(pairwise_time,2)) , (max(labels2)+2) .* ones(1,size(pairwise_phase,2)) ];

    PopAppData.settings.windows.time        = timeWin;
    PopAppData.settings.windows.phase       = phaseWin;
    PopAppData.settings.windows.baseline    = baseWin;
    PopAppData.session(aa).sessNum          = i;
    PopAppData.session(aa).pca.timeComps    = pcompsT;
    PopAppData.session(aa).pca.phaseComps   = pcompsP;
    PopAppData.session(aa).pca.meanProjT    = mappedA_t;
    PopAppData.session(aa).pca.meanProjP    = mappedA_w;
    PopAppData.session(aa).pwc.corrInds     = corrInds;
    PopAppData.session(aa).numUnits         = numUnits;    
    PopAppData.session(aa).numTrials        = numTrials;    
    PopAppData.session(aa).latencies        = latencies;
    PopAppData.session(aa).params           = params;
    PopAppData.session(aa).labels           = [labels1 labels2 labels3];    
 
    disp(['Computed ' num2str(size(params,2)) ' parameter values from ' num2str(numUnits) ' units recorded during ' num2str(size(params,1)) ' trials of behavior.']);    
    disp('%=====================================');
    disp(' ');

    %=====================================
    % Look at the feature matrix
    %=====================================
    forExamination = [ t_base_firing_rate, t_mean_firing_rate, t_max_firing_rate, t_min_firing_rate, t_slope_firing_rate, t_maxt_firing_rate, t_mint_firing_rate, t_dur_firing_rate, ...
                       t_pc1a_firing_rate, t_pc1b_firing_rate, t_pc1c_firing_rate, t_pc2a_firing_rate, t_pc2b_firing_rate, t_pc2c_firing_rate, t_pc3a_firing_rate, t_pc3b_firing_rate, t_pc3c_firing_rate, ...                
                       p_base_firing_rate, p_mean_firing_rate, p_max_firing_rate, p_min_firing_rate, p_slope_firing_rate, p_maxt_firing_rate, p_mint_firing_rate, p_dur_firing_rate, ...
                       p_pc1a_firing_rate, p_pc1b_firing_rate, p_pc1c_firing_rate, p_pc2a_firing_rate, p_pc2b_firing_rate, p_pc2c_firing_rate, p_pc3a_firing_rate, p_pc3b_firing_rate, p_pc3c_firing_rate];               
           
    figure(100);
    [newMap] = TNC_CreateRBColormap(1024,'mbr');
    imagesc(corr(forExamination),[-1 1]);
    colormap(newMap);
    
end
 
%% CODE TO GENERATE FIGURES

%% FIGURE 1 Example of the behavior, recording configuration, and dominant response types

% Get the 
% for i=1:NumSessions
count =0;
countN = 0;
countA = 0;

for aa = 1:numel(sessList)
  
        i = sessList(aa);
        numUnits = numel(PopData.session(i).unit);
        
        for j=1:numUnits
            
            countA = countA+1;
            
            if numel(PopData.session(i).unit(j).WfMean)>30
                disp('Waveform found.');                
                count = count+1;
                
                % get the spike width
                minV = min(PopData.session(i).unit(j).WfMean);
                fwhm = numel(find( abs(PopData.session(i).unit(j).WfMean) > (-minV/.5))) ./ 30;
                summaryData(count,1) = fwhm;
                
%                 summaryData(count,1) = trapz( abs(PopData.session(i).unit(j).WfMean)) ./ ( -min(PopData.session(i).unit(j).WfMean) );
                
                
                % get the average firing rate
                summaryData(count,2) = 1000/PopData.session(i).unit(j).isi.stats.mean;
                
                % get the fano factor
                summaryData(count,3) = PopData.session(i).unit(j).isi.stats.fano;
                
                % get daLogic
                summaryData(count,4) = PopData.session(i).unit(j).daLogic;
                
            else
                disp('Waveform not found.');
                countN = countN+1;
            end
        end
        
end

disp(['Found waveforms for ' num2str(count) ' units; Not found for ' num2str(countN) ' units.']);
figure(2);
DaInds = find(summaryData(:,4)==1);
plot(summaryData(:,1),summaryData(:,2),'ko',summaryData(DaInds,1),summaryData(DaInds,2),'r.');
xlabel('Spike FWHM (ms)'); ylabel('Mean Firing Rate (Hz)');

%% FIGURE 2 Warping / Tiling pseudopopulation response

% get the session mapping 
for i=1:69
    sessMap(i,1) = i;
    sessMap(i,2) = PopAppData.session(i).sessNum;
end

% 3 example cells (early middle and late) for the effect of warping (raster, psth, behavioral warping effects)
% good middle cells: 61 101
% candidate sessions to look at (manually annotated):
exampleCells = ...
[   71    12;
    71    13;
    73     6;
   101     4;
   101     8;
    89     2;
    89     6;
    89    16;
    78     9;
    78    16;
    85    11;
    85    12;
    87     8;
    87    16;
   181    46;
   107    12 ]


% for an example session get the time and phase rasters and compute the warped 
chosenOnes =     [1 3 6];

for q = 1:numel(chosenOnes)
% for q = 1:3

    [colors] = TNC_CreateRBColormap(numel(chosenOnes),'cpb');

    i = exampleCells(chosenOnes(q),1);
    j = exampleCells(chosenOnes(q),2);
    PopData.session(i).trials = size(PopData.session(i).events.CS.ts,1);

    numStamps = length(PopData.session(i).unit(j).ts);
    delta = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
    delta(1,round(PopData.session(i).unit(j).ts)+1) = 1;
    tmpSmooth = conv(delta,currParams.filter.kernel,'same');

    [respJCS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[1e3,6e3],1,1);
    
    [tN,tI] = sort(PopData.session(i).behavior.flCSon,'ascend');
    tmp = find(tN<5e3 & tN>750);
    avgCS2FL = mean(tN(tmp))
    validTrials = tI(tmp);
    
%     phaseDecay = floor( currParams.smthParams.decay .* (6280/avgCS2FL) ./ 3 );
    phaseDecay = floor( currParams.smthParams.decay .* (1570/avgCS2FL) ./ 3 );
    [phasekernel]  = TNC_CreateGaussian(phaseDecay.*15,phaseDecay,phaseDecay.*30,1);
    
    clear lat
        
    allPhaseStamps  = [];
    allPhaseRows    = [];
    allStamps   = [];
    allRows     = [];
    allCSpX     = [];
    allCSpY     = [];
    lat.flOn    = [];
    lat.flOff   = [];
    lat.preCS   = [];
%     deltaP      = zeros(1,6280+2e3);
%     imagePh     = zeros(numel(validTrials) , 6280+2e3);
    deltaP      = zeros(1,1570+2e3);
    imagePh     = zeros(numel(validTrials) , 1570+2e3);
    imageTm     = zeros(numel(validTrials) , numel(PopData.session(i).unit(j).respJCS.image(1,:)));

    for kk=1:numel(validTrials)

%         deltaP      = zeros(1,6280+2e3);
        deltaP      = zeros(1,1570+2e3);
        currTrial   = validTrials(kk);        
        theseStamps = PopData.session(i).unit(j).respJCS.raster.trial(currTrial).ts';
        theseRows   = kk + ones(1,numel(theseStamps));

        allStamps   = [ allStamps theseStamps ];
        allRows     = [ allRows theseRows ];
        allCSpX     = [ allCSpX PopData.session(i).behavior.flCSon(currTrial) ];
        allCSpY     = [ allCSpY , kk+1 ];

        % calculate warped stamps
        validStamps = find(theseStamps>0 & theseStamps<PopData.session(i).behavior.flCSon(currTrial));
        validStamps2= find(theseStamps>PopData.session(i).behavior.flCSon(currTrial) & theseStamps<PopData.session(i).behavior.flCSon(currTrial)+1000);
        validStamps3= find(theseStamps<0);

%         PopData.session(i).unit(j).respJCS.raster.trial(currTrial).ph = [   theseStamps(validStamps3) , ...                                                                                
%                                                                             (theseStamps(validStamps) .* (6280 ./ PopData.session(i).behavior.flCSon(currTrial))) , ...
%                                                                             (theseStamps(validStamps2)-PopData.session(i).behavior.flCSon(currTrial)) + 6280];

        PopData.session(i).unit(j).respJCS.raster.trial(currTrial).ph = [   theseStamps(validStamps3) , ...                                                                                
                                                                            (theseStamps(validStamps) .* (1570 ./ PopData.session(i).behavior.flCSon(currTrial))) , ...
                                                                            (theseStamps(validStamps2)-PopData.session(i).behavior.flCSon(currTrial)) + 1570];

        allPhaseStamps  = [ allPhaseStamps PopData.session(i).unit(j).respJCS.raster.trial(currTrial).ph ];                                                                        
        thesePhRows     = kk + ones(1,numel(PopData.session(i).unit(j).respJCS.raster.trial(currTrial).ph));
        allPhaseRows    = [ allPhaseRows thesePhRows ];
                                                                        
%         deltaP( ceil(theseStamps(validStamps3)) + 1e3 ) = 1e3 / 6280;
%         deltaP( ceil((theseStamps(validStamps) .* (6280 ./ PopData.session(i).behavior.flCSon(currTrial)))) + 1e3 ) = (avgCS2FL / PopData.session(i).behavior.flCSon(currTrial));
%         deltaP( ceil((theseStamps(validStamps2)-PopData.session(i).behavior.flCSon(currTrial)) + 7280) ) = 1e3 / 6280;

        deltaP( ceil(theseStamps(validStamps3)) + 1e3 ) = 1e3 / 1570;
        deltaP( ceil((theseStamps(validStamps) .* (1570 ./ PopData.session(i).behavior.flCSon(currTrial)))) + 1e3 ) = (avgCS2FL / PopData.session(i).behavior.flCSon(currTrial));
        deltaP( ceil((theseStamps(validStamps2)-PopData.session(i).behavior.flCSon(currTrial)) + 2570) ) = 1e3 / 1570;

        imagePh(kk,:)   = conv(deltaP,phasekernel,'same');
        imageTm(kk,:)   = PopData.session(i).unit(j).respJCS.image(currTrial,:);
        lat.flOn(kk)    = PopData.session(i).behavior.flCSon(currTrial);
        lat.flOff(kk)   = PopData.session(i).behavior.flCSoff(currTrial);
        lat.preCS(kk)   = PopData.session(i).timeSinceCS(currTrial);
        

    end
        WpsthAVG = mean(imagePh,1) - mean(mean(imagePh(:,1:1e3),1));
        WpsthSEM = std(imagePh,[],1) ./ sqrt(numel(validTrials)-1);

        TpsthAVG = mean(imageTm,1) - mean(mean(imageTm(:,1:1e3),1));
        TpsthSEM = std(imageTm,[],1) ./ sqrt(numel(validTrials)-1);

        figure(210+q); clf;
        subplot(2,5,[1 2 ]);
        plot(allStamps,allRows,'.','Color',colors(q,:),'MarkerSize',1); hold on; 
        plot(allCSpX,allCSpY,'k-');
        plot([0 0],[1 max(allCSpY)],'k-');
%         plot([0 0] , [0 PopData.session(i).trials] , 'k');
        axis([-1e3 3e3 1 numel(validTrials)]); set(gca,'TickDir','out'); box off;
        xlabel('Time from CS (ms)'); ylabel('Sorted trials'); title(['id = s' num2str(sessMap(find(sessMap(:,2)==i),1)) 'u' num2str(j)]);

        subplot(2,5,[6 7 ]);
        plot(allPhaseStamps,allPhaseRows,'.','Color',colors(q,:),'MarkerSize',1); hold on; 
%         plot((allCSpX.*0)+6280,allCSpY,'k-');
%         plot([0 0],[min(allCSpY) max(allCSpY)],'k-');
%         plot([0 0] , [0 PopData.session(i).trials] , 'k');
%         axis([-1e3 6e3 0 PopData.session(i).trials]);
        plot([1570 1570],[1 max(allCSpY)],'k-');
        plot([0 0],[1 max(allCSpY)],'k-');
        axis([-1e3 2570 1 numel(validTrials)]); set(gca,'TickDir','out'); box off; 
        xlabel('Phase from CS (mrad)'); ylabel('Sorted trials'); 

        subplot(2,5,[9 10]);
        shadedErrorBar(-999:2570,WpsthAVG.*1000,WpsthSEM.*1000,{'Color',colors(q,:)}); hold on;
%         plot([-1e3 7280] , [0 0] , 'k','Color',[0.75 0.75 0.75]);
        axis([-1e3 2570 min(TpsthAVG.*1000)-5 max(WpsthAVG.*1000)+10]);
        plot([-1e3 2570], [max(TpsthAVG.*1000) max(TpsthAVG.*1000)],'k--','color',[0.75 0.75 0.75]);
        ylabel('Norm. Firing Rate (sp/rad)'); set(gca,'TickDir','out'); box off;
        xlabel('Phase from CS (mrad)');
%         plot([6280 6280], [0 2],'k--','Color',[0.25 0.25 0.25]);
         
        subplot(2,5,[4 5]);
        shadedErrorBar(-1e3:6e3,TpsthAVG.*1000,TpsthSEM.*1000,{'Color',colors(q,:)}); hold on;
%         axis([-1e3 7280 0 2]);
        axis([-1e3 6e3 min(TpsthAVG.*1000)-5 max(WpsthAVG.*1000)+10]);
        ylabel('Norm. Firing Rate (sp/sec)'); set(gca,'TickDir','out'); box off;
        xlabel('Time from CS (ms)');
%         plot([0 0], [0 2],'k--');
%         plot([avgCS2FL avgCS2FL], [0 2],'k--','Color',[0.25 0.25 0.25]);
%     end
end

%         
%         PopData.session(i).unit(j).respJCS.imagePh  = imagePh;
%         PopData.session(i).unit(j).respJCS.imageTm  = imageTm;
%         PopData.session(i).unit(j).respJCS.lat      = lat;
%         
%         PopData.session(i).unit(j).respJCS.WpsthAVG = mean(imagePh,1) - mean(mean(imagePh(:,1:1e3),1));
%         PopData.session(i).unit(j).respJCS.WpsthSEM = std(imagePh,[],1) ./ sqrt(numel(validTrials)-1);
% 
%         PopData.session(i).unit(j).respJCS.TpsthAVG = mean(imageTm,1) - mean(mean(imageTm(:,1:1e3),1));
%         PopData.session(i).unit(j).respJCS.TpsthSEM = std(imageTm,[],1) ./ sqrt(numel(validTrials)-1);
%         
%         if abs(max(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280))) > abs(min(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)))
%             PopData.WpsthAVG(count,:) = PopData.session(i).unit(j).respJCS.WpsthAVG ./ abs(max(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)));
%             PopData.WpsthSEM(count,:) = PopData.session(i).unit(j).respJCS.WpsthSEM ./ abs(max(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)));
% 
%             PopData.TpsthAVG(count,:) = PopData.session(i).unit(j).respJCS.TpsthAVG ./ abs(max(PopData.session(i).unit(j).respJCS.TpsthAVG(1e3:1e3+avgCS2FL)));
%             PopData.TpsthSEM(count,:) = PopData.session(i).unit(j).respJCS.TpsthSEM ./ abs(max(PopData.session(i).unit(j).respJCS.TpsthAVG(1e3:1e3+avgCS2FL)));
%             
%             PopData.timeOfPeak(count,:) = find( PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)==max(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)) , 1 ) .* 1;
% 
%             PopData.ampOfPeak(count,1) = (max(mean(imagePh(:,1e3:7280),1)) - mean(mean(imagePh(:,1:1e3)))) ./ (max(mean(imagePh(:,1e3:7280),1)) + mean(mean(imagePh(:,1:1e3))));
%             PopData.ampOfPeak(count,2) = (max(mean(imageTm(:,1e3:3e3),1)) - mean(mean(imageTm(:,1:1e3)))) ./ (max(mean(imageTm(:,1e3:3e3),1)) + mean(mean(imageTm(:,1:1e3))));
% 
%         else
%             PopData.WpsthAVG(count,:) = PopData.session(i).unit(j).respJCS.WpsthAVG ./ abs(min(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)));
%             PopData.WpsthSEM(count,:) = PopData.session(i).unit(j).respJCS.WpsthSEM ./ abs(min(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)));
% 
%             PopData.TpsthAVG(count,:) = PopData.session(i).unit(j).respJCS.TpsthAVG ./ abs(min(PopData.session(i).unit(j).respJCS.TpsthAVG(1e3:1e3+avgCS2FL)));
%             PopData.TpsthSEM(count,:) = PopData.session(i).unit(j).respJCS.TpsthSEM ./ abs(min(PopData.session(i).unit(j).respJCS.TpsthAVG(1e3:1e3+avgCS2FL)));
% 
%             PopData.timeOfPeak(count,:) = find( PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)==min(PopData.session(i).unit(j).respJCS.WpsthAVG(1e3:7280)) , 1 ) .* -1;
%             
%             PopData.ampOfPeak(count,1) = (min(mean(imagePh(:,1e3:7280),1)) - mean(mean(imagePh(:,1:1e3)))) ./ (abs(min(mean(imagePh(:,1e3:7280),1))) + mean(mean(imagePh(:,1:1e3))));
%             PopData.ampOfPeak(count,2) = (min(mean(imageTm(:,1e3:3e3),1)) - mean(mean(imageTm(:,1:1e3)))) ./ (abs(min(mean(imageTm(:,1e3:3e3),1))) + mean(mean(imageTm(:,1:1e3))));
%         end
%         
%         count = count+1;
    
% figure(3); clf;
% plot(InitRespLatency,warpContrast,'ko'); hold on;
% tmp = find(warpContrast>25);
% for i=1:numel(tmp)
%     text(InitRespLatency(tmp(i))+5,warpContrast(tmp(i))+5,[num2str(tmp(i))]);
% end

% Population data showing the warping effect
% 1) Effect of warping on distribution of peak amplitudes
% 2) Effect of warping on the distribution of peak latencies
% 3) Precision effects of warping?

%% Code to export time raster, phase raster, and phase PSTH with sorted trials

export = 0; j=1;

        % unit
        m = 20;
        % session
        i = 83;

        disp(['Session ' num2str(i)]);

        validTrials = PopData.session(i).behavior.validPhaseTrials;
        numTrials = numel(validTrials);

        numUnits = 1;

        disp('Sorting by latency...');
        yLabelStr = '(sorted)';
        latencies = PopData.session(i).behavior.flCSon(validTrials);
        [vals,trialIndsTmp] = sort(latencies);
        trialInds = validTrials(trialIndsTmp);

            delta = zeros(1,7300);
            allStamps = []; allRows = []; allIDs = [];
            thisUnit = m;


            figure(1); clf;
            subplot(311);
            
            for k = 1:numTrials

                        thisTrialRow    = ((m-1).*numTrials) + (5.*m) + k;
                        theseTimeStamps = [PopData.session(i).unit(thisUnit).allPhase.seg(1).trial(trialIndsTmp(k)).ts,
                                           PopData.session(i).unit(thisUnit).allPhase.seg(2).trial(trialIndsTmp(k)).ts,
                                           PopData.session(i).unit(thisUnit).allPhase.seg(3).trial(trialIndsTmp(k)).ts,
                                           PopData.session(i).unit(thisUnit).allPhase.seg(4).trial(trialIndsTmp(k)).ts];

                       if numel(theseTimeStamps)>0
                           tmp = find(theseTimeStamps<6300 & theseTimeStamps>-1e3);

                           delta(ceil(theseTimeStamps(tmp)+1000)) = delta(ceil(theseTimeStamps(tmp)+1000)) + 1;

                           if rem(m,2)
                                subplot(5,1,2:5); hold on;
                                plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[0 0.67 1]);
                           else
                                subplot(5,1,2:5);  hold on;
                                plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[1 0.67 0]);            
                           end                            

                           plot(zeros(1,numTrials),((m-1).*numTrials) + (5.*m) + [1:numTrials],'b');
                           plot(PopData.session(i).behavior.flCSon(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');
                           plot(PopData.session(i).behavior.USon(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'r');
                           plot(PopData.session(i).behavior.flCSoff(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');

                           if export
                            allStamps   = [allStamps theseTimeStamps'];
                            allRows     = [allRows ones(1,numel(theseTimeStamps)).*thisTrialRow]; 
                            allIDs      = [allIDs ones(1,numel(theseTimeStamps)).*m./numUnits];
                           end

                       end
                    end

                    
                    if export
                        nameA = sprintf('/spikeTimesX%g',j-1);
                        nameB = sprintf('/spikeTimesY%g',j-1);
                        nameC = sprintf('/cellIDs%g',j-1);
                        hdf5write(['MatlabExport/egHahnStyleTime_' num2str(i) '_' num2str(m)  '.h5'],nameA,allStamps);
                        hdf5write(['MatlabExport/egHahnStyleTime_' num2str(i) '_' num2str(m)  '.h5'],nameB,allRows,'WriteMode','append');
                        hdf5write(['MatlabExport/egHahnStyleTime_' num2str(i) '_' num2str(m)  '.h5'],nameC,allIDs,'WriteMode','append');

                        nameD = sprintf('/psth%g',j-1);
                        hdf5write(['MatlabExport/egHahnStyleTime_' num2str(i) '_' num2str(m)  '.h5'],nameC,PopData.session(i).unit(thisUnit).allPhase.timePSTH.psth,'WriteMode','append');
                    end


            subplot(312);

                    allStamps = []; allRows = []; allIDs = [];                        
                    delta = zeros(1,7300);

                    for k = 1:numTrials

                        thisTrialRow    = ((m-1).*numTrials) + (5.*m) + k;
                        theseTimeStamps = [PopData.session(i).unit(thisUnit).allPhase.seg(1).trial(trialIndsTmp(k)).ph,
                                           PopData.session(i).unit(thisUnit).allPhase.seg(2).trial(trialIndsTmp(k)).ph,
                                           PopData.session(i).unit(thisUnit).allPhase.seg(3).trial(trialIndsTmp(k)).ph,
                                           PopData.session(i).unit(thisUnit).allPhase.seg(4).trial(trialIndsTmp(k)).ph];

                        tmp = find(theseTimeStamps>-1);          
                        delta(ceil((theseTimeStamps(tmp)+1).*1000)) = delta(ceil((theseTimeStamps(tmp)+1).*1000)) + 1;

                        if rem(m,2)
                            subplot(5,1,2:5); hold on;
                            plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[0 0.67 1]);
                        else
                            subplot(5,1,2:5);  hold on;
                            plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[1 0.67 0]);            
                        end

                        if export
                            allStamps   = [allStamps theseTimeStamps'];
                            allRows     = [allRows ones(1,numel(theseTimeStamps)).*thisTrialRow]; 
                            allIDs      = [allIDs ones(1,numel(theseTimeStamps)).*m./numUnits];
                        end

                        plot(zeros(1,numTrials),((m-1).*numTrials) + (5.*m) + [1:numTrials],'b');
                        plot(ones(1,numTrials).*0.718,((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');
                        plot(ones(1,numTrials).*1.795,((m-1).*numTrials) + (5.*m) + [1:numTrials],'r');
                        plot(ones(1,numTrials).*6.283,((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');

                    end

                        
                    if export
                        nameA = sprintf('/spikeTimesX%g',j-1);
                        nameB = sprintf('/spikeTimesY%g',j-1);
                        nameC = sprintf('/cellIDs%g',j-1);
                        hdf5write(['MatlabExport/egHahnStylePhase_' num2str(i) '_' num2str(m) '.h5'],nameA,allStamps);
                        hdf5write(['MatlabExport/egHahnStylePhase_' num2str(i) '_' num2str(m)  '.h5'],nameB,allRows,'WriteMode','append');
                        hdf5write(['MatlabExport/egHahnStylePhase_' num2str(i) '_' num2str(m)  '.h5'],nameC,allIDs,'WriteMode','append');
                        
                        nameD = sprintf('/psth%g',j-1);
                        hdf5write(['MatlabExport/egHahnStylePhase_' num2str(i) '_' num2str(m)  '.h5'],nameC,PopData.session(i).unit(thisUnit).allPhase.phasePSTH.psth,'WriteMode','append');
                    end
                    
                    
                    subplot(313);
                    plot(PopData.session(i).unit(thisUnit).allPhase.phasePSTH.psth,'r',PopData.session(i).unit(thisUnit).allPhase.timePSTH.psth,'k');

%% FIGURE 3 Schematic of alternative predictions of the gating and continuous control models

Predictions of a rate/phase code

sparsity vs. density in time
sparsity vs. density of neuronal responses

%% FIGURE 4 - FaRM works to explain behavior

% first, get the mapping of session number to raw session number
% 
% load('/Users/dudmanj/Dropbox/JD_Timing_Dataset/attempt4/ForParvez.mat');
% 
% for i=1:numel(PopAppData.session)
%     sessId(i,1) = i;
%     sessId(i,2) = PopAppData.session(i).sessNum;    
% end
% 
% clear PopAppData;
% 
% load('/Users/dudmanj/Dropbox/JD_Timing_Dataset/attempt4/analysisResults_att4_2014_06_05/PopAppData_processed.mat');

% Need a plot to show quality of prediction and convey the power of a trial by trial prediction
cMap = TNC_CreateRBColormap(7,'cpb');
sessNumber = 36
figure(300+sessNumber); clf;
for i=sessNumber
% for i=1:numel(PopAppData.session)


    subplot(5,1,1:3);
    numTrials = numel(PopAppData.session(i).latencies);
    
    % rough idea is a raster plot of measured latencies and overlaid predictions for each of the 7 levels
    for kk=1:7
        plot(PopAppData.session(i).latencies' - PopAppData.session(i).ErrorStats(:,kk), [1:numTrials].*1 , 'ko' , 'color' , cMap(kk,:) , 'MarkerSize' , kk.*2); hold on; 
%         plot(PopAppData.session(i).ErrorStats(:,kk), [1:numTrials].*1 , 'ko--' , 'color' , cMap(kk,:) , 'MarkerSize' , kk.*2); hold on; 
    end
    plot(PopAppData.session(i).latencies , [1:numTrials].*1 , 'k-', 'linewidth' , 2); 
    ylabel('Sorted Trials'); set(gca,'TickDir','out'); box off; axis([-100 5000 0 72]);
%     ylabel('Sorted Trials'); set(gca,'TickDir','out'); box off; axis([-2000 2000 0 72]);
    
    
    subplot(5,1,4:5);
    histXvals = 10.^[0:0.15:4];
    histXvals = 0:150:4050;
    latHist = hist(PopAppData.session(i).latencies,histXvals);
    bar(histXvals,latHist,'k'); hold on;
    for kk=1:7
        latHist = hist(PopAppData.session(i).latencies' - PopAppData.session(i).ErrorStats(:,kk),histXvals);
        plot(histXvals , latHist , 'ko' , 'color' , cMap(kk,:) , 'markersize' , kk*2 , 'linewidth' , 1);
    end    
    ylabel('Count'); xlabel('Latency to approach (ms)'); set(gca,'TickDir','out'); box off; axis([-100 5000 -1 max(latHist)+3]);
end

%% FIGURE 4 - Code to look at population prediction quality

ErrPerFeature = zeros(numel(PopAppData.session),7);
FeatureTested = zeros(numel(PopAppData.session),7);
FeaturePerNeuron = zeros(numel(PopAppData.session),2);

tmp = unique(CellCnt);
cntColors = TNC_CreateRBColormap(numel(tmp),'cpb');

figure(1); clf;
for i=1:numel(PopAppData.session)
    
   ErrPerFeature(i,1:numel(PopAppData.session(i).Error_MeanAbsError)) = PopAppData.session(i).Error_MeanAbsError;
   FeatureTested(i,1:numel(PopAppData.session(i).Error_MeanAbsError)) = 1;
   if numel(PopAppData.session(i).Error_MeanAbsError) < 7
       disp(['Session ' num2str(i) ' has only ' num2str(PopAppData.session(i).numElements(numel(PopAppData.session(i).Error_MeanAbsError))) ' features.']);
   end
   
   CellCnt(i) = PopAppData.session(i).numUnits;
   figure(1); plot(PopAppData.session(i).numElements , PopAppData.session(i).Error_MeanAbsError , 'ko' , 'markersize', find(tmp==CellCnt(i)), 'Color' , cntColors(find(tmp==CellCnt(i)),:)); hold on;
%    if PopAppData.session(i).numUnits<20
%    else
%        CellCnt(i) = 21;
%    end

    currFeatIdx = PopAppData.session(i).topKfeaturesIdx(1:60);
    counts      = hist(currFeatIdx,1:1:40);
    freqPerUnit = counts([1:8 18:25]) ./ PopAppData.session(i).numUnits;
    FeaturePerNeuron(i,1) = mean(freqPerUnit);   
    FeaturePerNeuron(i,2) = std(freqPerUnit); 
    
end

for i=1:21
        figure(1); plot(40 + sum([1:i])./2 , 1100 , 'ko' , 'markersize', i , 'Color' , cntColors(i,:));
end
text(68 , 1150 , 'Unit count (rank)');

for i=1:7
    tmp = find(FeatureTested(:,i)==1);
    medErrPerFeature(i) = median(ErrPerFeature(tmp,i));     
end

figure(1);
plot([2 10 20 40 60 80 120] , medErrPerFeature , 'k+-' , 'LineWidth' , 2 , 'markersize', 10);
plot([-10 130], [0 0], 'k--');
axis([-10 130 -100 1300]);
set(gca,'TickDir','out'); box off;
ylabel('Absolute Error (ms)');
xlabel('Number of Features');

figure(2)
boxplot(ErrPerFeature,'colors','k');
set(gca,'TickDir','out'); box off;

figure(3); clf;

subplot(131);
tmp = find(FeatureTested(:,4)==1);
plot(CellCnt(tmp),ErrPerFeature(tmp,4),'ko','color',[0.6 0.6 0.6]); hold on;
[binnedData] = TNC_BinAndMean(CellCnt(tmp),ErrPerFeature(tmp,4),12);
plot(binnedData.bins.center,binnedData.bins.avg,'k-+', 'LineWidth' , 2);
plot([0 25],[0 0],'k--');
title(['Features used: ' num2str(PopAppData.session(i).numElements(4))]);
xlabel('Number of Units'); ylabel('|error| (ms)');
axis([0 25 -50 600]);
set(gca,'TickDir','out'); box off;

subplot(132);
tmp = find(FeatureTested(:,5)==1);
plot(CellCnt(tmp),ErrPerFeature(tmp,5),'ko','color',[0.6 0.6 0.6]); hold on;
[binnedData] = TNC_BinAndMean(CellCnt(tmp),ErrPerFeature(tmp,5),12);
plot(binnedData.bins.center,binnedData.bins.avg,'k-+', 'LineWidth' , 2);
plot([0 25],[0 0],'k--');
title(['Features used: ' num2str(PopAppData.session(i).numElements(5))]);
xlabel('Number of Units'); ylabel('|error| (ms)');
axis([0 25 -34 400]);
set(gca,'TickDir','out'); box off;

subplot(133);
tmp = find(FeatureTested(:,6)==1);
plot(CellCnt(tmp),ErrPerFeature(tmp,6),'ko','color',[0.6 0.6 0.6]); hold on;
[binnedData] = TNC_BinAndMean(CellCnt(tmp),ErrPerFeature(tmp,6),12);
plot(binnedData.bins.center,binnedData.bins.avg,'k-+', 'LineWidth' , 2);
plot([0 25],[0 0],'k--');
title(['Features used: ' num2str(PopAppData.session(i).numElements(6))]);
xlabel('Number of Units'); ylabel('|error| (ms)');
axis([0 25 -25 300]);
set(gca,'TickDir','out'); box off;

% subplot(313);
% tmp = find(FeatureTested(:,7)==1);
% plot(CellCnt(tmp),ErrPerFeature(tmp,7),'ko','color',[0.6 0.6 0.6]); hold on;
% [binnedData] = TNC_BinAndMean(CellCnt(tmp),ErrPerFeature(tmp,7),12);
% plot(binnedData.bins.center,binnedData.bins.avg,'k+');
% axis([0 60 0 400]);

figure(4);
bar(0:0.025:0.5,hist(FeaturePerNeuron(:,1),0:0.025:0.5),'k');
axis([-0.05 0.4 -1 14]);
ylabel('Count'); xlabel('Frac. units / feature')
set(gca,'TickDir','out'); box off;

%% FIGURE 5 Distribution of features used (i.e. test of prediction 1 in fig. 3)

% pre-processing to get the summary distribution of feature frequencies
featVals = [1:17 19:37];
clear featureHist

for j=1:7
    featureHist.topK(j).vals=[];
end

for i=1:69
    for j=2:numel(PopAppData.session(i).numElements)
        cnts = hist(PopAppData.session(i).topKfeatureLabels(1:PopAppData.session(i).numElements(j)) , featVals);
        if size(featureHist.topK(j).vals,1)>0
            featureHist.topK(j).vals = [featureHist.topK(j).vals ; cnts ./ PopAppData.session(i).numUnits];
        else
            featureHist.topK(j).vals = cnts ./ PopAppData.session(i).numUnits;
        end
    end
end

% for a given session plot the distribution of feature types used to
% provide assymptotic descrpition
figure(501); clf;
assymFeats = 7
cntColors = TNC_CreateRBColormap(assymFeats-1,'cpb');
sumColors = TNC_CreateRBColormap(count,'cpb');
count = 0;

for i=36
    
    subplot(151);
    shadedErrorBar(PopAppData.session(i).numElements,mean(abs(PopAppData.session(i).ErrorStats),1),std(abs(PopAppData.session(i).ErrorStats),[],1)./(sqrt(size(PopAppData.session(i).ErrorStats,1))-1),'k');
    hold on;
    plot([0 120] , [0 0], 'k--');
    ylabel('|Error| (ms)'); xlabel('Training Features (count)'); title(['Example Session (id=' num2str(i) ')']);
    set(gca,'TickDir','out'); box off; axis([0 120 -50 500]);

    subplot(1,5,[2:3]);
    for j=assymFeats:-1:2
        cnts = hist(PopAppData.session(i).topKfeatureLabels(1:PopAppData.session(i).numElements(j)) , featVals);
        plot( featVals , cnts , 'k-' , 'linewidth' , 2 , 'Color' , cntColors(j-1,:)); hold on;
    end
    plot([1 numel(featVals)] , [PopAppData.session(i).numUnits , PopAppData.session(i).numUnits] , 'k--');
    axis([0 numel(featVals)+1 -1 PopAppData.session(i).numUnits+3]);
    ylabel('Training Features (count)'); xlabel('Feature Category'); title(['Example Session (id=' num2str(i) ')']);
    set(gca,'TickDir','out'); box off;
    
%     subplot(1,5,[4:5]);
%     for j=7:-1:2
%         cnts = mean(featureHist.topK(j).vals,1);
%         plot( featVals , cnts , 'k' , 'Color' , cntColors(j-1,:), 'linewidth' , 2 ); hold on;
%     end
%     plot([1 max(featVals)] , [1 1] , 'k--');
%     ylabel('Training Features / Unit'); xlabel('Feature Category'); title(['All Sessions']);
%     axis([0 max(featVals)+1 -0.1 1.5]);
%     set(gca,'TickDir','out'); box off;
    
    subplot(1,5,[4:5]);
    plot([0 122] , [1 1] , 'k--'); hold on;
    for kk=numel(featVals):-1:1
        k = featVals(kk);
        for j=2:7
            sums(k,j) = mean(featureHist.topK(j).vals(:,kk),1);
            errs(k,j) = std(featureHist.topK(j).vals(:,kk),[],1) ./ sqrt(size(featureHist.topK(j).vals,1)-1);
        end
        
        if max(sums(k,:)) > 0.3
            count = count+1;
%             shadedErrorBar( [2 10 20 40 60 80 120] , sums(k,:) , errs(k,:) , {'Color' , sumColors(count,:)});
            plot( [2 10 20 40 60 80 120] , sums(k,:) , 'k-' , 'Color' , sumColors(count,:)); % for creating a legend
            validFeats(count,1) = k;
        else
%             plot( [2 10 20 40 60 80 120] , sums(k,:) , 'k-' , 'color' , [0.8 0.8 0.8]);
        end
    end
    ylabel('Fraction of units contributing feature'); xlabel('Training Features (count)'); title(['All Sessions']);
    axis([0 122 -0.1 1.5]);
    set(gca,'TickDir','out'); box off;
    
end

%% FIGURE 6 Distribution of cells that donate a feature; number of cells required

%% FIGURE 7 All the way back to a raster plot of an example session and an explanation of the main result

%% Supp Figs

Recording locations
Distribution of firing rates and comparison with in vitro distribution

%% DOES INTRINSIC CONNECTIVITY PREDICT PHASE OFFSETS?

allCompares = [];

for zz=1:numel(sessList)

    thisSess = sessList(zz);
    plotOnline = 1;
    clear compareCorr;
    % how many units in this session?
    numUnits = numel(PopData.session(thisSess).unit)
    figure(1); clf;

    % Need to look at the pairwise xcorr (only use units recorded on different electrodes)
    for i=1:numUnits

        if numel(PopData.session(thisSess).unit(i).name)==7
            elecNums(i) = str2num(PopData.session(thisSess).unit(i).name(5:6));
        else
            elecNums(i) = str2num(PopData.session(thisSess).unit(i).name(5));
        end

    end

    count = 1;

    for i=1:numUnits

        thisElec = elecNums(i);

        for j=i:numUnits
           if thisElec~=elecNums(j)
               compareCorr.WPSTHcorr(count) = corr(PopData.session(thisSess).unit(i).respJCS.WpsthAVG(1000:7200)',PopData.session(thisSess).unit(j).respJCS.WpsthAVG(1000:7200)');
               compareCorr.WPSTHpkoff(count) = find(PopData.session(thisSess).unit(i).respJCS.WpsthAVG(1000:7200) == max(PopData.session(thisSess).unit(i).respJCS.WpsthAVG(1000:7200)) , 1) - find(PopData.session(thisSess).unit(j).respJCS.WpsthAVG(1000:7200) == max(PopData.session(thisSess).unit(j).respJCS.WpsthAVG(1000:7200)) , 1);

               if max(PopData.session(thisSess).unit(i).ts) >= max(PopData.session(thisSess).unit(j).ts)
                   delta1 = zeros(1,max(ceil(PopData.session(thisSess).unit(i).ts)));
                   delta2 = zeros(1,max(ceil(PopData.session(thisSess).unit(i).ts)));
                   deltaP = zeros(1,max(ceil(PopData.session(thisSess).unit(i).ts)));
               else
                   delta1 = zeros(1,max(ceil(PopData.session(thisSess).unit(j).ts)));
                   delta2 = zeros(1,max(ceil(PopData.session(thisSess).unit(j).ts)));
                   deltaP = zeros(1,max(ceil(PopData.session(thisSess).unit(j).ts)));
               end
               delta1(ceil(PopData.session(thisSess).unit(i).ts)) = 1;
               delta2(ceil(PopData.session(thisSess).unit(j).ts)) = 1;
    %                 p = randperm(numel(PopData.session(thisSess).unit(j).ts)-1);
    %                 p = p+1;
    %                 newArray(1) = PopData.session(thisSess).unit(j).ts(1);
    %                 for k=2:numel(PopData.session(thisSess).unit(j).ts)-1
    %                     newArray(k) = newArray(k-1) + ( PopData.session(thisSess).unit(j).ts(p(k)) - PopData.session(thisSess).unit(j).ts(p(k)-1) );
    %                 end
    %            deltaP(ceil(newArray)) = 1;

               compareCorr.rawXC(count,:) = xcorr(delta1,delta2,25,'unbiased');
               compareCorr.rawXC(count,:) = compareCorr.rawXC(count,:) - mean(compareCorr.rawXC(count,:));

                if plotOnline
                    figure(1); 
                    subplot(numUnits,numUnits,sub2ind([numUnits numUnits],i,j));
                    plot(-25:25,compareCorr.rawXC(count,:),'r',[-25 25], [0 0],'k--',[0 0],[-25 25], 'k--');
                    axis([-25 25 -5e-5 5e-5]); axis off; title(num2str(compareCorr.WPSTHpkoff(count)));
                end

                compareCorr.XCeffect(count) = min(compareCorr.rawXC(count,27:31)) - min(compareCorr.rawXC(count,21:25));
                count = count+1;
           end
        end

    end

    figure(2);
    subplot(121);
    plot(compareCorr.WPSTHpkoff,compareCorr.XCeffect,'ko');
    subplot(122);
    plot(compareCorr.WPSTHcorr,compareCorr.XCeffect,'ko');

    allCompares = [allCompares [compareCorr.WPSTHpkoff;compareCorr.WPSTHcorr;compareCorr.XCeffect] ];
    size(allCompares)
% Compare that to the offset in the phase of the peak response

end

%% FOR uSN PAPER

% idea 1: implement a model to show that distributed, weak inhibition is a unique way to achieve gain control
% idea 2: try to look for the hallmark of feedback gain control. does the 'state' of the network prior to stimulus predict the magnitude of the response 
%     -for each session of either acquisition or extinction use the smoothed rate psths to calculate the average activity before mean subtracted peak response

% testing idea 2:

clear gainEffect;

csTime  = 1e3;
baseWin = [csTime-250:csTime-50];
actWin  = [csTime:csTime+200];
gainEffect.completeX = [];
gainEffect.completeY = [];
disp(' ');disp(' ');disp(' ');


for aa=1:numel(sessList)
% for aa=1:1
    
    thisSession = sessList(aa);
%     thisSession = 83;

    numUnits = size(PopData.session(thisSession).unit,2);
    numTrials = size(PopData.session(thisSession).unit(1).respCSS.aligned,1);
    disp(['Session ' num2str(thisSession) ' ...units: ' num2str(numUnits) ' ...trials: ' num2str(numTrials)]);
    
    for k=1:numTrials
        for j = 1:numUnits
            gainEffect.session(aa).baseline(k,j) = sum(PopData.session(thisSession).unit(j).respCSS.aligned(k,baseWin).*1000);
            gainEffect.session(aa).activation(k,j) = sum(PopData.session(thisSession).unit(j).respCSS.aligned(k,actWin).*1000) - sum(PopData.session(thisSession).unit(j).respCSS.aligned(k,baseWin).*1000);
        end
        
        gainEffect.session(aa).response(k,1) = PopData.session(thisSession).behavior.flCSon(k);
        gainEffect.session(aa).response(k,2) = PopData.session(thisSession).behavior.anticipLicks(k);

    end
            
    count = 0;
    
    indArray = [1:numUnits];
    
    for j = 1:numUnits
        for k=1:numTrials
            count = count + 1;
            gainEffect.session(aa).y(1,count) = gainEffect.session(aa).activation(k,j) - mean(gainEffect.session(aa).activation(:,j));
            nonCellInds = find(indArray~=j);
            gainEffect.session(aa).x(1,count) = mean(gainEffect.session(aa).baseline(k,indArray(nonCellInds))) - mean(mean(gainEffect.session(aa).baseline(:,indArray(nonCellInds))));            
        end
    end
    
    p = polyfit(gainEffect.session(aa).x,gainEffect.session(aa).y,1);
    gainEffect.session(aa).xF = 0:0.25:15;
    gainEffect.session(aa).yF = polyval(p,[0:0.25:15]);
    gainEffect.session(aa).numUnits = numUnits;        

%     figure(1); plot(gainEffect.session(aa).x,gainEffect.session(aa).y,'k.',gainEffect.session(aa).xF,gainEffect.session(aa).yF,'r-'); drawnow;
    
    [gainEffect.rho(aa) , gainEffect.pval(aa)] = corr(gainEffect.session(aa).x' , gainEffect.session(aa).y');
    gainEffect.numU(aa) = numUnits;
    
    gainEffect.completeX = [gainEffect.completeX gainEffect.session(aa).x];
    gainEffect.completeY = [gainEffect.completeY gainEffect.session(aa).y];
    
    disp(['  ... rho=' num2str(gainEffect.rho(aa)) '... pval=' num2str(gainEffect.pval(aa))]);

end

[gainEffect.completeRho , gainEffect.completePval] = corr(gainEffect.completeX' , gainEffect.completeY');

p = polyfit(gainEffect.completeX , gainEffect.completeY , 1);
gainEffect.completeXF = -25:0.1:25;
gainEffect.completeYF = polyval(p,[-25:0.1:25]);

[vals,inds] = sort(gainEffect.completeX);
figure(1); plot(gainEffect.completeX, gainEffect.completeY, 'k.', gainEffect.completeX(inds), sgolayfilt(gainEffect.completeY(inds),3,75), 'b-', gainEffect.completeXF, gainEffect.completeYF , 'r-'); axis([-20 20 -30 30]);

figure(2); clf; hist(gainEffect.rho,-1:0.05:1); hold on; plot(gainEffect.completeRho,15,'rv','MarkerSize',10,'LineWidth',2); text(gainEffect.completeRho+0.25,15,[num2str(gainEffect.completeRho) '...' num2str(gainEffect.completePval)]); axis([-0.5 0.5 0 20]);

tmp = find(gainEffect.pval<0.01);
figure(3); plot(gainEffect.numU,gainEffect.rho,'ko',gainEffect.numU(tmp),gainEffect.rho(tmp),'k.'); axis([0 25 -0.5 0.5]);

%% Bin and mean the population data
clear binData;
dataY = gainEffect.completeY;
dataX = gainEffect.completeX;
binWindow = [-5000 , 5000];
binWidth = 500
binMin = binWindow(1)
binMax = binWindow(1)+binWidth
count = 1;

while binMax<=binWindow(2)

    if count==1
        valids                  = find(dataX<binMax);        
    else
        valids                  = find(dataX>=binMin & dataX<binMax);
    end
    binData.mean(count)     = mean(dataY(valids));
    binData.sem(count)      = std(dataY(valids)) ./ sqrt(numel(valids)-1);
    binData.cnt(count)      = numel(valids);
    binData.center(count)   = mean([binMin binMax]);
    
    binMin  = binWindow(1)+(count*binWidth);
    count   = count+1;
    binMax  = binWindow(1)+(count*binWidth);
    
end

figure(100);
plot(binData.center,binData.mean,'ko-');

%% perform a permutation test to make sure that this is not guaranteed to give this correlation

allSamples = numel(gainEffect.completeX);

for p=1:5000
    
    gainEffectPerm.y = gainEffect.completeY(1,:);
    gainEffectPerm.x = gainEffect.completeX(1,randperm(allSamples));

    [gainEffectPerm.rho(p) , gainEffectPerm.pval(p)] = corr(gainEffectPerm.x' , gainEffectPerm.y');

    
end

figure(4); hist(gainEffectPerm.rho,-0.5:0.01:0.5);

%% FILTER THIS SESSIONS CONTINUOUS DATA

numChan             = size(Ns4DATA.Data,1);

paramsC.tapers      = [3 5];
paramsC.pad         = 0;
paramsC.Fs          = 1000; % reflects the fact that the data was decimated down to 1 kHz sampling frequency
paramsC.fpass       = [0 100];
paramsC.err         = 0;
movingwin           = [0.5 0.05];

ConcatData.paramsC     = paramsC;
ConcatData.movingwin   = movingwin;

testDataRange = 1:7.5e6;
    
% cycle through all channels and filter/calculate power within each band
for i=1:numChan

    tstart = tic;

    % Filter the data to select for spikes
    disp(' ');
    disp(' ');
    disp(['Filtering data on channel ' num2str(i) ' ...']);
    [lowBandData,hiBandData] = TNC_FilterData(Ns4DATA.Data(i,testDataRange),Ns4DATA.MetaTags.SamplingFreq,Ns4DATA.MetaTags.SamplingFreq./1000,0);        

    data = rmlinesmovingwinc(lowBandData.values,movingwin,10,paramsC,0.05,'n',60); % filter out 60 Hz noise
    lowBandData.values=data';

    % Apply chronux multitaper method to extract power spectrum in time
    disp('Calculating the spectrogram of the lowBand data...');

    ConcatData.lfp.params      = paramsC;
    ConcatData.lfp.movingwin   = movingwin;
    [S,t,f] = mtspecgramc(lowBandData.values,movingwin,paramsC);

    ConcatData.elec(i).contData.t = t;
    ConcatData.elec(i).contData.f = f;
    ConcatData.elec(i).contData.S = S;
    
    ConcatData.bands(1).name = '2to5';
    ConcatData.bands(2).name = '5to12';
    ConcatData.bands(3).name = '15to30';
    ConcatData.bands(4).name = '30to50';
    ConcatData.bands(5).name = '70to80';

    ConcatData.bands(1).inds = find(f<2,1,'last')  : find(f>5,1,'first');
    ConcatData.bands(2).inds = find(f<5,1,'last')  : find(f>12,1,'first');
    ConcatData.bands(3).inds = find(f<12,1,'last') : find(f>24,1,'first');
    ConcatData.bands(4).inds = find(f<24,1,'last') : find(f>42,1,'first');
    ConcatData.bands(5).inds = find(f<65,1,'last') : find(f>85,1,'first');

    telapsed = toc(tstart);
    disp(['Filtered and extracted power spec from one channel in ' num2str(telapsed) ' seconds.']);
    
    figure(1); hold off;
    plot(f,std(S)./mean(S)); hold on;
    peak = max(std(S)./mean(S));
    plot([2 2],[0 peak],'r--');
    plot([5 5],[0 peak],'r--');
    plot([12 12],[0 peak],'r--');
    plot([24 24],[0 peak],'r--');
    plot([42 42],[0 peak],'r--');
    plot([65 65],[0 peak],'r--');
    plot([85 85],[0 peak],'r--');
    xlabel('Frequency (Hz)');
    ylabel('CV');
    title('Frequency bands selected for analysis'); 
    
end
  
%% EXPORT A SINGLE SESSION WITH UNITS AND LFPs CONCATENATED FOR ALL TRIALS
% cycle through each trial and grab continuous data and spike times
numTrials = 20;
expSess = 83;
winR = 8000;
winL = 1000;
tmpT = ConcatData.elec(1).contData.t.*1000; % get back into units of ms
numElecs = size(ConcatData.elec,2);

sortedByLat=1;

disp('Sorting...');
if sortedByLat==1
    [vals,trialSet] = sort(PopData.session(i).behavior.flCSon);
else
    trialSet = 1:numTrials;
end


for i=1:numTrials
    
%     currTrial = trialSet(i+10);
    currTrial = i;
    currTs = PopData.session(expSess).events.CS.ts(currTrial)

    window = find(tmpT<(currTs-winL),1,'last') : find(tmpT>(currTs+winR),1,'first');
    
    disp(['Trial ' num2str(currTrial) ' ...']);
    
    % cycle through all electrodes
    for j = 1:numElecs


        if i==1
            for k=1:5
                ConcatData.elec(j).bands(k).values = mean(ConcatData.elec(j).contData.S(window, ConcatData.bands(k).inds)');
            end
        else
            for k=1:5
                tmpValues = mean(ConcatData.elec(j).contData.S(window, ConcatData.bands(k).inds)');
                ConcatData.elec(j).bands(k).values = [ConcatData.elec(j).bands(k).values tmpValues];
            end
        end
    
        
    end
    
    % cycle through all units
    for j = 1:numUnits
            
            numStamps = length(PopData.session(expSess).unit(j).ts);
            delta = zeros(1,ceil(PopData.session(expSess).unit(j).ts(numStamps)));
            delta(1,round(PopData.session(expSess).unit(j).ts)+1) = 1;

            [respCS] = TNC_AlignRasters(delta,PopData.session(expSess).unit(j).ts,PopData.currParams.stableTrials,PopData.session(expSess).events.CS.ts,[winL,winR],1,1);
            PopData.session(expSess).unit(j).respCS.raster      = respCS.raster;

        if i==1
            if numel(PopData.session(expSess).unit(j).respCS.raster.trial(currTrial).ts) > 0
                ConcatData.unit(j).ts    =  PopData.session(expSess).unit(j).respCS.raster.trial(currTrial).ts' + ((winL+winR).*(i-1)) + winL;
            else
                ConcatData.unit(j).ts    = [];
            end
            ConcatData.unit(j).eName     =  PopData.session(expSess).unit(j).name;
        else
            clear tsTMP;
            tsTMP                        =  PopData.session(expSess).unit(j).respCS.raster.trial(currTrial).ts' + ((winL+winR).*(i-1)) + winL;
            if numel(tsTMP)>0
                ConcatData.unit(j).ts    =  [ConcatData.unit(j).ts tsTMP];
            end
            
        end
    
    end
    
end

for j = 1:numElecs

    for k=1:5
        miVal = min(ConcatData.elec(j).bands(k).values);
        maVal = max(ConcatData.elec(j).bands(k).values);
        if maVal>0
            ConcatData.elec(j).bands(k).values = (ConcatData.elec(j).bands(k).values-miVal) ./ (maVal-miVal);
        end
        miVal = min(ConcatData.elec(j).bands(k).values);
        maVal = max(ConcatData.elec(j).bands(k).values);
    end

    disp(['Electrode ' num2str(j) ' was normalized to the range 0...1']);

end

%% DISPLAY ALIGNED LFP RESPONSES FOR FIRST N TRIALS
% cycle through each trial and grab continuous data and spike times
numTrials = 10;
expSess = 83;
winR = 5000;
winL = 500;
tmpT = ConcatData.elec(1).contData.t.*1000; % get back into units of ms
numElecs = size(ConcatData.elec,2) - 1;
numUnits = size(PopData.session(expSess).unit,2);

sortedByLat=1;

disp('Sorting...');
if sortedByLat==1
    [vals,trialSet] = sort(PopData.session(i).behavior.flCSon);
else
    trialSet = 1:numTrials;
end


for i=1:numTrials
    currTrial = trialSet(i);
    currTs = PopData.session(expSess).events.CS.ts(currTrial);

    window = find(tmpT<(currTs-winL),1,'last') : find(tmpT>(currTs+winR),1,'first');
    
    disp(['Trial ' num2str(currTrial) ' ...']);
    
    % cycle through all electrodes
    for j = 1:numElecs

        ConcatData.raster.trial(i).elec(j).S = ConcatData.elec(j).contData.S(window, :);
        ConcatData.raster.trial(i).elec(j).t = -500:50:5050;
        ConcatData.raster.trial(i).elec(j).f = ConcatData.elec(j).contData.f;
        
        for k=1:numel(ConcatData.elec(j).contData.f)
            ConcatData.raster.trial(i).elec(j).S(:,k) = ConcatData.raster.trial(i).elec(j).S(:,k) ./ mean(ConcatData.raster.trial(i).elec(j).S(1:9,k));
        end
       
        % display the data for this trial
        figure(1); subplot(numTrials,numElecs,(i*numElecs)+j);
        imagesc(ConcatData.raster.trial(i).elec(j).t(1:50),ConcatData.raster.trial(i).elec(j).f,ConcatData.raster.trial(i).elec(j).S(1:50,:)',[0 5]);

    end
    
%     % cycle through all units
%     for j = 1:numUnits
%             
%         numStamps = length(PopData.session(expSess).unit(j).ts);
%         delta = zeros(1,ceil(PopData.session(expSess).unit(j).ts(numStamps)));
%         delta(1,round(PopData.session(expSess).unit(j).ts)+1) = 1;
% 
%         [respCS] = TNC_AlignRasters(delta,PopData.session(expSess).unit(j).ts,PopData.currParams.stableTrials,PopData.session(expSess).events.CS.ts,[winL,winR],1,1);
%         PopData.session(expSess).unit(j).respCS.raster      = respCS.raster;
% 
%         if numel(PopData.session(expSess).unit(j).respCS.raster.trial(i).ts) > 0
%             ConcatData.raster.trial(i).unit(j).ts    =  PopData.session(expSess).unit(j).respCS.raster.trial(i).ts';
%         else
%             ConcatData.raster.trial(i).unit(j).ts    = [];
%         end
%         ConcatData.raster.trial(i).unit(j).eName     =  PopData.session(expSess).unit(j).name;
%     
%     end    
    
end

%% EXAMINE CONCATENATED UNIT DATA TO CONFIRM THAT IT LOOKS CORRECT

% cycle through all units
for j = 1:numUnits

    delta = zeros(1,20.*8500);
    tDelta = zeros(1,20.*8500);
    tDelta((8500.*[0:19])+500) = 1;

    if numel(ConcatData.unit(j).ts)>0
        delta(ceil(ConcatData.unit(j).ts)) = 1;
        tmpSmooth = conv(delta,currParams.filter.kernel,'same');
    else
        tmpSmooth = delta;
    end
    
    % display the data for this trial
    figure(2); subplot(numUnits./3,3,j); hold off;
    plot(tmpSmooth,'k','LineWidth',2); hold on;
    plot(tDelta.*max(tmpSmooth),'r--');
    
end


%% SOME OPEN QUESTIONS

Does ability to predict relate to variability in "intervening behavioral states"?
Do delays propagate (i.e. a cascade model) or are they independent (i.e. an afferent drive for "sequences")?
Correlating activity trajectories with properties of the approach response

%% CLASSIFY UNITS

interested in:
delay period activity
burst UP lick transition
pause DOWN lick transition
sustained UP during lick
sustained DOWN during lick
nonzero slope during delay