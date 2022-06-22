% Approach behavior analysis script
%% STANDARD ANALYSIS PARAMETERS
% Define the smoothing to apply to the timeseries data
currParams.smthParams.rise     = 1;
% currParams.smthParams.decay    = 100;
% currParams.smthParams.decay    = 5;
% currParams.smthParams.decay    = 10;
currParams.smthParams.decay    = 25;
% currParams.smthParams.decay    = 60;
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
% for i=83;        
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

            [respCS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[1e3,1e4],1,1);
            PopData.session(i).unit(j).respCS.raster      = respCS.raster;

            tmpSmooth = conv(delta,currParams.filter.kernel,'same');
            [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[1e3,1e4],0,1);
            PopData.session(i).unit(j).respCSS.aligned    = respCSS.image.aligned;
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
k=1;

for i=1:NumSessions
% for i=83;        
    if strcmp(PopData.session(i).sessClass,'learning')

        numUnits = size(PopData.session(i).unit,2);
        countUnits = countUnits+numUnits;
        PopData.session(i).trials = size(PopData.session(i).events.CS.ts,1);
    
        if numUnits>8 & PopData.session(i).behavior.found==1
            sessList(k) = i;
            k=k+1;
            disp(['Session: ' num2str(i) ' is a candidate session | ' PopData.session(i).sessId])
        end
        

    end    
end
k = k-1;
disp(['___________________________________________________'])
disp(['COMPLETED ... Valid Sessions: ' num2str(k)]);

%% CHECK FIRST AND LAST LICK DETECTION

disp(['___________________________________________________'])
disp(['STARTED aligning behavior and updating the PopData structure...']);

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
validBehavSess = [];
k=1;
hahnloser = 1; 
numTrialsSpk = 40;
sortedByLat = 1;

% for i=1:NumSessions
% for aa = 1:numel(sessList)
for i=83 % specify the list of sessions to analyze
  
%         i = sessList(aa);
    
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
                    firstLickPostCS = find(theseTimeStamps>200,1);
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

% for i=1:NumSessions
for aa = 1:numel(sessList)
% for i=83 % specify the list of sessions to analyze
  
        i = sessList(aa);
        currFL = PopData.session(i).behavior.flCSon;
        PopDataSub.session(aa).firstLick = currFL;
        
        allFL = [allFL currFL];
        
end

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

% for aa = 1:numel(sessList)
for i = 83;
%     i=sessList(aa);
    
    if strcmp(PopData.session(i).sessClass,'learning')
        disp(['Session ' num2str(i)]);

        if isfield(PopData.session(i),'behavior')
            if isfield(PopData.session(i).behavior,'found')

                numUnits    = size(PopData.session(i).unit,2);
                
                figure(43); clf;
                figure(44); clf;
                figure(45); clf;
                figure(46); clf;
                
                validTrials = find(PopData.session(i).behavior.flCSon>600 & PopData.session(i).behavior.flCSon<1600 & PopData.session(i).behavior.flCSoff>(PopData.session(i).behavior.USon+500));
                numTrials   = numel(validTrials);
                PopData.session(i).behavior.validPhaseTrials = validTrials;

                for m=1:numUnits
                    thisUnit = m

                    for q=1:5
                        allPhase.seg(q).phases = [];
                        allPhase.seg(q).times = [];
                    end
                    allPhase.seg(1).pOff = -1;
                    allPhase.seg(2).pOff = 0;
                    allPhase.seg(3).pOff = 0.718;
                    allPhase.seg(4).pOff = 0.718 + 1.077;
                    
                    % find appropriate trials where the ts are in sequence
                    % CS on, CS off, FL, US, LL                    

                    
                    % cycle through all valid trials
                    for k = 1:numTrials
                        
                        thisTrial = validTrials(k);
                        theseTimeStamps = PopData.session(i).unit(thisUnit).respCS.raster.trial(thisTrial).ts;
                        
                        size(theseTimeStamps)
                        
                        % warp the baseline data (1 radian)
                        validStamps     = find( theseTimeStamps<0 );
                        thisSegDuration = 1000;
                        if numel(validStamps)>0
                            currPhases  = theseTimeStamps(validStamps) .* (1 ./ thisSegDuration); % multiply by rads / sec
                            currTimes   = theseTimeStamps(validStamps);
                            if abs(numel(currPhases)-numel(currTimes)) > 0
                                disp('An error has occurred.');
                            end
                            if k==1
                                allPhase.seg(1).phases = currPhases'; % multiply by rads / sec
                                allPhase.seg(1).times = currTimes'; % multiply by rads / sec
                                allPhase.seg(1).pOff = -1;
                            else
                                allPhase.seg(1).phases = [allPhase.seg(1).phases currPhases'];
                                allPhase.seg(1).times = [allPhase.seg(1).times currTimes'];
                            end
                            allPhase.seg(1).trial(k).ph = currPhases;
                            allPhase.seg(1).trial(k).ts = currTimes;
                        else
                            allPhase.seg(1).trial(k).ph = [];
                            allPhase.seg(1).trial(k).ts = [];                            
                        end                            
                        
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
                        else
                            allPhase.seg(2).trial(k).ph = [];
                            allPhase.seg(2).trial(k).ts = [];                            
                        end
                        
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
                        else
                            allPhase.seg(3).trial(k).ph = [];
                            allPhase.seg(3).trial(k).ts = [];                            
                        end

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
                        else
                            allPhase.seg(4).trial(k).ph = [];
                            allPhase.seg(4).trial(k).ts = [];                            
                        end
                        
                    end
                    
                    PopData.session(i).unit(thisUnit).allPhase = allPhase;
                    
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

                    drawnow;
                    
                end 
            end
        else
            disp('No good behavioral data... skipping');
        end
    
        drawnow;
    end
end

%% CREATE A HAHNLOSER TYPE RASTER PLOT FOR TIMES AND FOR PHASE
% Toggles between two styles. One is grouped by trial and one is grouped by unit

export =1;

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
% for aa = 1:numel(sessList)
for i=83 % specify the list of sessions to analyze
  
%         i = sessList(aa);
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
                        
                        delta = zeros(1,7300);

                        thisUnit = m;

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
                                    plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
                               else
                                    subplot(5,1,2:5);  hold on;
                                    plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);            
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
                    axis([-1000 8500 -1 ((numUnits).*numTrials) + (5.*(numUnits+1))]);
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
                        
                        delta = zeros(1,7300);
                        
                        thisUnit = m;

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
                    axis([-1 6.285 -1 ((numUnits).*numTrials) + (5.*(numUnits+1))]);
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
                        
                        delta = zeros(1,7300);
                        
                        thisUnit = sortedLats(m);

                        for k = 1:numTrials
                            
                            thisTrialRow    = ((numUnits+3).*k) + m;
                            theseTimeStamps = [PopData.session(i).unit(thisUnit).allPhase.seg(1).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(2).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(3).trial(trialIndsTmp(k)).ph,
                                               PopData.session(i).unit(thisUnit).allPhase.seg(4).trial(trialIndsTmp(k)).ph];
                                           
                            delta(ceil((theseTimeStamps+1).*1000)) = delta(ceil((theseTimeStamps+1).*1000)) + 1;

%                             if rem(m,2)
                                subplot(5,1,2:5); hold on;
                                plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',6,'Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
%                             else
%                                 subplot(5,1,2:5);  hold on;
%                                 plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',6,'Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);            
%                             end

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
% for i=85 % specify the list of sessions to analyze
  
% for aa=1:numel(sessList)
for i=83;    
    figure(39); clf;
%     i=sessList(aa);
    
    disp(['Session ' num2str(i)]);
    numUnits = size(PopData.session(i).unit,2);

    % sort by init response latency to make it easier to see
%     [vals,sortedLats] = sort(InitRespLatency); 
    [vals,sortedLats] = sort(1:numUnits ); 

                    
    PhaseMat = zeros(numUnits,7300);
    
    for m=1:numUnits
        
        thisUnit = sortedLats(m);
        
        tmpSmoothP = conv(PopData.session(i).unit(thisUnit).allPhase.phasePSTH.delta,currParams.filter.kernel,'same');
        PopData.session(i).unit(thisUnit).allPhase.phasePSTH.psth = tmpSmoothP;
        disp('Finished phase psth.');
        
        tmpSmoothT = conv(PopData.session(i).unit(thisUnit).allPhase.timePSTH.delta,currParams.filter.kernel,'same');
        PopData.session(i).unit(thisUnit).allPhase.timePSTH.psth = tmpSmoothT;
        disp('Finished time psth.');
        
        % compare to psth with no selection and no warping
        numTrials = size(PopData.session(i).unit(thisUnit).respCS.raster.trial,2);
        delta = zeros(1,7300);

        for k = 1:numTrials
           theseTimeStamps = PopData.session(i).unit(thisUnit).respCS.raster.trial(k).ts;
           tmp = find(theseTimeStamps<6300 & theseTimeStamps>-1000);
           delta(ceil(theseTimeStamps(tmp)+1000)) = delta(ceil(theseTimeStamps(tmp)+1000)) + 1;
        end
        
        delta = delta ./ numTrials;
        tmpSmoothU = conv(delta,currParams.filter.kernel,'same');
        disp('Finished U psth.');
        
        figure(39); hold on;
        plot(-1000:6299,(tmpSmoothP./max(tmpSmoothP))+m-1,'Color',[1 0 0]);
        plot(-1000:6299,(tmpSmoothT./max(tmpSmoothP))+m-1,'Color',[0 0.67 1]);            
        plot(-1000:6299,(tmpSmoothU./max(tmpSmoothP))+m-1,'Color',[0.5 0.5 0.5]);
        
        InitRespLatency(m) = find(tmpSmoothP==max(tmpSmoothP(1000:1718)),1);
        
        PhaseMat(thisUnit,:) = tmpSmoothP - mean(tmpSmoothP);

    end
    
    plot([0,0],[0,numUnits],'k--');
    plot([718,718],[0,numUnits],'k--');
    plot([1795,1795],[0,numUnits],'k--');
    plot([6283,6283],[0,numUnits],'k--');
    axis([-1100 6400 -1 numUnits+1]);
    
    PopData.session(i).PhaseMat = PhaseMat;

    figure(36);
    imagesc(corr(PhaseMat),[-1 1]);

    figure(37);
    imagesc(corr(PhaseMat'),[-1 1]);

end

%% DIMENSIONALITY REDUCTION OF POPULATION ACTIVITY - SINGLE SESSION
% This segment of the code is meant to calculate two kinds of reduced
% representations:
%   1) What are the minimal characteristic time courses of activation over the approach window?
%   2) What is the trajectory of the population vector over the approach window?

% Define the time window over which to reduce dimensions (window of
% approach seems most applicable)
window = 500:3500;
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
    binUnits = floor(numTrials./bins)
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

%% DIMENSIONALITY REDUCTION OF POPULATION ACTIVITY - PSEUDO-SIMULTANEOUS

%% LOADING OF CONTINUOUS AND BEHAVIORAL DATA
% hard coded at the moment
fileNameStr         = 'DA-Ma-05 100231 trace-d-001.ns4'
Ns4DATA             = openNSx('read',fileNameStr);
session.creation    = Ns4DATA.MetaTags.CreateDateTime;

%% WARPING OF BEHAVIORAL DATA
% uses interp and decimate to resample behavioral data into the warped coordinates

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
    
%% CLASSIFY UNITS

interested in:
delay period activity
burst UP lick transition
pause DOWN lick transition
sustained UP during lick
sustained DOWN during lick
nonzero slope during delay

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
numUnits = size(PopData.session(i).unit,2);
i=83;
clear PopVec
disp('Sorting...');

for aa=1:numel(sessList)

    i = sessList(aa);

    if sortedByLat==1
        [vals,trialInds] = sort(PopData.session(i).behavior.flCSon);
    else
        trialInds = 1:numTrials;
    end

    disp(['Collecting vectors for ' num2str(numUnits) ' units.']);
    for m=1:numUnits

        thisUnit = m;

        for k = 1:numTrialsSpk

            % create smoothed inst rate
            delta = zeros(1,10000);
            tmpTimeStamps = round(PopData.session(i).unit(thisUnit).respCS.raster.trial(trialInds(k)).ts);
            posStamps = find(tmpTimeStamps>-1000 & tmpTimeStamps<9000);
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

    [v,d] = eig(meanCovMat);
    PopVec.pca.vals = diag(d)./sum(diag(d));
    PopVec.pca.vecs(:,1) = v(:,numel(PopVec.pca.vals));
    PopVec.pca.vecs(:,2) = v(:,numel(PopVec.pca.vals)-1);
    PopVec.pca.vecs(:,3) = v(:,numel(PopVec.pca.vals)-2);
    PopVec.pca.varExp = sum(PopVec.pca.vals(numel(PopVec.pca.vals)-2:numel(PopVec.pca.vals)));
    disp(['Variance explained by first 3 components: ' num2str(sum(PopVec.pca.vals(numel(PopVec.pca.vals)-2:numel(PopVec.pca.vals))))]);

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
            
            % peak spikes relative to mean
            varSpks(k,p) = response - trapz(PopData.session(thisSession).unit(p).respCS.psthAVG());
            varSpksP(k,p) = response -


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

%% FIGURES

1. Example of the behavior and recording configuration
2. Microstructure of approach
3. Position data correlated with time to sample the port (i.e. can we use time to sample port as the read out?)
4. Compare population PSTHs aligned to the CS only and warped
5. Example activity trajectories in sample session and the pop activity in the approach behavior epoch
6. Comparing the number of neural states with the behavioral states. Possible? How to get 1st lick latency into behavioral data?
7. Ability to predict the approach time as a function of latency from the CS (or should i conceive of this as transitions through a series of states.)

Does ability to predict relate to variability in "intervening behavioral states"?
Do delays propagate (i.e. a cascade model) or are they independent (i.e. an afferent drive for "sequences")?

Potential
 Correlating activity trajectories with properties of the approach response
 Lag correlation structures (could i use the context-free grammar ideas?)
