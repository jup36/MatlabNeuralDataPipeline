% COMPLETE EXTINCTION PAPER ANALYSIS SEQUENCE
%
% Figure 1: Example of the training paradigm, recording schematic, behavioral control over learning and extinction, cs responses, latency determination, exemplar responses
% Figure 2: Extinction cs responses, contrast responses, example recording sessions, response amplitude scatter plot, response amplitude as a function of latency
% Figure 3: Example simultaneous recording, **comparison between extinction cells and dopamine cells**, comparison between extinction cells and dopamine cells and behavior
% __________________________________
% NOT INCLUDED IN THIS SCRIPT
% Figure 4: In vitro recording of monosynaptic inhibition, circuit mapping, functional role of inhibition

%% SWAP SPACE

% sessions with surprisingly long da latencies
% strange sessions: [15, 53, 57, 85, 87]

[15    53    57    59    85    87   109]

%% STANDARD ANALYSIS PARAMETERS
% Define the smoothing to apply to the timeseries data
currParams.smthParams.rise     = 1;
% currParams.smthParams.decay    = 4;
% currParams.smthParams.decay    = 10;
currParams.smthParams.decay    = 20;
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

% hand curated dopamine neuron list:
PopData.daList = [42,67,161,177,187,329,347,503,504,514,522,523,528,529,530,532,539,540,541,544,545,546];
PopData.daSessList = [15 25 53 57 59 85 87 109 111 113 115 117 119 122 123];

count=1;
for i=1:663
    if numel(find(PopData.daList==i))==0
        PopData.nonDaList(count) = i;
        count = count+1;
    end
end

count=1;
for i=1:663
    if numel(find(PopData.daSessList==allCSaligned.sess(i)))==0
        PopData.nonDaSessList(count) = i;
        count = count+1;
    end
end

PopData.currParams = currParams;
disp(' ');
disp(' ');
disp('________________________________');
disp(' ');
disp('Initialized analysis parameters.');
disp('________________________________');
disp(' ');
disp(' ');

%% Correct latency shift
NumSessions = size(PopData.session,2);
k=1;
undo = 0;

for i=1:NumSessions
    if undo
        if PopData.session(i).latencyShift==2
            PopData.session(i).events.CS.ts = PopData.session(i).events.CS.ts - 10;
            PopData.session(i).latencyShift = 1;
            disp([PopData.session(i).sessId ' >>> UN - SHIFT.']);
        else
            disp([PopData.session(i).sessId ' >>> skipping.']);
        end        
    else
        if PopData.session(i).latencyShift
            PopData.session(i).events.CS.ts = PopData.session(i).events.CS.ts + 10;
            PopData.session(i).latencyShift = 2;
            disp([PopData.session(i).sessId ' >>> SHIFT.']);
        else
            disp([PopData.session(i).sessId ' >>> skipping.']);
        end
    end
end

if undo
    disp('New latency shift undone.');
else
    disp('New corrected latency shift applied.');
end

%% LIST ALL UNITS AND THEIR INDEX

for i = 1:size(allCSaligned.psthZ,2)
    disp([num2str(i) ' >>> ' PopData.session(allCSaligned.sess(i)).sessId '   >>   ' PopData.session(allCSaligned.sess(i)).unit(allCSaligned.unit(i)).name]);
end

%% STORE DA LOGIC FOR EACH UNITS (i.e. CONVERT DALIST INTO SESSION AND UNIT NUMBERS)

NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
    
        numUnits = size(PopData.session(i).unit,2);

        for j = 1:numUnits
            PopData.session(i).unit(j).daLogic = 0;         
        end
end


for i=1:numel(PopData.daList)
    
    currIndex = PopData.daList(i);
    
    sessA = allCSaligned.sess(currIndex);
    unitA = allCSaligned.unit(currIndex);
    PopData.session(sessA).unit(unitA).daLogic = 1;

    sessE = allCSEXaligned.sess(currIndex);
    unitE = allCSEXaligned.unit(currIndex);
    PopData.session(sessE).unit(unitE).daLogic = 1;
    
end

%% EXTRACT THE ISI FROM ALL UNITS
NumSessions = size(PopData.session,2)
k=1;m=1;

for i=1:NumSessions
    % check if learning
    if i<10
        disp(['Session 00' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
    elseif i<100
        disp(['Session 0' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
    else
        disp(['Session ' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
    end

    if strcmp(PopData.session(i).sessClass,'learning')
            
        numUnits = size(PopData.session(i).unit,2);
        countUnits = countUnits+numUnits;

        for j = 1:numUnits
            if isnan(PopData.session(i).unit(j).isi.hist.logCount(1,1))
                disp('Invalid ISI distribution found');
            else
                allISI.acq.linear.times(1,:)    = PopData.session(i).unit(j).isi.hist.linTimes;
                allISI.acq.linear.counts(k,:)   = PopData.session(i).unit(j).isi.hist.linCount;
                
                allISI.acq.log.times(1,:)    = PopData.session(i).unit(j).isi.hist.logTimes;
                allISI.acq.log.counts(k,:)   = PopData.session(i).unit(j).isi.hist.logCount;

                k=k+1;
            end
        end
        
    else

        numUnits = size(PopData.session(i).unit,2);

        for l = 1:numUnits
            if isnan(PopData.session(i).unit(l).isi.hist.logCount(1,1))
                disp('Invalid ISI distribution found');
            else
                allISI.ext.linear.times(1,:)    = PopData.session(i).unit(l).isi.hist.linTimes;
                allISI.ext.linear.counts(m,:)   = PopData.session(i).unit(l).isi.hist.linCount;
                
                allISI.ext.log.times(1,:)    = PopData.session(i).unit(l).isi.hist.logTimes;
                allISI.ext.log.counts(m,:)   = PopData.session(i).unit(l).isi.hist.logCount;

                m=m+1;
            end
        end

        
    end
    
end

k = k-1;

disp(['Total units: ' num2str(k) ' | ' num2str(countUnits)]);

%% EXTRACT STABLE PSTHS FROM ALL SPIKE DATA
disp(['___________________________________________________'])
disp(['STARTED aligning all raster plots and updating the PopData structure...']);

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
    
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

        [respCS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[PopData.currParams.winParams.prior,PopData.currParams.winParams.after],1,1);
        PopData.session(i).unit(j).respCS.raster            = respCS.raster;
        PopData.session(i).unit(j).respCS.boxcar            = respCS.image.boxcar;

        tmpSmooth = conv(delta,currParams.filter.kernel,'same');
        [respCSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts,[PopData.currParams.winParams.prior,PopData.currParams.winParams.after],0,1);
        PopData.session(i).unit(j).respCSS.psthAVG    = respCSS.image.psthAVG;
        PopData.session(i).unit(j).respCSS.psthZ      = respCSS.image.psthZ;
        PopData.session(i).unit(j).respCSS.psthZe     = respCSS.image.psthZe;
        
        if strcmp(PopData.session(i).sessClass,'learning')
            if numel(PopData.session(i).events.US.ts)==1
                [respUS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts+2500,[PopData.currParams.winParams.prior,PopData.currParams.winParams.after],1,1);
                PopData.session(i).unit(j).respUS.raster    = respUS.raster;
                PopData.session(i).unit(j).respUS.boxcar    = respUS.image.boxcar;

                [respUSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.CS.ts+2500,[PopData.currParams.winParams.prior,PopData.currParams.winParams.after],0,1);
                PopData.session(i).unit(j).respUSS.psthAVG    = respUSS.image.psthAVG;
                PopData.session(i).unit(j).respUSS.psthZ      = respUSS.image.psthZ;
                PopData.session(i).unit(j).respUSS.psthZe     = respUSS.image.psthZe;
            else            
                [respUS] = TNC_AlignRasters(delta,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.US.ts,[PopData.currParams.winParams.prior,PopData.currParams.winParams.after],1,1);
                PopData.session(i).unit(j).respUS.raster    = respUS.raster;
                PopData.session(i).unit(j).respUS.boxcar    = respUS.image.boxcar;

                [respUSS] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,PopData.currParams.stableTrials,PopData.session(i).events.US.ts,[PopData.currParams.winParams.prior,PopData.currParams.winParams.after],0,1);
                PopData.session(i).unit(j).respUSS.psthAVG    = respUSS.image.psthAVG;
                PopData.session(i).unit(j).respUSS.psthZ      = respUSS.image.psthZ;
                PopData.session(i).unit(j).respUSS.psthZe     = respUSS.image.psthZe;
            end
        end
        
        k=k+1;

    end
    
    disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' of ' num2str(NumSessions) ' | ' PopData.session(i).sessId '_' PopData.session(i).sessClass])
    
end
k = k-1;
k = k./2;
disp(['___________________________________________________'])
disp(['COMPLETED ... Total units: ' num2str(k) ' | ' num2str(countUnits) ' | Pairwise: ' num2str(numPairwise)]);

%% ALIGN PSTHS FOR LEARNING, EXTINCTION, CONTRAST
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
                allCSaligned.psthS(:,k)     = PopData.session(i).unit(j).respCSS.psthAVG'.*1000;
                allCSaligned.psthZ(:,k)     = PopData.session(i).unit(j).respCSS.psthZ';
                allCSaligned.psthZe(:,k)    = PopData.session(i).unit(j).respCSS.psthZe';
                allCSaligned.boxcar(:,k)    = PopData.session(i).unit(j).respCS.boxcar'.*(1000./6);
                allCSaligned.sess(1,k)      = i;
                allCSaligned.unit(1,k)      = j;
                allCSaligned.id(1,k).sessID = PopData.session(i).sessId;
                allCSaligned.id(1,k).name   = PopData.session(i).unit(j).name;

                allUSaligned.psthZ(:,k)     = PopData.session(i).unit(j).respUSS.psthZ';
                allUSaligned.psthZe(:,k)    = PopData.session(i).unit(j).respUSS.psthZe';                
                allUSaligned.boxcar(:,k)    = PopData.session(i).unit(j).respUS.boxcar'.*(1000./6);
                
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

figure(5);
subplot(121);
[mapName] = TNC_CreateRBColormap(1024,'rb');
imagesc(allUSaligned.psthZ',[-15 15]);
colormap(mapName);
xlabel('Time (ms)');
ylabel('Cell index');
title('US Aligned ZDF');

subplot(122);
[mapName] = TNC_CreateRBColormap(1024,'rb');
imagesc(allCSaligned.psthZ' - allUSaligned.psthZ',[-15 15]);
colormap(mapName);
xlabel('Time (ms)');
ylabel('Cell index');
title('CS-US Aligned ZDF');

%__________________________________________________________________________________
%__________________________________________________________________________________

% NOW FOR EXTINCTION
%__________________________________________________________________________________
%__________________________________________________________________________________

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
                allCSEXaligned.psthS(:,k)   = PopData.session(i).unit(j).respCSS.psthAVG'.*1000;
                allCSEXaligned.psthZ(:,k)   = PopData.session(i).unit(j).respCSS.psthZ';
                allCSEXaligned.psthZe(:,k)  = PopData.session(i).unit(j).respCSS.psthZe';
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

%__________________________________________________________________________________
%__________________________________________________________________________________

% NOW OBTAIN THE CONTRAST
%__________________________________________________________________________________
%__________________________________________________________________________________

allCSContrast.psthZ = allCSEXaligned.psthZ - allCSaligned.psthZ;

%% GET PSTH RESPONSE LATENCIES v2
% switch to calculate psth for both CS and US 
% -> then use both CS and US responses to classify cells (i.e. long -> 2; short gets larger by including US)

winCS = 2001:2140;
winUS = winCS;

threshold = 3.8  
% 0.0001 error rate (i.e. ~100 comparisons at 0.01 error rate)

respQuant.cs.latency = zeros(size(allCSaligned.psthZ,2),1);
respQuant.us.latency = zeros(size(allCSaligned.psthZ,2),1);
respQuant.ex.latency = zeros(size(allCSaligned.psthZ,2),1);

respQuant.cs.peak = zeros(size(allCSaligned.psthZ,2),1);
respQuant.us.peak = zeros(size(allCSaligned.psthZ,2),1);
respQuant.ex.peak = zeros(size(allCSaligned.psthZ,2),1);

respQuant.cs.amp = zeros(size(allCSaligned.psthZ,2),1);
respQuant.us.amp = zeros(size(allCSaligned.psthZ,2),1);
respQuant.ex.amp = zeros(size(allCSaligned.psthZ,2),1);

respQuant.cs.max = zeros(size(allCSaligned.psthZ,2),1);
respQuant.us.max = zeros(size(allCSaligned.psthZ,2),1);
respQuant.ex.max = zeros(size(allCSaligned.psthZ,2),1);

respQuant.cs.restore = zeros(size(allCSaligned.psthZ,2),1);
respQuant.us.restore = zeros(size(allCSaligned.psthZ,2),1);
respQuant.ex.restore = zeros(size(allCSaligned.psthZ,2),1);


for i = 1:size(allCSaligned.psthZ,2)
    
    csData = allCSaligned.psthZ(winCS,i);
    csCheck  = find(csData > threshold       );
    csCheck2 = find(csData > (max(csData)./1.5));
    
    if numel(csCheck)>10 & numel(csCheck2)<150 % 2x the s.d. of the smoothing kernel 

        figure(50); clf; hold on;
        subplot(211);
        plot(1:numel(csData),csData,'k',csCheck,csData(csCheck),'ro',find(csData==max(csData)),max(csData),'b.');
        title(num2str(i));
        
        % go through csCheck and find phasic peaks
        peakStarts = find(diff(csCheck)>1);
        
        if numel(peakStarts)>0 % multiple peaks
            % find the largest peak
            candidatePeaks = zeros(1,numel(peakStarts)+1);
            for k = 1:numel(peakStarts)+1
                if k==1
                    candidatePeaks(1,k) = trapz(csData(csCheck(1:peakStarts(k)))) ./ numel(1:peakStarts(k));
                elseif k==numel(peakStarts)+1
                    candidatePeaks(1,k) = trapz(csData(csCheck(peakStarts(k-1)+1:numel(csCheck)))) ./ numel(peakStarts(k-1)+1:numel(csCheck));                                        
                else
                    candidatePeaks(1,k) = trapz(csData(csCheck(peakStarts(k-1)+1:peakStarts(k)))) ./ numel(peakStarts(k-1)+1:peakStarts(k));                    
                end
            end

            largestPeak = find(candidatePeaks==max(candidatePeaks));
            if largestPeak==1
                peakInds = csCheck(1:peakStarts(largestPeak));
            elseif largestPeak==numel(peakStarts)+1
                peakInds = csCheck(peakStarts(largestPeak-1)+1:numel(csCheck));                                        
            else
                peakInds = csCheck(peakStarts(largestPeak-1)+1:peakStarts(largestPeak));                   
            end
            subplot(212);
            plot(1:numel(csData),csData,'k',csCheck,csData(csCheck),'ro',peakInds,csData(peakInds),'c.');            
        else
            peakInds = csCheck;
%             subplot(212);
%             plot(1:numel(csData),csData,'k',csCheck,csData(csCheck),'ro',peakInds,csData(peakInds),'c.');
        end
                
        respQuant.cs.latency(i) = peakInds(1);
        respQuant.cs.peak(i)    = peakInds(1) + find(csData(peakInds)==max(csData(peakInds)));
        respQuant.cs.amp(i)     = trapz(csData(peakInds))./numel(peakInds);
        respQuant.cs.restore(i) = peakInds(numel(peakInds));
        respQuant.cs.max(i)     = max(csData(peakInds));
            
    else
        
        respQuant.cs.latency(i) = -1;
        respQuant.cs.peak(i)    = -1;
        respQuant.cs.amp(i)     = -1;
        respQuant.cs.restore(i) = -1;
        respQuant.cs.max(i)     = -1;
        
    end
    
    
end
disp('Completed mean peak detection routine | CS.');

for i = 1:size(allUSaligned.psthZ,2)
    
    usData = allUSaligned.psthZ(winUS,i);
    usCheck  = find(usData > threshold       );
    usCheck2 = find(usData > (max(usData)./1.5));
    
    if numel(usCheck)>10 & numel(usCheck2)<150 % 2x the s.d. of the smoothing kernel 

%         figure(50); clf; hold on;
%         subplot(211);
%         plot(1:numel(usData),usData,'k',usCheck,usData(usCheck),'ro',find(usData==max(usData)),max(usData),'b.');
%         title(num2str(i));
        
        % go through usCheck and find phasic peaks
        peakStarts = find(diff(usCheck)>1);
        
        if numel(peakStarts)>0 % multiple peaks
            % find the largest peak
            candidatePeaks = zeros(1,numel(peakStarts)+1);
            for k = 1:numel(peakStarts)+1
                if k==1
                    candidatePeaks(1,k) = trapz(usData(usCheck(1:peakStarts(k)))) ./ numel(1:peakStarts(k));
                elseif k==numel(peakStarts)+1
                    candidatePeaks(1,k) = trapz(usData(usCheck(peakStarts(k-1)+1:numel(usCheck)))) ./ numel(peakStarts(k-1)+1:numel(usCheck));                                        
                else
                    candidatePeaks(1,k) = trapz(usData(usCheck(peakStarts(k-1)+1:peakStarts(k)))) ./ numel(peakStarts(k-1)+1:peakStarts(k));                    
                end
            end

            largestPeak = find(candidatePeaks==max(candidatePeaks));
            if largestPeak==1
                peakInds = usCheck(1:peakStarts(largestPeak));
            elseif largestPeak==numel(peakStarts)+1
                peakInds = usCheck(peakStarts(largestPeak-1)+1:numel(usCheck));                                        
            else
                peakInds = usCheck(peakStarts(largestPeak-1)+1:peakStarts(largestPeak));                   
            end
%             subplot(212);
%             plot(1:numel(usData),usData,'k',usCheck,usData(usCheck),'ro',peakInds,usData(peakInds),'c.');            
        else
            peakInds = usCheck;
%             subplot(212);
%             plot(1:numel(usData),usData,'k',usCheck,usData(usCheck),'ro',peakInds,usData(peakInds),'c.');
        end
        
        respQuant.us.latency(i) = peakInds(1);
        respQuant.us.peak(i)    = peakInds(1) + find(usData(peakInds)==max(usData(peakInds)));
        respQuant.us.amp(i)     = trapz(usData(peakInds))./numel(peakInds);
        respQuant.us.restore(i) = peakInds(numel(peakInds));
        respQuant.us.max(i)     = max(usData(peakInds));
            
    else
        
        respQuant.us.latency(i) = -1;
        respQuant.us.peak(i)    = -1;
        respQuant.us.amp(i)     = -1;
        respQuant.us.restore(i) = -1;
        respQuant.us.max(i)     = -1;
        
    end
    
end
disp('Completed mean peak detection routine | US.');

for i = 1:size(allCSEXaligned.psthZ,2)
       
    csData = allCSaligned.psthZ(winCS,i);
    exData = allCSEXaligned.psthZ(winCS,i);
    exCheck  = find(exData > threshold       );
    exCheck2 = find(exData > (max(exData)./1.5));
    
    if numel(exCheck)>10 & numel(exCheck2)<150 % 2x the s.d. of the smoothing kernel 

%         figure(50); clf; hold on;
%         subplot(211);
%         plot(1:numel(exData),exData,'k',exCheck,exData(exCheck),'ro',find(exData==max(exData)),max(exData),'b.');
%         title(num2str(i));
        
        % go through exCheck and find phasic peaks
        peakStarts = find(diff(exCheck)>1);
        
        if numel(peakStarts)>0 % multiple peaks
            % find the largest peak
            candidatePeaks = zeros(1,numel(peakStarts)+1);
            for k = 1:numel(peakStarts)+1
                if k==1
                    candidatePeaks(1,k) = trapz(exData(exCheck(1:peakStarts(k)))) ./ numel(1:peakStarts(k));
                elseif k==numel(peakStarts)+1
                    candidatePeaks(1,k) = trapz(exData(exCheck(peakStarts(k-1)+1:numel(exCheck)))) ./ numel(peakStarts(k-1)+1:numel(exCheck));                                        
                else
                    candidatePeaks(1,k) = trapz(exData(exCheck(peakStarts(k-1)+1:peakStarts(k)))) ./ numel(peakStarts(k-1)+1:peakStarts(k));                    
                end
            end

            largestPeak = find(candidatePeaks==max(candidatePeaks));
            if largestPeak==1
                peakInds = exCheck(1:peakStarts(largestPeak));
            elseif largestPeak==numel(peakStarts)+1
                peakInds = exCheck(peakStarts(largestPeak-1)+1:numel(exCheck));                                        
            else
                peakInds = exCheck(peakStarts(largestPeak-1)+1:peakStarts(largestPeak));                   
            end
%             subplot(212);
%             plot(1:numel(exData),exData,'k',exCheck,exData(exCheck),'ro',peakInds,exData(peakInds),'c.');            
        else
            peakInds = exCheck;
%             subplot(212);
%             plot(1:numel(exData),exData,'k',exCheck,exData(exCheck),'ro',peakInds,exData(peakInds),'c.');
        end
        
        respQuant.ex.latency(i) = peakInds(1);
        respQuant.ex.peak(i)    = peakInds(1) + find(exData(peakInds)==max(exData(peakInds)));
        respQuant.ex.amp(i)     = trapz(exData(peakInds))./numel(peakInds);
        respQuant.ex.restore(i) = peakInds(numel(peakInds));
        respQuant.ex.max(i)     = max(exData(peakInds));
        respQuant.ex.contrast(i)= trapz(exData(peakInds)-csData(peakInds));
            
    else
        
        respQuant.ex.latency(i) = -1;
        respQuant.ex.peak(i)    = -1;
        respQuant.ex.amp(i)     = -1;
        respQuant.ex.restore(i) = -1;
        respQuant.ex.max(i)     = -1;
        respQuant.ex.contrast(i)= -1;
    end    
    
end

disp('Completed mean peak detection routine | EXTINCTION.');

%% SORT PSTH RESPONSES v2

%eval('home');
disp('BEGIN >>> sorting psth responses...');

% basic algorithm: 
% 1) for all cells look for a significant phasic response either to the CS or to the US
% 2) use the time to the latency to the peak for the larger of the response
% 3) try to classify according into separate groups: 
%       a) latency groups (3)
%           b) if there is both a cs and us response sort according to the difference in the amplitudes
%           a) CS only response
% US latencies are ~15ms delayed from CS responses

% find all units with a significant response during either the cs or us
totInds    = find( respQuant.cs.peak>0 | respQuant.us.peak>0 | respQuant.ex.peak>0 );
tonInds    = find( respQuant.cs.peak>0 | respQuant.ex.peak>0                       );
rewInds    = find( respQuant.us.peak>0                                             );
noRespInds = find( respQuant.cs.peak<0 & respQuant.us.peak<0 & respQuant.ex.peak<0 );

cellClasses.latGroup(4).inds = noRespInds;

for i = 1:numel(respQuant.ex.peak)

    if numel(find(totInds==i)) > 0
        
        testVals = [respQuant.cs.amp(i) respQuant.us.amp(i) respQuant.ex.amp(i)];
        
        biggest = find(testVals==max(testVals));
        
        switch biggest
            case 1
                bestPeakEst(i) = respQuant.cs.peak(i);
            case 2
                bestPeakEst(i) = respQuant.us.peak(i) - 15;
            case 3
                bestPeakEst(i) = respQuant.ex.peak(i);
        end
        
    else
        
        bestPeakEst(i)  = -1;
        classifier(i)   = 4;

    end

    % begin classifier
    if bestPeakEst(i) >   1 & bestPeakEst(i) < 20
        classifier(i) = 3;
    end
    
    if bestPeakEst(i) >= 20 & bestPeakEst(i) < 60
        classifier(i) = 2;
    end

    if bestPeakEst(i) >= 60 & bestPeakEst(i) < 122
        classifier(i) = 1;
    end
    
    if numel(find(PopData.daList==i))>0
        if classifier(i)==2
            disp('Changed a unit from 2 to 1');
            classifier(i) = 1;
        end
    end

end

% now sort indices within each classifier
displayIndexing=[];

for j = 1:4
    
    % get all indices of the classifier and create a new index
    groupIndsR  = find(classifier==j);
    groupIndsI  = 1:numel(groupIndsR);
%     groupAmps   = respQuant.cs.amp(groupIndsR);
%     groupAmps   = respQuant.cs.amp(groupIndsR) - respQuant.ex.amp(groupIndsR);
    groupAmps   = -respQuant.ex.contrast(groupIndsR);
%     groupAmps   = respQuant.cs.peak(groupIndsR) - respQuant.ex.peak(groupIndsR);
%     groupAmps   = respQuant.ex.amp(groupIndsR);

    % store the class index lists
    cellClasses.latGroup(j).inds    = groupIndsR;
    cellClasses.dopamine.inds       = PopData.daList;
    
    % generate sorted indices of the classifier and create a new index
    [sV,sI] = sort(groupAmps,'descend');
    groupIndsS  = groupIndsR(groupIndsI(sI));

    displayIndexing = [displayIndexing groupIndsS];
    
end

csRespIndex = (respQuant.cs.amp-respQuant.us.amp) ./ (abs(respQuant.cs.amp)+abs(respQuant.us.amp));
exRespIndex = (respQuant.ex.amp-respQuant.cs.amp) ./ (abs(respQuant.ex.amp)+abs(respQuant.cs.amp));

figure(23); clf;
[mapName] = TNC_CreateRBColormap(1024,'cb');

subplot(131);
imagesc(allCSaligned.psthZ(1900:2200,displayIndexing)',[-15 15])
colormap(mapName); title('CS'); set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
    % to export
    allCSaligned.sortedZ = allCSaligned.psthZ(:,displayIndexing);

subplot(132);
imagesc(allUSaligned.psthZ(1900:2200,displayIndexing)',[-15 15])
colormap(mapName); title('US'); set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
    % to export
    allUSaligned.sortedZ = allUSaligned.psthZ(:,displayIndexing);

subplot(133);
imagesc(allCSEXaligned.psthZ(1900:2200,displayIndexing)'-allCSaligned.psthZ(1900:2200,displayIndexing)',[-15 15])
colormap(mapName); title('CS_e_x - CS_a_q'); set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
    % to export
    allCSEXaligned.sortedZ = allCSEXaligned.psthZ(:,displayIndexing);


figure(24); clf;
subplot(211);
% subplot(5,5,[1:3 6:8 11:13]);
% plot3(respQuant.cs.peak(totInds),respQuant.cs.amp(totInds)-respQuant.us.amp(totInds),respQuant.ex.amp(totInds)-respQuant.cs.amp(totInds),'k.'); grid on;
% 
% subplot(5,5,[4:5 9:10 14:15]);
plot3(respQuant.cs.peak(totInds)./150,csRespIndex(totInds),exRespIndex(totInds),'k.'); grid on;
zlabel('EXT-ACQ'); ylabel('CS-US'); xlabel('Latency');

subplot(212);
plot(csRespIndex(totInds),exRespIndex(totInds),'k.'); grid on;
ylabel('EXT-ACQ'); xlabel('CS-US');

disp(['Found ' num2str(numel(totInds)) ' out of ' num2str(numel(respQuant.ex.peak)) ' units with a significant response to the stimulus or reward under acq or ex.']);
disp('COMPLETED. ');

%% RECALCULATE THE RESPONSE AMPLITUDES USING THE MAXIMAL RESPONSE

disp('begin.');

for i = 1:size(allCSaligned.psthZ,2)

    winCS = 2001:2140;
    winUS = winCS + 15;

    csData = allCSaligned.psthZ(winCS,i);
    usData = allUSaligned.psthZ(winUS,i);
    exData = allCSEXaligned.psthZ(winCS,i);

    csResp     = respQuant.cs.amp(i);
    usResp     = respQuant.us.amp(i);
    exResp     = respQuant.ex.amp(i);
%     csResp      = max(csData);
%     usResp      = max(usData);
%     exResp      = max(exData);

    testVals    = [csResp exResp];

    csRespL     = respQuant.cs.amp(i);
    usRespL     = respQuant.us.amp(i);
    exRespL     = respQuant.ex.amp(i);

    logicVals   = [csRespL exRespL usRespL];

    if max(logicVals)>0
        biggest = find(logicVals==max(logicVals));

        switch biggest
            case 1
                respWin = respQuant.cs.latency(i):respQuant.cs.restore(i);
                    fullRespQuant.cs.amp(i) = trapz(csData(respWin))./numel(respWin);
                    fullRespQuant.us.amp(i) = trapz(usData(respWin))./numel(respWin);
                    fullRespQuant.ex.amp(i) = trapz(exData(respWin))./numel(respWin);

                    fullRespQuant.cs.max(i) = max(csData(respWin));
                    fullRespQuant.us.max(i) = max(usData(respWin));
                    fullRespQuant.ex.max(i) = max(exData(respWin));

                    if respQuant.cs.peak(i)==0
                        disp('something is wrong');
                    end
                    
                    fullRespQuant.cs.peak(i) = respQuant.cs.peak(i);
                    fullRespQuant.us.peak(i) = respQuant.cs.peak(i);
                    fullRespQuant.ex.peak(i) = respQuant.cs.peak(i);

                    fullRespQuant.window(i).respWin = respWin;


            case 2
                respWin = respQuant.ex.latency(i):respQuant.ex.restore(i);
                    fullRespQuant.cs.amp(i) = trapz(csData(respWin))./numel(respWin);
                    fullRespQuant.us.amp(i) = trapz(usData(respWin))./numel(respWin);
                    fullRespQuant.ex.amp(i) = trapz(exData(respWin))./numel(respWin);

                    fullRespQuant.cs.max(i) = max(csData(respWin));
                    fullRespQuant.us.max(i) = max(usData(respWin));
                    fullRespQuant.ex.max(i) = max(exData(respWin));

                    if respQuant.ex.peak(i)==0
                        disp('something is wrong');
                    end

                    fullRespQuant.cs.peak(i) = respQuant.ex.peak(i);
                    fullRespQuant.us.peak(i) = respQuant.ex.peak(i);
                    fullRespQuant.ex.peak(i) = respQuant.ex.peak(i);

                    fullRespQuant.window(i).respWin = respWin;
                    
            case 3
                respWin = respQuant.us.latency(i):respQuant.us.restore(i);
                    fullRespQuant.cs.amp(i) = trapz(csData(respWin-13))./numel(respWin);
                    fullRespQuant.us.amp(i) = trapz(usData(respWin))./numel(respWin);
                    fullRespQuant.ex.amp(i) = trapz(exData(respWin-13))./numel(respWin);

                    fullRespQuant.cs.max(i) = max(csData(respWin));
                    fullRespQuant.us.max(i) = max(usData(respWin));
                    fullRespQuant.ex.max(i) = max(exData(respWin));

                    if respQuant.us.peak(i)==0
                        disp('something is wrong');
                    end
                    
                    fullRespQuant.cs.peak(i) = respQuant.us.peak(i);
                    fullRespQuant.us.peak(i) = respQuant.us.peak(i);
                    fullRespQuant.ex.peak(i) = respQuant.us.peak(i);

                    fullRespQuant.window(i).respWin = respWin;

            otherwise
                disp('this shouldn,t happen.')

        end

%                     figure(25); clf;
%                     plot([respQuant.cs.peak(i) respQuant.cs.peak(i)],[0 max(csData)],'g-',winCS-2000,csData,'k',winCS-2000,usData,'r',winCS-2000,exData,'b',[respWin(1,1) respWin(1,1)],[-10 30],'k--',[respWin(numel(respWin)) respWin(numel(respWin))],[-10 30],'k--');
%                     title(num2str(i));
%                     drawnow; pause(1);

    else
        
        fullRespQuant.cs.amp(i) = -1;
        fullRespQuant.us.amp(i) = -1;
        fullRespQuant.ex.amp(i) = -1;
        fullRespQuant.window(i).respWin = [];
        
    end
    
end

disp('completed full resp calculation.');

%% CLASSIFY CELLS INTO 7 MAJOR GROUPS v2

disp(' ');disp(' ');
disp('____________________________');
display('Begin final classifier...');

% Group 4 
% No response
currInds = cellClasses.latGroup(4).inds;
disp(['Total of Group 4: ' num2str(numel(currInds))]);

% Group 3
% Short, CS stable
% Short, Extinction response
currInds = cellClasses.latGroup(3).inds;
cellClasses.latGroup(3).sub(1).inds = [];
cellClasses.latGroup(3).sub(2).inds = [];
cellClasses.latGroup(3).sub(3).inds = [];
csAmps = fullRespQuant.cs.max(currInds);
usAmps = fullRespQuant.us.max(currInds);
exAmps = fullRespQuant.ex.max(currInds);
for j=1:numel(currInds)
    if csAmps(j)./usAmps(j) < 0.67 | csAmps(j)./usAmps(j) > 1.5
        cellClasses.latGroup(3).sub(3).inds = [cellClasses.latGroup(3).sub(3).inds currInds(j)]; % sensory
    else
        if csAmps(j) > usAmps(j)
            cellClasses.latGroup(3).sub(1).inds = [cellClasses.latGroup(3).sub(1).inds currInds(j)];  % cs
        else
            cellClasses.latGroup(3).sub(2).inds = [cellClasses.latGroup(3).sub(2).inds currInds(j)];  % us 
        end
    end
end
disp(['Total of Group 3: ' num2str(numel(currInds)) ' | CS: ' num2str(numel(cellClasses.latGroup(3).sub(1).inds)) ' US: ' num2str(numel(cellClasses.latGroup(3).sub(2).inds)) ' SENSORY: ' num2str(numel(cellClasses.latGroup(3).sub(3).inds))]);

% Group 2
% Med, Stable
% Med, Ext
currInds = cellClasses.latGroup(2).inds;
csAmps = fullRespQuant.cs.amp(currInds);
exAmps = fullRespQuant.ex.amp(currInds);
cellClasses.latGroup(2).sub(1).inds = [];
cellClasses.latGroup(2).sub(2).inds = [];

for j=1:numel(currInds)
    if csAmps(j) ./ exAmps(j) > 1
        cellClasses.latGroup(2).sub(1).inds = [cellClasses.latGroup(2).sub(1).inds currInds(j)];  % acq 
    else
        cellClasses.latGroup(2).sub(2).inds = [cellClasses.latGroup(2).sub(2).inds currInds(j)];  % ext
    end
end
disp(['Total of Group 2: ' num2str(numel(currInds)) ' | ACQ: ' num2str(numel(cellClasses.latGroup(2).sub(1).inds)) ' EXT: ' num2str(numel(cellClasses.latGroup(2).sub(2).inds))]);

% Group 1
% Long, DA
% Long, Non-DA
currInds = cellClasses.latGroup(1).inds;
cellClasses.latGroup(1).sub(1).inds = [];
cellClasses.latGroup(1).sub(2).inds = [];

for j=1:numel(currInds)
    if numel(find(PopData.daList == currInds(j)))>0
        cellClasses.latGroup(1).sub(1).inds = [cellClasses.latGroup(1).sub(1).inds currInds(j)];
    else
        cellClasses.latGroup(1).sub(2).inds = [cellClasses.latGroup(1).sub(2).inds currInds(j)];        
    end
end
disp(['Total of Group 1: ' num2str(numel(currInds)) ' | DA: ' num2str(numel(cellClasses.latGroup(1).sub(1).inds)) ' NDA: ' num2str(numel(cellClasses.latGroup(1).sub(2).inds))]);

display('Classifier completed.');
disp(' ');disp(' ');

%% FIND CLASS EXEMPLARS v2

winCS = 1950:2140;
winUS = winCS+10;

clear corrScores;

groupNum    = 2
subNum      = 2

currInds = cellClasses.latGroup(groupNum).sub(subNum).inds;

% class mean PSTH
meanResponse = [mean(allCSaligned.psthZ(winCS,currInds),2)' mean(allUSaligned.psthZ(winUS,currInds),2)'];
meanResponseE= [std(allCSaligned.psthZ(winCS,currInds),[],2)'./sqrt(numel(currInds)) , std(allUSaligned.psthZ(winUS,currInds),[],2)'./sqrt(numel(currInds))];

figure(51); clf; plot(meanResponse,'k'); hold on; plot(meanResponseE,'k--');
cellClasses.latGroup(groupNum).sub(subNum).meanResponse = meanResponse;
cellClasses.latGroup(groupNum).sub(subNum).meanResponseE = meanResponseE;

% unit-by-unit correlation
for i=1:numel(currInds)
    currResponse    = [allCSaligned.psthZ(winCS,currInds(i))' allUSaligned.psthZ(winUS,currInds(i))'];
    if groupNum == 2
        corrScores(i)   = max(currResponse(1:140));
    else
        corrScores(i)   = corr2(currResponse,meanResponse);
    end
    
    currSess = allCSaligned.sess(currInds(i));
    currUnit = allCSaligned.unit(currInds(i));
    rasterDS = PopData.session(currSess).unit(currUnit).respCS.raster;

    figure(51); clf;
    TNC_PlotRaster(51,rasterDS);
    title(num2str(currInds(i)));
    pause();
%     figure(51); plot(1:numel([winCS winUS]),meanResponse,'k',1:numel([winCS winUS]),currResponse,'r');
%     drawnow;
end

% max correlation index
max(corrScores)
find(corrScores==max(corrScores),1)
currInds(find(corrScores==max(corrScores),1))
cellClasses.latGroup(groupNum).exemplar(subNum).ind = currInds(find(corrScores==max(corrScores),1));

figure(51); plot(1:numel([winCS winUS]),meanResponse,'k',1:numel([winCS winUS]),[allCSaligned.psthZ(winCS,cellClasses.latGroup(groupNum).exemplar(subNum).ind)' allUSaligned.psthZ(winUS,cellClasses.latGroup(groupNum).exemplar(subNum).ind)'],'r');

%% EXPORT EXEMPLARS v2
m = 1;
saveFlag = 1;

% list of example indices:
% examples = [cellClasses.latGroup(3).exemplar(1).ind,
%     cellClasses.latGroup(3).exemplar(2).ind,
%     cellClasses.latGroup(3).exemplar(3).ind,
%     cellClasses.latGroup(2).exemplar(1).ind,
%     cellClasses.latGroup(2).exemplar(2).ind,
%     cellClasses.latGroup(1).exemplar(1).ind,
%     cellClasses.latGroup(1).exemplar(2).ind]
% numExamples = numel(examples);

examples = [cellClasses.latGroup(3).exemplar(1).ind,
    cellClasses.latGroup(3).exemplar(2).ind,
    cellClasses.latGroup(3).exemplar(3).ind,
    496,
    473,
    cellClasses.latGroup(1).exemplar(1).ind,
    cellClasses.latGroup(1).exemplar(2).ind]
numExamples = numel(examples);

exampMeans = [cellClasses.latGroup(3).sub(1).meanResponse;
    cellClasses.latGroup(3).sub(2).meanResponse;
    cellClasses.latGroup(3).sub(3).meanResponse;
    cellClasses.latGroup(2).sub(1).meanResponse;
    cellClasses.latGroup(2).sub(2).meanResponse;
    cellClasses.latGroup(1).sub(1).meanResponse;
    cellClasses.latGroup(1).sub(2).meanResponse]';

exampMeanE = [cellClasses.latGroup(3).sub(1).meanResponseE;
    cellClasses.latGroup(3).sub(2).meanResponseE;
    cellClasses.latGroup(3).sub(3).meanResponseE;
    cellClasses.latGroup(2).sub(1).meanResponseE;
    cellClasses.latGroup(2).sub(2).meanResponseE;
    cellClasses.latGroup(1).sub(1).meanResponseE;
    cellClasses.latGroup(1).sub(2).meanResponseE]';

for j=1:numExamples

    % grab the session and unit information
    currIndex = examples(j)
    currSess = allCSaligned.sess(currIndex);
    currUnit = allCSaligned.unit(currIndex);

    % index into the PopData and concatenate the timestamps
    numTrials = size(PopData.session(currSess).unit(currUnit).respCS.raster.trial,2)
    
    for k=1:numTrials
    
        % create a equal sized trial id vector
        thisTrialTS     = PopData.session(currSess).unit(currUnit).respCS.raster.trial(k).ts;
        numStamps       = size(thisTrialTS,1);
        thisTrialIDS    = ones(numStamps,1).*k;
        
        if k==1
            finalArrayTS = thisTrialTS;            
            finalArrayID = thisTrialIDS;            
        else
            finalArrayTS = [finalArrayTS;thisTrialTS];
            finalArrayID = [finalArrayID;thisTrialIDS];
        end

        % create a equal sized trial id vector
        thisTrialTSu     = PopData.session(currSess).unit(currUnit).respUS.raster.trial(k).ts;
        numStampsu       = size(thisTrialTSu,1);
        thisTrialIDSu    = ones(numStampsu,1).*k;
        
        if k==1
            finalArrayTSu = thisTrialTSu;            
            finalArrayIDu = thisTrialIDSu;            
        else
            finalArrayTSu = [finalArrayTSu;thisTrialTSu];
            finalArrayIDu = [finalArrayIDu;thisTrialIDSu];
        end

        
    end
    
    figure(1); subplot(4,1,1:3);
    plot(finalArrayTS,finalArrayID,'k.','MarkerSize',1);
    subplot(4,1,4);
    plot(-2000:3000, allCSaligned.psthZ(:,currIndex),'k');
    title(currIndex);
    drawnow; %pause(2);
    
    if saveFlag
        
        m = j;
        
        % save to an hdf5 structure for reading/plotting in igor
        nameA = sprintf('/spkXcs%g',0);
        nameB = sprintf('/spkYcs%g',0);
        nameAA = sprintf('/spkXus%g',0);
        nameBB = sprintf('/spkYus%g',0);
        nameC = sprintf('/boxCS%g',0);
        nameCC = sprintf('/boxUS%g',0);
        nameD = sprintf('/psthMcs%g',0);
        nameE = sprintf('/psthEcs%g',0);
        nameDD = sprintf('/psthMus%g',0);
        nameEE = sprintf('/psthEus%g',0);

%         if j==1
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameA,finalArrayTS);
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameB,finalArrayID,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameAA,finalArrayTSu,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameBB,finalArrayIDu,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameC,allCSaligned.boxcar(:,currIndex),'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameD,allCSaligned.psthZ(:,currIndex),'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameE,allCSaligned.psthZe(:,currIndex),'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameCC,allUSaligned.boxcar(:,currIndex),'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameDD,allUSaligned.psthZ(:,currIndex),'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameEE,allUSaligned.psthZe(:,currIndex),'WriteMode','append');
%         else
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameA,finalArrayTS,'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameB,finalArrayID,'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameAA,finalArrayTSu,'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameBB,finalArrayIDu,'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameC,allCSaligned.boxcar(:,currIndex),'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameD,allCSaligned.psthZ(:,currIndex),'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameE,allCSaligned.psthZe(:,currIndex),'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameCC,allUSaligned.boxcar(:,currIndex),'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameDD,allUSaligned.psthZ(:,currIndex),'WriteMode','append');
%             hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameEE,allUSaligned.psthZe(:,currIndex),'WriteMode','append');
%         end
    end
    
end

%% EXPORT RASTER, PSTH, AND BOXCAR FOR GIVEN INDEX

currIndex = 529;
% ex short:602
% da eg: 540
% create temporary waves for the RASTER CS and US

    % grab the session and unit information
    currSess = allCSaligned.sess(currIndex);
    currUnit = allCSaligned.unit(currIndex);

    % index into the PopData and concatenate the timestamps
    numTrialsACQ = size(PopData.session(currSess).unit(currUnit).respCS.raster.trial,2);
    
    for k=1:numTrialsACQ
    
        % create a equal sized trial id vector
        thisTrialTS     = PopData.session(currSess).unit(currUnit).respCS.raster.trial(k).ts;
        numStamps       = size(thisTrialTS,1);
        thisTrialIDS    = ones(numStamps,1).*k;
        
        if k==1
            finalArrayTS = thisTrialTS;            
            finalArrayID = thisTrialIDS;            
        else
            finalArrayTS = [finalArrayTS;thisTrialTS];
            finalArrayID = [finalArrayID;thisTrialIDS];
        end

        % create a equal sized trial id vector
        thisTrialTSu     = PopData.session(currSess).unit(currUnit).respUS.raster.trial(k).ts;
        numStampsu       = size(thisTrialTSu,1);
        thisTrialIDSu    = ones(numStampsu,1).*k;
        
        if k==1
            finalArrayTSu = thisTrialTSu;            
            finalArrayIDu = thisTrialIDSu;            
        else
            finalArrayTSu = [finalArrayTSu;thisTrialTSu];
            finalArrayIDu = [finalArrayIDu;thisTrialIDSu];
        end

        
    end
    
    psthZcs = allCSaligned.psthZ(:,currIndex);
    psthZus = allUSaligned.psthZ(:,currIndex);
    
    figure(1); subplot(4,3,[1 4 7]);
    plot(finalArrayTS,finalArrayID,'k.','MarkerSize',2);
    title(currIndex);
    axis([-50 250 0 numTrialsACQ]);
    subplot(4,3,[2 5 8]);
    plot(finalArrayTSu,finalArrayIDu,'k.','MarkerSize',2);
    axis([-50 250 0 numTrialsACQ]);
    subplot(4,3,10);
    plot(-50:250, allCSaligned.psthZ(1950:2250,currIndex),'k');
    axis([-50 250 -5 30]);
    subplot(4,3,11);
    plot(-50:250, allUSaligned.psthZ(1950:2250,currIndex),'k');
    axis([-50 250 -5 30]);
    
% create temporary waves for the RASTER EX

    % grab the session and unit information
    currSessE = allCSEXaligned.sess(currIndex);
    currUnitE = allCSEXaligned.unit(currIndex);

    % index into the PopData and concatenate the timestamps
    numTrialsEX = size(PopData.session(currSessE).unit(currUnitE).respCS.raster.trial,2);
    
    for k=1:numTrialsEX
    
        % create a equal sized trial id vector
        thisTrialTSe     = PopData.session(currSessE).unit(currUnitE).respCS.raster.trial(k).ts;
        numStampse       = size(thisTrialTSe,1);
        thisTrialIDSe    = ones(numStampse,1).*k;
        
        if k==1
            finalArrayTSe = thisTrialTSe;            
            finalArrayIDe = thisTrialIDSe;            
        else
            finalArrayTSe = [finalArrayTSe;thisTrialTSe];
            finalArrayIDe = [finalArrayIDe;thisTrialIDSe];
        end

        
    end
    
    psthZex = allCSEXaligned.psthZ(:,currIndex);
    
    figure(1);
     
    subplot(4,3,[3 6 9]);
    plot(finalArrayTSe,finalArrayIDe,'k.','MarkerSize',2);
    axis([-50 250 0 numTrialsEX]);
    title('Aligned to CS extinction');
    subplot(4,3,12);
    plot(-50:250, allCSEXaligned.psthZ(1950:2250,currIndex),'k');
    axis([-50 250 -5 30]);

% create h5 file
    name = '/rasterXcs';
    hdf5write('../../export_ExampCell.h5',name,finalArrayTS);
    name = '/rasterYcs';
    hdf5write('../../export_ExampCell.h5',name,finalArrayID,'WriteMode','append');
    name = '/rasterXus';
    hdf5write('../../export_ExampCell.h5',name,finalArrayTSu,'WriteMode','append');
    name = '/rasterYus';
    hdf5write('../../export_ExampCell.h5',name,finalArrayIDu,'WriteMode','append');
    name = '/rasterXex';
    hdf5write('../../export_ExampCell.h5',name,finalArrayTSe,'WriteMode','append');
    name = '/rasterYex';
    hdf5write('../../export_ExampCell.h5',name,finalArrayIDe,'WriteMode','append');

    name = '/psthZcs';
    hdf5write('../../export_ExampCell.h5',name,psthZcs,'WriteMode','append');
    name = '/psthZus';
    hdf5write('../../export_ExampCell.h5',name,psthZus,'WriteMode','append');
    name = '/psthZex';
    hdf5write('../../export_ExampCell.h5',name,psthZex,'WriteMode','append');

    name = '/matlabIndex';
    hdf5write('../../export_ExampCell.h5',name,currIndex,'WriteMode','append');

%% EXPORT ALIGNED PSTHS v2
clear currDat*
expWin = 1900:2200;

for i = 1:size(allCSaligned.sortedZ,2)
    
    currData    = allCSaligned.sortedZ(expWin,i);
    currData1   = allUSaligned.sortedZ(expWin,i);
    currData2   = allCSEXaligned.sortedZ(expWin,i);
    currData3   = allCSEXaligned.sortedZ(expWin,i) - allCSaligned.sortedZ(expWin,i);

    if i==1
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../exportAP_CS.h5',name,currData);
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../exportAP_US.h5',name,currData1);
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../exportAP_EX.h5',name,currData2);
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../exportAP_Diff.h5',name,currData3);
    else
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../exportAP_CS.h5',name,currData,'WriteMode','append');
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../exportAP_US.h5',name,currData1,'WriteMode','append');
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../exportAP_EX.h5',name,currData2,'WriteMode','append');
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../exportAP_Diff.h5',name,currData3,'WriteMode','append');
    end 
        
end

%% CALCULATE TRIALWISE CHANGES IN SPIKE COUNT v2

disp('>>> Build the respWin into the PopData structure...');
NumSessions = size(PopData.session,2);
k=1;

for i = 1:size(fullRespQuant.window,2)

    % Look up the session and unit information
    currSession     = allCSaligned.sess(i);
    currUnit        = allCSaligned.unit(i);
    
    disp([num2str(currSession) ' u' num2str(currUnit)]);
    
    % Store the phasic response window in the PopData structure
    if fullRespQuant.cs.amp(i) ~= -1
        PopData.session(currSession).unit(currUnit).respWin = fullRespQuant.window(i).respWin;
    else
        PopData.session(currSession).unit(currUnit).respWin = [];
    end        
    
    % Look up the session and unit information
    currSession     = allCSEXaligned.sess(i);
    currUnit        = allCSEXaligned.unit(i);
    
    % Store the phasic response window in the PopData structure
    if fullRespQuant.ex.amp(i) ~= -1
        PopData.session(currSession).unit(currUnit).respWin = fullRespQuant.window(i).respWin;
    else
        PopData.session(currSession).unit(currUnit).respWin = [];
    end        
    
    
end
disp('Completed. <<<');

for i=1:NumSessions
    
    numUnits = size(PopData.session(i).unit,2);

    for j = 1:numUnits

        numTrials   = size(PopData.session(i).unit(j).respCS.raster.trial,2);
        currWin     = PopData.session(i).unit(j).respWin;
        
        if numel(currWin)>0

            % need to expand this window to deal with the jittery trial by trial spiking.
            left        = 0; 
            right       = currWin(numel(currWin)) + ((currWin(numel(currWin))-currWin(1)) * 0.5);

            % Go through all trials and count the number of spikes or the summed response
            numStamps   = length(PopData.session(i).unit(j).ts);
            delta       = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
            delta(1,round(PopData.session(i).unit(j).ts)+1) = 1;
            tmpSmooth   = conv(delta,PopData.currParams.filter.kernel,'same');
            [spkResp] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,currParams.stableTrials,PopData.session(i).events.CS.ts,[PopData.currParams.winParams.prior,PopData.currParams.winParams.after],1,1);

            for m=1:numTrials

                thisTrialStamps = PopData.session(i).unit(j).respCS.raster.trial(m).ts;

                % count spikes in the phasic response window
                PopData.session(i).unit(j).phasicCS.trialwise.spkCnt(1,m)   = numel(find(thisTrialStamps > left & thisTrialStamps < right));

                % sum the smoothed spike data over the window
                PopData.session(i).unit(j).phasicCS.trialwise.smthCnt(1,m)  = sum(spkResp.image.aligned(m,currWin+PopData.currParams.winParams.prior),2);

            end        
        else
            % No phasic response was detected
            PopData.session(i).unit(j).phasicCS.trialwise.spkCnt = [];
            PopData.session(i).unit(j).phasicCS.trialwise.smthCnt= [];
        end
        
        k = k+1;
    
    end
    
    disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' | trials: ' num2str(numTrials)]);
    
end

k = k-1;
disp(['Total units: ' num2str(k)]);

%% COLLAPSE LEARNING - EXTINCTION PHASIC RESPONSE DATA INTO FLAT STRUCTURE 
% Collect all of the individual trialwise data into a large matrix that can be used to match units and do statistical
% testing.
% Also collect matrix for all of the stable session data for calculating contrast measures of significant, non-significant,
% and dopamine neurons.
clear trialWise;
NumSessions = size(PopData.session,2);
k=1;m=1;numValid = 0;

trialWise.acq.spkCnt        = zeros(663,200);
trialWise.acq.smthCnt       = zeros(663,200);

trialWise.ext.spkCnt        = zeros(663,200);
trialWise.ext.smthCnt       = zeros(663,200);


for i=1:NumSessions
    
    if strcmp(PopData.session(i).sessClass,'learning')
            
        numUnits = size(PopData.session(i).unit,2);

        for j = 1:numUnits

            % trialwise data
            if numel(PopData.session(i).unit(j).phasicCS.trialwise.spkCnt) > 0
                trialWise.acq.numTrials(k)                  = size(PopData.session(i).unit(j).phasicCS.trialwise.spkCnt,2);
                trialWise.acq.spkCnt(k,1:trialWise.acq.numTrials(k))         = PopData.session(i).unit(j).phasicCS.trialwise.spkCnt(1,1:trialWise.acq.numTrials(k));
                trialWise.acq.smthCnt(k,1:trialWise.acq.numTrials(k))        = PopData.session(i).unit(j).phasicCS.trialwise.smthCnt(1,1:trialWise.acq.numTrials(k));
            else
                trialWise.acq.numTrials(k)                  = -1;
                trialWise.acq.spkCnt(k,1)                   = -1;
                trialWise.acq.smthCnt(k,1)                  = -1;                
            end
            
            % increment the unit counter
            k = k+1;

        end
        
    elseif strcmp(PopData.session(i).sessClass,'extinction')
            
        numUnits = size(PopData.session(i).unit,2);

        for l = 1:numUnits
            
            % trialwise data
            if numel(PopData.session(i).unit(l).phasicCS.trialwise.spkCnt) > 0
                trialWise.ext.numTrials(m)                  = size(PopData.session(i).unit(l).phasicCS.trialwise.spkCnt,2);
                trialWise.ext.spkCnt(m,1:trialWise.ext.numTrials(m))         = PopData.session(i).unit(l).phasicCS.trialwise.spkCnt(1,1:trialWise.ext.numTrials(m));
                trialWise.ext.smthCnt(m,1:trialWise.ext.numTrials(m))        = PopData.session(i).unit(l).phasicCS.trialwise.smthCnt(1,1:trialWise.ext.numTrials(m));
            else
                trialWise.ext.numTrials(m)                  = -1;
                trialWise.ext.spkCnt(m,1)                   = -1;
                trialWise.ext.smthCnt(m,1)                  = -1;                
            end
            
            % increment the unit counter
            m = m+1;

        end
        
    else
        disp('Serious error: unrecognized session class');
    end
    
end

k = k-1; m = m-1;
disp('  ');
disp(['Total units: ' num2str(k) ' in acquisition | ' num2str(m) ' in extinction.']);
disp(['Valid extinction sessions: ' num2str(numValid)]);
disp('  ');

%% CALCULATE SIGNIFICANCE OF RESPONSES v2 

deltaExt.stats = zeros(663,4); % four columns corresponding to: kw test, rs test, significant logic, sign of change

for i = 1:663
    
%     % grab trial data v1
%     responseA = trialWise.acq.spkCnt(i,10:40);
%     responseE = trialWise.ext.spkCnt(i,20:50);

    % grab trial data v2
    responseA = trialWise.acq.spkCnt(i,:);
    responseE = trialWise.ext.spkCnt(i,:);

%     % grab trial data v3
%     responseA = trialWise.acq.spkCnt(i,:);
%     responseE = trialWise.ext.spkCnt(i,1:numel(responseA));

    % check significance
    deltaExt.stats(i,1) = kruskalwallis([responseA responseE], [zeros(1,numel(responseA)) ones(1,numel(responseE))], 'off');
    deltaExt.stats(i,2) = ranksum(responseA,responseE);

    if sessWise.stats(i,1)<0.05 | sessWise.stats(i,2)<0.05

        deltaExt.stats(i,3) = 1;
        deltaExt.stats(i,4) = mean(responseE) - mean(responseA);            
        
    else
        % everything should stay equal to 0
    end
    
    clear responseA responseE

end

tmpSig = find(deltaExt.stats(:,3)==1);
tmpPos = find(deltaExt.stats(:,4)>0);

disp(['Number of significant responses: ' num2str(numel(tmpSig)) ' | Number of increasers: ' num2str(numel(tmpPos))]);

%% PLOT EXTINCTION CHANGES AS A FUNCTION OF CS US CONTRAST

% allRespInds = find(respQuant.cs.latency>0);

allRespInds = cellClasses.latGroup(3).inds;

% allRespInds = PopData.daList;

figure(63); clf;
subplot(211);
plot(fullRespQuant.cs.amp(cellClasses.latGroup(2).inds),fullRespQuant.ex.amp(cellClasses.latGroup(2).inds)-fullRespQuant.cs.amp(cellClasses.latGroup(2).inds),'b.',fullRespQuant.cs.amp(allRespInds)-fullRespQuant.us.amp(allRespInds),fullRespQuant.ex.amp(allRespInds)-fullRespQuant.cs.amp(allRespInds),'k.',fullRespQuant.cs.amp(PopData.daList)-fullRespQuant.us.amp(PopData.daList),fullRespQuant.ex.amp(PopData.daList)-fullRespQuant.cs.amp(PopData.daList),'r.')
hold on;
plot([0 0],[-40 40],'k--',[-10 40],[0 0],'k--');
ylabel('Ext - Acq'); xlabel('CS - US');
% axis([-20 20 -20 20]);
subplot(212);
plot(fullRespQuant.cs.peak, fullRespQuant.ex.amp-fullRespQuant.cs.amp,'k.');

%% EXPORT EXAMPLES OF RESPONSE TYPES
% cs aligned timestamps for spikes are stored in:
% PopData.session(1).unit(1).respCS.raster.trial(1).ts

saveFlag = 1;

% iterate through each exemplar
numExamples = size(allCSaligned.exemplars,1);

% m = 1; % for the early latencies
% m = 2; % for the middle latencies
m = 3; % for the late latencies

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
    plot(-2000:3000, allCSaligned.psthZ(:,currIndex),'k');
    drawnow; pause(0.3);
    
    if saveFlag
        % save to an hdf5 structure for reading/plotting in igor
        nameA = sprintf('/spikeTimesX%g',j-1);
        nameB = sprintf('/spikeTimesY%g',j-1);
        nameC = sprintf('/psth%g',j-1);

        if j==1
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameA,finalArrayTS);
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameB,finalArrayID,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameC,allCSaligned.boxcar(:,currIndex),'WriteMode','append');
        else
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameA,finalArrayTS,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameB,finalArrayID,'WriteMode','append');
            hdf5write(['../../exemplarExport_' num2str(m) '.h5'],nameC,allCSaligned.boxcar(:,currIndex),'WriteMode','append');
        end
    end
    
end

%% FIND PEAK WINDOWS FOR EACH UNIT
% convert to store the data within the overall PopData structure

boxcarWindow = 201:1:225;
clear dataCS dataCSe baseAVGa baseSTDa baseAVGe baseSTDe

% cycle through all units
for i=1:663

    dataCS  = allCSaligned.boxcar(boxcarWindow,i);
    sessA = allCSaligned.sess(i);
    unitA = allCSaligned.unit(i);

    dataCSe = allCSEXaligned.boxcar(boxcarWindow,i);
    sessE = allCSEXaligned.sess(i);
    unitE = allCSEXaligned.unit(i);
    
    baseAVGa = mean(allCSaligned.boxcar(1:199,i));
    baseSTDa = std(allCSaligned.boxcar(1:199,i),0,1);

    baseAVGe = mean(allCSEXaligned.boxcar(1:199,i));
    baseSTDe = std(allCSEXaligned.boxcar(1:199,i),0,1);

    peakVal     = max(dataCS)-baseAVGa;
    peakVale    = max(dataCSe)-baseAVGe;

    % Determine the stronger of the two conditions
    if peakVal >= peakVale
        useLearn = 1;
    else
        useLearn = 0;
    end
    
    peakLoc    = find(dataCS(1:16) == max(dataCS(1:16)),1);
    peakLoce   = find(dataCSe(1:16) == max(dataCSe(1:16)),1);             

    if useLearn
        
        truePeakLoc = peakLoc;
        peakVal     = dataCS(peakLoc)-baseAVGa;
        peakVale    = dataCSe(peakLoc)-baseAVGe;
        
        if peakLoc > 13 | peakLoc==1
        
            disp(['No short latency peak; peak at ' num2str(peakLoc)]);
            peakLoc                = 3;
            peakLoce               = 3;
            truePeakLoc            = 3;
            respProps.fw(i).X(1,1) = 201; % this is the full window of the early latency peaks
            respProps.fw(i).X(1,2) = 206;
            thresh                 = baseAVGa + (baseSTDa.*3);
            respProps.fw(i).Y(1,1) = thresh;
            respProps.fw(i).Y(1,2) = thresh;
            peakVal         = dataCS(truePeakLoc)-baseAVGa;
            peakVale        = dataCSe(truePeakLoc)-baseAVGe;
            confirm = 1;
            
        else
            
            % Try to find where the peak emerges from the baseline
                
            thresh                  = baseAVGa + (baseSTDa.*3);
            respProps.fw(i).Y(1,1)  = thresh;
            respProps.fw(i).Y(1,2)  = thresh;

            % look around the peak to find where it falls back down to baseline
            tmpX = find(dataCS(1:peakLoc-1) < thresh,1,'last');

            if isempty(tmpX)
                respProps.fw(i).X(1,1) = 201;
            else
                respProps.fw(i).X(1,1) = tmpX + 200;
            end

            tmpX = find(dataCS(peakLoc+1:numel(boxcarWindow)) < thresh,1,'first');

            if isempty(tmpX)
                respProps.fw(i).X(1,2) = boxcarWindow(1,numel(boxcarWindow));
            else
                respProps.fw(i).X(1,2) = tmpX + peakLoc + 200;
            end

            confirm = 0;

        end
        % Put up a plot of the current unit to see how well i did to find the interval
%         figure(4); clf; plot(1:size(boxcarWindow,2),dataCS-baseAVGa,'k.-',[0 30],[thresh-baseAVGa thresh-baseAVGa],'r--',peakLoc,peakVal,'bo',respProps.fw(i).X(1,1)-boxcarWindow(1,1)+1,respProps.fw(i).Y(1,1)-baseAVGa,'g+',respProps.fw(i).X(1,2)-boxcarWindow(1,1)+1,respProps.fw(i).Y(1,2)-baseAVGa,'b+','LineWidth',2); axis([0 25 -50 150]);
%         title(num2str(i)); drawnow;
        
        if peakVal ~= dataCS(truePeakLoc)-baseAVGa
            disp('Peak value acq is wrong');
            pause;
        end
        
    else
                
        
        truePeakLoc = peakLoce;
        peakVal     = dataCS(peakLoce)-baseAVGa;
        peakVale    = dataCSe(peakLoce)-baseAVGe;

        if peakLoce > 13 | peakLoce==1
        
            disp(['No short latency peak; peak at ' num2str(peakLoce)]);
            peakLoc                = 3;
            peakLoce               = 3;
            truePeakLoc            = 3;
            respProps.fw(i).X(1,1) = 201; % this is the full window of the early latency peaks
            respProps.fw(i).X(1,2) = 206;
            thresh                 = baseAVGe + (baseSTDe.*3);
            respProps.fw(i).Y(1,1) = thresh;
            respProps.fw(i).Y(1,2) = thresh;
            peakVal         = dataCS(truePeakLoc)-baseAVGa;
            peakVale        = dataCSe(truePeakLoc)-baseAVGe;
            confirm = 1;
            
        else

            thresh                  = baseAVGe + (baseSTDe.*3);
            respProps.fw(i).Y(1,1)  = thresh;
            respProps.fw(i).Y(1,2)  = thresh;

            % look around the peak to find where it falls back down to baseline
            tmpX = find(dataCSe(1:peakLoce-1) < thresh,1,'last');
            if isempty(tmpX)
                respProps.fw(i).X(1,1) = 201;
            else
                respProps.fw(i).X(1,1) = tmpX + 200;
            end

            tmpX = find(dataCSe(peakLoce+1:numel(boxcarWindow)) < thresh,1,'first');

            if isempty(tmpX)
                respProps.fw(i).X(1,2) = boxcarWindow(1,numel(boxcarWindow));
            else
                respProps.fw(i).X(1,2) = tmpX + peakLoce + 200;
            end
            
            confirm = 0;
            
        end
        % Put up a plot of the current unit to see how well i did to find the interval
%         figure(4); clf; plot(1:size(boxcarWindow,2),dataCSe-baseAVGe,'k.-',[0 30],[thresh-baseAVGe thresh-baseAVGe],'r--',peakLoce,peakVale,'bo',respProps.fw(i).X(1,1)-200,respProps.fw(i).Y(1,1)-baseAVGe,'g+',respProps.fw(i).X(1,2)-200,respProps.fw(i).Y(1,2)-baseAVGe,'b+','LineWidth',2); axis([0 25 -50 150]);
%         title(num2str(i)); drawnow;

        if peakVale ~= dataCSe(truePeakLoc)-baseAVGe
            disp('Peak value ext is wrong');
            pause;
        end

    end

%     if confirm
%         pause;
%     else
%         pause(0.2);
%     end
    
    if respProps.fw(i).X(1,2)-respProps.fw(i).X(1,1) < 2
        disp(['ERROR 1: ' num2str(respProps.fw(i).X(1,1)) '__' num2str(respProps.fw(i).X(1,1)) '| unit: ' num2str(unitA) ' session: ' num2str(sessA)]);
        pause;
    end

    if respProps.fw(i).X(1,2) <= (truePeakLoc+200)
        disp(['ERROR 2: ' num2str(respProps.fw(i).X(1,2)) '__' num2str((truePeakLoc+boxcarWindow(1,1))) '| unit: ' num2str(unitA) ' session: ' num2str(sessA)]);
        pause;
    end

    
    PopData.session(sessA).unit(unitA).phasicCS.window.tenMs.left   = respProps.fw(i).X(1,1);    
    PopData.session(sessA).unit(unitA).phasicCS.window.tenMs.right  = respProps.fw(i).X(1,2);
    PopData.session(sessA).unit(unitA).phasicCS.window.ms.left      = (respProps.fw(i).X(1,1).*10)-9;    
    PopData.session(sessA).unit(unitA).phasicCS.window.ms.right     = (respProps.fw(i).X(1,2).*10)-9;
    PopData.session(sessA).unit(unitA).phasicCS.window.peakV        = peakVal;
    PopData.session(sessA).unit(unitA).phasicCS.window.peakL        = truePeakLoc;

    PopData.session(sessE).unit(unitE).phasicCS.window.tenMs.left   = respProps.fw(i).X(1,1);    
    PopData.session(sessE).unit(unitE).phasicCS.window.tenMs.right  = respProps.fw(i).X(1,2);
    PopData.session(sessE).unit(unitE).phasicCS.window.ms.left      = (respProps.fw(i).X(1,1).*10)-9;    
    PopData.session(sessE).unit(unitE).phasicCS.window.ms.right     = (respProps.fw(i).X(1,2).*10)-9;
    PopData.session(sessE).unit(unitE).phasicCS.window.peakV        = peakVale;
    PopData.session(sessE).unit(unitE).phasicCS.window.peakL        = truePeakLoc;
    
    if peakVal >= peakVale
        PopData.session(sessA).unit(unitA).phasicCS.learnLarger     = 1;
        PopData.session(sessE).unit(unitE).phasicCS.learnLarger     = 1;
    else
        PopData.session(sessA).unit(unitA).phasicCS.learnLarger     = 0;
        PopData.session(sessE).unit(unitE).phasicCS.learnLarger     = 0;
    end
    
    winL(i) = (respProps.fw(i).X(1,1).*10)-9;
    winR(i) = (respProps.fw(i).X(1,2).*10)-9;
end

figure(4); clf;
[mapName] = TNC_CreateRBColormap(1024,'rb');
imagesc(allCSaligned.psthZ'+allCSEXaligned.psthZ',[-10 10]);
colormap(mapName); hold on;
plot(winL,1:663,'b>',winR,1:663,'k<');
axis([2000 2250 0 665]);

%% QUANTIFY SESSION RESPONSE TO CS FOR ALL CELLS IN LEARNING AND EXTINCTION

NumSessions = size(PopData.session,2);
k=1;

latWindow.win(1).vals = 2003:1:2024;
latWindow.win(2).vals = 2025:1:2060;
latWindow.win(3).vals = 2061:1:2135;

for i=1:NumSessions
    
    numUnits = size(PopData.session(i).unit,2);

    for j = 1:numUnits

        baselineB    = mean(PopData.session(i).unit(j).respCS.boxcar(1,1:199));
        currWindowB  = PopData.session(i).unit(j).phasicCS.window.tenMs.left:1:PopData.session(i).unit(j).phasicCS.window.tenMs.right;
        PopData.session(i).unit(j).phasicCS.stable.ampb.total       = sum(PopData.session(i).unit(j).respCS.boxcar(currWindowB) - baselineB);
        PopData.session(i).unit(j).phasicCS.stable.ampb.peak        = max(PopData.session(i).unit(j).respCS.boxcar(currWindowB)) - baselineB;

        currWindowZ  = PopData.session(i).unit(j).phasicCS.window.ms.left:1:PopData.session(i).unit(j).phasicCS.window.ms.right;     
        PopData.session(i).unit(j).phasicCS.stable.ampz.int         = trapz(PopData.session(i).unit(j).respCSS.psthZ(1,currWindowZ));
        PopData.session(i).unit(j).phasicCS.stable.ampz.mean        = mean(PopData.session(i).unit(j).respCSS.psthZ(1,currWindowZ));
        PopData.session(i).unit(j).phasicCS.stable.ampz.peak        = max(PopData.session(i).unit(j).respCSS.psthZ(1,currWindowZ));

        % Calculate response magnitudes in the population-defined latency windows
        PopData.session(i).unit(j).phasicCS.stable.peak(1).int      = trapz(PopData.session(i).unit(j).respCSS.psthZ(1,latWindow.win(1).vals));
        PopData.session(i).unit(j).phasicCS.stable.peak(2).int      = trapz(PopData.session(i).unit(j).respCSS.psthZ(1,latWindow.win(1).vals));
        PopData.session(i).unit(j).phasicCS.stable.peak(3).int      = trapz(PopData.session(i).unit(j).respCSS.psthZ(1,latWindow.win(1).vals));
        
        % Reconfirm that the peak is valid.
        if PopData.session(i).unit(j).phasicCS.stable.ampz.mean > 3
            if PopData.session(i).unit(j).phasicCS.window.ms.left < 2140
                PopData.session(i).unit(j).phasicCS.stats.validPeak = 1;
            else
                PopData.session(i).unit(j).phasicCS.stats.validPeak = 0;
            end
        else
            PopData.session(i).unit(j).phasicCS.stats.validPeak = 0;
        end
        
%         numTrials = 1;
        
        % Go through all trials and count the number of spikes or the summed response
        numStamps   = length(PopData.session(i).unit(j).ts);
        delta       = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,round(PopData.session(i).unit(j).ts)+1) = 1;
        tmpSmooth   = conv(delta,PopData.currParams.filter.kernel,'same');

        [spkResp] = TNC_AlignRasters(tmpSmooth,PopData.session(i).unit(j).ts,currParams.stableTrials,PopData.session(i).events.CS.ts,[PopData.currParams.winParams.prior,PopData.currParams.winParams.after],1,1);
         
        numTrials   = size(spkResp.raster.trial,2);

        thisLeft    = PopData.session(i).unit(j).phasicCS.window.ms.left-2000;
        thisRight   = PopData.session(i).unit(j).phasicCS.window.ms.right-2000;

        if thisLeft>thisRight
            disp(['There is an error in the window calculation for unit ' num2str(j) ' in session ' num2str(i)]);
        end
        
        for m=1:numTrials

            % calculate baseline (spikes / ms)
            thisTrialStamps = spkResp.raster.trial(m).ts;
            numBaseSpikes   = numel(find(thisTrialStamps<0));
            baseSpkRate     = numBaseSpikes./2000;

            % count spikes in the latency windows
            PopData.session(i).unit(j).phasicCS.trialwise.latResp(1,m) = numel(find(thisTrialStamps> 3 & thisTrialStamps< 24)) - (baseSpkRate.*21);
            PopData.session(i).unit(j).phasicCS.trialwise.latResp(2,m) = numel(find(thisTrialStamps>24 & thisTrialStamps< 60)) - (baseSpkRate.*36);
            PopData.session(i).unit(j).phasicCS.trialwise.latResp(3,m) = numel(find(thisTrialStamps>60 & thisTrialStamps<140)) - (baseSpkRate.*80);
            
            % count spikes in the boxcar-derived window
            PopData.session(i).unit(j).phasicCS.trialwise.peakResp(1,m)= numel(find(thisTrialStamps > thisLeft & thisTrialStamps < thisRight)) - (baseSpkRate.*(thisRight-thisLeft));
            
            % sum the smoothed spike data over the window
            PopData.session(i).unit(j).phasicCS.trialwise.intResp(1,m) = (sum(spkResp.image.aligned(m,currWindowZ),2)./(thisRight-thisLeft)) - (baseSpkRate);
            
        end
        
 
        k = k+1;
    
    end
    
    disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' | trials: ' num2str(numTrials) ' | window: ' num2str(thisLeft) '_' num2str(thisRight)]);
    
end

k = k-1;
disp(['Total units: ' num2str(k)]);

%% CALCULATE TRIAL-BY-TRIAL ANTICIPATORY LICKS
NumSessions = size(PopData.session,2);
k=1; clear sumBehavior;
validIds = PopData.validBehavSess;

% licking is stored in PopData.session(i).events.EL.ts

for m=1:numel(validIds)

    i = validIds(m);
    
    if numel(PopData.session(i).events.EL.ts)>100
        % instantaneous ILI
        [isi] = TNC_QuantISI(PopData.session(i).events.EL.ts);
        sumBehavior.linTimes        = isi.hist.linTimes;
        sumBehavior.linCount(k,:)   = isi.hist.linCount;
        figure(11); plot(sumBehavior.linCount(k,:),'k'); drawnow;
        k = k+1;
    end

    disp(['Session number: ' num2str(i)]);

    PopData.session(i).behavior.numTrials = size(PopData.session(i).behavior.raster.trial,2);
    
    for j = 1:size(PopData.session(i).behavior.raster.trial,2)

        tmp                 = find(PopData.session(i).behavior.raster.trial(j).ts<2000);
        
        if numel(tmp) > 0
            
            licksToTest     = PopData.session(i).behavior.raster.trial(j).ts(tmp);

            allPostCSLicks  = find(PopData.session(i).behavior.raster.trial(j).ts>0);              
            allAntiCSLicks  = find(licksToTest>0);              
            allPreCSLicks   = find(PopData.session(i).behavior.raster.trial(j).ts<0);

            if numel(allAntiCSLicks) > numel(allPreCSLicks)
                PopData.session(i).behavior.anticipLicks(j) = numel(allAntiCSLicks) - numel(allPreCSLicks);
            else
                PopData.session(i).behavior.anticipLicks(j) = 0;
            end
            
            % Find the first lick after the cs
            if numel(allPostCSLicks) > 0
                PopData.session(i).behavior.flCSon(j)       = PopData.session(i).behavior.raster.trial(j).ts(allPostCSLicks(1,1));
            else
                PopData.session(i).behavior.flCSon(j)       = 0;
            end
        
            numPL           = numel(allPostCSLicks);

            if numel(allPostCSLicks) > 2

                % For all postLicks try looking for whether they were at a criterion rate (62.63 - 188.65 ms = +/-2sd)
                tmpPostLick = PopData.session(i).behavior.raster.trial(j).ts(allPostCSLicks);
                postLickIsi = tmpPostLick(2:numPL) - tmpPostLick(1:numPL-1);
                critLick    = find(postLickIsi<188.65,1);

                if numel(critLick) == 1
                    PopData.session(i).behavior.fCritLick(j) = PopData.session(i).behavior.raster.trial(j).ts(allPostCSLicks(critLick));
                else
                    PopData.session(i).behavior.fCritLick(j)= 0;                    
                end

            else

                PopData.session(i).behavior.flCSon(j)   = 0;
                PopData.session(i).behavior.fCritLick(j)= 0;

            end
            
        else

                PopData.session(i).behavior.anticipLicks(j) = 0;
                PopData.session(i).behavior.flCSon(j)       = 0;
                PopData.session(i).behavior.fCritLick(j)    = 0;

        end
        
        clear licksToTest allPostCSLicks allAntiCSLicks allPreCSLicks tmpPostLick postLickIsi;

    end
    
                
end

%% COLLAPSE LEARNING - EXTINCTION PHASIC RESPONSE DATA INTO FLAT STRUCTURE 
% Collect all of the individual trialwise data into a large matrix that can be used to match units and do statistical
% testing.
% Also collect matrix for all of the stable session data for calculating contrast measures of significant, non-significant,
% and dopamine neurons.

NumSessions = size(PopData.session,2);
k=1;m=1;numValid = 0;

sessWise.acq   = zeros(663,5);
sessWise.ext   = zeros(663,5);
sessWise.da    = zeros(663,1);
sessWise.behav = zeros(663,2);

trialWise.acq.trials        = zeros(663,1);
trialWise.acq.peak          = zeros(663,200);
trialWise.acq.int           = zeros(663,200);
trialWise.acq.lat(1).resp   = zeros(663,200);
trialWise.acq.lat(2).resp   = zeros(663,200);
trialWise.acq.lat(3).resp   = zeros(663,200);
trialWise.acq.antLicks      = zeros(663,200);

trialWise.ext.trials        = zeros(663,1);
trialWise.ext.peak          = zeros(663,200);
trialWise.ext.int           = zeros(663,200);
trialWise.ext.lat(1).resp   = zeros(663,200);
trialWise.ext.lat(2).resp   = zeros(663,200);
trialWise.ext.lat(3).resp   = zeros(663,200);
trialWise.ext.antLicks      = zeros(663,200);

for i=1:NumSessions
    
    if strcmp(PopData.session(i).sessClass,'learning')
            
        numUnits = size(PopData.session(i).unit,2);

        % valid behavior?
        if isfield(PopData.session(i).behavior,'anticipLicks')
            disp(['Session ' num2str(i) ' has valid behavior.']);
            validB = 1;
        else
            disp(['Session ' num2str(i) ' does not have valid behavior.']);
            validB = 0;
        end

        for j = 1:numUnits

            % store info about behavior
            if validB
                sessWise.behav(k,1) = 1;
            end
            
            % sessionwise data
            sessWise.acq(k,1) = PopData.session(i).unit(j).phasicCS.stable.ampb.total;
            sessWise.acq(k,2) = PopData.session(i).unit(j).phasicCS.stable.ampb.peak;
            sessWise.acq(k,3) = PopData.session(i).unit(j).phasicCS.stable.ampz.int;
            sessWise.acq(k,4) = PopData.session(i).unit(j).phasicCS.stable.ampz.mean;
            sessWise.acq(k,5) = PopData.session(i).unit(j).phasicCS.stable.ampz.peak;

            % is this a dopamine neuron?
            sessWise.da(k,1) = PopData.session(i).unit(j).daLogic; 

            % trialwise data
            numTrials                                   = size(PopData.session(i).unit(j).phasicCS.trialwise.peakResp,2);
            trialWise.acq.trials(k,1)                   = numTrials;
            trialWise.acq.peak(k,1:numTrials)           = PopData.session(i).unit(j).phasicCS.trialwise.peakResp(1,1:numTrials);
            trialWise.acq.int(k,1:numTrials)            = PopData.session(i).unit(j).phasicCS.trialwise.intResp(1,1:numTrials);
            trialWise.acq.lat(1).resp(k,1:numTrials)    = PopData.session(i).unit(j).phasicCS.trialwise.latResp(1,1:numTrials);
            trialWise.acq.lat(2).resp(k,1:numTrials)    = PopData.session(i).unit(j).phasicCS.trialwise.latResp(2,1:numTrials);
            trialWise.acq.lat(3).resp(k,1:numTrials)    = PopData.session(i).unit(j).phasicCS.trialwise.latResp(3,1:numTrials);

            % trialwise behavior
            if validB
                trialWise.acq.antLicks(k,1:numTrials)   = PopData.session(i).behavior.anticipLicks(1,1:numTrials);
            else
                trialWise.acq.antLicks(k,1:numTrials)   = -1.*ones(1,numTrials);
            end
            
            % increment the unit counter
            k = k+1;
        end
        
    elseif strcmp(PopData.session(i).sessClass,'extinction')
            
        numUnits = size(PopData.session(i).unit,2);

        % valid behavior?
        if isfield(PopData.session(i).behavior,'anticipLicks')
            disp(['Session ' num2str(i) ' has valid behavior.']);
            validB = 1;
            numValid = numValid +1;
        else
            disp(['Session ' num2str(i) ' does not have valid behavior.']);
            validB = 0;
        end

        for l = 1:numUnits
            
            % store info about behavior
            if validB
                sessWise.behav(m,2) = 1;
            end
            
            % sessionwise data
            sessWise.ext(m,1) = PopData.session(i).unit(l).phasicCS.stable.ampb.total;
            sessWise.ext(m,2) = PopData.session(i).unit(l).phasicCS.stable.ampb.peak;
            sessWise.ext(m,3) = PopData.session(i).unit(l).phasicCS.stable.ampz.int;
            sessWise.ext(m,4) = PopData.session(i).unit(l).phasicCS.stable.ampz.mean;
            sessWise.ext(m,5) = PopData.session(i).unit(l).phasicCS.stable.ampz.peak;
            
%             % is this a dopamine neuron?
%             sessWise.ids(m,1) = PopData.session(i).unit(l).daLogic; 

            % trialwise data
            numTrials                                   = size(PopData.session(i).unit(l).phasicCS.trialwise.peakResp,2);
            trialWise.ext.trials(m,1)                   = numTrials;
            trialWise.ext.peak(m,1:numTrials)           = PopData.session(i).unit(l).phasicCS.trialwise.peakResp(1,1:numTrials);
            trialWise.ext.int(m,1:numTrials)            = PopData.session(i).unit(l).phasicCS.trialwise.intResp(1,1:numTrials);
            trialWise.ext.lat(1).resp(m,1:numTrials)    = PopData.session(i).unit(l).phasicCS.trialwise.latResp(1,1:numTrials);
            trialWise.ext.lat(2).resp(m,1:numTrials)    = PopData.session(i).unit(l).phasicCS.trialwise.latResp(2,1:numTrials);
            trialWise.ext.lat(3).resp(m,1:numTrials)    = PopData.session(i).unit(l).phasicCS.trialwise.latResp(3,1:numTrials);

            % trialwise behavior
            if validB
                trialWise.ext.antLicks(m,1:numTrials)   = PopData.session(i).behavior.anticipLicks(1,1:numTrials);
            else
                trialWise.ext.antLicks(m,1:numTrials)   = -1.*ones(1,numTrials);
            end

            % increment the unit counter
            m = m+1;

        end
        
    else
        disp('Serious error: unrecognized session class');
    end
    
end

k = k-1; m = m-1;
disp('  ');
disp(['Total units: ' num2str(k) ' in acquisition | ' num2str(m) ' in extinction.']);
disp(['Valid extinction sessions: ' num2str(numValid)]);
disp('  ');
    
%% CALCULATE SIGNIFICANCE OF RESPONSE CONTRASTS 

sessWise.stats = zeros(663,4); % four columns corresponding to: kw test, rs test, significant logic, sign of change

for i = 1:663
    
    % grab trial data
    responseA = trialWise.acq.peak(i,15:40);
    responseE = trialWise.ext.peak(i,20:45);
    
    % check significance
    sessWise.stats(i,1) = kruskalwallis([responseA responseE], [zeros(1,numel(responseA)) ones(1,numel(responseE))], 'off');
    sessWise.stats(i,2) = ranksum(responseA,responseE);

    if sessWise.stats(i,1)<0.05 && sessWise.stats(i,2)<0.05

        sessWise.stats(i,3) = 1;
        sessWise.stats(i,4) = mean(responseE) - mean(responseA);            
        
    else
        % everything should stay equal to 0
    end
    
    clear responseA responseE

end

tmpSig = find(sessWise.stats(:,3)==1);
tmpPos = find(sessWise.stats(:,4)>0);

disp(['Number of significant responses: ' num2str(numel(tmpSig)) ' | Number of increasers: ' num2str(numel(tmpPos))]);

%% COLLECT SUMMARY RESPONSE DATA FOR SIGNIFICANT, NONSIGNIFICANT, AND DA UNITS

sessWise.sigInds    = [];
sessWise.nSigInds   = [];
sessWise.daInds     = [];

for i = 1:663

    if sessWise.da(i) == 1
        sessWise.daInds = [sessWise.daInds,i];
    end
    
    if sessWise.stats(i,3) == 1
        if sessWise.da(i) == 0
            sessWise.sigInds = [sessWise.sigInds,i];            
        end
    else
        if sessWise.da(i) == 0
            sessWise.nSigInds = [sessWise.nSigInds,i];
        end        
    end

end

%% COLLECT INDICES FOR DA CELL & SIGNIF EXTINCTION CELL IN VALID BEHAVIORAL SESSION

sessWise.sigIncrValidBehav      = [];
sessWise.daValidBehav           = [];
sessWise.sigIncrInvalidBehav    = [];
sessWise.daInvalidBehav         = [];

for i = 1:663

    if sessWise.da(i) == 1
        if sessWise.behav(i,2) == 1
            sessWise.daValidBehav = [sessWise.daValidBehav,i];
        else
            sessWise.daInvalidBehav = [sessWise.daInvalidBehav,i];            
        end
    end
    
    if sessWise.stats(i,3) == 1
        if sessWise.stats(i,4) > 0
            if sessWise.behav(i,2) == 1
                sessWise.sigIncrValidBehav = [sessWise.sigIncrValidBehav,i];
            else
                sessWise.sigIncrInvalidBehav = [sessWise.sigIncrInvalidBehav,i];
            end
        end
    end

end

disp('Completed');

%% COLLECT TRIAL-BY-TRIAL RESPONSE FOR ALL SIGNIFICANT NONDA CELLS IN NOT VALID EXTINCTION SESSIONS

%% CALCULATE TRIAL-BY-TRIAL CELL RESPONSES FOR DA CELLS

%% COLLECT TRIAL-BY-TRIAL BEHAVIOR DATA FOR VALID LEARNING-EXTINCTION SESSIONS

NumSessions = size(PopData.session,2);
k=1; m=1;
clear extLicks acqLicks

for i=1:NumSessions
    
    if strcmp(PopData.session(i).sessClass,'learning')

        % valid behavior?
        if isfield(PopData.session(i).behavior,'anticipLicks')
            disp(['Session ' num2str(i) ' has valid behavior.']);

            numTrials = numel(PopData.session(i).behavior.anticipLicks);
            acqLicks(k,1:numTrials) = PopData.session(i).behavior.anticipLicks;

            k = k+1;
        else
            disp(['Session ' num2str(i) ' does not have valid behavior.']);
        end

    else
        
        % valid behavior?
        if isfield(PopData.session(i).behavior,'anticipLicks')
            disp(['Session ' num2str(i) ' has valid behavior.']);

            numTrials = numel(PopData.session(i).behavior.anticipLicks);
            extLicks(m,1:numTrials) = PopData.session(i).behavior.anticipLicks;        

            m = m+1;
        else
            disp(['Session ' num2str(i) ' does not have valid behavior.']);
        end

    end
    
end

%% COLLECT ALL TRIAL-BY-TRIAL RESPONSES FOR SESSIONS INCLUDING DA NEURONS
NumSessions = size(PopData.session,2);
k=1;

sessionList.acq     = [];
sessionList.ext     = [];

disp(' ');
disp(' ');
disp(' ');

for i=1:NumSessions

    if isfield(PopData.session(i).behavior,'anticipLicks')
        if strcmp(PopData.session(i).sessClass,'learning')
            if i<85 | i>109

                numUnits = size(PopData.session(i).unit,2);

                for j = 1:numUnits

                    if PopData.session(i).unit(j).daLogic

                        disp(['Dopamine unit: ' num2str(j)]);

                        % Already in session list?
                        testSess = find(sessionList.acq==i);
                        if numel(testSess)==0
                            sessionList.acq = [sessionList.acq,i];
                            disp(['Added: ' num2str(i) ' >> acquisition session list.']);
                        end
                    end
                end
            end
        else

            numUnits = size(PopData.session(i).unit,2);

            for j = 1:numUnits
                if PopData.session(i).unit(j).daLogic

                    disp(['Dopamine unit: ' num2str(j)]);

                    % Already in session list?
                    testSess = find(sessionList.ext==i);
                    if numel(testSess)==0
                        sessionList.ext = [sessionList.ext,i];
                        disp(['Added: ' num2str(i) ' >> extinction session list.']);
                    end
                end
            end
        end
    end    
end

k = k-1;
disp(' ');
disp(' ');
sessionList.acq
sessionList.ext
disp(' ');

%% ACQ: FOR GOOD BEHAVIOR SESSIONS WITH DA NEURON RECORDED EXAMINE SIMULTANEOUS RECORDED 
k=1;
for i=1:numel(sessionList.acq)
    
    m = 1;
    n = 1;
    
    thisSess = sessionList.acq(i);
    
    if sessionList.acq(i)==53 % for 53 and 54; use 55 and 56 for simultaneous data

        daSessionData.acq.session(k).da.resp(m,:) = PopData.session(thisSess).unit(1).phasicCS.trialwise.intResp;

        thisSess = 55;

        numUnits = size(PopData.session(thisSess).unit,2);

        for j = 1:numUnits

            daSessionData.acq.session(k).nda.resp(n,:) = PopData.session(thisSess).unit(j).phasicCS.trialwise.intResp;                
            n = n+1;

        end

        k = k+1;
        
    else
                
        numUnits = size(PopData.session(sessionList.acq(i)).unit,2);

        for j = 1:numUnits
            
            if PopData.session(thisSess).unit(j).daLogic
                daSessionData.acq.session(k).da.resp(m,:) = PopData.session(thisSess).unit(j).phasicCS.trialwise.intResp;
                m = m+1;
            else
                daSessionData.acq.session(k).nda.resp(n,:) = PopData.session(thisSess).unit(j).phasicCS.trialwise.intResp;                
                n = n+1;
            end
        end
        
        k = k+1;
        
    end
    
end

%% ACQ: GET TRIAL BY TRIAL CORRELATIONS

numSessions = size(daSessionData.acq.session,2);
m = 1;
acqPlotDa = []; acqPlotNda = [];

for i=1:numSessions

    numDa       = size(daSessionData.acq.session(i).da.resp,1);
    numNonDa    = size(daSessionData.acq.session(i).nda.resp,1);

    for j=1:numDa
        
        for k=1:numNonDa
            
            [rho,pval]  = corr(daSessionData.acq.session(i).da.resp(j,:)',daSessionData.acq.session(i).nda.resp(k,:)');

            acqPlotDa   = [acqPlotDa,daSessionData.acq.session(i).da.resp(j,:)];
            acqPlotNda  = [acqPlotNda,daSessionData.acq.session(i).nda.resp(k,:)];

            daSessionData.acq.rhos(m)   = rho;
            daSessionData.acq.pvals(m)  = pval;
            
            m           = m+1;

        end
        
    end
    
end

%% EXT: FOR GOOD BEHAVIOR SESSIONS WITH DA NEURON RECORDED EXAMINE SIMULTANEOUS RECORDED 
k=1;
for i=1:numel(sessionList.ext)
    
    m = 1;
    n = 1;
    
    thisSess = sessionList.ext(i);
    
    if sessionList.ext(i)==54 % for 53 and 54; use 55 and 56 for simultaneous data

        daSessionData.ext.session(k).da.resp(m,:) = PopData.session(thisSess).unit(1).phasicCS.trialwise.intResp;

        thisSess = 56;

        numUnits = size(PopData.session(thisSess).unit,2);

        for j = 1:numUnits

            daSessionData.ext.session(k).nda.resp(n,:) = PopData.session(thisSess).unit(j).phasicCS.trialwise.intResp;                
            n = n+1;

        end

        k = k+1;
        
    else
                
        numUnits = size(PopData.session(sessionList.ext(i)).unit,2);

        for j = 1:numUnits
            
            if PopData.session(thisSess).unit(j).daLogic
                daSessionData.ext.session(k).da.resp(m,:) = PopData.session(thisSess).unit(j).phasicCS.trialwise.intResp;
                m = m+1;
            else
                daSessionData.ext.session(k).nda.resp(n,:) = PopData.session(thisSess).unit(j).phasicCS.trialwise.intResp;                
                n = n+1;
            end
        end
        
        k = k+1;
        
    end
    
end

%% EXT: GET TRIAL BY TRIAL CORRELATIONS

numSessions = size(daSessionData.ext.session,2);
m = 1;
extPlotDa = []; extPlotNda = [];

for i=1:numSessions

    numDa       = size(daSessionData.ext.session(i).da.resp,1);
    numNonDa    = size(daSessionData.ext.session(i).nda.resp,1);

    for j=1:numDa
        
        for k=1:numNonDa
            
            [rho,pval]  = corr(daSessionData.ext.session(i).da.resp(j,:)',daSessionData.ext.session(i).nda.resp(k,:)');

            extPlotDa   = [extPlotDa,daSessionData.ext.session(i).da.resp(j,:)];
            extPlotNda  = [extPlotNda,daSessionData.ext.session(i).nda.resp(k,:)];

            daSessionData.ext.rhos(m)   = rho;
            daSessionData.ext.pvals(m)  = pval;
            
            m           = m+1;

        end
        
    end
    
end

%% CHECK CROSS CORRELATIONS FOR SIMULTANEOUS UNITS
k=1;
for i=1:numel(sessionList.ext)
    
    
    m = 1;
    n = 1;
    
    thisSess = sessionList.ext(i)
    
    % find the timestamps for first and last trials and only calculate cross correlation within that fixed range
    trialCnt = size(PopData.session(thisSess).events.CS.ts,1);
    trialOne = PopData.session(thisSess).events.CS.ts(30);
    trialLst = PopData.session(thisSess).events.CS.ts(trialCnt)+1000
    
    
    if sessionList.ext(i)==54 % for 53 and 54; use 55 and 56 for simultaneous data
% 
%         daSessionData.ext.session(k).da.resp(m,:) = PopData.session(thisSess).unit(1).phasicCS.trialwise.intResp;
%             j                   = 1;
%             delta               = zeros(1,ceil(trialLst-trialOne));
%             validStampInds      = find(PopData.session(thisSess).unit(j).ts>trialOne & PopData.session(thisSess).unit(j).ts<trialLst);
%             validStamps         = PopData.session(thisSess).unit(j).ts(validStampInds);
%             delta(validStamps-trialOne)  = 1;
% 
%             daSessionData.ext.session(k).da(m).delta = delta;
%             m = m+1;
% 
%         thisSess = 56;
% 
%         numUnits = size(PopData.session(thisSess).unit,2);
% 
%         for j = 1:numUnits
%             
%             delta               = zeros(1,ceil(trialLst-trialOne));
%             validStampInds      = find(PopData.session(thisSess).unit(j).ts>trialOne & PopData.session(thisSess).unit(j).ts<trialLst);
%             validStamps         = round(PopData.session(thisSess).unit(j).ts(validStampInds)-trialOne);
%             delta(validStamps)  = 1;
% 
%             daSessionData.ext.session(k).nda(n).delta = delta;                
%             n = n+1;
% 
%         end
% 
%         k = k+1;
        
    else
                
        numUnits = size(PopData.session(sessionList.ext(i)).unit,2);

        for j = 1:numUnits
            
            delta               = zeros(1,ceil(trialLst-trialOne));
            validStampInds      = find(PopData.session(thisSess).unit(j).ts>trialOne & PopData.session(thisSess).unit(j).ts<trialLst);
            validStamps         = round(PopData.session(thisSess).unit(j).ts(validStampInds)-trialOne);
            delta(validStamps)  = 1;
            
            if PopData.session(thisSess).unit(j).daLogic==1
                disp('here.')
                daSessionData.ext.session(k).da(m).delta    = delta;
                daSessionData.ext.session(k).da(m).smth     = conv(delta,currParams.filter.kernel,'same');
                m = m+1;
            else
                disp('there.')
                daSessionData.ext.session(k).nda(n).delta   = delta;                
                daSessionData.ext.session(k).nda(n).smth    = conv(delta,currParams.filter.kernel,'same');                
                n = n+1;
            end
        end
        
        k = k+1;
        
    end
    
end

k=k-1;

%% FOR EACH SIMULTANEOUS SESSION CHECK THE XCORRS

count=1;
lags = 50;

clear crosscorr crosscorrS

for k = 1:size(daSessionData.ext.session,2)
% for k = 1:10

    k
    
    numDa = size(daSessionData.ext.session(k).da,2);
    
    for m=1:numDa
        
        
        numNDa = size(daSessionData.ext.session(k).nda,2);
        
        for n=1:numNDa
            
            crosscorr(count,:) = xcorr(daSessionData.ext.session(k).da(m).delta,daSessionData.ext.session(k).nda(n).delta,lags);
            crosscorrS(count,:) = xcorr(daSessionData.ext.session(k).da(m).smth,daSessionData.ext.session(k).nda(n).smth,lags);
            count = count+1;
        end
    end
   
end


figure(1)
for p=1:81
    subplot(9,9,p);
    bar(-lags:lags,crosscorrS(p,:))
end

%% GET POPULATION MEANS FOR SIMULTANEOUS SESSIONS

k=1;
m = 1;
n = 1;
    
for i=1:numel(sessionList.acq)
    

    
    thisSess = sessionList.acq(i);
    
    if sessionList.acq(i)==53 % for 53 and 54; use 55 and 56 for simultaneous data

        allDaSessData.acq.da.psth(m,:) = PopData.session(thisSess).unit(1).respCSS.psthZ;

        thisSess = 55;

        numUnits = size(PopData.session(thisSess).unit,2);

        for j = 1:numUnits

            allDaSessData.acq.nda.psth(n,:) = PopData.session(thisSess).unit(j).respCSS.psthZ;
            n = n+1;

        end

        k = k+1;
        
    else
                
        numUnits = size(PopData.session(sessionList.acq(i)).unit,2);

        for j = 1:numUnits
            
            if PopData.session(thisSess).unit(j).daLogic
                allDaSessData.acq.da.psth(m,:) = PopData.session(thisSess).unit(j).respCSS.psthZ;
                m = m+1;
            else
                allDaSessData.acq.nda.psth(n,:) = PopData.session(thisSess).unit(j).respCSS.psthZ;
                n = n+1;
            end
        end
        
        k = k+1;
        
    end
    
end

k=1;
m = 1;
n = 1;
    
for i=1:numel(sessionList.ext)
    

    
    thisSess = sessionList.ext(i);
    
    if sessionList.acq(i)==54 % for 53 and 54; use 55 and 56 for simultaneous data

        allDaSessData.ext.da.psth(m,:) = PopData.session(thisSess).unit(1).respCSS.psthZ;

        thisSess = 56;

        numUnits = size(PopData.session(thisSess).unit,2);

        for j = 1:numUnits

            allDaSessData.ext.nda.psth(n,:) = PopData.session(thisSess).unit(j).respCSS.psthZ;
            n = n+1;

        end

        k = k+1;
        
    else
                
        numUnits = size(PopData.session(sessionList.acq(i)).unit,2);

        for j = 1:numUnits
            
            if PopData.session(thisSess).unit(j).daLogic
                allDaSessData.ext.da.psth(m,:) = PopData.session(thisSess).unit(j).respCSS.psthZ;
                m = m+1;
            else
                allDaSessData.ext.nda.psth(n,:) = PopData.session(thisSess).unit(j).respCSS.psthZ;
                n = n+1;
            end
        end
        
        k = k+1;
        
    end
    
end

exportSimulMean = [mean(allDaSessData.acq.da.psth,1)',mean(allDaSessData.ext.da.psth,1)',mean(allDaSessData.acq.nda.psth,1)',mean(allDaSessData.ext.nda.psth,1)'];
exportSimulERR = [std(allDaSessData.acq.da.psth,0,1)'./sqrt(17),std(allDaSessData.ext.da.psth,0,1)'./sqrt(17),std(allDaSessData.acq.nda.psth,0,1)'./sqrt(64),std(allDaSessData.ext.nda.psth,0,1)'./sqrt(64)];

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

%% ANALYZE TRIAL-BY-TRIAL RESPONSES FOR SESSIONS WITH MEASUREABLE BEHAVIORAL EXTINCTION

%% COLLECT TRIAL-BY-TRIAL RESPONSES FOR SESSIONS WITH 'SPONTANEOUS' EXTINCTION

%% CORE ITERATOR CODE...
NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
    
    if strcmp(PopData.session(i).sessClass,'learning')
            
        numUnits = size(PopData.session(i).unit,2);

        for j = 1:numUnits
            
            PopData.session(i).unit(j).
            
            k = k+1;
        end
        
    end
    
end

k = k-1;
disp(['Total units: ' num2str(k)]);

%% SOME ADDITIONAL VISUALIZATION CODE

%-------------------------------------------------------------------------------------------
% ADD IN DISPLAY FOR VALIDATION OF THE ALL THE MEASUREMENTS    
%-------------------------------------------------------------------------------------------

left    = (respProps.fw(currTotIndex).X(1,1).*10) - 9;
right   = (respProps.fw(currTotIndex).X(1,2).*10) - 9;

figure(1); clf; 

% plot boxcar and indicate the window used
subplot(411);
plot(1:500,allCSaligned.boxcar(:,currTotIndex),'r.-',1:500,respCSS.image.boxcar.*100,'k.-',[left./10 left./10],[0 400],'g',[right./10 right./10],[0 400],'b','LineWidth',2);
title(currTotIndex);
axis([200 225 0 150]);

% plot the raster plot over the window from boxcar
subplot(412);
plot(1:5001,respCSS.image.psthAVG.*1000,'k',[left left],[0 500],'g',[right right],[0 500],'b','LineWidth',2);
axis([2000 2250 0 150]);

% show the smoothed trial-by-trial estimates and indicate the window used
subplot(413);
plot(1:5001,respCSS.image.aligned(10:25,1:5001).*1000,'--'); hold on;
plot([left left],[0 500],'g',[right right],[0 500],'b','LineWidth',3);
axis([2000 2250 0 500]);

% plot all three of the output measures and check how well correlated they are
subplot(414);

plot(1:numTrials,respVectorA,'k.-',1:numTrials,spkRespVectorA,'ro-'); axis([0 numTrials -3 10]);
drawnow;

%-------------------------------------------------------------------------------------------

%% LOOK AT DA NEURON CS-EVOKED PAUSES IN EXTINCTION
clear allDAEXisis;
NumSessions = size(PopData.session,2);
k=1;
allDAEXisis = [];
allDAACisis = [];
latToFirstA = [];
latToFirstE = [];
countA=0;
countE=0;

for i=1:NumSessions

    if strcmp(PopData.session(i).sessClass,'extinction')
                
        numUnits = size(PopData.session(i).unit,2);

        for j = 1:numUnits

            if PopData.session(i).unit(j).daLogic

                numTrials = size(PopData.session(i).unit(j).respCS.raster.trial,2);
                if numTrials > 35
                    numTrials=35;
                end
                for m=5:numTrials
                    posTs   = find(PopData.session(i).unit(j).respCS.raster.trial(m).ts>0 & PopData.session(i).unit(j).respCS.raster.trial(m).ts<1000);
                    if numel(posTs)>1
                        currTs  = PopData.session(i).unit(j).respCS.raster.trial(m).ts(posTs);
                        numTs   = numel(currTs);
                        currISI = currTs(2:numTs) - currTs(1:numTs-1);
                        allDAEXisis = [allDAEXisis,currISI'];
                        if currTs(find(currISI==min(currISI),1,'first'))<200
                            latToFirstE = [latToFirstE currTs(find(currISI==min(currISI),1,'first'))];
                        end
                    end
                end
                
                countE = countE+1;
                latency(countE,2) = mean(latToFirstE);
                latToFirstE = [];

            end
        end
        
    else
        
        numUnits = size(PopData.session(i).unit,2);

        for j = 1:numUnits

            if PopData.session(i).unit(j).daLogic

                numTrials = size(PopData.session(i).unit(j).respCS.raster.trial,2);
                if numTrials > 40
                    numTrials=40;
                end                
                for m=10:numTrials
                    posTs   = find(PopData.session(i).unit(j).respCS.raster.trial(m).ts>0 & PopData.session(i).unit(j).respCS.raster.trial(m).ts<1000);
                    if numel(posTs)>1
                        currTs  = PopData.session(i).unit(j).respCS.raster.trial(m).ts(posTs);
                        numTs   = numel(currTs);
                        currISI = currTs(2:numTs) - currTs(1:numTs-1);
                        allDAACisis = [allDAACisis,currISI'];
                        if currTs(find(currISI==min(currISI),1,'first'))<200
                            latToFirstA = [latToFirstA currTs(find(currISI==min(currISI),1,'first'))];
                        end
                    end
                end
                
                countA = countA+1;
                latency(countA,1) = mean(latToFirstA);
                latToFirstA = [];
    
            end
        end
        
    end
disp(i);    
end

xVals = 10.^[0:0.05:3];
ISIcompare.acHist = hist(allDAACisis,xVals)./numel(allDAACisis);
ISIcompare.exHist = hist(allDAEXisis,xVals)./numel(allDAEXisis);
figure(1); clf; semilogx(xVals,ISIcompare.acHist,'k',xVals,ISIcompare.exHist,'r',xVals,ISIcompare.exHist-ISIcompare.acHist,'b');
figure(2); clf; plot(0:5:1000,hist(latToFirstA,0:5:1000),'k',0:5:1000,hist(latToFirstE,0:5:1000),'r');

%% DO I HAVE US TIMINGS?

NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
    
    if strcmp(PopData.session(i).sessClass,'learning')

        if numel(PopData.session(i).events.US.ts)==1
            disp(PopData.session(i).sessId);
        end
        
    end
    
end

k = k-1;
disp(['Total units: ' num2str(k)]);

%% DEPRECATED CODE
% GET PSTH RESPONSE LATENCIES

window  = 2000:2150;
window2 = 2500:4000;

responseLatency = zeros(size(allCSaligned.psthZ,2),1);
responsePeak    = zeros(size(allCSaligned.psthZ,2),1);
responseWidth   = zeros(size(allCSaligned.psthZ,2),1);
responseAmp     = zeros(size(allCSaligned.psthZ,2),1);
responseMid     = zeros(size(allCSaligned.psthZ,2),1);

for i = 1:size(allCSaligned.psthZ,2)
    
    tmpCheck = allCSaligned.psthZ(2001:2125,i);
    check = find(tmpCheck>3);

    tmpDATA = allCSaligned.psthZ(window,i);
    
    if numel(check)>8

        peakLoc     = find(tmpDATA(2:149)==max(tmpDATA(2:149)),1);
        beginPeak   = find(tmpDATA(1:peakLoc)<3,1,'last');
        if numel(beginPeak)==0
            beginPeak = peakLoc-1;
        end
        endPeak     = peakLoc+find(tmpDATA(peakLoc+1:numel(tmpDATA))<3,1,'first');
        if numel(endPeak)==0
            endPeak = peakLoc+1;
        end
            
        peakInds = beginPeak:1:endPeak;
        
        responseAmp(i)      = trapz(tmpDATA(peakInds));
        responseLatency(i)  = beginPeak;
        responsePeak(i)     = peakLoc;
        responseMid(i,1)    = mean(peakInds);
        responseWidth(i)    = endPeak - beginPeak;


%     elseif min(tmpDATA)>-2
%         
%         peakLoc     = find(tmpDATA(2:149)==min(tmpDATA(2:149)),1);
%         beginPeak   = find(tmpDATA(1:peakLoc)>-2,1,'last');
%         if numel(beginPeak)==0
%             beginPeak = peakLoc-1;
%         end
%         endPeak     = peakLoc+find(tmpDATA(peakLoc+1:numel(tmpDATA))>-2,1,'first');
%         if numel(endPeak)==0
%             endPeak = peakLoc+1;
%         end
%         beginPeak
%         endPeak
% 
%         peakInds = beginPeak:1:endPeak;
% 
%         responseAmp(i) = trapz(tmpDATA(peakInds))./numel(peakInds);
% 
%         responseLatency(i)  = beginPeak;
%         responsePeak(i)     = peakLoc;
%         responseMid(i,1)    = mean(peakInds);

    else
        
        responseLatency(i)  = 0;
        responsePeak(i)     = 0;
        responseAmp(i)      = 0;
        responseMid(i,1)    = 0;
        responseWidth(i)    = 0;

    end
    
    
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
imagesc(allCSaligned.sort.amp.sortedZ(window,:)',[-10 10]);
xlabel('Time (ms)');
ylabel('Cell index');
title('CS Aligned PSTH');

% subplot(1,5,3);
% plot(allCSaligned.sort.km.ids(allCSaligned.sort.amp.inds3),701:-1:1,'k.'); axis([0 numClusts+1 1 701]); axis off;

subplot(1,5,4:5);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
% imagesc(corr(allCSaligned.sort.amp.sorted),[-1,1]);
imagesc(allCSaligned.sort.amp.sortedZ2(window,:)',[-10 10]);
ylabel('Cell index');
xlabel('Time (ms)');
title('CS Aligned PSTH');

% Random Code

validIndices = find(validList);
disp(['Number of valid responses: ' num2str(length(validIndices))]);

daList = PopData.daList;

% sort the data by the size of the response and make sure to keep track of the correct indices
responseMagnitude   = respProps.ampB.ext.integral(validIndices) - respProps.ampB.acq.integral(validIndices);
responseMagnitudeZ  = respProps.ampZ.ext.integral(validIndices) - respProps.ampZ.acq.integral(validIndices);
clear sortedResponses;
[values,sortInds]   = sort(responseMagnitude,2,'ascend');
sortedResponses(1,:)= responseMagnitude(sortInds);
sortedResponses(2,:)= validIndices(sortInds);

figure(4);
plot([0 0],[-10 70],'k-',[-100 700],[0 0],'k-',[-100 700],[-100 700],'r',[-50 350],[-50 350].*2,'r--',[-50 350].*2,[-50 350],'r--',respProps.ampZ.acq.integral(sigMod.nonSigInd),respProps.ampZ.ext.integral(sigMod.nonSigInd),'k+',respProps.ampZ.acq.integral(sigMod.sigInd),respProps.ampZ.ext.integral(sigMod.sigInd),'r.',respProps.ampZ.acq.integral(daList),respProps.ampZ.ext.integral(daList),'b.')

figure(5);
% plot([0 150],[0 0],'r',respProps.ampZ.learn.center(validIndices),respProps.ampZ.ext.integral(validIndices)-respProps.ampZ.learn.integral(validIndices),'k.');
plot([0 150],[0 0],'k--',(respProps.winProps.center(sigMod.nonSigInd).*10)-50,respProps.ampZ.ext.integral(sigMod.nonSigInd) - respProps.ampZ.acq.integral(sigMod.nonSigInd),'k+',(respProps.winProps.center(daList).*10)-50,respProps.ampZ.ext.integral(daList) - respProps.ampZ.acq.integral(daList),'b.',(respProps.winProps.center(sigMod.sigInd).*10)-50,respProps.ampZ.ext.integral(sigMod.sigInd) - respProps.ampZ.acq.integral(sigMod.sigInd),'ro');
% (respProps.winProps.center(validIndices).*10)-50,responseMagnitudeZ,'ko',
figure(10);
plot(1:length(validIndices),sortedResponses(1,:),'ko',[0 length(validIndices)],[0 0],'r-',[length(validIndices)./2 length(validIndices)./2],[-50 50],'k--');

% SORT ALIGNED PSTHS BY LATENCY AND AMPLITUDE

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
peaks       = [0,12.4,45.5,71.1];
fwhmsH      = [1,13.8,16.8,108.0];
thresholds  = peaks+fwhmsH

% figure(6); plot(responsePeak,responseAmp,'k.',peaks,amps,'r.',peaks+fwhmsH,amps,'ro',peaks-fwhmsL,amps,'ro'); drawnow;

% Classify the neurons based upon responseLatency
numCells    = size(responseMid,1);
[y,inds]    = sort(responseMid,1,'descend');

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

%__________________________________________________________________________________
%__________________________________________________________________________________

% NOW SORT BY AMPLITUDE WITHIN THE LATENCY PEAKS
%__________________________________________________________________________________
%__________________________________________________________________________________


totalGroups = max(latencyGroups);

for i=1:totalGroups

    % find the cells matching the current latency peak
    indsToGroup     = trackSort.real(i).inds;
    tmpAmpData      = responseAmp(indsToGroup);
    tmpSubData      = allCSaligned.psthZ(:,indsToGroup);
    tmpSubDataUS    = allUSaligned.psthZ(:,indsToGroup);
    tmpSubDataEX    = allCSEXaligned.psthZ(:,indsToGroup);
    [y,inds]        = sort(tmpAmpData,1,'descend');

    tmpSubDataSort      = tmpSubData(:,inds);
    tmpSubDataUSSort    = tmpSubDataUS(:,inds);
    tmpSubDataEXSort    = tmpSubDataEX(:,inds);
    
    allCSaligned.sort.latency.group(i).indAMP = inds;
    allUSaligned.sort.latency.group(i).indAMP = inds;
    
    % within those cells sort by response amplitude
    if i>1
        allCSaligned.sort.ampPeak.sortedZ = [tmpSubDataSort,allCSaligned.sort.ampPeak.sortedZ];
        allUSaligned.sort.ampPeak.sortedZ = [tmpSubDataUSSort,allUSaligned.sort.ampPeak.sortedZ];
        allCSEXaligned.sort.ampPeak.sortedZ = [tmpSubDataEXSort,allCSEXaligned.sort.ampPeak.sortedZ];
    else
        allCSaligned.sort.ampPeak.sortedZ = tmpSubDataSort;
        allUSaligned.sort.ampPeak.sortedZ = tmpSubDataUSSort;
        allCSEXaligned.sort.ampPeak.sortedZ = tmpSubDataEXSort;
    end

    figure(3); subplot(1,7,1:3); colormap(mapName);
    imagesc(allCSaligned.sort.ampPeak.sortedZ(1900:2200,:)',[-15 15]); drawnow; pause(1);

    figure(3); subplot(1,7,5:7); colormap(mapName);
    imagesc(allUSaligned.sort.ampPeak.sortedZ(1900:2200,:)',[-15 15]); drawnow; pause(1);

    figure(7); subplot(1,7,1:3); colormap(mapName);
    imagesc(allCSEXaligned.sort.ampPeak.sortedZ(1900:2200,:)',[-15 15]); drawnow; pause(1);

%     figure(8); 
%     imagesc(allCSEXaligned.sort.ampPeak.sortedZ(1500:4500,:)'-allCSaligned.sort.ampPeak.sortedZ(1500:2500,:)',[-15 15]); drawnow; pause(1);
%     colormap(mapName);
end

% FIND EXEMPLARS OF LATENCY GROUPS

% get indices of 5 largest responses for each group
numExamples = 15;
allCSaligned.exemplars = zeros(numExamples,3);

for m = 2:4
    tmpGroupData                        = responseAmp(trackSort.real(m).inds);
    [sortedGroupData,sortedInds]        = sort(tmpGroupData,1,'descend');

    allCSaligned.exemplars(:,m-1)         = trackSort.real(m).inds(sortedInds(1:numExamples));
end

%% COLLECT ALL AVERAGE FIRING RATES
NumSessions = size(PopData.session,2)
k=1;m=1;

for i=1:NumSessions
    % check if learning
    if i<10
        disp(['Session 00' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
    elseif i<100
        disp(['Session 0' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
    else
        disp(['Session ' num2str(i) ' | ' num2str(k) ' - ' PopData.session(i).sessId]);
    end

    if strcmp(PopData.session(i).sessClass,'learning')

        numUnits = size(PopData.session(i).unit,2);

        for j = 1:numUnits
            
            meanISI(k) = PopData.session(i).unit(j).isi.hist.logTimes(find(PopData.session(i).unit(j).isi.hist.logCount==max(PopData.session(i).unit(j).isi.hist.logCount),1,'first'));
            k = k + 1;
            
        end
        
    else
        
        numUnits = size(PopData.session(i).unit,2);

        for j = 1:numUnits
            
            meanISIe(m) = PopData.session(i).unit(j).isi.hist.logTimes(find(PopData.session(i).unit(j).isi.hist.logCount==max(PopData.session(i).unit(j).isi.hist.logCount),1,'first'));
            m = m + 1;
            
        end
    end
    
end
            
