% SCRIPT OVERVIEW: Analysis of photostimulation during extracellular recording
% _________________________________________________________________________
% NOTES:
%
% _________________________________________________________________________
% PART OF THE TONIC SOFTWARE PACKAGE
%   developed by JOSHUA DUDMAN
%   begun in the Siegelbaum/Kandel Labs at Columbia University
%   primary development at HHMI / JFRC
%   requires MATLAB 2011a; Tested on MacOS (10.6.8)
% _________________________________________________________________________
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: http://dudmanlab.org/html/projects.html

%% STANDARD ANALYSIS PARAMETERS
% Define the smoothing to apply to the timeseries data
currParams.smthParams.rise     = 0.01;
currParams.smthParams.decay    = 5;
% currParams.smthParams.decay    = 10;
% currParams.smthParams.decay    = 25;
% currParams.smthParams.decay    = 50;
currParams.filter.causal = 0;

if currParams.filter.causal
    [currParams.filter.kernel]  = TNC_CreateCausalKernel(currParams.smthParams.rise,currParams.smthParams.decay,1);
else
    [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
end

currParams.winParams.prior     = 0.2e3;
currParams.winParams.after     = 0.5e3; 

PopData.currParams = currParams;

disp(' ');
disp(' ');
disp('________________________________');
disp(' ');
disp('Initialized analysis parameters.');
disp('________________________________');
disp(' ');
disp(' ');

%% LIST ALL UNITS AND SESSION NAMES
numSessions = size(PopData.session,2);
count = 0;

for j=1:numSessions
    numUnits = size(PopData.session(j).unit,2);

    for i=1:numUnits
        count = count+1;
        disp([num2str(count) ' >>> ' PopData.session(j).sessId ' >>> '  PopData.session(j).unit(i).name ]);

    end
end

% Dopamine neuron list:
% 88 >>> G03_110928_007nnDAch53bcde>>6 >>> elec53b
% 89 >>> G03_110928_007nnDAch53bcde>>6 >>> elec53c
% 90 >>> G03_110928_007nnDAch53bcde>>6 >>> elec53d
% 91 >>> G03_110928_007nnDAch53bcde>>6 >>> elec53e
DAlist = 88:91;

%% EXTRACT THE ISI FROM ALL UNITS
NumSessions = size(PopData.session,2)
k=1;m=1;countUnits=0;
eval('home');
for i=1:NumSessions
    % check if learning
    if i<10
        disp(['Session 00' num2str(i) ' | ' num2str(k) ' >>> ' PopData.session(i).sessId]);
    elseif i<100
        disp(['Session 0' num2str(i) ' | ' num2str(k) ' >>> ' PopData.session(i).sessId]);
    else
        disp(['Session ' num2str(i) ' | ' num2str(k) ' >>> ' PopData.session(i).sessId]);
    end

    numUnits = size(PopData.session(i).unit,2);
    countUnits = countUnits+numUnits;

    for j = 1:numUnits
        if isnan(PopData.session(i).unit(j).isi.hist.logCount(1,1))
            disp('Invalid ISI distribution found');
        else
            allISI.stim.linear.times(1,:)    = PopData.session(i).unit(j).isi.hist.linTimes;
            allISI.stim.linear.counts(k,:)   = PopData.session(i).unit(j).isi.hist.linCount;

            allISI.stim.log.times(1,:)    = PopData.session(i).unit(j).isi.hist.logTimes;
            allISI.stim.log.counts(k,:)   = PopData.session(i).unit(j).isi.hist.logCount;

            k=k+1;
        end
    end

    
end

k = k-1;

disp(['COMPLETED ISI EXTRACTION >>> Total units: ' num2str(k) ' | ' num2str(countUnits)]);

%% EXTRACT STIMULUS PETHS FROM ALL UNITS FOR LITE AND TONE
eval('home');
disp(['___________________________________________________']);
disp(['STARTED aligning all raster plots and updating the PopData structure...']);

PopData.currParams = currParams;

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);
k=1;

for i=1:NumSessions
    
    disp('>>>>>>>>>>');
    disp(['Begin session: ' num2str(i) ' of ' num2str(NumSessions) ' | ' PopData.session(i).sessId])

    numUnits = size(PopData.session(i).unit,2);
    countUnits = countUnits+numUnits;


    PopData.session(i).trials = size(PopData.session(i).events.LITE.ts,1);
    
    for j = 1:numUnits

        numStamps   = length(PopData.session(i).unit(j).ts);
        delta       = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
        delta(1,round(PopData.session(i).unit(j).ts)+1) = 1;
        tmpSmooth   = conv(delta,currParams.filter.kernel,'same');

        % for reference
        % [response] = TNC_AlignRasters(delta,spkStamps,stableTrials,alignStamps,window,rasterFlag,boxcar)
        [rLITE] = TNC_AlignRasters(delta , PopData.session(i).unit(j).ts , -1 , PopData.session(i).events.LITE.ts , [PopData.currParams.winParams.prior,PopData.currParams.winParams.after],1,1);
        PopData.session(i).unit(j).rLITE.raster      = rLITE.raster;
        PopData.session(i).unit(j).rLITE.psthAVG     = rLITE.image.psthAVG;
        PopData.session(i).unit(j).rLITE.boxcar      = rLITE.image.boxcar;

        [rsLITE] = TNC_AlignRasters(tmpSmooth , PopData.session(i).unit(j).ts , -1 , PopData.session(i).events.LITE.ts , [PopData.currParams.winParams.prior ,PopData.currParams.winParams.after],0,1);
        PopData.session(i).unit(j).rsLITE.psthAVG    = rsLITE.image.psthAVG;
        PopData.session(i).unit(j).rsLITE.psthZ      = rsLITE.image.psthZ;
        PopData.session(i).unit(j).rsLITE.psthZe     = rsLITE.image.psthZe;
                
        [rTONE] = TNC_AlignRasters(delta , PopData.session(i).unit(j).ts , -1 , PopData.session(i).events.TONE.ts , [PopData.currParams.winParams.prior , PopData.currParams.winParams.after],1,1);
        PopData.session(i).unit(j).rTONE.raster      = rTONE.raster;
        PopData.session(i).unit(j).rTONE.psthAVG     = rTONE.image.psthAVG;
        PopData.session(i).unit(j).rTONE.boxcar      = rTONE.image.boxcar;

        [rsTONE] = TNC_AlignRasters(tmpSmooth , PopData.session(i).unit(j).ts , -1 , PopData.session(i).events.TONE.ts , [PopData.currParams.winParams.prior,PopData.currParams.winParams.after],0,1);
        PopData.session(i).unit(j).rsTONE.psthAVG    = rsTONE.image.psthAVG;
        PopData.session(i).unit(j).rsTONE.psthZ      = rsTONE.image.psthZ;
        PopData.session(i).unit(j).rsTONE.psthZe     = rsTONE.image.psthZe;

        k=k+1;

        disp(['Completed unit: ' num2str(j) ' of ' num2str(numUnits) ' ... session: ' num2str(i) ' of ' num2str(NumSessions) ' | ' PopData.session(i).sessId])
    end

    disp(' ');
    
end
k = k-1;
disp(['___________________________________________________'])
disp(['COMPLETED ... Total units: ' num2str(k) ' | ' num2str(countUnits)]);

%% ALIGN PSTHS FOR LITE AND TONE STIMULATION
clear lightALIGN toneALIGN
countUnits = 0;
NumSessions = size(PopData.session,2)
k=1;

figure(6); clf;
figure(7); clf;

for i=1:NumSessions
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
            lightALIGN.psthA(:,k)     = PopData.session(i).unit(j).rLITE.psthAVG'.*1000;
            lightALIGN.psthZ(:,k)     = PopData.session(i).unit(j).rsLITE.psthZ';
            lightALIGN.psthZe(:,k)    = PopData.session(i).unit(j).rsLITE.psthZe';
            lightALIGN.boxcar(:,k)    = PopData.session(i).unit(j).rLITE.boxcar'.*167;
            lightALIGN.sess(1,k)      = i;
            lightALIGN.unit(1,k)      = j;
            lightALIGN.id(1,k).sessID = PopData.session(i).sessId;
            lightALIGN.id(1,k).name   = PopData.session(i).unit(j).name;

            toneALIGN.psthA(:,k)     = PopData.session(i).unit(j).rTONE.psthAVG'.*1000;
            toneALIGN.psthZ(:,k)     = PopData.session(i).unit(j).rsTONE.psthZ';
            toneALIGN.psthZe(:,k)    = PopData.session(i).unit(j).rsTONE.psthZe';
            toneALIGN.boxcar(:,k)    = PopData.session(i).unit(j).rTONE.boxcar'.*167;
            toneALIGN.sess(1,k)      = i;
            toneALIGN.unit(1,k)      = j;
            toneALIGN.id(1,k).sessID = PopData.session(i).sessId;
            toneALIGN.id(1,k).name   = PopData.session(i).unit(j).name;

            k=k+1;
        end
    end    
end

k = k-1;

disp(['Total units: ' num2str(k) ' | ' num2str(countUnits)])

numUnits = countUnits;

[mapName] = TNC_CreateRBColormap(1024,'rb');

figure(6);
% subplot(4,2,[1,3,5]);
% [mapName] = TNC_CreateRBColormap(1024,'rb');
% imagesc(lightALIGN.boxcar');
% ylabel('Cell index');
% title('LITE Aligned PSTH');
% 
% subplot(4,2,[2,4,6]);
% imagesc(toneALIGN.boxcar');
% colormap(mapName);
% ylabel('Cell index');
% title('TONE Aligned PSTH');
% 
% subplot(4,2,7);
% plot(mean(lightALIGN.boxcar,2),'k');
% colormap(mapName);
% axis([0 91 -10 40]);
% xlabel('Window Number (5 ms)');
% ylabel('Response (spk / 5ms)');
% title('LITE Aligned PSTH');
% 
% subplot(4,2,8);
% bar(mean(toneALIGN.boxcar,2),'k');
% axis([0 91 0 40]);
% xlabel('Window Number (5 ms)');
% ylabel('Response (spk / 5ms)');
% title('TONE Aligned PSTH');

% figure(8);
subplot(4,2,[1,3,5]);
imagesc(lightALIGN.psthZ',[-5 5]);
axis([0 1000 0 275]);
colormap(mapName);
ylabel('Cell index');
title('LITE Aligned PSTH');
box off;
set(gca,'TickDir','out');

subplot(4,2,[2,4,6]);
imagesc(toneALIGN.psthZ',[-3.5 3.5]);
colormap(mapName);
axis([0 1000 0 275]);
ylabel('Cell index');
title('TONE Aligned PSTH');
box off;
set(gca,'TickDir','out');

subplot(4,2,7);
bar(mean(lightALIGN.psthZ,2),'k');
axis([0 1000 -0.75 2]);
xlabel('Time (ms)');
ylabel('Z-score');
title('LITE Aligned PSTH');
box off;
set(gca,'TickDir','out');

subplot(4,2,8); 
bar(mean(toneALIGN.psthZ,2),'k'); hold on;
plot(mean(lightALIGN.psthZ,2)./2.5,'r','LineWidth',2);
axis([0 1000 -0.25 0.8]);
xlabel('Time (ms)');
ylabel('Z-score');
title('TONE Aligned PSTH');
box off;
set(gca,'TickDir','out');

%% GET TONE PSTH RESPONSE LATENCIES

window  = 145:245;
window2 = 250:350;

threshold = 4.42;
minPeakWidth = 5;

responseLatency = zeros(size(toneALIGN.psthZ,2),1);
responsePeak    = zeros(size(toneALIGN.psthZ,2),1);
responseWidth   = zeros(size(toneALIGN.psthZ,2),1);
responseAmp     = zeros(size(toneALIGN.psthZ,2),1);
responseMid     = zeros(size(toneALIGN.psthZ,2),1);

for i = 1:size(toneALIGN.psthZ,2)
    
    tmpCheck = toneALIGN.psthZ(window,i);
    check = find(tmpCheck>threshold);

    tmpDATA = toneALIGN.psthZ(window,i);
    
    if numel(check)>minPeakWidth

        peakLoc     = find(tmpDATA(2:99)==max(tmpDATA(2:99)),1);

        beginPeak   = find(tmpDATA(1:peakLoc)<3,1,'last');
        if numel(beginPeak)==0
            beginPeak = peakLoc-1;
        end
        endPeak     = peakLoc+find(tmpDATA(peakLoc+1:numel(tmpDATA))<3,1,'first');
        if numel(endPeak)==0
            endPeak = peakLoc+1;
        end
            
        peakInds = beginPeak:1:endPeak;

%         % check for an earlier peak to make sure nothing is missed
%         firstSupra          = find(tmpDATA(2:99)>threshold,1);
%         
%         if firstSupra < beginPeak
%             responseLatency(i)  = firstSupra;
%         else
            responseLatency(i)  = beginPeak;
%         end
        
        responseAmp(i)      = trapz(tmpDATA(peakInds));
        responsePeak(i)     = peakLoc;
        responseMid(i,1)    = mean(peakInds);
        responseWidth(i)    = endPeak - beginPeak;

    else
        
        responseLatency(i)  = 0;
        responsePeak(i)     = 0;
        responseAmp(i)      = 0;
        responseMid(i,1)    = 0;
        responseWidth(i)    = 0;

    end
    
    toneRESP.latency = responseLatency(i);    
    toneRESP.amp    = responseAmp(i);
    toneRESP.peak   = responsePeak(i);
    toneRESP.mid    = responseMid(i);
    toneRESP.width  = responseWidth(i);

    % store the response parameters back in the complete data structure
    PopData.session(lightALIGN.sess(i)).unit(lightALIGN.unit(i)).toneRESP = toneRESP;
 
end


% store the response parameters back in the complete data structure

[y,toneALIGN.sort.amp.inds]  = sort(responseAmp,1,'descend');
toneALIGN.sort.amp.sortedZ   = toneALIGN.psthZ(:,toneALIGN.sort.amp.inds);
toneALIGN.sort.amp.sortedA   = toneALIGN.psthA(:,toneALIGN.sort.amp.inds);

[y,toneALIGN.sort.peak.inds2] = sort(responsePeak,1,'descend');
toneALIGN.sort.peak.sortedZ  = toneALIGN.psthZ(:,toneALIGN.sort.peak.inds2);
toneALIGN.sort.peak.sortedA  = toneALIGN.psthA(:,toneALIGN.sort.peak.inds2);

% toneALIGN.sort.km.ids        = kmeans([responseAmp,responseLatency],numClusts);
% [y,toneALIGN.sort.amp.inds3] = sort(toneALIGN.sort.km.ids,1,'ascend');
% toneALIGN.sort.amp.sortedZ3  = toneALIGN.psthZ(:,toneALIGN.sort.amp.inds3);
% toneALIGN.sort.amp.sortedS3  = toneALIGN.psthS(:,toneALIGN.sort.amp.inds3);


figure(6);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
subplot(1,5,1:2);
imagesc(toneALIGN.sort.amp.sortedZ(window,:)',[-10 10]);
xlabel('Time (ms)');
ylabel('Cell index');
title('TONE Aligned PSTH');

% subplot(1,5,3);
% plot(toneALIGN.sort.km.ids(toneALIGN.sort.amp.inds3),701:-1:1,'k.'); axis([0 numClusts+1 1 701]); axis off;

subplot(1,5,4:5);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
% imagesc(corr(toneALIGN.sort.amp.sorted),[-1,1]);
imagesc(toneALIGN.sort.peak.sortedZ(window,:)',[-10 10]);
ylabel('Cell index');
xlabel('Time (ms)');
title('TONE Aligned PSTH');

%% GET LIGHT PSTH RESPONSE LATENCIES

window  = 145:295;
window2 = 250:350;

lightRESP.latency = zeros(size(lightALIGN.psthZ,2),1);
lightRESP.peak    = zeros(size(lightALIGN.psthZ,2),1);
lightRESP.width   = zeros(size(lightALIGN.psthZ,2),1);
lightRESP.amp     = zeros(size(lightALIGN.psthZ,2),1);
lightRESP.mid     = zeros(size(lightALIGN.psthZ,2),1);

for i = 1:size(lightALIGN.psthZ,2)
    
    tmpCheck = lightALIGN.psthA(window,i);
    check = find(tmpCheck>4.42);

    tmpDATA = lightALIGN.psthA(window,i);
    
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
        
        lightRESP.amp(i)      = trapz(tmpDATA(peakInds));
        lightRESP.latency(i)  = beginPeak;
        lightRESP.peak(i)     = peakLoc;
        lightRESP.mid(i,1)    = mean(peakInds);
        lightRESP.width(i)    = endPeak - beginPeak;

    else
        
        lightRESP.latency(i)  = 0;
        lightRESP.peak(i)     = 0;
        lightRESP.amp(i)      = 0;
        lightRESP.mid(i,1)    = 0;
        lightRESP.width(i)    = 0;

    end

    lightRESPl.latency = lightRESP.latency(i);    
    lightRESPl.amp    = lightRESP.amp(i);
    lightRESPl.peak   = lightRESP.peak(i);
    lightRESPl.mid    = lightRESP.mid(i);
    lightRESPl.width  = lightRESP.width(i);

    % store the response parameters back in the complete data structure
    PopData.session(lightALIGN.sess(i)).unit(lightALIGN.unit(i)).lightRESP = lightRESPl;
    
end

[y,lightALIGN.sort.amp.inds]  = sort(lightRESP.amp,1,'descend');
lightALIGN.sort.amp.sortedByLight   = lightALIGN.psthZ(:,lightALIGN.sort.amp.inds);
lightALIGN.sort.amp.sortedByTone    = lightALIGN.psthZ(:,toneALIGN.sort.amp.inds);

[y,lightALIGN.sort.peak.inds2] = sort(lightRESP.peak,1,'descend');
lightALIGN.sort.peak.sortedByLight   = lightALIGN.psthZ(:,lightALIGN.sort.peak.inds2);
lightALIGN.sort.peak.sortedByTone    = lightALIGN.psthZ(:,toneALIGN.sort.peak.inds2);

% lightALIGN.sort.km.ids        = kmeans([responseAmp,lightRESP.latency],numClusts);
% [y,lightALIGN.sort.amp.inds3] = sort(lightALIGN.sort.km.ids,1,'ascend');
% lightALIGN.sort.amp.sortedZ3  = lightALIGN.psthZ(:,lightALIGN.sort.amp.inds3);
% lightALIGN.sort.amp.sortedS3  = lightALIGN.psthS(:,lightALIGN.sort.amp.inds3);


figure(7);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
subplot(1,5,1:2);
imagesc(toneALIGN.sort.peak.sortedZ(window,:)',[-30 30]);
xlabel('Time (ms)');
ylabel('Cell index');
title('TONE sorted by TONE PEAK');

% subplot(1,5,3);
% plot(lightALIGN.sort.km.ids(lightALIGN.sort.amp.inds3),701:-1:1,'k.'); axis([0 numClusts+1 1 701]); axis off;

subplot(1,5,4:5);
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);
% imagesc(corr(lightALIGN.sort.amp.sorted),[-1,1]);
imagesc(lightALIGN.sort.peak.sortedByTone(window,:)',[-30 30]);
ylabel('Cell index');
xlabel('Time (ms)');
title('LIGHT sorted by TONE PEAK');

%% SORT BOXCAR RESPONSES BY VALUE OF FIRST POST STIMULUS BIN, THEN BY SECOND POST-STIMULUS BIN
clear exampleInhib forDisp inhibTest


% already normalized to baseline
numCells = size(lightALIGN.boxcar,2);
count = 0; count2 = 0; valid = 1;
inhibLogic = [];
directEx = [];
    
for i = 1:numCells
    if mean(lightALIGN.boxcar(1:24,i)) >= 8
        if mean(lightALIGN.boxcar(1:24,i)) > (mean(lightALIGN.boxcar(25:33,i)).*1.25)% & (lightALIGN.boxcar(26,i)<(mean(lightALIGN.boxcar(1:24,i)) ))
            inhibTest(valid) = mean(lightALIGN.boxcar(26:29,i)) - mean(lightALIGN.boxcar(1:24,i)) +  (max(lightALIGN.boxcar(26:29,i))-min(lightALIGN.boxcar(26:29,i)));
            count = count+1;
            inhibLogic = [inhibLogic i];            
        elseif (mean(lightALIGN.boxcar(26:27,i))>mean(lightALIGN.boxcar(1:24,i).*3))
            inhibTest(valid) = mean(lightALIGN.boxcar(26:29,i)) - mean(lightALIGN.boxcar(1:24,i)) + (max(lightALIGN.boxcar(26:29,i))-min(lightALIGN.boxcar(26:29,i))) + 10000;
            count2 = count2+1;
            directEx = [directEx i];
        else
            inhibTest(valid) = mean(lightALIGN.boxcar(26:29,i)) - mean(lightALIGN.boxcar(1:24,i)) + (max(lightALIGN.boxcar(26:29,i))-min(lightALIGN.boxcar(26:29,i))) + 5000;
        end
        
%         forDisp(:,valid) = (lightALIGN.boxcar(:,i) - mean(lightALIGN.boxcar(1:24,i))) ./ std(lightALIGN.boxcar(1:24,i));
        forDisp(:,valid) = (lightALIGN.psthZ(:,i) - mean(lightALIGN.psthZ(1:150,i))) ./ std(lightALIGN.psthZ(1:150,i));
    
        valid = valid + 1;
    else
        
%         inhibTest(i) = mean(lightALIGN.boxcar(26:40,i)) - mean(lightALIGN.boxcar(1:24,i)) + 5000;
%         forDisp(:,i) = (lightALIGN.boxcar(:,i) - mean(lightALIGN.boxcar(1:24,i)));

    end
    

end
[values,indices] = sort(inhibTest,'ascend');

figure(1);
imagesc(forDisp(:,indices)',[-7 7]);
colormap(mapName);

figure(5);
plot(1:91,mean(lightALIGN.boxcar(:,inhibLogic)'),'ko-',1:91,mean(lightALIGN.boxcar(:,directEx)'),'bo-','LineWidth',2);

disp(' ');
disp(' ');
disp(['Total units matching criterion: ' num2str(valid)]);
disp(['Fraction of population with suppression: N=' num2str(count) ' ... ' num2str(count./valid.*100) '%']);
disp(['Fraction of population with activation: N=' num2str(count2)  ' ... ' num2str(count2./valid.*100) '%']);
disp(' ');
disp(' ');
disp(' ');
exampleInhib    = inhibLogic;
exampleExcite   = directEx;

% Now plot the example cells
forExport = forDisp(:,indices)';
% normalize negative and positive directions
negInds = find(forExport<0);
forExport(negInds) = -1 .* (forExport(negInds) ./ min(min(forExport)));
posInds = find(forExport>0);
forExport(posInds) = forExport(posInds) ./ max(max(forExport));

TNC_ExportMatToIgor(forExport,[1:size(forExport,2)],'AllInVivoStim');

%% RETURN THE ELECTRODE LOCATION OF INHIBITED UNITS AND DIRECT EXCITED UNITS

exNum = numel(directEx);
exTrodes = [];

inNum = numel(inhibLogic);
inTrodes = [];

for i=1:exNum
    
    % current index
    exIn = directEx(i);
    
    % retrieve the session and unit number
    currSess = lightALIGN.sess(exIn);
    currUnit = lightALIGN.unit(exIn);
    
    % get the electrode number
    currName = PopData.session(currSess).unit(currUnit).name;
    
    % append to the list
    exTrodes = [exTrodes , str2num(currName(5:6))];
    
end

for i =1:inNum
    % current index
    inIn = inhibLogic(i);
    
    % retrieve the session and unit number
    currSess = lightALIGN.sess(inIn);
    currUnit = lightALIGN.unit(inIn);
    
    % get the electrode number
    currName = PopData.session(currSess).unit(currUnit).name;
    
    % append to the list
    inTrodes = [inTrodes , str2num(currName(5:6))];
    
end

figure(10); clf;
exFreq = hist(exTrodes,32:1:64)./exNum;
inFreq = hist(inTrodes,32:1:64)./inNum;
corr(exFreq',inFreq')
plot(exFreq,inFreq,'ko'); xlabel('Frequency (excitation)'); ylabel('Frequency (inhibition)');
axis([-0.05 0.25 -0.05 0.25]);

%% GO THROUGH SESSIONS AND PUT UNITS IN ORDER ACCORDING TO ELECTRODE NUMBER
initialize = 0;
numSessions = size(PopData.session,2);
allEs = [];
clear DataByTrode

for m = 2:-1:1
    initialize = m-1;
    
%     for j=1:numSessions
    for j=13
        numUnits = size(PopData.session(j).unit,2);
        if initialize
            disp([num2str(j) ' >>> ' num2str(numUnits)]);
        end

        for i =1:numUnits
            eNumDigits = numel(PopData.session(j).unit(i).name)-1;
%             disp(['Electrode number: ' PopData.session(j).unit(i).name(5:eNumDigits)]);
            currElec = str2num(PopData.session(j).unit(i).name(5:eNumDigits)) - 32;
            allEs = [allEs currElec];

            if initialize
                DataByTrode.elec(currElec).lite=[];
                DataByTrode.elec(currElec).tone=[];
            else
                % append current data to a matrix associated with that electrode
                if numel(DataByTrode.elec(currElec).lite) > 0
                    DataByTrode.elec(currElec).lite = [DataByTrode.elec(currElec).lite PopData.session(j).unit(i).rLITE.psthAVG'];
                    DataByTrode.elec(currElec).tone = [DataByTrode.elec(currElec).tone PopData.session(j).unit(i).rTONE.psthAVG'];
                else
                    DataByTrode.elec(currElec).lite = PopData.session(j).unit(i).rLITE.psthAVG';
                    DataByTrode.elec(currElec).tone = PopData.session(j).unit(i).rTONE.psthAVG';
                end
            end
        end        
    end
end

numElecMats = size(DataByTrode.elec,2);

figure(13); clf; hold on;
for k=1:numElecMats
    subplot(16,2,k);
    plot(mean(DataByTrode.elec(k).lite,2),'LineWidth', 2, 'Color', [k./numElecMats 0 1-(k./numElecMats)]);
%     imagesc(DataByTrode.elec(k).lite',[0 1]);
    axis([140 180 -0.025 0.6]); 
    axis off; 
    drawnow;
end

figure(14); clf;  hold on;
for k=1:numElecMats
    subplot(16,2,k);
    plot(mean(DataByTrode.elec(k).tone,2),'LineWidth', 2, 'Color', [k./numElecMats 0 1-(k./numElecMats)]);
%     imagesc(DataByTrode.elec(k).lite',[0 1]);
    axis([140 180 -0.025 0.4]); 
    axis off; 
    drawnow;
end

%% FIND TRULY 'LIGHT-ACTIVATED' UNITS DEFINED BY LOW JITTER RESPONSES

numSessions = size(PopData.session,2);
allEs = []; count = 0; clear possDirectActivation
plotOn = 1;
window = 150:160;

for j=1:numSessions
    numUnits = size(PopData.session(j).unit,2);
    
    for i=1:numUnits
        % if there is a significant modulation by light and a short latency peak examine the raster
        if PopData.session(j).unit(i).lightRESP.amp > 0

            % look for a modulation of the psthAVG (i.e. low jitter) within 15 ms of light
            % stimulation
            stdOfPavg = std(PopData.session(j).unit(i).rLITE.psthAVG(1:149));            
            if max(PopData.session(j).unit(i).rLITE.psthAVG(window)) > (stdOfPavg*5)              
                posTrials   = 0;
                preciseSpk  = [];
                count       = count+1;
                tmpLat      = find(PopData.session(j).unit(i).rLITE.psthAVG(window) == max(PopData.session(j).unit(i).rLITE.psthAVG(window)),1)-2;
                
                if  plotOn
                    % display the raster plot with psthA to examine jitter manually
                    figure(31); clf;
                    subplot(5,1,1:3); hold on;
                        numTrials = size(PopData.session(j).unit(i).rLITE.raster.trial,2);
                        for k = 1:numTrials
                            numSpks = size(PopData.session(j).unit(i).rLITE.raster.trial(k).ts,1);
                            plot(PopData.session(j).unit(i).rLITE.raster.trial(k).ts,ones(numSpks,1).*k,'k.');
                            numSpkInWin = find(PopData.session(j).unit(i).rLITE.raster.trial(k).ts > (tmpLat-2) & PopData.session(j).unit(i).rLITE.raster.trial(k).ts < (tmpLat+1),1);
                            if numel(numSpkInWin)>0
                                plot(PopData.session(j).unit(i).rLITE.raster.trial(k).ts(numSpkInWin),k,'ro');
                                preciseSpk = [preciseSpk PopData.session(j).unit(i).rLITE.raster.trial(k).ts(numSpkInWin)];
                            end
                        end
                        axis([-24 60 0 numTrials]);
                    title([PopData.session(j).sessId ' >>> ' PopData.session(j).unit(i).name ' >>> ' num2str(count)]);
                    subplot(514);
                       bar(-150:400,PopData.session(j).unit(i).rLITE.psthAVG,'k'); hold on;
                       plot(tmpLat,0,'r^','MarkerSize',10);
                       axis([-24 60 0 max(PopData.session(j).unit(i).rLITE.psthAVG).*1.1]);
                    subplot(515);
                       bar(-24.5:65.5,PopData.session(j).unit(i).rLITE.boxcar,'k');
                       axis([-4.5 10.5 0 max(PopData.session(j).unit(i).rLITE.boxcar).*1.1]);

                    % wait for the user to clear
                    % pause();
                end
                
                % update the candidate list
                possDirectActivation.count(count)  = count;
                possDirectActivation.sessId(count) = j;
                possDirectActivation.unitId(count) = i;
            
                eNumDigits = numel(PopData.session(j).unit(i).name)-1;
                currElec = str2num(PopData.session(j).unit(i).name(5:eNumDigits)) - 32;
                possDirectActivation.elecNum(count) = currElec;

                possDirectActivation.latency(count) = tmpLat;
                numTrials = size(PopData.session(j).unit(i).rLITE.raster.trial,2);
                for k = 1:numTrials
                    numSpkInWin = find(PopData.session(j).unit(i).rLITE.raster.trial(k).ts > (tmpLat-2) & PopData.session(j).unit(i).rLITE.raster.trial(k).ts < (tmpLat+2),1);
                    if numel(numSpkInWin)>0
                        posTrials = posTrials+1;
                    end
                end
                possDirectActivation.reliability(count) = posTrials./numTrials;

                possDirectActivation.jitter(count) = var(preciseSpk);

            end
        end
    end
end

disp(['___________________________________________________'])
disp(['Found ' num2str(count) ' putative direct light responses.']);

%% DEFINE DIRECT ACTIVATION

% Imperfect definition but go from reliability>0.33, jitter<1, latency<8ms
lightEvoked = [];
numCandidates = numel(possDirectActivation.reliability);

for i = 1:numCandidates
    if possDirectActivation.reliability(i) > 0.4
        if possDirectActivation.jitter(i) < 1
            if possDirectActivation.latency < 10
                lightEvoked = [lightEvoked i];
            end
        end
    end    
end

disp(['___________________________________________________'])
disp(['Found ' num2str(numel(lightEvoked)) ' valid light responses.']);

figure(3); plot(possDirectActivation.unitId(lightEvoked),possDirectActivation.latency(lightEvoked),'k.'); grid on;
axis([0 32 0 10]);

%% GET AVERAGE FIRING RATE OF CONFIRMED DIRECT ACTIVATION CELLS

for i =1:numel(lightEvoked)

    currSess = lightALIGN.sess(lightEvoked(i));
    currUnit = lightALIGN.unit(lightEvoked(i));
    
    rate(i) = 1000./mean(PopData.session(currSess).unit(currUnit).isi.instant);
    
end

mean(rate)

%% LOOK FOR A TONE RESPONSE
clear validLight

for i = 1:numel(lightEvoked)

    % plot the light and tone responses
    validLight.toneLAT(i) = PopData.session(possDirectActivation.sessId(i)).unit(possDirectActivation.unitId(i)).toneRESP.latency;
    validLight.toneAMP(i) = PopData.session(possDirectActivation.sessId(i)).unit(possDirectActivation.unitId(i)).toneRESP.amp;    

end

posTone = find(validLight.toneLAT>0);
figure(14); subplot(311);
hist(validLight.toneLAT(posTone),0:5:100);
figure(14); subplot(312);
hist(validLight.toneAMP(posTone));
figure(14); subplot(313);
plot(validLight.toneLAT(posTone),validLight.toneAMP(posTone),'k.');

disp(['___________________________________________________'])
disp(['Found ' num2str(numel(posTone)) ' units positive for both direct light response and tone response.']);

%% FOR A GIVEN SESSION PLOT THE TEMPORAL PATTERN OF CELL ACTIVATION
numSessions = size(PopData.session,2);
% numSessions = 10;
figure(17); clf;

for j=1:numSessions
    numUnits = size(PopData.session(j).unit,2);

    for i=1:numUnits
        if j <= numSessions./2
            subplot(4,numSessions./2,j); hold on;
        else
            subplot(4,numSessions./2,j+(numSessions./2)); hold on;            
        end
        plot(20:80,PopData.session(j).unit(i).rTONE.boxcar(20:80),'Color',[i./numUnits 0 1-(i./numUnits)],'LineWidth',1);
%         axis([140 250 -20 20]); axis off;
        axis([20 80 0 2]); axis off;

        if j <= numSessions./2
            subplot(4,numSessions./2,j+(numSessions./2)); hold on;
        else
            subplot(4,numSessions./2,j+numSessions); hold on;            
        end
        plot(20:80,PopData.session(j).unit(i).rLITE.boxcar(20:80),'Color',[i./numUnits 0 1-(i./numUnits)],'LineWidth',1);
%         axis([140 250 -20 80]); axis off;        
        axis([20 80 0 4]); axis off;
    end

end

%% FOR A GIVEN SESSION PLOT THE TEMPORAL PATTERN OF CELL ACTIVATION
numSessions = size(PopData.session,2);
figure(18); clf;
allSessPlot = 0;

for j=1:numSessions
    numUnits = size(PopData.session(j).unit,2);

    toneMatForCorr=[];
    liteMatForCorr=[];
    
    figure(19); clf;

    for i=1:numUnits
        if allSessPlot
            if j <= numSessions./2
                subplot(4,numSessions./2,j); hold on;
            else
                subplot(4,numSessions./2,j+(numSessions./2)); hold on;            
            end
    %         plot(-10:130,PopData.session(j).unit(i).rsTONE.psthZ(140:280),'Color',[i./numUnits 0 1-(i./numUnits)],'LineWidth',1);
    % %         axis([140 250 -20 20]); axis off;
    %         axis([-10 130 -20 20]); axis off;
            toneMatForCorr = [toneMatForCorr PopData.session(j).unit(i).rsTONE.psthZ(140:280)'];
            imagesc(corr(toneMatForCorr),[-1 1]); axis tight; axis off;

            if j <= numSessions./2
                subplot(4,numSessions./2,j+(numSessions./2)); hold on;
            else
                subplot(4,numSessions./2,j+numSessions); hold on;            
            end
    %         plot(-10:130,PopData.session(j).unit(i).rsLITE.psthZ(140:280),'Color',[i./numUnits 0 1-(i./numUnits)],'LineWidth',1);
    % %         axis([140 250 -20 80]); axis off;        
    %         axis([-10 130 -20 50]); axis off;
            liteMatForCorr = [liteMatForCorr PopData.session(j).unit(i).rsLITE.psthZ(140:280)'];
            imagesc(corr(liteMatForCorr),[-1 1]); axis tight; axis off;
        else
            figure(18);
            subplot(numUnits,1,i);
            bar(PopData.session(j).unit(i).rLITE.psthAVG(140:280),'k');
            axis([-10 130 -0.01 0.35]); axis off;
            drawnow;
%             figure(19); hold on;
%             plot(-10:130,PopData.session(j).unit(i).rsLITE.psthZ(140:280),'LineWidth',2,'Color',[i./numUnits 0 1-(i./numUnits)]);

            figure(19);
            subplot(numUnits,1,i);
            bar(PopData.session(j).unit(i).rTONE.psthAVG(140:280),'k');
            axis([-10 130 -0.01 0.35]); axis off;
            drawnow;
%             figure(19); hold on;
%             plot(-10:130,PopData.session(j).unit(i).rsTONE.psthZ(140:280),'LineWidth',2,'Color',[i./numUnits 0 1-(i./numUnits)]);

        end
    end
    
    if allSessPlot==0         
       pause(2);
    end
end

%% PLOT HAHNLOSER STYLE RASTER PLOTS OF REPEATED TRIALS
window          = 140:190;
sessionToPlot   = 15;
j               = sessionToPlot;
numUnits        = size(PopData.session(j).unit,2);
numTrials       = 30;
% numTrials       = size(PopData.session(j).events.LITE.ts,1);
numRows         = (numUnits.*numTrials) + ((numUnits+1).*20);
firstPeakLats   = zeros(1,numUnits);
avgPSTHz        = [];
offset          = 150-window(1)+3;

timesExport = [];
trialExport = [];
colorExport = [];

for i=1:numUnits

    % look for the order at which the units peak in their psthA
    currPSTHavg         = PopData.session(j).unit(i).rLITE.psthAVG(window);
    tmpIndex            = find(currPSTHavg==max(currPSTHavg),1);
    if tmpIndex>offset
        firstPeakLats(i)    = tmpIndex;
    else
        firstPeakLats(i)    = numel(window);        
    end
    currPSTHz           = PopData.session(j).unit(i).rsLITE.psthZ(window);
    avgPSTHz            = [avgPSTHz currPSTHz'];
end

figure(22); clf;
subplot(811); bar(window-150,mean(avgPSTHz,2),'k'); axis([window(1)-150 window(numel(window))-150 0 max(mean(avgPSTHz,2)).*1.1])

[y,unitsByLat] = sort(firstPeakLats,'ascend');

subplot(8,1,2:8); hold on;
for m=1:numUnits
    thisUnit = unitsByLat(m);

    for k = 1:numTrials
        thisTrialRow    = ((m-1).*numTrials) + (10.*m) + k;
        theseTimeStamps = PopData.session(j).unit(thisUnit).rLITE.raster.trial(k).ts;
        validStamps     = find(theseTimeStamps>window(1)-150 & theseTimeStamps<window(numel(window))-150);

        if rem(m,2)
            colorSpecM = [m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)];
        else
            colorSpecM = [1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits];
        end
        
        plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*thisTrialRow,'.','MarkerSize',6,'Color',colorSpecM);
        if k==1
            plot([y(m)-offset y(m)-offset],[thisTrialRow thisTrialRow+numTrials],'-','LineWidth',2,'Color',colorSpecM);
        end
        
        timesExport = [timesExport ; theseTimeStamps(validStamps)];
        trialExport = [trialExport ; (ones(numel(validStamps),1).*thisTrialRow)];
        colorExport = [colorExport ; (ones(numel(validStamps),1) * colorSpecM)];
        
    end
    
    axis([window(1)-150 window(numel(window))-150 0 numRows]);
    set(gca,'ydir','reverse','TickDir','out');

end

plot([0 0],[0 numRows],'-','Color',[0 0 0]);
% plot([5 5],[0 numRows],'-','Color',[0.25 0.25 0.25]);
% plot([10 10],[0 numRows],'-','Color',[0 0 0]+0.5);

%% PLOT HAHNLOSER STYLE RASTER PLOTS OF REPEATED TRIALS
window          = 100:300;
sessionToPlot   = 10;
j               = sessionToPlot;
numUnits        = size(PopData.session(j).unit,2);
numTrials       = 20;
numTrials       = size(PopData.session(j).events.TONE.ts,1);
numRows         = numUnits.*numTrials;
firstPeakLats   = zeros(1,numUnits);
avgPSTHz        = [];

for i=1:numUnits

    % look for the order at which the units peak in their psthA
    currPSTHavg         = PopData.session(j).unit(i).rTONE.psthAVG(window);
    tmpIndex            = find(currPSTHavg==max(currPSTHavg),1);
    if tmpIndex>50
        firstPeakLats(i)    = tmpIndex;
    else
        firstPeakLats(i)    = numel(window);        
    end
    currPSTHz           = PopData.session(j).unit(i).rsTONE.psthZ(window);
    avgPSTHz            = [avgPSTHz currPSTHz'];
end

figure(23); clf;
subplot(811); bar(window-150,mean(avgPSTHz,2),'k'); axis([window(1)-150 window(numel(window))-150 0 max(mean(avgPSTHz,2)).*1.1])

[y,unitsByLat] = sort(firstPeakLats,'ascend');

subplot(8,1,2:8); hold on;
for m=1:numUnits
    thisUnit = unitsByLat(m);

    for k = 1:numTrials
        thisTrialRow    = ((m-1).*numTrials) + (10.*m) + k;
        theseTimeStamps = PopData.session(j).unit(thisUnit).rTONE.raster.trial(k).ts;
        validStamps     = find(theseTimeStamps>window(1)-150 & theseTimeStamps<window(numel(window))-150);
        if rem(m,2)
            plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*thisTrialRow,'.','Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
            if k==1
                plot([y(m)-52 y(m)-52],[thisTrialRow thisTrialRow+numTrials],'-','LineWidth',2,'Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
            end
        else
            plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*thisTrialRow,'.','Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);            
            if k==1
                plot([y(m)-52 y(m)-52],[thisTrialRow thisTrialRow+numTrials],'-','LineWidth',2,'Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);
            end
        end
    end
    
    axis([window(1)-150 window(numel(window))-150 0 numRows]);
    set(gca,'ydir','reverse','TickDir','out');

end

plot([0 0],[0 numRows],'k--');

%% PLOT HAHNLOSER STYLE RASTER PLOTS OF SPECIFIED UNITS INDICES
window          = 90:210;
listOfInds      = exampleInhib([2:4,6:7,9:18]);

numTrials       = 40;
% numTrials       = size(PopData.session(j).events.LITE.ts,1);
firstPeakLats   = zeros(1,numUnits);
avgPSTHz        = [];
offset          = 150-window(1)+3;

timesExport = [];
trialExport = [];
colorExport = [];

numExamples     = numel(listOfInds);
numUnits        = size(PopData.session(j).unit,2);
numRows         = (numTrials.*numExamples) + ((numExamples+1).*20);

figure(9); clf; hold on;

for m=1:numExamples
    thisIndex   = listOfInds(m);
    thisSession = lightALIGN.sess(thisIndex)
    thisUnit    = lightALIGN.unit(thisIndex)

    for k = 1:numTrials
        thisTrialRow    = ((m-1).*numTrials) + (10.*m) + k;
        theseTimeStamps = PopData.session(thisSession).unit(thisUnit).rLITE.raster.trial(k).ts;
        validStamps     = find(theseTimeStamps>window(1)-150 & theseTimeStamps<window(numel(window))-150);

        if rem(m,2)
            colorSpecM = [m./numExamples 0.67-(0.67.*(m./numExamples)) 1-(m./numExamples)];
        else
            colorSpecM = [1-(m./numExamples) 0.67.*(m./numExamples) m./numExamples];
        end
        
        plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*thisTrialRow,'.','MarkerSize',6,'Color',colorSpecM);
%         if k==1
%             plot([y(m)-offset y(m)-offset],[thisTrialRow thisTrialRow+numTrials],'-','LineWidth',2,'Color',colorSpecM);
%         end
        
        timesExport = [timesExport ; theseTimeStamps(validStamps)];
        trialExport = [trialExport ; (ones(numel(validStamps),1).*thisTrialRow)];
        colorExport = [colorExport ; (ones(numel(validStamps),1) * colorSpecM)];
        
    end
    
%     axis([window(1)-150 window(numel(window))-150 0 numRows]);
    set(gca,'ydir','reverse','TickDir','out');

end

plot([0 0],[0 numRows],'-','Color',[0 0 0]);
% plot([5 5],[0 numRows],'-','Color',[0.25 0.25 0.25]);
% plot([10 10],[0 numRows],'-','Color',[0 0 0]+0.5);

%% SORT ALIGNED PSTHS BY LATENCY AND AMPLITUDE

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

%% SEARCH FOR UNITS WITH LARGE TONE AND LARGE LIGHT RESPONSES
NumSessions = size(PopData.session,2);
k=0;

for i=1:NumSessions
    numUnits = size(PopData.session(i).unit,2);

    for j = 1:numUnits
        k=k+1;
        summaryData.toneAmp(k,1)    = PopData.session(i).unit(j).toneRESP.amp;
        summaryData.liteAmp(k,1)    = PopData.session(i).unit(j).lightRESP.amp;
        summaryData.sessInd(k,1)    = i;
        summaryData.unitInd(k,1)    = j;
    end
    
end

summaryData.sumAmp     = (summaryData.toneAmp(:,1)./max(summaryData.toneAmp(:,1))) + (summaryData.liteAmp(:,1)./max(summaryData.liteAmp(:,1)));
summaryData.difAmp     = (summaryData.toneAmp(:,1)./max(summaryData.toneAmp(:,1))) - (summaryData.liteAmp(:,1)./max(summaryData.liteAmp(:,1)));
summaryData.respInd    = summaryData.sumAmp-summaryData.difAmp;

figure(41);
loglog(summaryData.toneAmp,summaryData.liteAmp,'k.');
axis([1 max(summaryData.toneAmp(:,1)).*10 1 max(summaryData.liteAmp(:,1)).*10]);

figure(42); subplot(5,1,1:3);
plot(summaryData.toneAmp(:,1)./max(summaryData.toneAmp(:,1)),summaryData.liteAmp(:,1)./max(summaryData.liteAmp(:,1)),'k.');
axis([0.01 1 0.01 1]);
subplot(5,1,4);
plot(summaryData.respInd,'k.'); axis([0 k 0 2]);

[summaryData.sorted.value,summaryData.sorted.inds] = sort(summaryData.respInd,'descend')
subplot(5,1,5);
plot(summaryData.sorted.value,'k.'); axis([0 k 0 2]);

disp('done.');

%% FIND THE PSTH BEST MATCHED TO THE MEAN EXCITE AND MEAN INHIB

exNum = numel(directEx);
exCorrs = [];

inNum = numel(inhibLogic);
inCorrs = [];

exMean = mean(lightALIGN.psthA(:,directEx)');
inMean = mean(lightALIGN.psthA(:,inhibLogic)');

for i=1:exNum
   
    % current index
    exIn = directEx(i);
    
    % retrieve the session and unit number
    currSess = lightALIGN.sess(exIn);
    currUnit = lightALIGN.unit(exIn);
    
    % get the electrode number
    currCorr = corr(PopData.session(currSess).unit(currUnit).rLITE.psthAVG(1:551)',exMean');
    
    % append to the list
    exCorrs = [exCorrs , currCorr];
    
end

for i =1:inNum
    % current index
    inIn = inhibLogic(i);
    
    % retrieve the session and unit number
    currSess = lightALIGN.sess(inIn);
    currUnit = lightALIGN.unit(inIn);
    
    % get the electrode number
    currCorr = corr(PopData.session(currSess).unit(currUnit).rLITE.psthAVG(1:551)',inMean');
    
    % append to the list
    inCorrs = [inCorrs , currCorr];
    
end

bestExcite = directEx(find(exCorrs==max(exCorrs)))
bestInhib = inhibLogic(find(inCorrs==max(inCorrs)))

%% EXPORT A GOOD EXAMPLE UNIT

% possUnits = find(summaryData.toneAmp(:,1)./max(summaryData.toneAmp(:,1))>0.5 & summaryData.liteAmp(:,1)./max(summaryData.liteAmp(:,1))>0.2);
% egIndex = possUnits(1)

% egIndex = 89; % DA
egIndex = bestExcite;
% egIndex = bestInhib;

sessInd = lightALIGN.sess(egIndex);
unitInd = lightALIGN.unit(egIndex);
spkSTAMPS = [];
trialLIST = [];

numTrials = size(PopData.session(sessInd).unit(unitInd).rLITE.raster.trial,2);
% numTrials = size(PopData.session(sessInd).unit(unitInd).rTONE.raster.trial,2);

for n = 1:numTrials
    thisStamps = PopData.session(sessInd).unit(unitInd).rLITE.raster.trial(n).ts;
%     thisStamps = PopData.session(sessInd).unit(unitInd).rTONE.raster.trial(n).ts;

    spkSTAMPS = [spkSTAMPS thisStamps'];
    trialLIST = [trialLIST ones(1,numel(thisStamps)).*n];
    
end

exportToIgor1 = [trialLIST' spkSTAMPS'];
exportToIgor2 = [PopData.session(sessInd).unit(unitInd).rLITE.boxcar'];
% exportToIgor2 = [PopData.session(sessInd).unit(unitInd).rTONE.boxcar'];
exportToIgor3 = [PopData.session(sessInd).unit(unitInd).rsLITE.psthZ' PopData.session(sessInd).unit(unitInd).rsLITE.psthZe'];
% exportToIgor3 = [PopData.session(sessInd).unit(unitInd).rsTONE.psthZ' PopData.session(sessInd).unit(unitInd).rsTONE.psthZe'];

figure(22); clf;
plot(PopData.session(sessInd).unit(unitInd).rsLITE.psthZ)
% plot(PopData.session(sessInd).unit(unitInd).rsTONE.psthZ)

%% LOAD THE IN VITRO DATA
figure(19); clf;
figure(20); clf;
report=[];
n=0;o=0;
% Load Jenny's data into a local structure
allFiles = dir('Spk_*');

for i = 1:size(allFiles,1)

    disp(' ');
    report = 'DETAILS: ';

    disp('__________________________________________________');
    disp(allFiles(i).name);   
    disp(' ');    

    filename = allFiles(i).name;

    info = h5info(filename,'/spkTimes');

    data.pairs(i).name = allFiles(i).name;
    
    for j = 1:size(info.Groups.Groups,1)

        data.pairs(i).pharm(j).numGroups = size(info.Groups.Groups,1);
        
        for k=1:size(info.Groups.Groups(j).Datasets,1)

            nameToLoad  = [info.Groups.Groups(j).Name '/' info.Groups.Groups(j).Datasets(k).Name];
            tmpTimes    = h5read(filename,nameToLoad);

            msTimes = ceil(double(tmpTimes).*1000);
            
            if max(msTimes)<5001
                msTimes = msTimes+2500;
                data.pairs(i).pharm(j).stim(k).times = msTimes;
                disp('Shifting stim time.');
                totWin=[2500:3500];
            else
                totWin = 1:10000;
                data.pairs(i).pharm(j).stim(k).times = msTimes;
            end

            if currParams.filter.causal==1
                msTimes = msTimes + (numel(currParams.filter.kernel)./2);
            end
            
            delta = zeros(1,10000 + (numel(currParams.filter.kernel)./2));            
            for m=1:numel(msTimes)                
                delta(msTimes(m)) = delta(msTimes(m)) + 1;
            end
            
            data.pairs(i).pharm(j).stim(k).delta = delta;            
            tmpSmooth = conv(delta,currParams.filter.kernel,'same');

%              data.pairs(i).pharm(j).stim(k).psthZ = (tmpSmooth-mean(tmpSmooth(totWin))) ./ std(tmpSmooth(totWin));
             data.pairs(i).pharm(j).stim(k).psthZ = (tmpSmooth-mean(tmpSmooth(totWin)));
            
            if numel(strfind(info.Groups.Groups(j).Name,'GABA'))>0
                data.pairs(i).pharm(j).id = 'gabazine';
            else
                data.pairs(i).pharm(j).id = 'control';
            end
            
            if numel(strfind(info.Groups.Groups(j).Datasets(k).Name,'G1'))>0 | numel(strfind(info.Groups.Groups(j).Datasets(k).Name,'C1'))>0
                data.pairs(i).pharm(j).stim(k).cellId = 1;
            else
                data.pairs(i).pharm(j).stim(k).cellId = 2;                
            end
            
            lch = numel(info.Groups.Groups(j).Datasets(k).Name);
            uin = strfind(info.Groups.Groups(j).Datasets(k).Name,'_');            
            data.pairs(i).pharm(j).stim(k).duration = str2num(info.Groups.Groups(j).Datasets(k).Name(uin(2)+1:lch-2));


            if numel(strfind(data.pairs(i).pharm(j).id,'control'))>0
                n=n+1;
                allIVpsth.cntrl.psthZ(:,n) = data.pairs(i).pharm(j).stim(k).psthZ;
                figure(19); hold on;
                plot(data.pairs(i).pharm(j).stim(k).psthZ); hold on;
                title(nameToLoad); %axis([4800 5200 -1 2]); drawnow;
            else
                o=o+1;
                allIVpsth.gabaz.psthZ(:,o) = data.pairs(i).pharm(j).stim(k).psthZ;                
                figure(20); hold on;
                plot(data.pairs(i).pharm(j).stim(k).psthZ); hold on;
                title(nameToLoad); %axis([4800 5200 -1 2]); drawnow;
            end
            
        end
        
    end
    
end

inVitro.data = data;
disp(report);
disp('__________________________________________________');
disp(' ');

disp(' ');
disp('COMPLETED.');
disp(' ');

%% Calculat traditional FI curves for light stimulation

%% CALCULATE THE GABAZINE to CONTROL CONTRAST

numPairs = size(data.pairs,2);
count = 1;
cntrlVec = [];
gabaVec = [];
stimVecC = [];
stimVecG = [];
cellIDC = [];
cellIDG = [];
cntrlBG = [];
gabaBG  = [];

for i=1:numPairs
    for j=1:size(data.pairs(i).pharm,2)        
        for k=1:size(data.pairs(i).pharm(j).stim,2)

            win = 4000:5000;
%             win = 2600:3400;
            bgR = sum(data.pairs(i).pharm(j).stim(k).delta(win));% - (sum(data.pairs(i).pharm(j).stim(k).delta(cntrl))./numel(cntrl).*numel(win));
            
            win = 5000:5150;
%             win = 3500:4600;
            tmp = trapz(data.pairs(i).pharm(j).stim(k).psthZ(win));            
            
                if j==1
                    cntrlVec    = [cntrlVec tmp./10];
                    cntrlBG     = [cntrlBG bgR./10];
                    stimVecC    = [stimVecC data.pairs(i).pharm(j).stim(k).duration];
                    cellIDC     = [cellIDC i+(9.*(data.pairs(i).pharm(j).stim(k).cellId-1))];
                else
                    gabaBG      = [gabaBG bgR./10];
                    gabaVec     = [gabaVec tmp./10];
                    stimVecG    = [stimVecG data.pairs(i).pharm(j).stim(k).duration];
                    cellIDG     = [cellIDG i+(9.*(data.pairs(i).pharm(j).stim(k).cellId-1))];
                end                       
        end
    end
end

allCells = unique(cellIDC);
cntrlVecN = cntrlVec;
gabaVecN = gabaVec;

% toChange = find(stimVecC==16);
% stimVecC(toChange) = 12;
% toChange = find(stimVecG==16);
% stimVecG(toChange) = 12;
% 
% toChange = find(stimVecC==25);
% stimVecC(toChange) = 20;
% toChange = find(stimVecG==25);
% stimVecG(toChange) = 20;

for k=1:numel(allCells)
    
    inds = find(cellIDC == allCells(k));

    % Normalize by the max amplitude of GABA
        peakResponse    = max(gabaVec(inds));
        gabaVecN(inds)  = gabaVec(inds) ./ peakResponse;
        cntrlVecN(inds) = cntrlVec(inds) ./ peakResponse;
  
end

% Allow negative numbers?
% negs = find(cntrlVecN<0);
% cntrlVecN(negs) = 0;

% figure(5);
% p = anovan([cntrlVec,gabaVec]',{[stimVecC,stimVecG] [zeros(1,numel(cntrlVec)),ones(1,numel(gabaVec))]})
% ranksum(cntrlVec,gabaVec)

stimBins = unique([stimVecC stimVecG])

cntrlBinMean = zeros(1,numel(stimBins));
cntrlBinSEM = zeros(1,numel(stimBins));
gbzBinMean = zeros(1,numel(stimBins));
gbzBinSEM = zeros(1,numel(stimBins));

for m=1:numel(stimBins)
    theseInds       = find(stimVecC==stimBins(m));
    cntrlBinMean(m) = mean(cntrlVecN(theseInds));
    cntrlBinSEM(m)  = std(cntrlVecN(theseInds)) ./ sqrt(numel(theseInds));
    
    theseInds       = find(stimVecG==stimBins(m));
    gbzBinMean(m)   = mean(gabaVecN(theseInds));
    gbzBinSEM(m)    = std(gabaVecN(theseInds)) ./ sqrt(numel(theseInds));
end


% figure(1); clf;
% scatter(gabaVec,cntrlVec,20,cellIDC);
% hold on; plot([-2 12],[-2 12],'r--');
% % axis([-2 12 -2 12]);
% xlabel('GABAzine');
% ylabel('Control');

figure(2); clf;
scatter(stimVecC,cntrlVecN,30,[0 0.67 1]); hold on;
scatter(stimVecG,gabaVecN,30,[1 0 0]);
% axis([0 21 -0.25 1.25]);
xlabel('Stim');
ylabel('Normalized Response');

figure(3);
plot(stimBins,cntrlBinMean,'co-',stimBins,gbzBinMean,'ro-');

forIgor = [cntrlVecN',gabaVecN',stimVecG'];
forIgor2 = [cntrlBinMean',cntrlBinSEM',gbzBinMean',gbzBinSEM',stimBins'];

cntMax = find(stimVecC==max(stimBins));
gbzMax = find(stimVecG==max(stimBins));

forIgor3 = [cntrlVec(cntMax)',cellIDC(cntMax)',gabaVec(gbzMax)',cellIDG(gbzMax)'];

%% GRAB IN VITRO PHOTOSTIM DATA ACCORDING TO STIM STRENGTH
aa=1; bb=1; cc=1; dd=1; ee=1;
founda=0;foundb=0;foundc=0;foundd=0;founde=0;

disp(' ');disp(' ');disp(' ');disp(' ');
disp(['Analyzing all paired spiking experiments in memory']);
disp('__________________________________________________');

numPairs = size(inVitro.data.pairs,2)

for i = 1:numPairs

    if founda
        aa = aa+1;
        founda=0;
    end

    if foundb
        bb = bb+1;
        foundb=0;
    end

    if foundc
        cc = cc+1;
        foundc=0;
    end

    if foundd
        dd = dd+1;
        foundd=0;
    end

    if founde
        ee = ee+1;
        founde=0;
    end

    numConditions = size(inVitro.data.pairs(i).pharm,2);
    for j=1:numConditions
        
        numStims = size(inVitro.data.pairs(i).pharm(j).stim,2);
        for k = 1:numStims
            
            disp(['Found stim duration of: ' num2str(inVitro.data.pairs(i).pharm(j).stim(k).duration) ' ms; cell ' num2str(inVitro.data.pairs(i).pharm(j).stim(k).cellId) ' in pair ' num2str(i)]);
            disp(inVitro.data.pairs(i).pharm(j).id)
            switch inVitro.data.pairs(i).pharm(j).stim(k).duration

                case 4
                    if aa==5
                        inVitro.data.pairs(i)
                    end
                    if numel(strfind(inVitro.data.pairs(i).pharm(j).id,'gaba'))>0
                        inVitro.stimGroup(1).gabazine(aa).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;
                    else
                        inVitro.stimGroup(1).control(aa).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;               
                    end
                    founda = 1;

                case 8
                    if numel(strfind(inVitro.data.pairs(i).pharm(j).id,'gaba'))>0
                        inVitro.stimGroup(2).gabazine(bb).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;
                    else
                        inVitro.stimGroup(2).control(bb).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;               
                    end
                    foundb = 1;

                case 12
                    if numel(strfind(inVitro.data.pairs(i).pharm(j).id,'gaba'))>0
                        inVitro.stimGroup(3).gabazine(cc).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;
                    else
                        inVitro.stimGroup(3).control(cc).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;               
                    end
                    foundc = 1;

                case 16
                    if numel(strfind(inVitro.data.pairs(i).pharm(j).id,'gaba'))>0
                        inVitro.stimGroup(4).gabazine(dd).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;
                    else
                        inVitro.stimGroup(4).control(dd).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;               
                    end
                    foundd = 1;

                case 20
                    if numel(strfind(inVitro.data.pairs(i).pharm(j).id,'gaba'))>0
                        inVitro.stimGroup(5).gabazine(ee).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;
                    else
                        inVitro.stimGroup(5).control(ee).psth(inVitro.data.pairs(i).pharm(j).stim(k).cellId,:) = data.pairs(i).pharm(j).stim(k).psthZ;               
                    end
                    founde = 1;

                otherwise
                    disp('Oops. found an unknown stim condition.');
                    inVitro.data.pairs(i).pharm(j).stim(k).duration
            end
        end    
    end
end

disp(['Finished building stim database.']);
disp('__________________________________________________');
disp(' ');disp(' ');disp(' ');disp(' ');

%% GET THE GABAZ EFFECT SIZE FOR EACH PAIR AND RECORD FOR X vs. Y PLOT
win = 4995:5020;
winD= 4990:5100;
count = 1;
countD= 1;
numGroups = size(inVitro.stimGroup,2);

for i = 1:numGroups

    if i ~= 4
        numExps = size(inVitro.stimGroup(i).control,2);

        for j=1:numExps
            pair = size(inVitro.stimGroup(i).control(j).psth,1)
            if pair==2
                pairDeltaGaba(1,count) = (sum(inVitro.stimGroup(i).gabazine(j).psth(1,win)) - sum(inVitro.stimGroup(i).control(j).psth(1,win))) ./ sum(inVitro.stimGroup(i).gabazine(j).psth(1,win));
                pairDeltaGaba(2,count) = (sum(inVitro.stimGroup(i).gabazine(j).psth(2,win)) - sum(inVitro.stimGroup(i).control(j).psth(2,win))) ./ sum(inVitro.stimGroup(i).gabazine(j).psth(2,win));
                count = count+1;
                
                % get the duration effect
                tmpC1 = inVitro.stimGroup(i).control(j).psth(1,winD);
                tmpG1 = inVitro.stimGroup(i).gabazine(j).psth(1,winD);
                tmpC2 = inVitro.stimGroup(i).control(j).psth(2,winD);
                tmpG2 = inVitro.stimGroup(i).gabazine(j).psth(2,winD);
                
                respWidth(2,countD) = find(tmpG1>max(tmpG1)/2,1,'last') - find(tmpG1>max(tmpG1)/2,1,'first');
                countD = countD+1;
                if max(tmpC1)>1.5
                    respWidth(1,countD) = find(tmpC1>max(tmpC1)/2,1,'last') - find(tmpC1>max(tmpC1)/2,1,'first');
                else
                    respWidth(1,countD) = 0;
                end

                
                respWidth(2,countD) = find(tmpG2>max(tmpG2)/2,1,'last') - find(tmpG2>max(tmpG2)/2,1,'first');
                countD = countD+1;
                if max(tmpC1)>1.5
                    respWidth(1,countD) = find(tmpC2>max(tmpC2)/2,1,'last') - find(tmpC2>max(tmpC2)/2,1,'first');
                else
                    respWidth(1,countD) = 0;
                end
                
            end
        end
    end
end

plot(pairDeltaGaba(1,:),pairDeltaGaba(2,:),'ko')

%% AVERAGE EACH STIMULUS CONDITION FOR CONTROL, GABA, CONTRAST
win = 4950:5350;
numGroups = size(inVitro.stimGroup,2);

for i = 1:numGroups

    numExps = size(inVitro.stimGroup(i).control,2);
    figure(10); clf;
    
    for j=1:numExps
        subplot(3,numExps,j);
        if size(inVitro.stimGroup(i).control(j).psth,1)>1
            bar(mean(lightALIGN.psthZ(:,inhibLogic),2),'c'); hold on;
            plot(win-4845,inVitro.stimGroup(i).control(j).psth(1,win),'b',win-4845,inVitro.stimGroup(i).control(j).psth(2,win),'r');            
            axis([100 300 -10 10]);
        else
            bar(mean(lightALIGN.psthZ(:,inhibLogic),2),'c'); hold on;
            plot(win-4845,inVitro.stimGroup(i).control(j).psth(1,win),'b');
            axis([100 300 -10 10]);
        end
        
        subplot(3,numExps,j+numExps);
        if size(inVitro.stimGroup(i).control(j).psth,1)>1
            hold on;
            plot(win-4845,inVitro.stimGroup(i).control(j).psth(1,win),'b',win-4845,inVitro.stimGroup(i).control(j).psth(2,win),'r');            
            plot(win-4845,inVitro.stimGroup(i).gabazine(j).psth(1,win),'k',win-4845,inVitro.stimGroup(i).gabazine(j).psth(2,win),'c');            
            axis([100 300 -10 10]);
        else
            hold on;
            plot(win-4845,inVitro.stimGroup(i).control(j).psth(1,win),'b');
            plot(win-4845,inVitro.stimGroup(i).gabazine(j).psth(1,win),'k');
            axis([100 300 -10 10]);
        end
        
        subplot(3,numExps,j+(2*numExps));
        if size(inVitro.stimGroup(i).control(j).psth,1)>1
            bar(mean(lightALIGN.psthZ(:,inhibLogic),2),'c'); hold on;
            plot(win-4845,inVitro.stimGroup(i).control(j).psth(1,win)-inVitro.stimGroup(i).gabazine(j).psth(1,win),'b',win-4845,inVitro.stimGroup(i).control(j).psth(2,win)-inVitro.stimGroup(i).gabazine(j).psth(2,win),'r');            
            axis([100 300 -10 10]);
        else
            bar(mean(lightALIGN.psthZ(:,inhibLogic),2),'c'); hold on;
            plot(win-4845,inVitro.stimGroup(i).control(j).psth(1,win)-inVitro.stimGroup(i).gabazine(j).psth(1,win),'b');
            axis([100 300 -10 10]);
        end
        
    end
    
    drawnow; pause();
    
end

%% Create IO Curves for control and GABA

%% AVERAGE ALL CONTROL IN VITRO STIM PSTHS
numGroups = size(inVitro.stimGroup,2);
count = 1;
clear sorter

for i=1:numGroups
    
    if i~=4
        numControls = size(inVitro.stimGroup(i).control,2);
        for j =1:numControls

            numPsths = size(inVitro.stimGroup(i).control(j).psth,1);
            for k=1:numPsths

                allInVitro(count,:) = inVitro.stimGroup(i).control(j).psth(k,:);
                allInVitroG(count,:) = inVitro.stimGroup(i).gabazine(j).psth(k,:);
                count = count+1;

            end
        end
    end
end

count = count-1;

for p=1:count
    sorter(p) = mean(allInVitro(p,5000:5025)) - mean(allInVitro(p,1:4500));    
%     sorter(p) = mean(allInVitro(p,5000:5025)) - mean(allInVitroG(p,5000:5025));    
end

[vals,inds] = sort(sorter,'ascend');
figure(11); clf;

% Now export the data
forExport = allInVitro(inds,5000-currParams.winParams.prior:5000+currParams.winParams.after);
% % normalize negative and positive directions
% negInds = find(forExport<0);
% forExport(negInds) = -1 .* (forExport(negInds) ./ min(min(forExport)));
% posInds = find(forExport>0);
% forExport(posInds) = forExport(posInds) ./ max(max(forExport));
% TNC_ExportMatToIgor(forExport,[1:651],'AllInVitroControls');

subplot(311);imagesc(forExport,[-1 1]);
colormap(mapName);

% Now export the data
forExport = allInVitroG(inds,5000-currParams.winParams.prior:5000+currParams.winParams.after);
% % normalize negative and positive directions
% negInds = find(forExport<0);
% forExport(negInds) = -1 .* (forExport(negInds) ./ min(min(forExport)));
% posInds = find(forExport>0);
% forExport(posInds) = forExport(posInds) ./ max(max(forExport));
TNC_ExportMatToIgor(forExport,[1:651],'AllInVitroGabazine');

subplot(312);imagesc(forExport,[-1 1]);
colormap(mapName);

% Now export the data
forExport = allInVitro(inds,5000-currParams.winParams.prior:5000+currParams.winParams.after)-allInVitroG(inds,5000-currParams.winParams.prior:5000+currParams.winParams.after);
% % normalize negative and positive directions
% negInds = find(forExport<0);
% forExport(negInds) = -1 .* (forExport(negInds) ./ min(min(forExport)));
% posInds = find(forExport>0);
% forExport(posInds) = forExport(posInds) ./ max(max(forExport));
TNC_ExportMatToIgor(forExport,[1:651],'AllInVitroContrast');
    
subplot(313);imagesc(forExport,[-1 1]);
colormap(mapName);

%% HAHNLOSER STYLE PLOT FOR EACH STIMULUS CONDITION
m=0;
stimDur = 20; group = 5;
% stimDur = 16; group = 4;
% stimDur = 12; group = 3;
% stimDur = 8; group = 2;
% stimDur = 4; group = 1;
pharmTest = 'gabazine';
% pharmTest = 'control';
daCell = 1;

timesExport = [];
trialExport = [];
colorExport = [];

numGroups = size(inVitro.data.pairs,2);

for i = 1:numGroups

    numExps = size(inVitro.data.pairs(i).pharm,2);
    for j=1:numExps

        if strcmp(inVitro.data.pairs(i).pharm(j).id,pharmTest)
            
            numStim = size(inVitro.data.pairs(i).pharm(j).stim,2);
            for k = 1:numStim

                if inVitro.data.pairs(i).pharm(j).stim(k).duration == stimDur

                    m = m+1;
                    inVitro.stimGroup(group).raster(m).times          = inVitro.data.pairs(i).pharm(j).stim(k).times;
                    inVitro.stimGroup(group).raster(m).trialLabels    = inVitro.data.pairs(i).pharm(j).stim(k).trialLabels;
                    inVitro.stimGroup(group).psthData2(m,:)           = inVitro.data.pairs(i).pharm(j).stim(k).psthZ;
                    inVitro.stimGroup(group).delta(m,:)               = inVitro.data.pairs(i).pharm(j).stim(k).delta;
                    
                    if strcmp('gabazine',pharmTest)
                        if daCell
                            numPost = numel(find(inVitro.stimGroup(group).raster(m).times>=5000 & inVitro.stimGroup(group).raster(m).times<5100));
                            numPre  = numel(find(inVitro.stimGroup(group).raster(m).times>0 & inVitro.stimGroup(group).raster(m).times<5000));                            
                            numPre = numPre ./ 50;
                        else
                            numPost = numel(find(inVitro.stimGroup(group).raster(m).times>=5000 & inVitro.stimGroup(group).raster(m).times<5000+(stimDur*2)));
                            numPre  = numel(find(inVitro.stimGroup(group).raster(m).times>(5000-(stimDur*2)) & inVitro.stimGroup(group).raster(m).times<5000));
                        end
                        evokedGabaSpks(group,m) = (numPost - numPre) ./ max(inVitro.data.pairs(i).pharm(j).stim(k).trialLabels);
                    else
                        if daCell
                            numPost = numel(find(inVitro.stimGroup(group).raster(m).times>=5000 & inVitro.stimGroup(group).raster(m).times<5100));
                            numPre  = numel(find(inVitro.stimGroup(group).raster(m).times>0 & inVitro.stimGroup(group).raster(m).times<5000)); 
                            numPre = numPre ./ 50;
                        else
                            numPost = numel(find(inVitro.stimGroup(group).raster(m).times>=5000 & inVitro.stimGroup(group).raster(m).times<5000+(stimDur*2)));
                            numPre  = numel(find(inVitro.stimGroup(group).raster(m).times>(5000-(stimDur*2)) & inVitro.stimGroup(group).raster(m).times<5000));
                        end
                        evokedCntrlSpks(group,m) = (numPost - numPre) ./ max(inVitro.data.pairs(i).pharm(j).stim(k).trialLabels);
                    end
                    
                end
            end        
        end
    end
    
end

% if strcmp('control',pharmTest)
%    tmpCnt = mean(inVitro.stimGroup(group).delta,1);
% else
%    tmpGaba = mean(inVitro.stimGroup(group).delta,1);    
% end
% 
% % now the hahnloser style plot
% numRasters = size(inVitro.stimGroup(group).raster,2)
% if strcmp('gabazine',pharmTest)
%     figure(14); clf; hold on;
% else
%     figure(13); clf; hold on;
% end
% subplot(8,1,1);
% % bar(mean(inVitro.stimGroup(group).psthData2,1),'k');
% plot(1:10000,inVitro.stimGroup(group).psthData2,'LineWidth',2);
% axis([4950 5100 -15 70]);
% clear peakInd
% 
% for n=1:numRasters
%     peakInd(n) = find(inVitro.stimGroup(group).psthData2(n,5000:5150) == max(inVitro.stimGroup(group).psthData2(n,5000:5150)));
% end
% 
% % if strcmp('control',pharmTest)
% %     clear ys
% %     clear inds
% %     clear p
% %     [ys,inds] = sort(peakInd+5000,'descend');
% % end
% 
% subplot(8,1,2:8); hold on;
% for n=1:numRasters
% 
% %     if strcmp('control',pharmTest)
% %         p = inds(n);        
% %     else
%         p=n;
% %     end
% 
%     if rem(n,2)
%         colorSpecM = [n./numRasters 0.67-(0.67.*(n./numRasters)) 1-(n./numRasters)];
%     else
%         colorSpecM = [1-(n./numRasters) 0.67.*(n./numRasters) n./numRasters];
%     end
%     
%     plot([ys(n) , ys(n)],[((n-1).*10)+(5*n) , ((n).*10)+(5*n)],'-','LineWidth',2,'Color',colorSpecM);
%     plot(inVitro.stimGroup(group).raster(p).times,inVitro.stimGroup(group).raster(p).trialLabels+((n-1).*10)+(5*n),'.','MarkerSize',6,'Color',colorSpecM)
%     
%     % export for igor
%     timesExport = [timesExport ; inVitro.stimGroup(group).raster(p).times];
%     trialExport = [trialExport ; inVitro.stimGroup(group).raster(p).trialLabels+((n-1).*10)+(2*n)];
%     colorExport = [colorExport ; (ones(numel(inVitro.stimGroup(group).raster(p).times),1) * colorSpecM)];
% 
% end
% 
% axis([4950 5100 0 (numRasters*10) + (5*n)]);

%% LOOK AT ISI DISTS +/- GBZ
m=0;n=0;
numPairs = size(inVitro.data.pairs,2);

for i = 1:numPairs

    numExps = size(inVitro.data.pairs(i).pharm, 2);
    
    for j=1:numExps
        
        numStim = size(inVitro.data.pairs(i).pharm(j).stim,2);

            if numel(strfind(inVitro.data.pairs(i).pharm(j).id,'gaba'))>0
                for k=1:numStim
                    
                    m = m+1;

                    currTimes = inVitro.data.pairs(i).pharm(j).stim(k).times;
                    baseTimes = find(currTimes<5000);
                    isiTimes = diff(currTimes);
                    trialBounds = find(isiTimes<0);
                    trialLabels = zeros(numel(currTimes),1);
                    
                    % put in trial markers
                    for p =1:numel(trialBounds)
                        if p==1
                            trialLabels(1:trialBounds(p)) = p.*ones(trialBounds(p),1);
                        else
                            trialLabels(trialBounds(p-1):trialBounds(p)) = p.*ones(trialBounds(p)-trialBounds(p-1)+1,1);                            
                        end
                    end
                                        
                    % get baseline isis
                    baseIsi = isiTimes(baseTimes(1:numel(baseTimes))); 
%                     baseIsi = isiTimes(find(isiTimes>0));
                    
                    % store hist
                    baseIsiH = hist(baseIsi,0:1:200);
                    
                    % store cv
                    cv = std(baseIsi)./median(baseIsi);

                    figure(3); clf;
%                     plot(currTimes,trialLabels,'k.','MarkerSize',2); drawnow;
                    bar(0:1:200,baseIsiH,'k'); title(num2str(cv));

                    inVitro.data.pairs(i).pharm(j).stim(k).isi = baseIsi;
                    inVitro.data.pairs(i).pharm(j).stim(k).trialLabels = trialLabels;

                    if cv>1
%                         inVitro.data.pairs(i).name
                        plot(currTimes,trialLabels,'k.','MarkerSize',2); drawnow;
                        pause();
                        allGabaCV(m) = 0;                        
                        allGabaMean(m) = mean(baseIsi);
                    else
                        % compile all gaba
                        allGabaCV(m) = cv;
                        allGabaMean(m) = mean(baseIsi);
                    end
                    
%                     disp([num2str(m) ',' num2str(cv) ',' num2str(inVitro.data.pairs(i).pharm(j).stim(k).duration)]);
                    
                end        
            else
                for k=1:numStim
                    
                    n = n+1;

                    currTimes = inVitro.data.pairs(i).pharm(j).stim(k).times;
                    baseTimes = find(currTimes<5000);
                    isiTimes = diff(currTimes);
                    trialBounds = find(isiTimes<0);
                    trialLabels = zeros(numel(currTimes),1);
                    
                    % put in trial markers
                    for p =1:numel(trialBounds)
                        if p==1
                            trialLabels(1:trialBounds(p)) = p.*ones(trialBounds(p),1);
                        else
                            trialLabels(trialBounds(p-1):trialBounds(p)) = p.*ones(trialBounds(p)-trialBounds(p-1)+1,1);                            
                        end
                    end

                    % get baseline isis
                    baseIsi = isiTimes(baseTimes(1:numel(baseTimes))); 
%                     baseIsi = isiTimes(find(isiTimes>0));
                    
                    % store hist
                    baseIsiH = hist(baseIsi,0:1:200);

                    % store cv
                    cv = std(baseIsi)./median(baseIsi);

                    figure(3); clf;
%                     plot(currTimes,trialLabels,'k.','MarkerSize',2); drawnow;
                    bar(0:1:200,baseIsiH,'k'); title(num2str(cv));
                                        
                    if cv>1
                        inVitro.data.pairs(i).name
                    end
                    
                    inVitro.data.pairs(i).pharm(j).stim(k).isi = baseIsi;
                    inVitro.data.pairs(i).pharm(j).stim(k).trialLabels = trialLabels;

                    
                    % compile all gaba
                    allCntrCV(n) = cv;
                    allCntrMean(n) = mean(baseIsi);

%                     disp([num2str(n) ',' num2str(cv) ',' num2str(inVitro.data.pairs(i).pharm(j).stim(k).duration)]);
                    disp([num2str(n) ',' inVitro.data.pairs(i).name ',' num2str(inVitro.data.pairs(i).pharm(j).stim(k).cellId)]);

                end        
            end
    end
    
end

figure(3); clf;
plot(allGabaCV,allCntrCV,'k.',[0 1], [0 1],'k--')
xlabel('+GBZ'); ylabel('Control');

%% average in vitro psths

figure(8); clf;
subplot(4,2,[1,3,5]);
[mapName] = TNC_CreateRBColormap(1024,'rb');
imagesc(lightALIGN.psthZ',[-10 10]);
colormap(mapName);
ylabel('Cell index');
title('LITE Aligned PSTH');

subplot(4,2,[2,4,6]);
imagesc((allIVpsth.cntrl.psthZ(4850:5400,:)-allIVpsth.gabaz.psthZ(4850:5400,:))',[-10 10]);
colormap(mapName);
ylabel('Cell index');
title('IVlite Aligned PSTH');

subplot(4,2,7);
bar(mean(lightALIGN.psthZ,2),'k');
axis([0 550 -5 15]);
xlabel('Time (ms)');
ylabel('Z-score');
title('LITE Aligned PSTH');

subplot(4,2,8); 
bar(mean(allIVpsth.cntrl.psthZ(4850:5400,:),2),'k'); hold on;
plot(mean(allIVpsth.gabaz.psthZ(4850:5400,:),2),'b');
plot(mean(lightALIGN.psthZ,2),'r');

axis([0 550 -5 15]);
xlabel('Time (ms)');
ylabel('Z-score');
title('IVlite Aligned PSTH');

avgCntrl = mean(allIVpsth.cntrl.psthZ(4850:5400,:),2);
avgGabaz = mean(allIVpsth.gabaz.psthZ(4850:5400,:),2);
avgInViv = mean(lightALIGN.psthZ,2);
avgDiff  = mean(allIVpsth.cntrl.psthZ(4850:5400,:)-allIVpsth.gabaz.psthZ(4850:5400,:),2);

figure(9);
plot(-146:404,avgCntrl./max(avgCntrl),'k',-146:404,avgGabaz./max(avgGabaz),'b',-150:400,avgInViv./max(avgInViv),'r',-146:404,avgDiff,'k--');

%% COMPARE IN VITRO INHIB PSTHS WITH IN VIVO

[vals,inds] = sort(sorter,'ascend');
maxVal = 7;
stimVITROc = mean(allInVitro(inds(1:maxVal),5000-31:5000+54));
stimVITROg = mean(allInVitroG(inds(18:35),5000-31:5000+54));
stimVIVOi = mean(lightALIGN.psthZ(150-25:150+60,inhibLogic)');
stimVIVOe = mean(lightALIGN.psthZ(150-25:150+60,directEx)');

stimVITROcE = std(allInVitro(inds(1:maxVal),5000-31:5000+54)) ./ sqrt(maxVal);
stimVITROgE = std(allInVitroG(inds(18:35),5000-31:5000+54)) ./ sqrt(maxVal);
stimVIVOiE  = std(lightALIGN.psthZ(150-25:150+60,inhibLogic)') ./ sqrt(numel(inhibLogic));
stimVIVOeE  = std(lightALIGN.psthZ(150-25:150+60,directEx)') ./ sqrt(numel(directEx));

forIgor = [stimVITROc'./max(stimVITROc),stimVITROg'./max(stimVITROg),stimVIVOi'./max(stimVIVOi),stimVIVOe'./max(stimVIVOe),stimVITROcE'./max(stimVITROc),stimVITROgE'./max(stimVITROg),stimVIVOiE'./max(stimVIVOi),stimVIVOeE'./max(stimVIVOe)];

figure(100);
clf;
subplot(212); plot(-25:60,stimVITROc./max(stimVITROc),'k',-25:60,stimVIVOi./max(stimVIVOi),'b','LineWidth',2);
subplot(211); plot(-25:60,stimVITROg./max(stimVITROg),'k',-25:60,stimVIVOe./max(stimVIVOe),'b','LineWidth',2);

%% IN VIVO SI PROBE >>>> LOAD ALL DATA FROM SORTED NEV FILES

numChannels = 64;
count = 1;

eval('home');

% Get names of all files
fileStruct = dir('*.nev');
numFiles = size(fileStruct,1);

disp('___________________________________________________________');
disp('___________________________________________________________');
disp(['Begin extracting data for the IN VIVO SI PROBE PLUS STIM dataset.']);
disp(' ');    disp(' ');


% index through files 
for i=1:numFiles
    
    nevData = TNC_LoadData(0, 0, fileStruct(i).name);
    count = 1;
    [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,'NN_w64');
    PopData.session(i).name         = fileStruct(i).name;
    PopData.session(i).trode.elMat  = electrodeMatrix;
    PopData.session(i).trode.rSp    = rSp;
    PopData.session(i).trode.cSp    = cSp;    
    
    for j=1:numChannels
        
        for k = 1:6
        
            allTs = find(nevData.Data.Spikes.Electrode==j & nevData.Data.Spikes.Unit==k);
            
            if numel(allTs)>0

                
                
                % index through channels looking for sorted units
                PopData.session(i).unit(count).ts   = round(double(nevData.Data.Spikes.Timestamps(allTs))./30);
                PopData.session(i).unit(count).el   = j;
                PopData.session(i).unit(count).un   = k;
                
                % Get the electrode map
                [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(j,'NN_w64');
                PopData.session(i).unit(count).row = row;
                PopData.session(i).unit(count).col = col;
                
                count = count+1;
                
            end
        end
    end
    count = count-1;
    
    % get stim times
    tempSamps = double(nevData.Data.Spikes.Timestamps(find(nevData.Data.Spikes.Electrode==144 & nevData.Data.Spikes.Unit==1)));
    numStim = numel(tempSamps);
    
    % get stim durations
    endSamps = double(nevData.Data.Spikes.Timestamps(find(nevData.Data.Spikes.Electrode==144 & nevData.Data.Spikes.Unit==0)));

    for p=1:numStim
        
        currSamp    = tempSamps(p);
        currEnd     = endSamps(find(endSamps>currSamp,1));
    
        PopData.session(i).events.ls.dur(p) = floor((currEnd-currSamp) ./ 30);
        if PopData.session(i).events.ls.dur(p) > 75
            PopData.session(i).events.ls.dur(p) = 1;            
        end

    end
    
    PopData.session(i).events.ls.ts = round(tempSamps./30);

    disp(['File: ' fileStruct(i).name ' completed.']);
    disp(['Completed session: ' num2str(i) ' of ' num2str(numFiles) ' ... units: ' num2str(count) ' | stimuli: ' num2str(numStim)])
    disp(['Stimuli used: ' num2str(unique(PopData.session(i).events.ls.dur))])
    disp(' ');    disp(' ');

end

% save popdata structure

for i=1:numFiles
    
    numUnits = size(PopData.session(i).unit,2);
    elList = [];
    
    for j=1:numUnits
        elList = [elList PopData.session(i).unit(j).el];
    end

    elList

end

%% IN VIVO SI PROBE >>>> STIMULUS PETHS FOR IN VIVO SI PROBE LIGHT STIM

eval('home');
disp('___________________________________________________');
disp('Get aligned raster plots and update the PopData structure...');

PopData.currParams = currParams;

countUnits  = 0; numPairwise=0;
NumSessions = size(PopData.session,2);

possStim = [1 2 5 10 20 50]
numStim = numel(possStim);

for i=1:NumSessions
    
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    disp(['Begin session: ' num2str(i) ' of ' num2str(NumSessions)])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');

    numUnits    = size(PopData.session(i).unit,2);
    countUnits  = countUnits+numUnits;

    for p=1:numStim
        
        allStimInds = find(PopData.session(i).events.ls.dur == possStim(p));                
        PopData.session(i).stim(p).numTrials = numel(allStimInds);

        if numel(allStimInds)==0
            
            PopData.session(i).stim(p).ts = [];          

        else
            
            % check to be sure that there were spikes before and after the
            % stimulus.... (some units may have drifted away)
            
            PopData.session(i).stim(p).ts = PopData.session(i).events.ls.ts(allStimInds);
            
            for j = 1:numUnits

                numStamps   = length(PopData.session(i).unit(j).ts);              
                
                if PopData.session(i).unit(j).ts(1) < PopData.session(i).stim(p).ts(1) & PopData.session(i).unit(j).ts(numStamps)>PopData.session(i).stim(p).ts(numel(allStimInds))
                
                    PopData.session(i).unit(j).stim(p).valid = 1;

                    delta       = zeros(1,ceil(PopData.session(i).unit(j).ts(numStamps)));
                    delta(PopData.session(i).unit(j).ts+1) = 1;
                    tmpSmooth   = conv(delta,currParams.filter.kernel,'same');

                    % for reference
                    % [response] = TNC_AlignRasters(delta,spkStamps,stableTrials,alignStamps,window,rasterFlag,boxcar)
                    [rLITE] = TNC_AlignRasters(delta , PopData.session(i).unit(j).ts , -1 , PopData.session(i).stim(p).ts , [PopData.currParams.winParams.prior,PopData.currParams.winParams.after],1,1);
                    PopData.session(i).unit(j).stim(p).raster     = rLITE.raster;
                    PopData.session(i).unit(j).stim(p).boxcar     = rLITE.image.boxcar;

                    [rsLITE] = TNC_AlignRasters(tmpSmooth , PopData.session(i).unit(j).ts , -1 , PopData.session(i).stim(p).ts - (numel(currParams.filter.kernel)./2), [PopData.currParams.winParams.prior ,PopData.currParams.winParams.after],0,1);
                    PopData.session(i).unit(j).stim(p).psthAVG    = rsLITE.image.psthAVG;
                    PopData.session(i).unit(j).stim(p).psthSEM    = rsLITE.image.psthSEM;
                    PopData.session(i).unit(j).stim(p).psthZ      = rsLITE.image.psthZ;
                    PopData.session(i).unit(j).stim(p).psthZe     = rsLITE.image.psthZe;
                    
                else
                    
                    disp('Found a unit that was not stable over stim period.');
                    PopData.session(i).unit(j).stim(p).valid = 0;

                end

            end 
            
            disp(['Completed all units for stim: ' num2str(p) ' of ' num2str(numStim)]);
            
        end   
    end
    
    disp(' ');
    
end

disp(['___________________________________________________'])
disp(['...COMPLETE...']);

%% IN VIVO SI PROBE >>>> GET BASELINE FIRING RATES

eval('home');
disp('___________________________________________________');
disp('Get background firing rates...');

PopData.currParams = currParams;

countUnits  = 0; 
count = 1;
NumSessions = size(PopData.session,2);

possStim = [1 2 5 10 20 50]
numStim = numel(possStim);
allRates = [];            
allRates2 = [];            

for i=1:NumSessions
    
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    disp(['Begin session: ' num2str(i) ' of ' num2str(NumSessions)])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');

    numUnits    = size(PopData.session(i).unit,2);
    countUnits  = countUnits+numUnits;
    numStamps   = length(PopData.session(i).events.ls.ts);
            
    for j = 1:numUnits

        allISIs = [];

        for q=1:numStamps

            thisStamp   = PopData.session(i).events.ls.ts(q);
            candStamps  = find(PopData.session(i).unit(j).ts>thisStamp-750 & PopData.session(i).unit(j).ts<thisStamp-50);
            allISIs     = [allISIs diff(PopData.session(i).unit(j).ts(candStamps))'];

        end

        PopData.session(i).unit(j).meanRate = mean(1./allISIs);
        allRates(count,1) = mean(1./allISIs);
        allRates2(count,1) = mean(1./diff(PopData.session(i).unit(j).ts));
        count = count+1;

    end
    
    disp(' ');
    
end

disp(['___________________________________________________'])
disp(['...COMPLETE...']);

%% IN VIVO SI PROBE >>>> SORT LIGHT PETHS BY ELECTRODE LOCATION AND STIM DURATION
clear lightALIGN
countUnits = 0;
NumSessions = size(PopData.session,2)

for p = 1:numel(possStim)
    
    k=1;

    for i=1:NumSessions
        
        if PopData.session(i).stim(p).numTrials > 0 % no stim of this duration presented

            numUnits = size(PopData.session(i).unit,2);
            
            for j = 1:numUnits

                if PopData.session(i).unit(j).stim(p).valid==1
                    lightALIGN.stim(p).dur            = possStim(p);

                    lightALIGN.stim(p).psthA(:,k)     = (PopData.session(i).unit(j).stim(p).psthAVG'-mean(PopData.session(i).unit(j).stim(p).psthAVG(1:currParams.winParams.prior))).*1000;
                    lightALIGN.stim(p).psthAe(:,k)    = (PopData.session(i).unit(j).stim(p).psthSEM').*1000;
                    lightALIGN.stim(p).psthZ(:,k)     = PopData.session(i).unit(j).stim(p).psthZ';
                    lightALIGN.stim(p).psthZe(:,k)    = PopData.session(i).unit(j).stim(p).psthZe';
                    lightALIGN.stim(p).boxcar(:,k)    = PopData.session(i).unit(j).stim(p).boxcar'-mean(PopData.session(i).unit(j).stim(p).boxcar(1:80));
                    lightALIGN.stim(p).sess(1,k)      = i;
                    lightALIGN.stim(p).unit(1,k)      = j;
                    lightALIGN.stim(p).row(1,k)       = PopData.session(i).unit(j).row;
                    lightALIGN.stim(p).col(1,k)       = PopData.session(i).unit(j).col;

                    k=k+1;
                end
            end
        end
    end
end

%% IN VIVO SI PROBE >>>> SORT RESPONSES BY MAGNITUDE

for p=1:6
    
    numUnits = size(lightALIGN.stim(p).psthZ,2);
    
    for k=1:numUnits
        if lightALIGN.stim(p).dur>10
            lightALIGN.stim(p).respSyn(k) = trapz(lightALIGN.stim(p).psthZ(currParams.winParams.prior+5:currParams.winParams.prior+lightALIGN.stim(p).dur+16,k));
        else
            lightALIGN.stim(p).respSyn(k) = trapz(lightALIGN.stim(p).psthZ(currParams.winParams.prior+5:currParams.winParams.prior+16,k));
        end
        lightALIGN.stim(p).respDir(k) = trapz(lightALIGN.stim(p).psthZ(currParams.winParams.prior-5:currParams.winParams.prior+lightALIGN.stim(p).dur,k));
    end
    
    
%     figure(1); subplot(6,1,p); hist(lightALIGN.stim(p).respDir,-40:1:40);
    
    [y,lightALIGN.stim(p).sortDir] = sort(lightALIGN.stim(p).respDir,'ascend');
    [y,lightALIGN.stim(p).sortSyn] = sort(lightALIGN.stim(p).respSyn,'ascend');
    directEx = find(lightALIGN.stim(p).respDir>2.5);
    directInhib = find(lightALIGN.stim(p).respSyn<-1);
    quartInhib = 1:round(numUnits./4);
    
    [mapName] = TNC_CreateRBColormap(500,'rb');
    
    figure(200);
    if p==1
        clf;
    end
    subplot(1,6,p); imagesc(-500:500,1:numUnits,lightALIGN.stim(p).psthZ(:,lightALIGN.stim(p).sortSyn)',[-1.5 1.5]);
    set(gca,'TickDir','out'); box off;
    colormap(mapName);
    axis([-100 100 1 numUnits]);
    
    figure(201);
    if p==1
        clf;
    end
    subplot(1,6,p); plot(-500:500,mean(lightALIGN.stim(p).psthZ(:,lightALIGN.stim(p).sortSyn(quartInhib)),2),'Color',[0 0.67 1]);
    hold on;
    subplot(1,6,p); plot(-500:500,mean(lightALIGN.stim(p).psthZ(:,directEx),2),'Color',[1 0 0]);
    axis([-150 250 -1.2 1]); box off; 
    title([num2str(numel(quartInhib)) ' | ' num2str(numel(directEx))]);
    if p==1
        ylabel('Z-score');
        xlabel('Time (ms)');
    end
    drawnow;

    figure(202);
    if p==1
        clf;
    end
    if p>4
        subplot(1,2,p-4); plot([1:8 1:8 1:8 1:8 1:8 1:8 1:8 1:8],[ones(1,8) ones(1,8).*2 ones(1,8).*3 ones(1,8).*4 ones(1,8).*5 ones(1,8).*6 ones(1,8).*7 ones(1,8).*8],'k.','Color',[0.5 0.5 0.5]); hold on;
        subplot(1,2,p-4); scatter(lightALIGN.stim(p).row(directInhib),lightALIGN.stim(p).col(directInhib),(lightALIGN.stim(p).respSyn(directInhib).*-40),[0 0.67 1],'LineWidth',2); hold on;
        subplot(1,2,p-4); scatter(lightALIGN.stim(p).row(directEx),lightALIGN.stim(p).col(directEx),(lightALIGN.stim(p).respDir(directEx).*40),[1 0 0],'^','LineWidth',2);
        axis([1 8 1 8]); axis off;
    end
    
    % get response magnitude as a function of distance (two ways: cum hist and image);
    if p>1
        SmoothedMaps.inhib.stim(p).map = zeros(8,8);
        SmoothedMaps.inhib.stim(p).cnt = ones(8,8);
        SmoothedMaps.excite.stim(p).map = zeros(8,8);
        SmoothedMaps.excite.stim(p).cnt = ones(8,8);
        for pp = 1:numel(directInhib)
            qq = directInhib(pp);
            if SmoothedMaps.inhib.stim(p).map(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) > lightALIGN.stim(p).respSyn(qq) 
%             SmoothedMaps.inhib.stim(p).map(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) = SmoothedMaps.inhib.stim(p).map(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) + lightALIGN.stim(p).respSyn(qq);
                SmoothedMaps.inhib.stim(p).map(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) = lightALIGN.stim(p).respSyn(qq);
%             SmoothedMaps.inhib.stim(p).cnt(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) = SmoothedMaps.inhib.stim(p).cnt(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) + 1;
            end
        end

        for pp = 1:numel(directEx)
            qq = directEx(pp);
            if SmoothedMaps.excite.stim(p).map(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) < lightALIGN.stim(p).respDir(qq)
%                 SmoothedMaps.excite.stim(p).map(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) = SmoothedMaps.excite.stim(p).map(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) + lightALIGN.stim(p).respDir(qq);
                SmoothedMaps.excite.stim(p).map(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) = lightALIGN.stim(p).respDir(qq);
%             SmoothedMaps.excite.stim(p).cnt(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) = SmoothedMaps.excite.stim(p).cnt(lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)) + 1;
            end
        end

        figure(210);
        subplot(2,5,p+4);        
        imagesc(SmoothedMaps.inhib.stim(p).map./SmoothedMaps.inhib.stim(p).cnt,[min(min(SmoothedMaps.inhib.stim(p).map./SmoothedMaps.inhib.stim(p).cnt)) -min(min(SmoothedMaps.inhib.stim(p).map./SmoothedMaps.inhib.stim(p).cnt))]);
        colormap(mapName); title(num2str(min(min(SmoothedMaps.inhib.stim(p).map./SmoothedMaps.inhib.stim(p).cnt))));
        subplot(2,5,p-1);        
        imagesc(SmoothedMaps.excite.stim(p).map./SmoothedMaps.excite.stim(p).cnt,[-max(max(SmoothedMaps.excite.stim(p).map./SmoothedMaps.excite.stim(p).cnt)) max(max(SmoothedMaps.excite.stim(p).map./SmoothedMaps.excite.stim(p).cnt))]);
        colormap(mapName); title(num2str(max(max(SmoothedMaps.excite.stim(p).map./SmoothedMaps.excite.stim(p).cnt))));

        [inCnt.row,inCnt.col] = find(  SmoothedMaps.inhib.stim(p).map == min(min(SmoothedMaps.inhib.stim(p).map)) );
        [exCnt.row,exCnt.col] = find(  SmoothedMaps.excite.stim(p).map == max(max(SmoothedMaps.excite.stim(p).map)) );
        
        for pp = 1:numel(directEx)
            qq = directEx(pp);
%             SmoothedMaps.excite.stim(p).dis(pp) = pdist2([lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)],[exCnt.row,exCnt.col]) .* 200;
            SmoothedMaps.excite.stim(p).dis(pp) = pdist2([lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)],[4,4]) .* 200;
            SmoothedMaps.excite.stim(p).mag(pp) = lightALIGN.stim(p).respDir(qq);
        end
        
        for pp = 1:numel(directInhib)
            qq = directInhib(pp);
%             SmoothedMaps.inhib.stim(p).dis(pp) = pdist2([lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)],[inCnt.row,inCnt.col]) .* 200;
            SmoothedMaps.inhib.stim(p).dis(pp) = pdist2([lightALIGN.stim(p).row(qq),lightALIGN.stim(p).col(qq)],[4,4]) .* 200;
            SmoothedMaps.inhib.stim(p).mag(pp) = lightALIGN.stim(p).respSyn(qq);
        end
        
        [values,indices] = sort(SmoothedMaps.inhib.stim(p).dis,'ascend');
        [valuesE,indicesE] = sort(SmoothedMaps.excite.stim(p).dis,'ascend');
        
        figure(211);
        if p==2
            clf;
        end
        subplot(2,5,p-1);
%         plot(values,cumsum(-SmoothedMaps.inhib.stim(p).mag(indices)),'b','LineWidth',2);
        plot(values,cumsum(SmoothedMaps.inhib.stim(p).mag(indices))./sum(SmoothedMaps.inhib.stim(p).mag(indices)),'b','LineWidth',2);
        axis([0 1200 0 1]); hold on;
%         plot(valuesE,cumsum(SmoothedMaps.excite.stim(p).mag(indicesE)),'r','LineWidth',2);
        plot(valuesE,cumsum(SmoothedMaps.excite.stim(p).mag(indicesE))./sum(SmoothedMaps.excite.stim(p).mag(indicesE)),'r','LineWidth',2);        
        axis([0 1200 0 1]);
        xlabel('Distance (um)'); ylabel('Cumulative Response'); title(num2str(lightALIGN.stim(p).dur));

    end
    
    inMag(p,1) = possStim(p);
    if numel(directInhib)>0
        inMag(p,2) = mean(lightALIGN.stim(p).respSyn(lightALIGN.stim(p).sortSyn(quartInhib)));
        inMag(p,3) = std(lightALIGN.stim(p).respSyn(lightALIGN.stim(p).sortSyn(quartInhib))) ./ sqrt(numel(quartInhib));
    else
        inMag(p,2) = 0;
        inMag(p,3) = 0;
    end
    
    deMag(p,1) = possStim(p);
    deMag(p,2) = mean(lightALIGN.stim(p).respDir(directEx));
    deMag(p,3) = std(lightALIGN.stim(p).respDir(directEx)) ./ sqrt(numel(directEx));

    exMag(p,1) = possStim(p);
    exMag(p,2) = mean(lightALIGN.stim(p).respDir(directEx)+ lightALIGN.stim(p).respSyn(directEx));
    exMag(p,3) = std(lightALIGN.stim(p).respDir(directEx)+lightALIGN.stim(p).respSyn(directEx)) ./ sqrt(numel(directEx));

end

%% IN VIVO SI PROBE >>>> EXPORT RASTER DATA FOR STIM CONDITIONS
sN = 6
uN = 1
clear psthTMP*

numStim = size(PopData.session(sN).unit(uN).stim,2);
x = [];
y = [];
s = [];
count=0;

for i=2:numStim
    numTrials = size(PopData.session(sN).unit(uN).stim(i).raster.trial,2)
    for j=1:numTrials
        x = [ x ; PopData.session(sN).unit(uN).stim(i).raster.trial(j).ts ];
        y = [ y ; count .* ones(numel(PopData.session(sN).unit(uN).stim(i).raster.trial(j).ts),1) ];
        s = [ s ; i     .* ones(numel(PopData.session(sN).unit(uN).stim(i).raster.trial(j).ts),1) ];        
        count = count+1;
    end
    psthTMP1(:,i) = PopData.session(sN).unit(uN).stim(i).boxcar';
end

exportToIgor = [x y s];
exportToIgor2 = [psthTMP1 [-500:6:494]'];

save forIgor exportToIgor -ascii
save forIgor2 exportToIgor2 -ascii

%% COMAPRE PSTHS BEFORE AND AFTER CNO

% unitElList = [1,2,6,7,10,39,47,48,49,64]
unitElList = [1,2,6,64]

figure(1); clf;

whichStim = 5;

for i=1:3 % look at control, and 10,20,30 min after CNO
   
    switch i
        case 1
            colorspec = [1 0 0];
        case 2
            colorspec = [0 0.67 1];
        case 3
            colorspec = [0 0.45 .67];
        case 4
            colorspec = [0 0.22 0.33];
    end
    
    for j = 1:numel(unitElList)
        
        % get the correct index
        numUnits = numel(PopData.session(i).unit);
        for k=1:numUnits
            if unitElList(j)==PopData.session(i).unit(k).el
                uI = k;
            end
        end
        
        figure(1);
        subplot(2,2,j); 
        title(PopData.session(i).unit(uI).el);
        plot(-currParams.winParams.prior:currParams.winParams.after,PopData.session(i).unit(uI).stim(whichStim).psthAVG.*1000,'Color',colorspec,'LineWidth',2); hold on;
        ylabel('Firing Rate (Hz)');
        xlabel('Time from Stimulus (ms)');
    end    
    
end

%% OLD PLOTTING ROUTINES
% disp(['Total units: ' num2str(k) ' | ' num2str(countUnits)])
% 
% numUnits = countUnits;
% 
% [mapName] = TNC_CreateRBColormap(1024,'rb');
% 
% figure(6);
% subplot(4,2,[1,3,5]);
% [mapName] = TNC_CreateRBColormap(1024,'rb');
% imagesc(lightALIGN.boxcar',[-50 50]);
% 
% ylabel('Cell index');
% title('LITE Aligned PSTH');
% 
% subplot(4,2,[2,4,6]);
% imagesc(toneALIGN.boxcar',[0 100]);
% colormap(mapName);
% ylabel('Cell index');
% title('TONE Aligned PSTH');
% 
% subplot(4,2,7);
% plot(mean(lightALIGN.boxcar,2),'k');
% colormap(mapName);
% axis([0 91 -10 40]);
% xlabel('Window Number (5 ms)');
% ylabel('Response (spk / 5ms)');
% title('LITE Aligned PSTH');
% 
% subplot(4,2,8);
% bar(mean(toneALIGN.boxcar,2),'k');
% axis([0 91 0 40]);
% xlabel('Window Number (5 ms)');
% ylabel('Response (spk / 5ms)');
% title('TONE Aligned PSTH');
% 
% figure(7);
% subplot(4,2,[1,3,5]);
% imagesc(lightALIGN.psthZ',[-10 10]);
% colormap(mapName);
% ylabel('Cell index');
% title('LITE Aligned PSTH');
% 
% subplot(4,2,[2,4,6]);
% imagesc(toneALIGN.psthZ',[-10 10]);
% colormap(mapName);
% ylabel('Cell index');
% title('TONE Aligned PSTH');
% 
% subplot(4,2,7);
% bar(mean(lightALIGN.psthZ,2),'k');
% axis([0 550 -5 15]);
% xlabel('Time (ms)');
% ylabel('Z-score');
% title('LITE Aligned PSTH');
% 
% subplot(4,2,8); 
% bar(mean(toneALIGN.psthZ,2),'k'); hold on;
% plot(mean(lightALIGN.psthZ,2)./2.9,'r');
% 
% axis([0 550 -5 15]);
% xlabel('Time (ms)');
% ylabel('Z-score');
% title('TONE Aligned PSTH');
