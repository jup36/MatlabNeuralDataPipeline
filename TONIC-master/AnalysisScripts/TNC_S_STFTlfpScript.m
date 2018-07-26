%% WALK THROUGH DIRECTORY STRUCTURE

paramsC.tapers      = [3 5];
paramsC.pad         = 0;
paramsC.Fs          = 1000; % reflects the fact that the data was decimated down to 1 kHz sampling frequency
paramsC.fpass       = [2 150];
paramsC.err         = 0;
movingwin           = [0.3 0.01];

ContData.paramsC    = paramsC;
ContData.movingwin  = movingwin;

filenamestr = 'WT1RD1LSTR004_2012_01_06.ns5';
% filenamestr = 'MP2RD2LSTR002_2011_12_02.ns5';

elecType = 'short';
justMiddle = 1;     % just use the middle electrodes of the array for LFP analysis
onlyOut = 1;        % select for only outward movements

disp(' ');
disp(' ');
disp('Initialized the parameters for the analysis of LFP data...');
disp(' ');
disp(' ');

%% FILTER THIS SESSIONS CONTINUOUS DATA
% NOTE: THIS STORES CONTINUOUS POWER SPEC DATA AT 50 ms RESOLUTION

justMiddle = 1;

switch elecType
    
    case 'short'
        if justMiddle==1
            chanList = [   20 12 5 61 58 50   ];            
        else
            chanList = [28 20 12 5 61 58 50 42];
        end
        ContData.chanList   = chanList;
        
    case 'long'
        
        if justMiddle==1
            chanList = [   24 16 8 33 64 54   ];     
        else
            chanList = [27 24 16 8 33 64 54 46];
        end
        ContData.chanList   = chanList;
        
end

for j=1:numel(chanList)
    ContData.shank(j).repChannel = chanList(j);
end

% cycle through all channels and filter/calculate power within each band
for j=1:numel(chanList)

    tstart = tic;

    i = chanList(j);
    
    electrode = ['e:' num2str(i)];
    Ns5DATA = openNSx('report','read',filenamestr, electrode);    
    
    rawData = decimate(Ns5DATA.Data(1,:),30);    
    
    % Filter the data to select for spikes
    disp(' ');
    disp(' ');
    disp(['Filtering data on channel ' num2str(i) ' ...']);

    % design the bandpass filter to use for field potentials
    dLo = fdesign.bandpass('n,fc1,fc2', 2000, 1, 150, 1000);
    HdLo = design(dLo);
    lowBandData.lowCutOff = 1;
    lowBandData.hiCutOff = 150;    
    lowBandDataTMP = filtfilt(HdLo.Numerator,1,rawData); % zero-phase filtering

    decimateCheck = numel(Ns5DATA.Data(1,:)) / numel(lowBandData.values);
    disp(['Number of samples in lowBand: ' num2str(numel(lowBandDataTMP)) ' ... from ' num2str(numel(Ns5DATA.Data(1,:))) ' original samples. Downsampled by ' num2str(decimateCheck)]);
    
    finalData = rmlinesmovingwinc(lowBandDataTMP,movingwin,10,paramsC,0.05,'n',60); % filter out 60 Hz noise
    lowBandData.values=finalData';

    % Apply chronux multitaper method to extract power spectrum in time
    disp('Calculating the spectrogram of the lowBand data...');

    ContData.lfp.params      = paramsC;
    ContData.lfp.movingwin   = movingwin;
    [S,t,f] = mtspecgramc(lowBandData.values,movingwin,paramsC);

    ContData.shank(j).contData.t = t;
    ContData.shank(j).contData.f = f;
    ContData.shank(j).contData.S = S;
    
    ContData.bands(1).name = '2to5';
    ContData.bands(2).name = '5to12';
    ContData.bands(3).name = '14to34';
    ContData.bands(4).name = '40to55';
    ContData.bands(5).name = '65to130';

    ContData.bands(1).inds = find(f<2,1,'last')  : find(f>5,1,'first');
    ContData.bands(2).inds = find(f<5,1,'last')  : find(f>12,1,'first');
    ContData.bands(3).inds = find(f<14,1,'last') : find(f>34,1,'first');
    ContData.bands(4).inds = find(f<40,1,'last') : find(f>55,1,'first');
    ContData.bands(5).inds = find(f<65,1,'last') : find(f>130,1,'first');


    for k=1:5
        ContData.shank(j).bands(k).values = mean(ContData.shank(j).contData.S(:, ContData.bands(k).inds)');
    end
    
    telapsed = toc(tstart);
    disp(['Filtered and extracted power spec from one channel in ' num2str(telapsed) ' seconds.']);
    
    figure(1);
    subplot(numel(chanList),1,j);
    hold off;
    plot(f,std(S)./mean(S)); hold on;
    peak = max(std(S)./mean(S));
    plot([2 2],[0 peak],'r--');
    plot([5 5],[0 peak],'r--');
    plot([14 14],[0 peak],'r--');
    plot([34 34],[0 peak],'r--');
    plot([40 40],[0 peak],'r--');
    plot([55 55],[0 peak],'r--');
    plot([65 65],[0 peak],'r--');
    plot([130 130],[0 peak],'r--');
    xlabel('Frequency (Hz)');
    ylabel('CV');
    title('Frequency bands selected for analysis'); 
    drawnow;
    
end
  
%% LOAD EVENT DATA
% NOTE: THIS STORES EVENT DATA AT 1 MS RESOLUTION
filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'nev']
dataEvents = openNEV(filenamestrE,'read','nomat','nosave');

nevRes = dataEvents.MetaTags.SampleRes;

rewardInds = find(dataEvents.Data.Spikes.Electrode==142); % reward
ContData.behavior.rewardInds = round(dataEvents.Data.Spikes.Timestamps(rewardInds)./30);

threshInds = find(dataEvents.Data.Spikes.Electrode==143); % threshold crossing
ContData.behavior.threshInds = round(dataEvents.Data.Spikes.Timestamps(threshInds)./30);

lickInds = find(dataEvents.Data.Spikes.Electrode==139); % licking
ContData.behavior.lickInds = round(dataEvents.Data.Spikes.Timestamps(lickInds)./30);

%% LOAD CONTINUOUS LEVER DATA
% NOTE: THIS STORES LEVER DATA AT 1 MS RESOLUTION

filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'ns4']
Ns4DATA = openNSx('report','read',filenamestrE);
clear leverData sLeverData tmpLeverData 

leverData(1,:) = decimate(Ns4DATA.Data(1,:),10);
leverData(2,:) = decimate(Ns4DATA.Data(2,:),10);

sLeverData(1,:) = sgolayfilt(leverData(1,:),5,101);
sLeverData(2,:) = sgolayfilt(leverData(2,:),5,101);

ContData.behavior.sLeverData = sLeverData;

if strcmp(filenamestr(1),'M')
    figure(55); clf;
else
    figure(5); clf;
end
plotArray = 1:250000;
plot3(sLeverData(1,plotArray),sLeverData(2,plotArray),plotArray,'k');
    
rawLick = decimate(Ns4DATA.Data(3,:),10);
difLick = diff(rawLick);
evLick=zeros(1,numel(rawLick));

ContData.behavior.evLick    = evLick;
ContData.behavior.rawLick   = rawLick;

% figure(6); clf;
% plot(rawLick-mean(rawLick),'k'); hold on; plot(difLick,'b');
% tmp = find(difLick<-400);
% evLick(tmp) = difLick(tmp);
% plot(evLick,'r','LineWidth',2);

%% EXTRACT VELOCITY DATA FROM CONT LEVER DATA
% NOTE: THIS STORES LEVER (VELOCITY) DATA AT 1 MS RESOLUTION

numSamples = size(ContData.behavior.sLeverData,2);
% numSamples = 1000;
tmpLeverData(1,:) = sgolayfilt(ContData.behavior.sLeverData(1,:),3,21);
tmpLeverData(2,:) = sgolayfilt(ContData.behavior.sLeverData(2,:),3,21);

% tmpLeverData(1,:) = medfilt1(ContData.behavior.sLeverData(1,:),91);
% tmpLeverData(2,:) = medfilt1(ContData.behavior.sLeverData(2,:),91);

% ContData.behavior.sLeverV = zeros(1,numSamples);
sLeverV = zeros(1,numSamples);

disp(' ');disp(' ');disp('Extracting velocity...');

dX = diff(tmpLeverData(1,:));
dY = diff(tmpLeverData(2,:));
    
sLeverV = sqrt( dX.^2 + dY.^2 );

disp(' ');disp(' Complete. ');disp(' ');

ContData.behavior.sLeverV = sLeverV;
ContData.behavior.sLeverVm = medfilt1(sLeverV,11);
clear sLeverV;

figure(2); clf;
subplot(211);
plot(ContData.behavior.sLeverData(1,:),'k'); hold on;
subplot(212);
plot(ContData.behavior.sLeverV,'k'); hold on;
plot(ContData.behavior.sLeverVm,'r');

%% FIND MVMTS
% NOTE: THIS STORES REACH EVENT DATA AT 1 MS RESOLUTION

% method 1 is based upon 'center changes'
% method 2 is based upon thresholding the velocities

clear reach progS*

method = 'vel';
numC = 5;

currX       = ContData.behavior.sLeverData(1,:);
currY       = ContData.behavior.sLeverData(2,:);
currV       = ContData.behavior.sLeverVm;

switch method

    case 'vel'
        
        % threshold the velocities
        allValidSamps   = find(currV>4);
        
        pre     = 10;
        post    = 10;       
        minSpace = 100;
        count = 1;
        
%         progStartTMP(count,1)  = currStamps(allValidSamps(1));
        progStartTMP(count,2)  = currX(allValidSamps(1));
        progStartTMP(count,3)  = currY(allValidSamps(1));
        progStartTMP(count,4)  = allValidSamps(1);
        
        for j=2:numel(allValidSamps)
            if allValidSamps(j)>allValidSamps(j-1)+minSpace
%                 progStopTMP(count,1)   = currStamps(allValidSamps(j-1)+post);
                progStopTMP(count,2)   = currX(allValidSamps(j-1)+post);
                progStopTMP(count,3)   = currY(allValidSamps(j-1)+post);
                progStopTMP(count,4)   = allValidSamps(j-1)+post;

                count                  = count+1;

%                 progStartTMP(count,1)  = currStamps(allValidSamps(j)-pre);
                progStartTMP(count,2)  = currX(allValidSamps(j)-pre);
                progStartTMP(count,3)  = currY(allValidSamps(j)-pre);
                progStartTMP(count,4)  = allValidSamps(j)-pre;                
            end
            
            if j==numel(allValidSamps)
%                 progStopTMP(count,1)   = currStamps(allValidSamps(j)+post);
                progStopTMP(count,2)   = currX(allValidSamps(j)+post);
                progStopTMP(count,3)   = currY(allValidSamps(j)+post);
                progStopTMP(count,4)   = allValidSamps(j)+post;
            end
            
        end
        
        count = 1;
        for k = 1:size(progStartTMP,1)
            
            % reaches must be at least 50 ms long
            if progStopTMP(k,4)-progStartTMP(k,4)>=90            

                trajAngle   = atan2(progStopTMP(k,3)-progStartTMP(k,3),progStopTMP(k,2)-progStartTMP(k,2));
                
                if (pdist2([progStopTMP(k,2),progStopTMP(k,3)],[mean(currX),mean(currY)]) > pdist2([progStartTMP(k,2),progStartTMP(k,3)],[mean(currX),mean(currY)]))
                    reach.out(count) = 1;
                else
                    reach.out(count) = 0;
                end
                velTraj = ContData.behavior.sLeverV(progStartTMP(k,4) : progStopTMP(k,4));
                reach.start(count,:)  = progStartTMP(k,:);
                reach.stop(count,:)   = progStopTMP(k,:);
                reach.angle(count,1)  = trajAngle;
                reach.dist(count,:)   = sum(velTraj);        
                tmp = findpeaks(velTraj);
                reach.numpks(count,1) = numel(tmp.loc);
                reach.dur(count,1)    = progStopTMP(k,4) - progStartTMP(k,4);
                reach.maxV(count,1)   = max(velTraj);
                reach.maxV(count,2)   = progStartTMP(k,4) + find(velTraj==max(velTraj),1);
                reach.maxV(count,3)   = mean(velTraj);
                count                 = count+1;
            
            end            
        end

        
end

level = ones(1,size(progStartTMP,1));
figure(3); clf;
plot(ContData.behavior.sLeverV,'k'); hold on;
plot(progStartTMP(:,4),level,'r^');
plot(progStopTMP(:,4),level,'bo');

numReaches = size(reach.start,1);
if strcmp(filenamestr(1),'M')
    figure(55); clf;
else
    figure(5); clf;
end
subplot(5,1,1:4);

for l = 1:numReaches
%   offset = reach.start(l,4).*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
    offset = l.*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
    plot3(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4)), currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4)) , offset , 'r','LineWidth',1);
%   plot3(currX(reach.start(l,4):reach.stop(l,4)), currY(reach.start(l,4):reach.stop(l,4)) , offset);
    hold on;
end

plot3([0 0],[0 0],[-5 numReaches+5],'k-','LineWidth',2);
grid on; box on;
xlabel('X (a.u.)');
ylabel('Y (a.u.)');
zlabel('Time (samples)');
title(['Number of reaches: ' num2str(numReaches)]);
view([-45 10]);
axis tight;

subplot(5,1,5);
rose(reach.angle);

reach.numReaches = numReaches;

ContData.behavior.reach = reach;

%% SAVE this file

disp(' ');
disp(' ');
disp('--------------------------------------------------------------------------');
disp(['saving structure to file: ' filenamestr(1,1:length(filenamestr)-4) '.tcd']);
disp('--------------------------------------------------------------------------');
disp(' ');
disp(' ');

eval(['save ' filenamestr(1,1:length(filenamestr)-4) '.tcd ContData']);

%% DISPLAY LFP POWER SPECTRA ALIGNED TO REACHES

for j = 1:numReaches
    
    % did this reach trigger a reward?
    tmp = find(ContData.behavior.threshInds > reach.start(j,4) & ContData.behavior.threshInds < reach.stop(j,4));
    threshLogic(j) = numel(tmp);
    
end

tmp2 = find(threshLogic);
tmp3 = 1:numReaches;
tmp4 = find(reach.out);
tmp5 = find(reach.out==0)
figure(7); clf;
plot(tmp3(tmp4), reach.angle(tmp4),'ro','MarkerSize',10,'LineWidth',2); hold on;
plot(tmp3(tmp5), reach.angle(tmp5),'bo','MarkerSize',10,'LineWidth',2); hold on;
% plot(tmp3, reach.angle,'ro-','MarkerSize',10); hold on;
plot(tmp3(tmp2), reach.angle(tmp2),'ko-','MarkerSize',6,'MarkerFace',[0 0 0],'LineWidth',2); 

reach.threshLogic = threshLogic;

numReward = numel(find(threshLogic))

%% REACH TRIGGERED LFP AVERAGES

win1ms  = [-500:1:500];
win50ms = [-500./50:1:500./50];
clear reachTrig*;

% Pick a band and a shank to view
vBand   = 3;
vBand2  = 5;
vShank  = 2;
%     For reference:
%     ContData.bands(1).name = '2to5';
%     ContData.bands(2).name = '5to12';
%     ContData.bands(3).name = '18to32';
%     ContData.bands(4).name = '34to50';
%     ContData.bands(5).name = '65to120';

[vals,inds] = sort(reach.threshLogic ,'descend');
% [vals,inds] = sort(reach.dist ,'descend');
% [vals,inds] = sort(reach.maxV(:,1) ,'descend');

% for j = 1:numel(inds)
for j = 1:numReward
    
    sInd = inds(j);
    zeroInd = reach.start(sInd,4)./50;
%     zeroInd = round(reach.maxV(sInd,2)./50);

    reachTrigLFP1(j,1:71) = ContData.shank(vShank).bands(vBand).values(zeroInd-25:zeroInd+45) - mean(ContData.shank(vShank).bands(vBand).values(zeroInd-25:zeroInd-20));
    reachTrigLFP2(j,1:71) = ContData.shank(vShank).bands(vBand2).values(zeroInd-25:zeroInd+45) - mean(ContData.shank(vShank).bands(vBand2).values(zeroInd-25:zeroInd-20));

%     reachTrigVEL(j,1:3501) = ContData.behavior.sLeverV(reach.maxV(sInd,2)-1250:reach.maxV(sInd,2)+2250);
%     reachTrigPOS(j,1:3501) = abs(ContData.behavior.sLeverData(1,reach.maxV(sInd,2)-1250:reach.maxV(sInd,2)+2250) - mean(ContData.behavior.sLeverData(1,reach.maxV(sInd,2)-1250:reach.maxV(sInd,2)-1000)) ) + abs(ContData.behavior.sLeverData(2,reach.maxV(sInd,2)-1250:reach.maxV(sInd,2)+2250) - mean(ContData.behavior.sLeverData(2,reach.maxV(sInd,2)-1250:reach.maxV(sInd,2)-1000)) );

    reachTrigVEL(j,1:2501) = ContData.behavior.sLeverV(reach.start(sInd,4)-1250:reach.start(sInd,4)+1250);
    reachTrigPOS(j,1:2501) = abs(ContData.behavior.sLeverData(1,reach.start(sInd,4)-1250:reach.start(sInd,4)+1250) - mean(ContData.behavior.sLeverData(1,reach.start(sInd,4)-1250:reach.start(sInd,4)-1000)) ) + abs(ContData.behavior.sLeverData(2,reach.start(sInd,4)-1250:reach.start(sInd,4)+1250) - mean(ContData.behavior.sLeverData(2,reach.start(sInd,4)-1250:reach.start(sInd,4)-1000)) );

end

figure(9); clf;
subplot(131);
imagesc(reachTrigPOS);
subplot(132);
imagesc(reachTrigLFP1);
subplot(133);
imagesc(reachTrigLFP2);


figure(8); clf;
subplot(7,1,1);
imagesc(reachTrigVEL);
subplot(7,1,2:3);
imagesc(reachTrigPOS);
subplot(7,1,4);
imagesc(reachTrigLFP1);
subplot(7,1,5);
imagesc(reachTrigLFP2);
subplot(7,1,6:7);
plot(([1:71]-25).*50,mean(reachTrigLFP1)-mean(mean(reachTrigLFP1(:,1:10))),'LineWidth',2,'Color',[0 0 0]); hold on;
plot(([1:71]-25).*50,mean(reachTrigLFP2)-mean(mean(reachTrigLFP2(:,1:10))),'LineWidth',2,'Color',[1 0 0]); 
plot([0 0],[-20 10],'k--');
plot([-1000 2000],[0 0],'k--'); 
axis tight

%% EXTRACT SPIKE EVENTS

eventWindow = [-15:45]; % in samples (corresponds to -0.5 ms to 1.5 ms)
snrTh = [4,30];
justMiddle = 0;

switch elecType
    
    case 'short'
        if justMiddle==1
            chanList = [   20 12 5 61 58 50   ];            
        else
            chanList = [28 20 12 5 61 58 50 42];
        end
        ContData.chanList   = chanList;
        
    case 'long'
        
        if justMiddle==1
            chanList = [   24 16 8 33 64 54   ];     
        else
            chanList = [27 24 16 8 33 64 54 46];
        end
        ContData.chanList   = chanList;
        
end

for j=1:numel(chanList)
    ContData.shank(j).repChannel = chanList(j);
end

% cycle through all channels and filter/calculate power within each band
% for j=1:numel(chanList)
for j=7:7

    tstart = tic;
    i = chanList(j);
    
    electrode = ['e:' num2str(i)];
    Ns5DATA = openNSx('report','read',filenamestr, electrode);    
    
    % Filter the data to select for spikes
    disp(' ');
    disp(' ');
    disp(['Filtering data on channel ' num2str(i) ' ...']);
    
    rawData = sgolayfilt(Ns5DATA.Data(1,1:1000000),7,25);

    [lowBandData,hiBandData] = TNC_FilterData(rawData,Ns5DATA.MetaTags.SamplingFreq,0,0);        

    % do event detection
    disp('Looking for events...');
    [events] = TNC_EventDetect(hiBandData.values,30,2.5);

    % do event extraction
    disp('Extract spike waveforms...');
    [spikes] = TNC_EventExtract(hiBandData.values,rawData,events.inds,1,[12,25]);
    spikes.inds = events.inds;
    
    % pull out a bandpass filtered beta
    dLo  = fdesign.bandpass('n,fc1,fc2', 2000, 20, 29, Ns5DATA.MetaTags.SamplingFreq);
    HdLo = design(dLo);
    beta = filtfilt(HdLo.Numerator,1,rawData); %zero-phase filtering

    % quick pass event quantitation to clean out some apparent noise
    [spikes] = TNC_EventHeuristic(spikes,0.01);
    
    % quantify each spike to look for easy sorting
    [spikes] = TNC_EventQuant(spikes,'pca','dot',1);
    [spikes] = TNC_EventQuant(spikes,'scalar','dot',1);
   
    % display the MUA rate as a function of time
    delta = zeros(1,numel(rawData));
    delta(spikes.inds) = 1;
    
    currParams.smthParams.rise      = 1;
    currParams.smthParams.decay     = 75;
    [currParams.filter.kernel]      = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
    MUAdata.shank(j).endSite.smth   = conv(delta,currParams.filter.kernel,'same');

%     % get spike triggered average
%     figure(2); clf;
%     subplot(6,1,j);
%     
%     plot(MUAtrig.range,(MUAtrig.avg-mean(MUAtrig.avg))./max(MUAtrig.avg-mean(MUAtrig.avg)),'k','LineWidth',2);
%     hold on; 
%     % plot(MUAtrig.range,MUAtrig.avg-mean(MUAdata.smth)-MUAtrig.err,'k--',MUAtrig.range,MUAtrig.avg-mean(MUAdata.smth)+MUAtrig.err,'k--')
%     [MVtrig] = TNC_ExtTrigWins(ContData.behavior.sLeverV,ContData.behavior.threshInds(tmpEV),[1500,2000]);
%     plot(MVtrig.range,(MVtrig.avg-mean(MVtrig.avg))./-max((MVtrig.avg-mean(MVtrig.avg))),'r','LineWidth',2);

end

%% LOAD CONTINUOUS DATA FOR A CHOSEN SHANK

currShank = ;
eventWindow = [-15:45]; % in samples (corresponds to -0.5 ms to 1.5 ms)
snrTh = [4,30];

switch elecType
    
    case 'short'
        if justMiddle==1
            chanList = [   20 12 5 61 58 50   ];            
        else
            chanList = [28 20 12 5 61 58 50 42];
        end
        ContData.chanList   = chanList;
        
    case 'long'
        
        if justMiddle==1
            chanList = [   24 16 8 33 64 54   ];     
        else
            chanList = [27 24 16 8 33 64 54 46];
        end
        ContData.chanList   = chanList;
        
end



%% PLOT THE EVENT TRIGGERED DATA

figure(1); clf;
figure(2); clf;

for j=1:numel(chanList)
    figure(1);
    subplot(numel(chanList)+1,1,j);
        plot(MUAdata.shank(j).endSite.smth,'k');
        tmpEV = find(ContData.behavior.threshInds<1000000);
        hold on; 
        plot(ContData.behavior.threshInds(tmpEV) ,zeros(1,numel(tmpEV)),'b.');
        title(chanList(j));

    if j==1
        figure(1);
        subplot(numel(chanList)+1,1,numel(chanList)+1);
        plot(ContData.behavior.sLeverV(1,1:1000000)./-1000,'r'); 
        figure(2);
        subplot(numel(chanList)+1,1,numel(chanList)+1);
        [MVtrig] = TNC_ExtTrigWins(ContData.behavior.sLeverV,ContData.behavior.threshInds(tmpEV),[1500,2000]);
        plot(MVtrig.range,(MVtrig.avg-mean(MVtrig.avg))./-max((MVtrig.avg-mean(MVtrig.avg))),'r','LineWidth',2);
    end
    
    [MUAtrig] = TNC_ExtTrigWins(MUAdata.shank(j).endSite.smth,ContData.behavior.threshInds(tmpEV),[1500,2000]);
    figure(2);
    subplot(numel(chanList)+1,1,j);
    plot(MUAtrig.range,MUAtrig.avg-mean(MUAdata.shank(j).endSite.smth),'LineWidth',2);
    title(chanList(j));
    
end    
drawnow;