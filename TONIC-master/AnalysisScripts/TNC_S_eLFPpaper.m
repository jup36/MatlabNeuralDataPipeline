%% TNC_S_eLFPpaper

% Active directory
% /Volumes/huxley/_MANUS_ACTIVE/2013-Pan-DAFIELD/LFP/LFPdata

% Outline of figures / arguments in the paper
Spatial organization corresponds to the borders of midbrain nuclei
Timecourse and amplitude match the latencies of midbrain neuron populations
Learning related changes in LFP mirror learning related changes in single unit responses
Shared pharmacology
Trial by trial correlations - evidence for subpopulations?

%% STANDARD ANALYSIS PARAMETERS
% Define the smoothing to apply to the timeseries data
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 50;
currParams.filter.causal       = 1;

if currParams.filter.causal
    [currParams.filter.kernel]  = TNC_CreateCausalKernel(currParams.smthParams.rise,currParams.smthParams.decay,1);
else
    [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
end

currParams.winParams.prior     = 0.5e3;
currParams.winParams.after     = 4e3; 

PopData.currParams = currParams;

disp(' ');
disp(' ');
disp('________________________________');
disp(' ');
disp('Initialized analysis parameters.');
disp('________________________________');
disp(' ');
disp(' ');

%% OUTLINE OF STANDARD SESSION STRUCTURE
eLFPpopData.session(n).unit(j).ts % time stamps at 1 kHz res for all spikes
eLFPpopData.session(n).unit(j).el % electrode number on which unit was recorded
eLFPpopData.session(n).unit(j).wf % waveforms for all spikes
eLFPpopData.session(n).unit(j).da % logical value for whether called da or not

eLFPpopData.session(n).cont.trode(k).onek % continuous filtered data (bandwidth: )

eLFPpopData.session(n).learning % logical value to class the session as either learning or extinction

eLFPpopData.session(n).events.cs % time stamps for the tone
eLFPpopData.session(n).events.us % time stamps for the reward
eLFPpopData.session(n).events.fl % time stamps for the programmatically recorded first lick

eLFPpopData.session(n).behavior.trl   % time stamps for each trials
eLFPpopData.session(n).behavior.rew   % whether the trial was rewarded
eLFPpopData.session(n).behavior.lat   % time to sample the port for every trial

eLFPpopData.session(n).behavior.lick.ts   % time stamps for every lick
eLFPpopData.session(n).behavior.lick.fl   % first lick by trial
eLFPpopData.session(n).behavior.lick.ll   % last lick by trial

%% LOAD DATA (ACQ & EXT) INTO POPULATION DATA STRUCTURE
currSes = 0;
allFiles = dir('*.nex');
firstTrode  = 33
numTrodes   = 32
extForCont  = 'ns3'

for j = 1:size(allFiles,1)
    disp(allFiles(j).name);   

    currSes = currSes+1;
    disp(['Current session: ' num2str(currSes)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD SORTED SPIKE DATA FROM NEUROEXPLORER     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
    fileNameStr = [allFiles(j).name];                    
    RD_Active = TNC_LoadData(0, 0, fileNameStr);

    propUA = []; propCS = 0; propEL = 0; clear propNames; p = 1;

    for m = 1:size(RD_Active.waves,1)    
        
        if numel(findstr(RD_Active.waves{m}.name,'Chan'))>0 % electrode
            newtype=0;

            testi = findstr(RD_Active.waves{m}.name,'u');
            if size(testi,1) == 0
                propUA = [propUA,m];
                propNames.names(p).str = RD_Active.waves{m}.name;
                disp(['Found unit ' num2str(p) ': ' RD_Active.waves{m}.name]);
                p = p+1;
            end
            
        elseif numel(findstr(RD_Active.waves{m}.name,'elec'))>0

            newtype=1;

%             if numel(findstr(RD_Active.waves{m}.name,'_ts'))>0 && numel(findstr(RD_Active.waves{m}.name,'_u'))>0
            tmpInd = findstr(RD_Active.waves{m}.name,'_u');
            if numel(tmpInd)>0 && str2num(RD_Active.waves{m}.name(tmpInd+2))~=0
                propUA = [propUA,m];
                propNames.names(p).str = RD_Active.waves{m}.name;
                disp(['Found unit ' num2str(p) ': ' RD_Active.waves{m}.name]);
                p = p+1;                
            end
        end
    end

    p=p-1;

    % Classify the session
    if size(findstr(allFiles(j).name,'-ex')) > 0
        disp([fileNameStr ' is an extinction session.']);
        PopData.session(currSes).learning = 0;
    else
        disp([fileNameStr ' is a learning session.']);
        PopData.session(currSes).learning = 1;
    end

    if newtype
        % Write into the pop data structure
        for j=1:p
            i = propUA(j);
            eLFPpopData.session(currSes).unit(j).ts = RD_Active.waves{i}.timestamps.*1000;
            eLFPpopData.session(currSes).unit(j).el = str2double(RD_Active.waves{i}.name(5:6));        
            eLFPpopData.session(currSes).unit(j).wf = RD_Active.waves{i}.waveforms;
            eLFPpopData.session(currSes).unit(j).ui = str2double( RD_Active.waves{i}.name(findstr(RD_Active.waves{i}.name,'_u')+2) );
            eLFPpopData.session(currSes).unit(j).da = 0; % to be set later...

            [isi] = TNC_QuantISI(eLFPpopData.session(currSes).unit(j).ts);
            eLFPpopData.session(currSes).unit(j).isi = isi; 

            disp(['Unit ' num2str(j) ': elec' num2str(eLFPpopData.session(currSes).unit(j).el) '_unit' num2str(eLFPpopData.session(currSes).unit(j).ui)]);
        end
    else
        % Write into the pop data structure
        for j=1:p
            i = propUA(j);
            eLFPpopData.session(currSes).unit(j).ts = RD_Active.waves{i}.timestamps.*1000;
            if str2num(propNames.names(j).str(6))==0
                eLFPpopData.session(currSes).unit(j).el = str2num(propNames.names(j).str(7));        
            else
                eLFPpopData.session(currSes).unit(j).el = str2num(propNames.names(j).str(6:7));        
            end
            eLFPpopData.session(currSes).unit(j).wf = RD_Active.waves{i}.waveforms;
            eLFPpopData.session(currSes).unit(j).ui = RD_Active.waves{i}.unitNumber;
            eLFPpopData.session(currSes).unit(j).da = 0; % to be set later...

            [isi] = TNC_QuantISI(eLFPpopData.session(currSes).unit(j).ts);
            eLFPpopData.session(currSes).unit(j).isi = isi; 

            disp(propNames.names(j).str);
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD EVENT DATA FROM NEV FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fNs = [fileNameStr(1:numel(fileNameStr)-3) 'nev'];
    [NEV] = TNC_LoadData(0, 0, fNs);
    
    CSraw = NEV.Data.Spikes.Timestamps(find(NEV.Data.Spikes.Electrode==137 & NEV.Data.Spikes.Unit>0));
    eLFPpopData.session(currSes).events.cs = CSraw ./ 30;
    eLFPpopData.session(currSes).events.us = eLFPpopData.session(currSes).events.cs + 2000;
    ALraw = NEV.Data.Spikes.Timestamps(find(NEV.Data.Spikes.Electrode==138));
    eLFPpopData.session(currSes).events.al = ALraw ./ 30;
    FLraw = NEV.Data.Spikes.Timestamps(find(NEV.Data.Spikes.Electrode==139 & NEV.Data.Spikes.Unit>0));
    eLFPpopData.session(currSes).events.fl = FLraw ./ 30;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD CONTINUOUS DATA FROM NS3-5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fNsC = [fileNameStr(1:numel(fileNameStr)-3) extForCont];

    % load data
    NsDATA = openNSx(fNsC,'read','report');
    numTrodes = NsDATA.MetaTags.ChannelCount;

    for k=1:numTrodes
%         if k>firstTrode
%             % load new channel data
%             NsDATA = openNSx(fNsC,'read','report','p: double',['e: ' num2str(k)]);
%         end
        
        % decimate data
        scaling = NsDATA.MetaTags.SamplingFreq./1000;
        rawData = decimate(NsDATA.Data(k,:),scaling);
        dataPnts = numel(rawData);
        
        % filter using Savitzky-Golay smoothing
        eLFPpopData.session(currSes).cont.trode(k).onek = sgolayfilt(rawData,5,31);
        eLFPpopData.session(currSes).cont.trode(k).id   = NsDATA.MetaTags.ChannelID(k);
        eLFPpopData.session(currSes).cont.trode(k).hasDA   = 0;
        
        % show the effect of filtering for a sanity check
        figure(1); clf;
        plot(1:dataPnts,rawData,'k',1:dataPnts,eLFPpopData.session(currSes).cont.trode(k).onek,'r');
        title(['Channel ' num2str(k)]);
        drawnow;
    
    end
    
end

%% ANALYSIS OF LFP TRIAL BY TRIAL

% window          = [100 2700];
window          = [100 500];
windows(1,:)    = [5 20];
windows(2,:)    = [20 40];
windows(3,:)    = [55 150];
offset          = window(1);
showEachTrial   = 0; 
spacing         = 150;
decay           = 6;
kernel          = TNC_CreateGaussian(decay.*15,decay,decay.*30,1);
kernel(1:decay.*15) = 0;
    figure(507); plot(kernel);
% kernel          = TNC_CreateCausalKernel(0.5,decay,1);
numChannels     = 32;

% Calculate the trialwise windows of spike density functions
numUnits = size(eLFPpopData.session(currSes).unit,2);
numSamp1k = numel(eLFPpopData.session(currSes).cont.trode(1).onek);

for p=1:2

    currSes = p
    
    if currSes == 2
        eLFPpopData.session(currSes).events.cs = eLFPpopData.session(currSes).events.cs(1:60);
        disp('Using only the first 60 trials of extinction.');
    end
    
    for j=1:numUnits

        disp(['Unit ' num2str(j) '...']);
        
        delta = zeros(1,numSamp1k);
        delta(round(eLFPpopData.session(currSes).unit(j).ts+1)) = 1;
        eLFPpopData.session(currSes).unit(j).sdf = conv(delta,kernel,'same');
        [rCS] = TNC_ExtTrigWins(eLFPpopData.session(currSes).unit(j).sdf,eLFPpopData.session(currSes).events.cs,window);

        numTrials = size(rCS.wins,1);
        for i=1:numTrials            
            eLFPpopData.session(currSes).unit(j).rAmp(1,i) = trapz( rCS.wins(i,offset+[1:windows(3,2)]) - mean(rCS.wins(i,[1:windows(3,2)]+offset-100))  );
        end
        
        unitData.corrMat(currSes).matVals(:,j) = eLFPpopData.session(currSes).unit(j).rAmp';

        figure(100); clf; subplot(5,1,1:4);
        imagesc(rCS.wins); axis tight;
        [mapName] = TNC_CreateRBColormap(1024,'mbr');
        colormap(mapName);
        title(['Unit: ' num2str(j) ' ... ' num2str(eLFPpopData.session(currSes).unit(j).el)]);
        eLFPpopData.session(currSes).unit(j).rCS = rCS;
        subplot(5,1,5);
        plot(abs(mean(rCS.wins)),'LineWidth',2); axis tight; drawnow;

    end

    % Take the integral of the absolute power for each trial over a peristimulus window
    tmpChan = size(eLFPpopData.session(currSes).cont.trode,2)
    if tmpChan > 32
        numChannels = 32;
    else
        numChannels = 16;
    end

    for j=1:numChannels

        disp(['Continuous channel ' num2str(j) '...']);

        [sink] = TNC_ExtTrigWins(eLFPpopData.session(currSes).cont.trode(j).onek,eLFPpopData.session(currSes).events.cs,window);
        if p==1
            [sink2] = TNC_ExtTrigWins(eLFPpopData.session(currSes).cont.trode(j).onek,eLFPpopData.session(currSes).events.cs+2450,window);
        end
        
        figure(200);
        subplot(5,1,1:4);
        imagesc(sink.wins ./ std(mean(sink.wins)),[-6 6]); axis tight;
        [mapName] = TNC_CreateRBColormap(1024,'mbr');
        colormap(mapName); 
        ylabel('Trials');
        title(j+32);
        subplot(5,1,5);
        plot(-window(1):window(2), abs( (mean(sink.wins) - mean(mean(sink.wins))) ./ std(mean(sink.wins)) ),'k','LineWidth',2); hold off;
        axis([-window(1) window(2) -6 6]);
        xlabel('Time (ms)');
        drawnow;
        

        for i=1:numTrials
            eLFPpopData.session(currSes).cont.trode(j).rAmp(1,i) = trapz( abs(sink.wins(i,[windows(1,1):windows(1,2)]+offset)) );
            eLFPpopData.session(currSes).cont.trode(j).rAmp(2,i) = trapz( abs(sink.wins(i,[windows(2,1):windows(2,2)]+offset)) );
            eLFPpopData.session(currSes).cont.trode(j).rAmp(3,i) = trapz( abs(sink.wins(i,[windows(3,1):windows(3,2)]+offset)) );
            if p==1
                eLFPpopData.session(currSes).cont.trode(j).rAmp(4,i) = trapz(abs(sink2.wins(i,[windows(3,1):windows(3,2)]+offset)) );
            end
        end

    end
    
%     for j=1:numUnits
%         tmpCS = eLFPpopData.session(currSes).unit(j).rCS;
%         uAmp(1,i) = trapz((tmpCS.wins(i,[windows(1,1):windows(1,2)]+offset)));
%         eLFPpopData.session(currSes).unit(j).uAmp = uAmp;
%     end

end

%% List all units so that I can annotate DA subtype

    for j=1:32
        eLFPpopData.session(currSes).cont.trode(j).hasDA = 0;
        eLFPpopData.session(currSes).cont.trode(j).hasNDA = 0;
    end

    daUnits = [5 12 13 14 16 17 21 22 23 24 26 34 35 38 39 40]

    for j=1:numUnits
        disp(['Unit ' num2str(j) ': elec' num2str(eLFPpopData.session(currSes).unit(j).el) '_unit' num2str(eLFPpopData.session(currSes).unit(j).ui)]);

        eLFPpopData.session(currSes).cont.trode(eLFPpopData.session(currSes).unit(j).el-32).hasNDA   = -1;

        if numel( find(daUnits==j) )>0
            disp(['                       Unit ' num2str(j) ' is a putative dopamine unit.']);
            eLFPpopData.session(currSes).unit(j).da = 1;
            eLFPpopData.session(currSes).cont.trode(eLFPpopData.session(currSes).unit(j).el-32).hasDA   = 1;
        end
    end
        
%% Find the number of clusters

Y = 2:0.2:4;

for i = 1:numel(Y);

   for k = 1:20; %define possible sizes of k

       Color = {'m','b','g','r','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.


       %Kmeans requires rows to contain observations (so cells), columns to contain variables

       FullChar_2 = [trialWise1 ; trialWise2 ; trialWise3];

       [idx,C,sumd,D] = kmeans(FullChar_2,k,'replicates',5);


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

%% SPATIAL ORGANIZATION OF ELECTRODES


numClust = 12;
extinctExists = 1;
clustered=1;

currSes = 1;
elecMap = [  1  3  5  7  9 11 13 15 ; 
            17 19 21 23 25 27 29 31 ;
             2  4  6  8 10 12 14 16 ;
            18 20 22 24 26 28 30 32 ];
        
% colorMap                = [ 14, 83,167;
%                             92,200,200;
%                            230, 60,137;
%                            255,116,  0;
%                              2,146,146;
%                            247,254,  0;
%                            204,  3, 95;
%                            186,239,110;                               
%                            230, 60,137;
%                              0,153,153;
%                            113,  9,170;
%                            230, 60,137;
%                            134,223,  4;
%                            186,239,110;                               
%                            204,  3, 95;
%                            247,254,  0;
%                            249,175,114 ];
%                            
%                        colorMap = colorMap ./ 255;
% % [colorMap] = TNC_CreateRBColormap(numClust,'bo');
% % colorMap = colorMap ./ 1024;
% 
% colorMap = [ [0:11]' , ones(12,1).*9 , [11:-1:0]' ] ./ 12;

     colorMap   = [ 166,206,227;
                     31,120,180;
                    178,223,138;
                     51,160, 44;
                    251,154,153;
                    227, 26, 28;
                    253,191,111;
                    255,127,  0;
                    202,178,214;
                    106, 61,154;
                    255,255,153;
                    177, 89, 40] ./ 256;
xSpc = 700;
ySpc = 150;
numTrials = numel(eLFPpopData.session(currSes).cont.trode(1).rAmp(3,:));
trialWise1 = zeros(numTrials,32);
trialWise2 = zeros(numTrials,32);
trialWise3 = zeros(numTrials,32);
trialWise4 = zeros(numTrials,32);
trialWise1s = zeros(numTrials,32);
trialWise2s = zeros(numTrials,32);
trialWise3s = zeros(numTrials,32);
trialWise4s = zeros(numTrials,32);


if clustered==1 && currSes==1
    figure(400); clf;
end

clear forAveraging; count = 1;
for i=1:32
    
    if clustered==1 && currSes==1
        [r,c] = find(elecMap==i);
        
%         c=-c

        [sink1] = TNC_ExtTrigWins(eLFPpopData.session(1).cont.trode(i).onek,eLFPpopData.session(1).events.cs,window);
        [sink2] = TNC_ExtTrigWins(eLFPpopData.session(1).cont.trode(i).onek,eLFPpopData.session(1).events.cs+2450,window);
        [sink3] = TNC_ExtTrigWins(eLFPpopData.session(2).cont.trode(i).onek,eLFPpopData.session(2).events.cs,window);

        norm1   = (mean(sink1.wins(1:numTrials,:))- mean(mean(sink1.wins(1:numTrials,1:window(1)))) ) ./ max(abs(mean(sink1.wins(1:numTrials,:))));        
        norm2   = (mean(sink2.wins(1:numTrials,:))- mean(mean(sink2.wins(1:numTrials,1:window(1))))  ) ./ max(abs(mean(sink1.wins(1:numTrials,:))));        
        norm3   = (mean(sink3.wins(30:numel(eLFPpopData.session(2).events.cs),:))- mean(mean(sink3.wins(30:numel(eLFPpopData.session(2).events.cs),1:window(1))))  ) ./ max(abs(mean(sink1.wins(1:numTrials,:))));        

%         if any(i~=[18,20,22,24,26,28,30,32])
        if any(i==[1 3 5 7 9 11 13 19 21 23 25 27 29 6 8])
            forAveraging(count,:) = abs(norm1);
            count = count+1;
        end
        
        
%         plot([-window(1):window(2)]+(xSpc*c),(norm1*50)+(r*ySpc),'Color',colorMap(round(indx(i)./numClust.*1024),:),'LineWidth',2); hold on;
%         plot([-window(1):window(2)]+(xSpc*c),(norm2*50)+(r*ySpc),'Color',colorMap(round(indx(i)./numClust.*1024),:)./1.5,'LineWidth',2); hold on;
%         plot([-window(1):window(2)]+(xSpc*c),(norm3*50)+(r*ySpc),'Color',colorMap(round(indx(i)./numClust.*1024),:)./4,'LineWidth',2); hold on;

        plot([-window(1):window(2)]+(xSpc*c),(norm1*50)+(r*ySpc),'Color',colorMap(indx(i),:),'LineWidth',2); hold on;
        plot([-window(1):window(2)]+(xSpc*c),(norm2*50)+(r*ySpc),'Color',colorMap(indx(i),:)./1.5,'LineWidth',2); hold on;
        plot([-window(1):window(2)]+(xSpc*c),(norm3*50)+(r*ySpc),'Color',colorMap(indx(i),:)./4,'LineWidth',2); hold on;


        ind = 8*(r-1) + c; set(gcf,'Color','white'); set(gca,'TickDir','out','TickLength',[0.005 0]); axis off;
        
        text((xSpc*c)+185,(r*ySpc+58),[num2str(indx(i)) '.' num2str(i+32) '.' num2str(eLFPpopData.session(1).cont.trode(i).hasDA) '.' num2str(eLFPpopData.session(1).cont.trode(i).hasNDA)]);            
        
%         if corrCalc==1
%             text((xSpc*c)-100,(r*ySpc+80),[num2str(eLFPpopData.session(j).events.corrScores(i,2),'%0.1g') '...' num2str(eLFPpopData.session(j).events.corrScores(i,3),'%0.1g')]);
%         end
    end
    
    trialWise1(:,i) = (eLFPpopData.session(currSes).cont.trode(i).rAmp(1,:)' - mean(eLFPpopData.session(currSes).cont.trode(i).rAmp(1,:))) ./ std(eLFPpopData.session(currSes).cont.trode(i).rAmp(1,:));
    trialWise2(:,i) = (eLFPpopData.session(currSes).cont.trode(i).rAmp(2,:)' - mean(eLFPpopData.session(currSes).cont.trode(i).rAmp(2,:))) ./ std(eLFPpopData.session(currSes).cont.trode(i).rAmp(2,:));
    trialWise3(:,i) = (eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:)' - mean(eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:))) ./ std(eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:));
    trialWise4(:,i) = (eLFPpopData.session(currSes).cont.trode(i).rAmp(4,:)' - mean(eLFPpopData.session(currSes).cont.trode(i).rAmp(4,:))) ./ std(eLFPpopData.session(currSes).cont.trode(i).rAmp(4,:));

    if extinctExists==1
        if i==1
            disp('Using extinction session data in the clustering');
        end
        currSes = 2;
        trialWise5(:,i) = (eLFPpopData.session(currSes).cont.trode(i).rAmp(1,:)' - mean(eLFPpopData.session(currSes).cont.trode(i).rAmp(1,:))) ./ std(eLFPpopData.session(currSes).cont.trode(i).rAmp(1,:));
        trialWise6(:,i) = (eLFPpopData.session(currSes).cont.trode(i).rAmp(2,:)' - mean(eLFPpopData.session(currSes).cont.trode(i).rAmp(2,:))) ./ std(eLFPpopData.session(currSes).cont.trode(i).rAmp(2,:));
        trialWise7(:,i) = (eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:)' - mean(eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:))) ./ std(eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:));
        currSes = 1;
    end
    %     trialWise1(:,i) = eLFPpopData.session(currSes).cont.trode(i).rAmp(1,:)' ./ max(abs(eLFPpopData.session(currSes).cont.trode(i).rAmp(1,:)));
%     trialWise2(:,i) = eLFPpopData.session(currSes).cont.trode(i).rAmp(2,:)' ./ max(abs(eLFPpopData.session(currSes).cont.trode(i).rAmp(2,:)));
%     trialWise3(:,i) = eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:)' ./ max(abs(eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:)));
end

if currSes==1
    figure(401);
    if clustered==0
        
        [indx,c] = kmeans([trialWise1 ; trialWise2 ; trialWise3],numClust);
        
        if extinctExists
%             [indx, clustCenters, SUMD, D] = kmeans([trialWise1 ; trialWise2 ; trialWise3 ; trialWise4 ; trialWise5 ; trialWise6 ; trialWise7]', numClust);
            [indx, clustCenters, SUMD, D] = kmeans([trialWise1 ; trialWise2 ; trialWise3 ; trialWise4]',numClust);
            [ys,indxS] = sort(indx);
            clustered = 1;
            corrCalc = 0;
        else
            [indx, clustCenters, SUMD, D] = kmeans([trialWise1 ; trialWise2 ; trialWise3 ; trialWise4]',numClust);
            [ys,indxS] = sort(indx);
            clustered = 1;
            corrCalc = 0;      
        end
    end
else
    figure(402);
end

for i=1:numel(indxS)
    trialWise1s(:,i) = trialWise1(:,indxS(i));
    trialWise2s(:,i) = trialWise2(:,indxS(i));
    trialWise3s(:,i) = trialWise3(:,indxS(i));
    trialWise4s(:,i) = trialWise4(:,indxS(i));
    if extinctExists==1
        trialWise5s(:,i) = trialWise5(:,indxS(i));
        trialWise6s(:,i) = trialWise6(:,indxS(i));
        trialWise7s(:,i) = trialWise7(:,indxS(i));
    end
end

eLFPpopData.session(currSes).trialWise1     = trialWise1;
eLFPpopData.session(currSes).trialWise2     = trialWise2;
eLFPpopData.session(currSes).trialWise3     = trialWise3;
eLFPpopData.session(currSes).trialWise4     = trialWise4;

eLFPpopData.session(currSes).trialWise1s    = trialWise1s;
eLFPpopData.session(currSes).trialWise2s    = trialWise2s;
eLFPpopData.session(currSes).trialWise3s    = trialWise3s;
eLFPpopData.session(currSes).trialWise4s    = trialWise4s;

[colorMap] = TNC_CreateRBColormap(1024,'mbr');

subplot(241);
imagesc((trialWise1s)); title('C15 | CS');
colormap(colorMap);

subplot(245);
imagesc(corr(abs(trialWise1s)),[0 1]);
colormap(colorMap);

subplot(242);
colormap(colorMap);
imagesc((trialWise2s)); title('C30 | CS');

subplot(246);
imagesc(corr(abs(trialWise2s)),[0 1]);
colormap(colorMap);

subplot(243);
imagesc((trialWise3s)); title('C90 | CS');
colormap(colorMap);

subplot(247);
imagesc(corr(abs(trialWise3s)),[0 1]);
colormap(colorMap);

subplot(244);
imagesc((trialWise4s));  title('C90 | US');
colormap(colorMap);

subplot(248);
imagesc(corr(abs(trialWise4s)),[0 1]);
colormap(colorMap);

figure(403);
subplot(241);
imagesc((trialWise5s)); title('C15 | eCS');
colormap(colorMap);

subplot(245);
imagesc(corr(abs(trialWise5s)),[0 1]);
colormap(colorMap);

subplot(242);
imagesc((trialWise6s)); title('C30 | eCS');
colormap(colorMap);

subplot(246);
imagesc(corr(abs(trialWise6s)),[0 1]);
colormap(colorMap);

subplot(243);
imagesc((trialWise7s)); title('C90 | eCS');
colormap(colorMap);

subplot(247);
imagesc(corr(abs(trialWise7s)),[0 1]);
colormap(colorMap);

figure(402);
plot(trialWise3(:,11),trialWise4(:,11),'k.',trialWise3(:,17),trialWise4(:,17),'r.');
title(corr(trialWise3(:,11),trialWise4(:,11)));

%% Create matrices to examine the trialwise change in LFP magnitude according to cluster
currSes = 1; i=1;

trialMag1=zeros(numel(eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(1,:)),numel(indxS));
trialMag2=zeros(numel(eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(1,:)),numel(indxS));
trialMag3=zeros(numel(eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(1,:)),numel(indxS));
trialMag4=zeros(numel(eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(1,:)),numel(indxS));
trialMag5=zeros(numel(eLFPpopData.session(2).cont.trode(indxS(i)).rAmp(1,:)),numel(indxS));
trialMag6=zeros(numel(eLFPpopData.session(2).cont.trode(indxS(i)).rAmp(1,:)),numel(indxS));
trialMag7=zeros(numel(eLFPpopData.session(2).cont.trode(indxS(i)).rAmp(1,:)),numel(indxS));

for i=1:numel(indxS)

    currSes=1;        
    trialMag1(:,i) = eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(1,:)' ./ max(eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(1,:));
    trialMag2(:,i) = eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(2,:)' ./ max(eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(2,:));
    trialMag3(:,i) = eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(3,:)' ./ max(eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(3,:));
    trialMag4(:,i) = eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(4,:)' ./ max(eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(4,:));

    currSes=2;
    trialMag5(:,i) = eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(1,:)' ./ max(eLFPpopData.session(1).cont.trode(indxS(i)).rAmp(1,:));
    trialMag6(:,i) = eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(2,:)' ./ max(eLFPpopData.session(1).cont.trode(indxS(i)).rAmp(2,:));
    trialMag7(:,i) = eLFPpopData.session(currSes).cont.trode(indxS(i)).rAmp(3,:)' ./ max(eLFPpopData.session(1).cont.trode(indxS(i)).rAmp(3,:));

    figure(900);
    subplot(6,6,i);
    channelResp = [eLFPpopData.session(1).cont.trode(indxS(i)).rAmp(3,:) eLFPpopData.session(2).cont.trode(indxS(i)).rAmp(3,:)];
    plot(1:148,channelResp./max(eLFPpopData.session(1).cont.trode(indxS(i)).rAmp(3,:)).*8000,'ko-', 1:148, 8000-[eLFPpopData.session(1).events.fl eLFPpopData.session(2).events.fl],'b-');
    title([num2str(indxS(i)+32) '...' num2str(corr(channelResp',[eLFPpopData.session(1).events.fl eLFPpopData.session(2).events.fl]'))]);
    
    disp(indxS(i)+32);
    [rho, pval] = corr(channelResp',[eLFPpopData.session(1).events.fl eLFPpopData.session(2).events.fl]')

    figure(901);
    subplot(6,6,i);
    channelResp = [eLFPpopData.session(1).cont.trode(indxS(i)).rAmp(2,:) eLFPpopData.session(2).cont.trode(indxS(i)).rAmp(2,:)];
    plot(1:148,channelResp./max(eLFPpopData.session(1).cont.trode(indxS(i)).rAmp(2,:)).*8000,'ko-', 1:148, 8000-[eLFPpopData.session(1).events.fl eLFPpopData.session(2).events.fl],'b-');
    title(indxS(i)+32);
    
end

for i=1:7
    figure(801);
    subplot(1,8,i);
    eval(['imagesc(trialMag' num2str(i) ',[0 1])']);
    colormap(mapName);
end

% average class 6 electrodes and class 10 electrodes to compare their
% correlations with the latency to sample variable

behaviorData = (8000-[eLFPpopData.session(1).events.fl eLFPpopData.session(2).events.fl]) ./ 8000;

indxSsel = [45 57 59] - 32

for i=1:numel(indxSsel)
    i
    
    indexToUse = indxSsel(i)
    channelResp6(i,:) = [eLFPpopData.session(1).cont.trode(indexToUse).rAmp(3,:) eLFPpopData.session(2).cont.trode(indexToUse).rAmp(3,:)];
    channelResp6(i,:) = channelResp6(i,:) ./ max(channelResp6(i,:));
[rho6(i), pval6(i)] = corr(channelResp6(i,:)',behaviorData');
    
end

classResp6avg = mean(channelResp6);
classResp6sem = std(channelResp6,[],1) ./ sqrt(numel(indxSsel)-1);

indxSsel = [37 41 43 61] - 32

for i=1:numel(indxSsel)
    
    indexToUse = indxSsel(i);
    channelResp10(i,:) = [eLFPpopData.session(1).cont.trode(indexToUse).rAmp(3,:) eLFPpopData.session(2).cont.trode(indexToUse).rAmp(3,:)];
    channelResp10(i,:) = channelResp10(i,:) ./ max(channelResp10(i,:));
[rho10(i), pvall0(i)] = corr(channelResp10(i,:)',behaviorData');
    
end

classResp10avg = mean(channelResp10);
classResp10sem = std(channelResp10,[],1) ./ sqrt(numel(indxSsel)-1);



figure(920); clf;
subplot(121); 
plot(1:148,behaviorData,'ko','Color',[0.73 0.73 0.73]);hold on;
shadedErrorBar(1:148,classResp6avg,classResp6sem,'r'); 
plot(1:148,sgolayfilt(behaviorData,3,9),'k-','LineWidth',0.5);

subplot(122); 
shadedErrorBar(1:148,classResp10avg,classResp10sem,'r');
 hold on;
plot(1:148,behaviorData,'k.');


exportCorrs = [rho6 rho10 ; pval6 pvall0]'

%% The two DA clusters ([37, 41, 43, 61] & [45,57,59]) look distinctive in their ACQ - EXT transition and correlation with latency
    % To examine this in more detail look at 
%         the aligned raster plots of DA units on these channels separately
%         the correlation of mEP components and FL (across learning and extinction)
        
% sort acquisition phase by first lick latency
[lickTimes , trialInds] = sort( double(eLFPpopData.session(1).events.fl) );
% trialInds=1:88;
clear rhoc* pvalc*
% Find the relevant unit numbers.
class1.trodes   = [37, 41, 43, 61];
class2.trodes   = [45, 57, 59];
class1.uID      = [];
class2.uID      = [];

numUnits = numel(eLFPpopData.session(1).unit);
for i=1:numUnits
    if numel(find(class1.trodes == eLFPpopData.session(1).unit(i).el))>0
        if eLFPpopData.session(1).unit(i).da==1
            class1.uID = [class1.uID i]
        end
    elseif numel(find(class2.trodes == eLFPpopData.session(1).unit(i).el))>0
        if eLFPpopData.session(1).unit(i).da==1
            class2.uID = [class2.uID i]            
        end        
    end        
end

winSamples = numel([-window(1):window(2)]);
avgResponse1 = zeros(148,winSamples);
avgResponse2 = zeros(148,winSamples);

figure(1001); clf;
for j=1:numel(class1.uID)    
    subplot(1,numel(class1.uID),j);
    currResponse = [eLFPpopData.session(1).unit(class1.uID(j)).rCS.wins ; eLFPpopData.session(2).unit(class1.uID(j)).rCS.wins];
    avgResponse1 = avgResponse1 + (currResponse.*1000);
    imagesc(currResponse); hold on;
    colormap(mapName);
    plot([-window(1) window(2)], [88.5 88.5], 'k-');
    title(['u' num2str(eLFPpopData.session(1).unit(class1.uID(j)).ui) '-e' num2str(eLFPpopData.session(1).unit(class1.uID(j)).el)]);
    
    respMag(j,:) = [eLFPpopData.session(1).unit(class1.uID(j)).rAmp eLFPpopData.session(2).unit(class1.uID(j)).rAmp];
    [rhoc1(j,1), pvalc1(j,1)] = corr(respMag(j,:)',behaviorData');

end

figure(1002); clf;
for j=1:numel(class2.uID)    
    subplot(1,numel(class2.uID),j);
    currResponse = [eLFPpopData.session(1).unit(class2.uID(j)).rCS.wins ; eLFPpopData.session(2).unit(class2.uID(j)).rCS.wins];
    avgResponse2 = avgResponse2 + (currResponse.*1000);
    imagesc(currResponse); hold on;
    colormap(mapName);
    plot([-window(1) window(2)], [88.5 88.5], 'k-');
    title(['u' num2str(eLFPpopData.session(1).unit(class2.uID(j)).ui) '-e' num2str(eLFPpopData.session(1).unit(class2.uID(j)).el)]);

    respMag2(j,:) = [eLFPpopData.session(1).unit(class2.uID(j)).rAmp eLFPpopData.session(2).unit(class2.uID(j)).rAmp];
    [rhoc2(j,1), pvalc2(j,1)] = corr(respMag2(j,:)',behaviorData');
end

figure(1003); clf;
avgResponse1 = avgResponse1 ./ numel(class1.uID);
avgResponse2 = avgResponse2 ./ numel(class2.uID);
subplot(2,3,1:2); imagesc(avgResponse1([1:88 89:148],:), [-70 70]); hold on; plot([1 winSamples], [88.5 88.5], 'k-'); title(class1.trodes); colormap(mapName);
subplot(2,3,4:5); imagesc(avgResponse2([1:88 89:148],:), [-70 70]); hold on; plot([1 winSamples], [88.5 88.5], 'k-'); title(class2.trodes); colormap(mapName);
subplot(2,3,[3]); 
plot([-window(1):window(2)],mean(avgResponse1(1:88,1:winSamples)),'k-',[-window(1):window(2)],mean(avgResponse2(1:88,1:winSamples)),'r-'); 
axis([-window(1) window(2) 0 70]);
subplot(2,3,[6]); plot([-window(1):window(2)],mean(avgResponse1(108:148,1:winSamples)),'k-',[-window(1):window(2)],mean(avgResponse2(118:148,1:winSamples)),'r-'); 
axis([-window(1) window(2) 0 70]);

%% Also need to look at the pairwise correlations within, and across functional groups
% DA unit list ordered by class membership

DAinds = [class1.uID class2.uID];
DAindsT = [class1.trodes class2.trodes];

DAcorrMat.rho = zeros(numel(DAinds),numel(DAinds));
DAcorrMat.rhoE = zeros(numel(DAinds),numel(DAinds));

for i=1:numel(DAinds)
    for j=1:numel(DAinds)
        [rho , pval] = corr(eLFPpopData.session(1).unit(DAinds(i)).rAmp(1,:)' , eLFPpopData.session(1).unit(DAinds(j)).rAmp(1,:)');
%         if pval < 0.05
            DAcorrMat.rho(i,j) = rho;
            DAcorrMat.pval(i,j) = pval;
%         else
%             DAcorrMat.rho(i,j) = 0;
%             DAcorrMat.pval(i,j) = 1;            
%         end 

        [rho , pval] = corr(eLFPpopData.session(2).unit(DAinds(i)).rAmp(1,:)' , eLFPpopData.session(2).unit(DAinds(j)).rAmp(1,:)');
%         if pval < 0.05
            DAcorrMat.rhoE(i,j) = rho;
            DAcorrMat.pvalE(i,j) = pval;
%         else
%             DAcorrMat.rhoE(i,j) = 0;
%             DAcorrMat.pvalE(i,j) = 1;            
%         end
    end
    
%     [rho , pval] = corr([eLFPpopData.session(1).unit(DAinds(i)).rAmp(1,:) eLFPpopData.session(2).unit(DAinds(i)).rAmp(1,:)]' , [double(eLFPpopData.session(1).events.fl) double(eLFPpopData.session(2).events.fl)]');
%     [rho , pval] = corr(eLFPpopData.session(1).unit(DAinds(i)).rAmp(1,:)' , double(eLFPpopData.session(1).events.fl)');
    [rho , pval] = corr(eLFPpopData.session(2).unit(DAinds(i)).rAmp(1,:)' , double(eLFPpopData.session(2).events.fl)');
    if pval < 0.05
        DAcorrMat.rhoB(i) = rho;
        DAcorrMat.pvalB(i) = pval;
    else
        DAcorrMat.rhoB(i) = 0;
        DAcorrMat.pvalB(i) = 1;            
    end
    
    [rho , pval] = corr(eLFPpopData.session(1).unit(DAinds(i)).rAmp(1,:)' , eLFPpopData.session(1).cont.trode(eLFPpopData.session(1).unit(DAinds(i)).el-32).rAmp(3,:)');
    if pval < 0.05
        DAcorrMat.rhoL(i) = rho;
        DAcorrMat.pvalL(i) = pval;
    else
        DAcorrMat.rhoL(i) = 0;
        DAcorrMat.pvalL(i) = 1;            
    end    
end

[mapName] = TNC_CreateRBColormap(1024,'mbr');

figure(37); clf; 
subplot(121); imagesc(DAcorrMat.rho,[-1 1]); 
    hold on; plot([6.5 6.5] , [0.5 15.5] , 'k', [0.5 15.5] , [6.5 6.5] ,  'k'); title('Acquisition'); xlabel('Unit ID'); set(gca, 'TickDir', 'out'); box off;
subplot(122); imagesc(DAcorrMat.rhoE,[-1 1]); colormap(mapName);
    hold on; plot([6.5 6.5] , [0.5 15.5] , 'k', [0.5 15.5] , [6.5 6.5] ,  'k'); title('Extinction'); xlabel('Unit ID'); set(gca, 'TickDir', 'out'); box off;

%% EXTRACT THE LATENCY TO FIRST LICK



for j = 1:2

    allCS = numel(eLFPpopData.session(j).events.cs)
    
    eLFPpopData.session(j).events.fl = zeros(1, allCS);
    
    for i=1:allCS
        
        currCSstamp = eLFPpopData.session(j).events.cs(i);
        flInd = find(eLFPpopData.session(j).events.al>currCSstamp+500 & eLFPpopData.session(j).events.al<currCSstamp+8000);

        if numel(flInd)>0
            latency = eLFPpopData.session(j).events.al(flInd(1)) - currCSstamp;
            eLFPpopData.session(j).events.fl(i) = latency;
        else
            eLFPpopData.session(j).events.fl(i) = 8000;
        end
        
    end

end

%% CORRELTAION OF FL AND LFP DATA

currSes = 1;
numChan = size(eLFPpopData.session(currSes).trialWise3,2);
eLFPpopData.session(currSes).events.sig = zeros(numChan,3);
eLFPpopData.session(currSes).events.corrScores = zeros(numChan,4);

for p=1:numChan
    eLFPpopData.session(currSes).events.corrScores(p,4) = p+32;
    [eLFPpopData.session(currSes).events.corrScores(p,1),eLFPpopData.session(currSes).events.sig(p,1)] = corr( eLFPpopData.session(currSes).trialWise1(:,p),double(eLFPpopData.session(currSes).events.fl') );
    if eLFPpopData.session(currSes).events.sig(p,1) > 0.05
        eLFPpopData.session(currSes).events.corrScores(p,1) = 0;
    end
    [eLFPpopData.session(currSes).events.corrScores(p,2),eLFPpopData.session(currSes).events.sig(p,2)] = corr( eLFPpopData.session(currSes).trialWise2(:,p),double(eLFPpopData.session(currSes).events.fl') );
    if eLFPpopData.session(currSes).events.sig(p,2) > 0.05
        eLFPpopData.session(currSes).events.corrScores(p,2) = 0;
    end
    [eLFPpopData.session(currSes).events.corrScores(p,3),eLFPpopData.session(currSes).events.sig(p,3)] = corr( eLFPpopData.session(currSes).trialWise3(:,p),double(eLFPpopData.session(currSes).events.fl') );
    eLFPpopData.session(currSes).events.sig(p,3)
    if eLFPpopData.session(currSes).events.sig(p,3) > 0.05
        eLFPpopData.session(currSes).events.corrScores(p,3) = 0;
    end
    
    [rhoA(p), pA(p)] = corr( eLFPpopData.session(1).trialWise1(:,p), double(eLFPpopData.session(1).events.fl'));
    [rhoB(p), pB(p)] = corr( eLFPpopData.session(1).trialWise2(:,p), double(eLFPpopData.session(1).events.fl'));
    [rhoC(p), pC(p)] = corr( eLFPpopData.session(1).trialWise3(:,p), double(eLFPpopData.session(1).events.fl'));
        
    [rhoD(p), pD(p)] = corr( eLFPpopData.session(currSes).cont.trode(p).rAmp(3,:)', double(eLFPpopData.session(currSes).events.fl)');
    
end

figure(702);
    semilogy([1:32]+32,pA,'bo-',[1:32]+32,pB,'ro-',[1:32]+32,pC,'ko-',[1:32]+32,ones(1,32).*0.05,'k--');
    
    
%%

if currSes==1
    figure(500);
else
    figure(501);
end
plot(indx,eLFPpopData.session(currSes).events.corrScores(:,1),'ko',indx+0.2,eLFPpopData.session(currSes).events.corrScores(:,2),'r+',indx+0.4,eLFPpopData.session(currSes).events.corrScores(:,3),'bs');
axis([0 numClust+1 -1 1]);


% Calculate trial-based lick rasters
%     if currParams.filter.causal==1
%         stimTimes   = double(eLFPpopData.session(currSes).events.cs) - (numel(currParams.filter.kernel)./2);
%     else
        stimTimes   = double(eLFPpopData.session(currSes).events.cs);
%     end
    numStamps   = numel(eLFPpopData.session(currSes).events.al);
    delta       = zeros(1,ceil(eLFPpopData.session(currSes).events.al(numStamps)));
    [rLITE]     = TNC_AlignRasters(delta , double(eLFPpopData.session(currSes).events.al) , -1 , stimTimes , [1000,currParams.winParams.after],1,1);
    eLFPpopData.session(currSes).licksCS.raster      = rLITE.raster;

    exportToIgor = [];
        eLFPpopData.session(j).events.fl = zeros(numel(eLFPpopData.session(currSes).licksCS.raster.trial),1);
    for i=1:numel(eLFPpopData.session(currSes).licksCS.raster.trial)
        numLicks = numel(eLFPpopData.session(currSes).licksCS.raster.trial(i).ts);
        firstLick = find(eLFPpopData.session(currSes).licksCS.raster.trial(i).ts>0,1);

        
        if numel(firstLick)>0
            eLFPpopData.session(j).events.fl(i) = eLFPpopData.session(currSes).licksCS.raster.trial(i).ts(firstLick(1));
        else
            eLFPpopData.session(j).events.fl(i) = 5000;            
        end
        
        if i==1
            exportToIgor = [ones(numLicks,1).*i eLFPpopData.session(currSes).licksCS.raster.trial(i).ts];
        else
            exportToIgorTmp = [ones(numLicks,1).*i eLFPpopData.session(currSes).licksCS.raster.trial(i).ts];
            exportToIgor = [exportToIgor ; exportToIgorTmp];
        end
    end
corrCalc=1;

%% Calculate the smoothed psth for each unit   
numUnits = size(eLFPpopData.session(currSes).unit,2);
clear physCorr
unitToLFPcorrs = []; 
elList = [];

currSes = 2;

for k=1:numUnits

    numStamps   = length(eLFPpopData.session(currSes).unit(k).ts);
    delta       = zeros(1,ceil(eLFPpopData.session(currSes).unit(k).ts(numStamps)));
    delta(1,round(eLFPpopData.session(currSes).unit(k).ts)+1) = 1;
    tmpSmooth   = conv(delta,currParams.filter.kernel,'same');
%     if currParams.filter.causal==1
%         stimTimes   = double(eLFPpopData.session(currSes).events.cs) - (numel(currParams.filter.kernel)./2);
%     else
        stimTimes   = double(eLFPpopData.session(currSes).events.cs);
%     end

    if currSes == 2
        stimTimes = stimTimes(20:60);
    end
    
    % for reference
    % [response] = TNC_AlignRasters(delta,spkStamps,stableTrials,alignStamps,window,rasterFlag,boxcar)
    [rLITE] = TNC_AlignRasters(delta , double(eLFPpopData.session(currSes).unit(k).ts+1) , -1 , stimTimes , [currParams.winParams.prior,currParams.winParams.after],1,1);
    eLFPpopData.session(currSes).unit(k).respCS.raster      = rLITE.raster;
    eLFPpopData.session(currSes).unit(k).respCS.boxcar      = rLITE.image.boxcar;

    [rsLITE] = TNC_AlignRasters(tmpSmooth , double(eLFPpopData.session(currSes).unit(k).ts+1) , -1 , stimTimes , [currParams.winParams.prior,currParams.winParams.after],0,1);
    eLFPpopData.session(currSes).unit(k).respCS.psthAVG    = rsLITE.image.psthAVG;
    eLFPpopData.session(currSes).unit(k).respCS.psthSEM    = rsLITE.image.psthSEM;
    eLFPpopData.session(currSes).unit(k).respCS.psthZ      = rsLITE.image.psthZ;
    eLFPpopData.session(currSes).unit(k).respCS.psthZe     = rsLITE.image.psthZe;
    eLFPpopData.session(currSes).unit(k).respCS.trialwise  = rsLITE.image.aligned;
  
    for n = 1:size(eLFPpopData.session(currSes).unit(k).respCS.trialwise,1)
        eLFPpopData.session(currSes).unit(k).respCS.intResp(n) = trapz(rsLITE.image.aligned(n,500:650)) - trapz(rsLITE.image.aligned(n,300:450));
        if currSes==1
            eLFPpopData.session(currSes).unit(k).respUS.intResp(n) = trapz(rsLITE.image.aligned(n,3000:3150)) - trapz(rsLITE.image.aligned(n,300:450));
        end
    end
    
    if currSes == 1
        alignCS.psthzMat(:,k) = eLFPpopData.session(currSes).unit(k).respCS.psthZ;
    else
        alignCSe.psthzMat(:,k) = eLFPpopData.session(currSes).unit(k).respCS.psthZ;
    end
    
end

figure(2);
subplot(2,1,currSes);
    if currSes == 1
        imagesc(alignCS.psthzMat',[-max(max(alignCS.psthzMat))./5 max(max(alignCS.psthzMat))./5]);
    else
        imagesc(alignCSe.psthzMat',[-max(max(alignCS.psthzMat))./5 max(max(alignCS.psthzMat))./5]);
    end
colormap(colorMap);
set(gcf,'Color','white');
set(gca,'TickDir','out','TickLength',[0.005 0]);

%%
for q = 1:2
    unitsALL = [];
    unitsALLlist.da = [];
    unitsALLlist.el = [];

    for k=1:numUnits
        if abs(mean(eLFPpopData.session(q).unit(k).respCS.intResp))>0.5
            unitsALLlist.da = [unitsALLlist.da eLFPpopData.session(q).unit(k).da];
            unitsALLlist.el = [unitsALLlist.el eLFPpopData.session(q).unit(k).el];
            unitsALL = [unitsALL eLFPpopData.session(q).unit(k).respCS.intResp'];
        end
    end

    figure(3);
    subplot(2,1,q);
    [corrMat,pMat] = corr(unitsALL);
    nSigInds = find(pMat>0.05);
    corrMat(nSigInds) = 0;
    imagesc(tril(corrMat,-1),[-1 1]);
    colormap(colorMap); 
    colorbar;
    set(gcf,'Color','white');
    set(gca,'TickDir','out','TickLength',[0.005 0]);
    
%     if q==1
%         priorMat = tril(corrMat,-1);
%     else
%         subplot(313);
%         imagesc(priorMat-tril(corrMat,-1),[-1 1]);
%         colorbar;
%         set(gcf,'Color','white');
%         set(gca,'TickDir','out','TickLength',[0.005 0]);        
%     end
    
    eLFPpopData.session(q).unitsALL = unitsALL;
end

%% FOR EACH DOPAMINE UNIT EXAMINE THE CS AND US RESPONSE MAG AS FUNCTION OF FIRST LICK

figure(800); clf;

[lickTimes , trialInds] = sort( double(eLFPpopData.session(currSes).events.fl) );

lickTimes = lickTimes(1:86);
trialInds = trialInds(1:86);
    
for i=1:numel(daUnits)
    
    currUnit = daUnits(i);
    figure(800);
    subplot(4,4,i);
%     plot(abs(eLFPpopData.session(currSes).events.fl-2500) , eLFPpopData.session(currSes).unit(currUnit).respCS.intResp , 'ko'); hold on;
    plot( lickTimes , eLFPpopData.session(currSes).unit(currUnit).respCS.intResp(trialInds) , 'k-');
    plot( lickTimes , eLFPpopData.session(currSes).unit(currUnit).respUS.intResp(trialInds) , 'b-');

    if i==1
        intRespUS = eLFPpopData.session(currSes).unit(currUnit).respUS.intResp(trialInds);
    else
        intRespUS = intRespUS + eLFPpopData.session(currSes).unit(currUnit).respUS.intResp(trialInds);        
    end    
    
end

intRespUS = intRespUS ./ numel(daUnits);
figure(801); clf;
p = polyfit(lickTimes-2500,intRespUS,2);
fitData = polyval(p,lickTimes-2500);
plot(lickTimes-2500 , intRespUS , 'k.' , 'Color', [0.5 0.5 0.5]); hold on;
plot(lickTimes-2500 , fitData , 'r-', 'LineWidth' , 2);
plot([-4000 4000] , [0 0], 'k--' , [0 0] , [-1 4], 'k--');
axis([-2500 4000 -1 4]);
title('Summary data for 16 Da units');
xlabel('First Lick relative to reward (ms)');
ylabel('DA Neuron Response Amplitude (z)');

%% FOR EACH ELECTRODE LOOK AT CORRELATION BETWEEN UNITS AND LFP

for j=1:numUnits
    elList = [elList eLFPpopData.session(currSes).unit(j).el];
end

UelList = unique(elList);

for m=1:numel(UelList)
    
    currChan = UelList(m);

    units = [];

    for k=1:numUnits
        if eLFPpopData.session(currSes).unit(k).el == currChan
            units = [units eLFPpopData.session(currSes).unit(k).respCS.intResp'];
        end
    end

    physCorr.el(m).numU = numel(find(elList==currChan));
    physCorr.el(m).lfp  = eLFPpopData.session(currSes).trialWise3(:,currChan-32);
    physCorr.el(m).el   = currChan;
    physCorr.el(m).units= units;
    
end

for k=1:numUnits
    [rho,pval] = corr(eLFPpopData.session(currSes).unit(k).respCS.intResp',double(eLFPpopData.session(currSes).events.fl)')
    if pval<0.05
        physCorr.uCorr.withFL(k,1) = rho;
        physCorr.uCorr.withFL(k,2) = pval;   
        physCorr.uCorr.withFL(k,3) = eLFPpopData.session(currSes).unit(k).el  
    else
        physCorr.uCorr.withFL(k,1) = NaN;
        physCorr.uCorr.withFL(k,2) = 1;           
        physCorr.uCorr.withFL(k,3) = eLFPpopData.session(currSes).unit(k).el  
    end
end

for k=1:numUnits
    [rho,pval] = corr(eLFPpopData.session(currSes).unit(k).respCS.intResp',eLFPpopData.session(currSes).trialWise1(:,physCorr.uCorr.withFL(k,3)-32));
    if pval<0.05
        physCorr.uCorr.withFL(k,4) = rho;
        physCorr.uCorr.withFL(k,5) = pval;   
    else
        physCorr.uCorr.withFL(k,4) = NaN;
        physCorr.uCorr.withFL(k,5) = 1;           
    end
end

for k=1:numUnits
    [rho,pval] = corr(eLFPpopData.session(currSes).unit(k).respCS.intResp',eLFPpopData.session(currSes).trialWise2(:,physCorr.uCorr.withFL(k,3)-32));
    if pval<0.05
        physCorr.uCorr.withFL(k,6) = rho;
        physCorr.uCorr.withFL(k,7) = pval;   
    else
        physCorr.uCorr.withFL(k,6) = NaN;
        physCorr.uCorr.withFL(k,7) = 1;           
    end
end

for k=1:numUnits
    [rho,pval] = corr(eLFPpopData.session(currSes).unit(k).respCS.intResp',eLFPpopData.session(currSes).trialWise3(:,physCorr.uCorr.withFL(k,3)-32));
    if pval<0.05
        physCorr.uCorr.withFL(k,8) = rho;
        physCorr.uCorr.withFL(k,9) = pval;   
    else
        physCorr.uCorr.withFL(k,8) = NaN;
        physCorr.uCorr.withFL(k,9) = 1;           
    end
end

currSes
physCorr.uCorr.withFL

eLFPpopData.session(currSes).physCorr = physCorr;
eLFPpopData.session(currSes).alignCS = alignCS;

%% PLOT POSITIONS WHERE THERE WERE SIGNIFICANT CORRELATIONS BETWEEN LFP+BEH, UNIT+BEH, LFP+UNIT 


for currSes=1:1
    figure(499+currSes); clf;
    figure(599+currSes); clf;

    figure(499+currSes);
    for i=1:32
        [r,c] = find(elecMap==i);
        plot3([c c],[r r],[1 3],'-','Color',[0.7 0.7 0.7],'LineWidth',numel(find(i==eLFPpopData.session(currSes).physCorr.uCorr.withFL(:,3)-32))+1); hold on;
%         if numel(find(i==eLFPpopData.session(currSes).physCorr.uCorr.withFL(:,3)-32)) > 0
%             plot3(c,r,1,'s','Color',[0.7 0.7 0.7],'MarkerSize',numel(find(i==eLFPpopData.session(currSes).physCorr.uCorr.withFL(:,3)-32)),'LineWidth',2);
%             plot3(c,r,3,'s','Color',[0.7 0.7 0.7],'MarkerSize',numel(find(i==eLFPpopData.session(currSes).physCorr.uCorr.withFL(:,3)-32)),'LineWidth',2);
%         end
    end

    for i=1:size(eLFPpopData.session(currSes).physCorr.uCorr.withFL,1)

        [r,c] = find(elecMap==eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,3)-32);

        if eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,5)<1
            plot3(c,r,3,'o','Color',[1 0.5 0],'LineWidth',2,'MarkerSize',round(abs(eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,4)).*75));
        end
        if eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,7)<1
            plot3(c,r,3,'o','Color',[0.5 0.5 0.5],'LineWidth',2,'MarkerSize',round(abs(eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,6)).*75));
        end
        if eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,9)<1
            plot3(c,r,3,'o','Color',[0 0.5 1],'LineWidth',2,'MarkerSize',round(abs(eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,8)).*75));
        end
        
        axis([0 9 0 5 1 3]);

        if eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,2)<1
            plot3(c,r,1,'o','Color',[1 0 0],'LineWidth',2,'MarkerSize',round(abs(eLFPpopData.session(currSes).physCorr.uCorr.withFL(i,1)).*75));
            axis([0 9 0 5 1 3]); 
        end
    end
    
    for i=1:32

        [r,c] = find(elecMap==i);
        
        if sum(eLFPpopData.session(currSes).events.corrScores(i,:))>0
            plot3(c,r,2,'o','Color',[0 0 0],'LineWidth',2,'MarkerSize',round(abs(max(eLFPpopData.session(currSes).events.corrScores(i,:))).*75));
            axis([0 9 0 5 1 3]); 
            view([-28 28]);
        end    

    end
    
    set(gcf,'Color','white');
    axis off;
end

%%

    figure(599+currSes); clf;
%     for i=1:32
%         [r,c] = find(elecMap==i);
%         subplot(121);
%         plot3(c,r,1,'o','Color',[0.75 0.75 0.75], 'MarkerSize',2,'LineWidth',2); hold on;
%         subplot(122);
%         plot3(c,r,1,'o','Color',[0.75 0.75 0.75], 'MarkerSize',2,'LineWidth',2); hold on;
%     end

%     for i=1:size(eLFPpopData.session(currSes).unit,2)
%         [r,c] = find(elecMap==eLFPpopData.session(currSes).unit(i).el-32);
% %         if eLFPpopData.session(currSes).unit(i).ui>1
% %             subplot(121);
% %             plot3([c,c],[r,r],[eLFPpopData.session(currSes).unit(i).ui,eLFPpopData.session(currSes).unit(i).ui-1],'-','Color',[0.75 0.75 0.75], 'LineWidth',1); hold on;
% %             subplot(122);
% %             plot3([c,c],[r,r],[eLFPpopData.session(currSes).unit(i).ui,eLFPpopData.session(currSes).unit(i).ui-1],'-','Color',[0.75 0.75 0.75], 'LineWidth',1); hold on;            
% %         end
%         subplot(121);
%         plot3(c,r,eLFPpopData.session(currSes).unit(i).ui+1,'o','Color',[0.5 0.5 0.5], 'MarkerSize',6,'LineWidth',2); hold on;
%         subplot(122);
%         plot3(c,r,eLFPpopData.session(currSes).unit(i).ui+1,'o','Color',[0.5 0.5 0.5], 'MarkerSize',6,'LineWidth',2); hold on;
% 
%     end
    
    [corrMat,pMat] = corr(eLFPpopData.session(currSes).unitsALL);
    nSigInds = find(pMat>0.05);
    corrMat(nSigInds) = 0;
    forPlotMat = tril(corrMat,-1);
    [toPlotR, toPlotC] = find(abs(forPlotMat)>0);
    for j=1:numel(toPlotC)
        [r1,c1] = find(elecMap==eLFPpopData.session(currSes).unit(toPlotR(j)).el-32);
        [r2,c2] = find(elecMap==eLFPpopData.session(currSes).unit(toPlotC(j)).el-32);
        if eLFPpopData.session(currSes).unit(toPlotR(j)).el ~= eLFPpopData.session(currSes).unit(toPlotC(j)).el
            if forPlotMat(toPlotR(j),toPlotC(j))>0
                subplot(131);
%                 plot3([c1,c2],[r1,r2],[eLFPpopData.session(currSes).unit(toPlotR(j)).ui+1,eLFPpopData.session(currSes).unit(toPlotC(j)).ui+1],'-','Color',[1 0 0],'LineWidth',round( abs(forPlotMat(toPlotR(j),toPlotC(j)))*5 )); hold on;
                plot([c1,c2],[r1,r2],'-','Color',[1 0 0],'LineWidth',round( abs(forPlotMat(toPlotR(j),toPlotC(j)))*5 )); hold on;
                axis([0 9 0 5]);
                set(gcf,'Color','white');
                axis off;
%                 view([-40 50]);
            else
                subplot(132);
%                 plot3([c1,c2],[r1,r2],[eLFPpopData.session(currSes).unit(toPlotR(j)).ui+1,eLFPpopData.session(currSes).unit(toPlotC(j)).ui+1],'-','Color',[0 0.67 1],'LineWidth',round( abs(forPlotMat(toPlotR(j),toPlotC(j)))*5 )); hold on;
                plot([c1,c2],[r1,r2],'-','Color',[0 0.67 1],'LineWidth',round( abs(forPlotMat(toPlotR(j),toPlotC(j)))*5 )); hold on;
                axis([0 9 0 5]);
                set(gcf,'Color','white');
                axis off;
%                 view([-40 50]);
            end
        end
    end
    
    subplot(133)
    % compare this to the pairwise noise correlations between channels for
    % the mep components

%% MAKE A HIVE PLOT TO EXAMINE CORRELATION STRUCTURE

% The electrode on which the unit was recorded gives the dimension (1-4)
    % clust 1 = count(1),0
    % clust 1 = 0,count(2)
    % clust 1 = -count(3)
    % clust 1 = 0,-count(4)
    
    count= zeros(5,1);
    figure(700); clf;
    
    for i=1:numUnits
        thisUnitClass = indx(eLFPpopData.session(currSes).unit(i).el-32);
        count(thisUnitClass) = count(thisUnitClass)+1;
        
        switch thisUnitClass
            
            case 1
                eLFPpopData.session(currSes).unit(i).hive.node.x = 0;
                eLFPpopData.session(currSes).unit(i).hive.node.y = count(1);
                
            case 2
                eLFPpopData.session(currSes).unit(i).hive.node.x = count(2);
                eLFPpopData.session(currSes).unit(i).hive.node.y = -count(2);
                
            case 3
                eLFPpopData.session(currSes).unit(i).hive.node.x = count(3);
                eLFPpopData.session(currSes).unit(i).hive.node.y = 0;
                
            case 4
                eLFPpopData.session(currSes).unit(i).hive.node.x = -count(4);
                eLFPpopData.session(currSes).unit(i).hive.node.y = -count(4);
        end

    end
    

    [corrMat,pMat] = corr(eLFPpopData.session(currSes).unitsALL);
    nSigInds = find(pMat>0.05);
    corrMat(nSigInds) = 0;
    forPlotMat = tril(corrMat,-1);
    [toPlotR, toPlotC] = find(abs(forPlotMat)>0);
    
    for j=1:numel(toPlotC)
        
        x1 = eLFPpopData.session(currSes).unit(toPlotR(j)).hive.node.x;
        y1 = eLFPpopData.session(currSes).unit(toPlotR(j)).hive.node.y;

        x2 = eLFPpopData.session(currSes).unit(toPlotC(j)).hive.node.x;
        y2 = eLFPpopData.session(currSes).unit(toPlotC(j)).hive.node.y;
        
        if eLFPpopData.session(currSes).unit(toPlotR(j)).el ~= eLFPpopData.session(currSes).unit(toPlotC(j)).el
            if forPlotMat(toPlotR(j),toPlotC(j))>0
                subplot(121);
                plot([x1,x2],[y1,y2],'-','Color',[1 0 0],'LineWidth',round( abs(forPlotMat(toPlotR(j),toPlotC(j)))*10 )); hold on;               
            else
                subplot(122);
                plot([x1,x2],[y1,y2],'-','Color',[0 0.67 1],'LineWidth',round( abs(forPlotMat(toPlotR(j),toPlotC(j)))*10 )); hold on;
            end
        end
    end
    
    for i=1:numUnits
        figure(700);
        subplot(121);
        plot(eLFPpopData.session(currSes).unit(i).hive.node.x,eLFPpopData.session(currSes).unit(i).hive.node.y,'k.'); hold on;            
        text(eLFPpopData.session(currSes).unit(i).hive.node.x+1,eLFPpopData.session(currSes).unit(i).hive.node.y,num2str(unitsALLlist.el(i))); hold on;
        if eLFPpopData.session(currSes).unit(i).da == 1
            plot(eLFPpopData.session(currSes).unit(i).hive.node.x,eLFPpopData.session(currSes).unit(i).hive.node.y,'ko','MarkerSize',8); hold on;
        end
        axis([-25 25 -25 25 -25 25]); axis off;
        set(gcf,'Color','white');
        subplot(122);
        plot(eLFPpopData.session(currSes).unit(i).hive.node.x,eLFPpopData.session(currSes).unit(i).hive.node.y,'k.');
        text(eLFPpopData.session(currSes).unit(i).hive.node.x+1,eLFPpopData.session(currSes).unit(i).hive.node.y,num2str(unitsALLlist.el(i))); hold on;
        if eLFPpopData.session(currSes).unit(i).da == 1
            plot(eLFPpopData.session(currSes).unit(i).hive.node.x,eLFPpopData.session(currSes).unit(i).hive.node.y,'ko','MarkerSize',8); hold on;
        end
        axis([-25 25 -25 25 -25 25]); axis off;
        set(gcf,'Color','white');
    end
    
%% Examine mEP and Single Unit correlations for all dopamine units
    
    figure(701);
    
    zScoredBehav = (eLFPpopData.session(1).events.fl - mean(eLFPpopData.session(1).events.fl)) ./ std(eLFPpopData.session(1).events.fl);
    
    for i=1:16
        subplot(4,4,i);
        zScoredData = (eLFPpopData.session(1).unit(daUnits(i)).respCS.intResp - mean(eLFPpopData.session(1).unit(daUnits(i)).respCS.intResp)) ./ std(eLFPpopData.session(1).unit(daUnits(i)).respCS.intResp);
        
        daPopZScored(:,i) = zScoredData';
        
%         plot(eLFPpopData.session(1).trialWise3(:,eLFPpopData.session(1).unit(daUnits(i)).el-32), zScoredData , 'k.')
%         [rho2, p2] = corr( eLFPpopData.session(1).trialWise3(:,eLFPpopData.session(1).unit(daUnits(i)).el-32), eLFPpopData.session(1).events.fl);
        
        [rho, p] = corr( eLFPpopData.session(1).trialWise3(:,eLFPpopData.session(1).unit(daUnits(i)).el-32), zScoredData');

        plot(zScoredBehav , zScoredData, 'k.')
        [rho2, p2] = corr( zScoredData' , zScoredBehav);    
        
        title([num2str( rho2 )  ' | ' num2str(eLFPpopData.session(1).unit(daUnits(i)).el) ' | ' num2str( rho ) ] );
    end
    
    axis([-5 5 -5 5]);        
    
    figure(702);
    for i=1:32
        [rhoA(i), pA(i)] = corr( eLFPpopData.session(1).trialWise1(:,i), double(eLFPpopData.session(1).events.fl));
        [rhoB(i), pB(i)] = corr( eLFPpopData.session(1).trialWise2(:,i), double(eLFPpopData.session(1).events.fl));
        [rhoC(i), pC(i)] = corr( eLFPpopData.session(1).trialWise3(:,i), double(eLFPpopData.session(1).events.fl));
    end
    
    pC
    
    semilogy([1:32]+32,pA,'bo-',[1:32]+32,pB,'ro-',[1:32]+32,pC,'ko-',[1:32]+32,ones(1,32).*0.05,'k--');
    
%% Show all PSTHs from a given electrode

figNum = 405
elNum = 41
sessNum = 1
clear inds
count = 0;

for i=1:numUnits
    if eLFPpopData.session(sessNum).unit(i).el == elNum
        count=count+1;
        inds(count) = i;
    end    
end

figure(elNum); clf;

for j=1:numel(inds)
    subplot(1,numel(inds),j)
    plot(eLFPpopData.session(sessNum).unit(inds(j)).respCS.boxcar);
    title([num2str(elNum) ' - u' num2str(inds(j))]);
end


%% Population analysis of the field potential / behavior correlations (trial by trial)

%% DAma05

clear DAma05

% parameter settings
window          = [100 2500];
windows(1,:)    = [5 20];
windows(2,:)    = [20 40];
windows(3,:)    = [55 150];

% Load field potential data
RD_Active = TNC_LoadData(0, 0, 'DA-Ma-05 100520trace-001001.ns4');

% Load behavioral data and calculate latency to sample the port
EV_Active = TNC_LoadData(0, 0, 'DA-Ma-05 100520trace-001001.nev');

figure(1); clf;
% Calculate the windowed LFP data from the CS event times
DAma05.events.cs = EV_Active.Data.Spikes.Timestamps(find(EV_Active.Data.Spikes.Electrode==137));
% DAma05.events.al = ;

for i=1:16
    rawLFPdata = decimate(sgolayfilt(RD_Active.Data(i,:),3,21),10);
    [rCS] = TNC_ExtTrigWins(rawLFPdata , round(DAma05.events.cs./30) , window);
    DAma05.trode(i).wins = rCS.wins;
    DAma05.trode(i).avg = rCS.avg;
    DAma05.trode(i).err = rCS.err;
end

for i=1:16
    figure(1); subplot(4,4,i);
    shadedErrorBar(-100:2500,DAma05.trode(i).avg,DAma05.trode(i).err,'-r',1);
    hold on;
    plot([-100 2500],[0 0],'k--',[0 0],[-100 25],'k--'); axis([-100 2500 -100 25]); title(i);        
end

% Get the amplitude of the LFP over the event windows

% Calculate the latency to first lick
allCS = numel(DAma05.events.cs)

DAma05.events.fl = zeros(1, allCS);

for i=1:allCS

    currCSstamp = DAma05.events.cs(i);
    flInd = find(DAma05.events.al>currCSstamp+500 & DAma05.events.al<currCSstamp+8000);

    if numel(flInd)>0
        latency = eLFPpopData.session(j).events.al(flInd(1)) - currCSstamp;
        eLFPpopData.session(j).events.fl(i) = latency;
    else
        eLFPpopData.session(j).events.fl(i) = 8000;
    end

end

% Compute the correlations between LFP amplitude and behavior

% Compute the pairwise noise correlations between lfp channels to look for clustering

%% DaMa10

clear DaMa10

egTrode = 5;
[mapName] = TNC_CreateRBColormap(20,'gp');

% parameter settings
window          = [100 2500];
windows(1,:)    = [5 20];
windows(2,:)    = [20 40];
windows(3,:)    = [55 150];

windowSlope     = [55 90];

% Load field potential data
RD_Active = TNC_LoadData(0, 0, 'DaMa10-110329tracet5-005.ns5');
numChan = numel(RD_Active.MetaTags.ChannelID);

% Load behavioral data and calculate latency to sample the port
EV_Active = TNC_LoadData(0, 0, 'DaMa10-110329tracet5-005.nev');

figure(1); clf;
% Calculate the windowed LFP data from the CS event times
DaMa10.events.cs = round(EV_Active.Data.Spikes.Timestamps(find(EV_Active.Data.Spikes.Electrode==137))./30);
DaMa10.events.al = round(EV_Active.Data.Spikes.Timestamps(find(EV_Active.Data.Spikes.Electrode==138))./30);

DaMa10.events

for i=1:numChan
    rawLFPdata = decimate(sgolayfilt(RD_Active.Data(i,:),3,21),30);
    [rCS] = TNC_ExtTrigWins(rawLFPdata , DaMa10.events.cs , window);
    DaMa10.trode(i).wins = rCS.wins;
    DaMa10.trode(i).avg = rCS.avg;
    DaMa10.trode(i).err = rCS.err;
end

for i=1:numChan
    figure(1); subplot(4,ceil(numChan./4),i);
    shadedErrorBar(-100:window(2),DaMa10.trode(i).avg,DaMa10.trode(i).err,'-r',1);
    hold on;
    plot([-100 window(2)],[0 0],'k--',[0 0],[-100 25],'k--'); axis([-100 window(2) -100 25]); title(RD_Active.MetaTags.ChannelID(i));        
end

% Calculate the latency to first lick
allCS = numel(DaMa10.events.cs)
DaMa10.events.fl = zeros(1, allCS);

for i=1:allCS
    currCSstamp = DaMa10.events.cs(i);
    flInd = find(DaMa10.events.al>currCSstamp+500 & DaMa10.events.al<currCSstamp+8000);

    if numel(flInd)>0
        latency = DaMa10.events.al(flInd(1)) - currCSstamp;
        DaMa10.events.fl(i) = latency;
    else
        DaMa10.events.fl(i) = 8000;
    end
end
figure(2); hist(DaMa10.events.fl,0:250:8000);

[vals, inds] = sort(DaMa10.events.fl,'ascend');
normTrials = find(vals<2250);
figure(3); clf; imagesc(DaMa10.trode(egTrode).wins(inds,:)); hold on; plot(vals,1:allCS,'k-','LineWidth',2);

% Get the amplitude of the LFP over the event windows
for i=1:8
    for j=1:allCS
        for k=1:3
            DaMa10.trode(i).amp(k,j) = trapz(DaMa10.trode(i).wins(j , window(1)+windows(k,1):window(1)+windows(k,2))) - trapz(DaMa10.trode(i).wins(j , 1:window(1)));        
        end 
        
        p = polyfit(windowSlope(1):windowSlope(2) , DaMa10.trode(i).wins(j , [windowSlope(1):windowSlope(2)] + window(1) ) , 1);
        DaMa10.trode(i).slope(1,j) = p(1);
                
        DaMa10.trode(i).winsBS(j,:) = DaMa10.trode(i).wins(j , :) - mean(DaMa10.trode(i).wins(j , 1:window(1)));
    end
    
    DaMa10.compiledSlopes(i,:) = DaMa10.trode(i).slope;
    
    % Compute the correlations between LFP amplitude and behavior
    [rho,pval] = corr( DaMa10.events.fl(inds(normTrials))' , DaMa10.trode(i).slope(1,inds(normTrials))' );
    DaMa10.trode(i).analysis.correl.pval = pval;
    DaMa10.trode(i).analysis.correl.rho = rho;
    DaMa10.trode(i).analysis.correl.p = polyfit( DaMa10.events.fl(inds(normTrials))' , DaMa10.trode(i).slope(1,inds(normTrials))' , 1 );
    
    if pval < 0.05
       disp(['The initial slope of the C90 component on electrode ' num2str(i) ' is correlated with approach latency (' num2str(pval) ')' ]); 
    else
       disp(['The initial slope of the C90 component on electrode ' num2str(i) ' is not correlated with approach latency' ]); 
    end
    
end

figure(3); clf; subplot(2,2,[1 3]); imagesc(DaMa10.trode(egTrode).winsBS(inds,:) , [-400 400]); hold on; plot(vals,1:allCS,'k-','LineWidth',2); title(['Example data for electrode ' num2str(egTrode)]);
colorbar; colormap(mapName);
% subplot(222); semilogx(DaMa10.events.fl , DaMa10.trode(egTrode).amp(3,:) , 'ko'); 
subplot(2,2,[2 4]); plot(DaMa10.events.fl , DaMa10.trode(egTrode).slope(1,:) , 'ko' , [250:2250] , polyval(DaMa10.trode(i).analysis.correl.p,[250:2250])); axis([0 2500 -15 5]);


% Compute the pairwise noise correlations between lfp channels to look for clustering

figure(5); imagesc(    corr(DaMa10.compiledSlopes') , [-1 1] ); colorbar;
colormap(mapName);

%% DA02

clear

egTrode = 5;
[mapName] = TNC_CreateRBColormap(20,'gp');

% parameter settings
window          = [100 2500];
windows(1,:)    = [5 20];
windows(2,:)    = [20 40];
windows(3,:)    = [55 150];

windowSlope     = [55 90];

% Load field potential data
RD_Active = TNC_LoadData(0, 0, 'DA02-090610-trace-ref14-001.ns4');
numChan = numel(RD_Active.MetaTags.ChannelID);

% Load behavioral data and calculate latency to sample the port
EV_Active = TNC_LoadData(0, 0, 'DA02-090610-trace-ref14-001.nev');

figure(1); clf;
% Calculate the windowed LFP data from the CS event times
DA02.events.cs = round(EV_Active.Data.Spikes.Timestamps(find(EV_Active.Data.Spikes.Electrode==137))./30);
DA02.events.al = round(EV_Active.Data.Spikes.Timestamps(find(EV_Active.Data.Spikes.Electrode==139))./30);

DA02.events

for i=1:numChan
    rawLFPdata = decimate(sgolayfilt(RD_Active.Data(i,:),3,21),30);
    [rCS] = TNC_ExtTrigWins(rawLFPdata , DaMa10.events.cs , window);
    DA02.trode(i).wins = rCS.wins;
    DA02.trode(i).avg = rCS.avg;
    DA02.trode(i).err = rCS.err;
end

for i=1:numChan
    figure(1); subplot(4,ceil(numChan./4),i);
    shadedErrorBar(-100:window(2),DA02.trode(i).avg,DA02.trode(i).err,'-r',1);
    hold on;
    plot([-100 window(2)],[0 0],'k--',[0 0],[-100 25],'k--'); axis([-100 window(2) -100 25]); title(RD_Active.MetaTags.ChannelID(i));        
end


%% deprecated
        
        
%         matForComp = zeros(numTrials,numel(unitList)+1);
%         matForComp(:,1) = abs( eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:)' ./ max(eLFPpopData.session(currSes).cont.trode(i).rAmp(3,:)) );
%     
%         for j=1:numel(unitList)
%             matForComp(:,j+1) = eLFPpopData.session(currSes).unit(unitList(j)).rAmp(1,:)' ./ max(eLFPpopData.session(currSes).unit(unitList(j)).rAmp(1,:));
%         end
%         
%         corrStruc       = corr(matForComp);
%         unitToLFPcorrs  = [unitToLFPcorrs corrStruc(2:numel(unitList)+1,1)'];
%         
%         figure(302);
%         i+32
%         unitList
%         corrStruc(2:numel(unitList)+1,1)
%         subplot(211);
%         imagesc(matForComp);
%         title(['Channel ' num2str(i+32)]);
%         subplot(212);
%         imagesc(corr(matForComp),[-1 1]);
%         drawnow; 
%         pause();


% numEvents = size(sink.wins,1);
% % Calculate the trialwise windows of the lfp
% channel = 9;
% [sink] = TNC_ExtTrigWins(eLFPpopData.session(currSes).cont.trode(channel).onek,eLFPpopData.session(currSes).events.cs,window);
% 
% figure(channel); clf;
% if showEachTrial
%     for i=1:numEvents
%         plot(-window(1):window(2),sink.wins(i,:)+(i.*spacing),'k');
%         hold on;
%     end
% else
%     subplot(5,1,1:4);
%     imagesc(sink.wins); axis tight;
%     subplot(5,1,5);
%     plot(abs(mean(sink.wins)),'LineWidth',2); axis tight;
% end
% drawnow;