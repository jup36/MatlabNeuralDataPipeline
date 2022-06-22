% _________RUNTONIC_________

%% INITIALIZATION 
disp('  ');
disp('  ');
disp('  ');
disp('----------------------------------------------');
disp('  ');
disp(['RunTONIC initialized at ' datestr(now)]);
disp('  ');
disp(pwd);
disp('  ');
disp('----------------------------------------------');
disp('  ');

% clear memory
clear;
% close all figures

RunTonicInfo.loadedFiles = 0;

%% LOAD DATA FROM OVATION DATABASE INTO MEMORY

%% LOAD DATA FROM FILE INTO MEMORY

loadWholeDir = 0; % USE THIS FLAG TO SWITCH MODES BETWEEN LOADING ALL FILES IN DIRECTORY (=1) AND ONLY LOADING SINGLE FILES (=0)
[fname,pname]= uigetfile('*.*');

if fname == 0
    return
end

if loadWholeDir == 1
    allFiles    = dir('*.*'); % creates a structure, use the "name" field
    totalFiles  = size(allFiles,1);
else
    totalFiles = 1;
end

RunTonicInfo.activeSet = RunTonicInfo.loadedFiles;
RunTonicInfo.loadedFiles = RunTonicInfo.loadedFiles+totalFiles;
RunTonicInfo.path = pname;

for i=1:totalFiles
    
    if loadWholeDir == 1
        fileNameStr = [pname allFiles(i).name];
        tmpName = fname(1:strfind(allFiles(i).name,'.')-1);
    else
        fileNameStr = [pname fname];
        tmpName = fname(1:strfind(fname,'.')-1);
    end

    % remove '-'s from the name
    nonDashInds = strfind(tmpName,'-');
    tmpName(nonDashInds) = '_';
    nameInMem = tmpName;
    
    RunTonicInfo.activeSet = RunTonicInfo.activeSet+1;
    RunTonicInfo.fileNames(RunTonicInfo.activeSet).name     = nameInMem;
    RunTonicInfo.fileNames(RunTonicInfo.activeSet).loadTime = datestr(now);
    eval([nameInMem ' = TNC_LoadData(0, 0, fileNameStr);']);
    
end

%% CONVERT DATA STRUCTURES TO STANDARD TONIC FORMATING

%% FILTER, DOWNSAMPLE AND EXTRACT LFP RESPONSES

%% SORT SPIKES AND CREATE UNITS FROM SINGLE CHANNEL

%% SORT SPIKE IMAGES FROM NS5 DATA

chanArray = [];
candidateStamps = [];

% detect events (retrieve timestamps and the event image)
    for i=1:length(chanArray)
        [events] = TNC_EventDetect(dataSource,snrThresh);
        candidateStamps = [candidateStamps,events.ts];
    end

% reconcile events (remove overlaps)
    [recStamps] = TNC_EventReconcile(candidateStamps,minEventTime);

% extract event images
    [eventImages] = TNC_EventExtractImage(dataSource,recStamps,chanArray,windowPnts);
    
% produce scalar quantification of the event images
    [chVector,chScalar] = TNC_EventQuant_MCV(image,imgCompress,projection);

% cluster the quantified events
    [events] = TNC_EventCluster(events,dMethod,cMethod);


%% CREATE A STANDARD NAME FOR THE ACTIVE DATASTRUCTURE TO USE AND REPORT THE FIELDS PRESENT

RunTonicInfo.activeSet = 1; % User should select the data set to examine with this variable

propUA = []; propCS = 0; propEL = 0;
eval(['RD_Active =' RunTonicInfo.fileNames(RunTonicInfo.activeSet).name ';']);
disp(' ');
disp('----------------------------------------------');
disp(['ACTIVE STRUCTURE :: ' RunTonicInfo.fileNames(RunTonicInfo.activeSet).name])
disp(['SOURCE FILE      :: ' RD_Active.fileName])
disp('----------------------------------------------');
for i = 1:size(RD_Active.neurons,1)
    
    if findstr(RD_Active.waves{i}.name,'elec') == 1% electrode
        testi = findstr(RD_Active.waves{i}.name,'U');
        if size(testi,1) == 0
            propUA = [propUA,i];
        end
    else % analog input
        if findstr(RD_Active.waves{i}.name,'ainp9a')==1
            propCS = i;
        end
        if findstr(RD_Active.waves{i}.name,'ainp10a')==1
            propUS = i;
        end
    end
    
    if i < 10
        disp([num2str(i) '....' RD_Active.waves{i}.name]); 
    else
        disp([num2str(i) '...' RD_Active.waves{i}.name]);         
    end
end
disp('----------------------------------------------');
disp('Proposed event and unit pins');
disp('----------------------------------------------');
disp(['Unit Array (non-U) : ' num2str(propUA)]);
disp(['Cond Stim (ainp9)  : ' num2str(propCS)]);
disp(['Every Lick (ainp10): ' num2str(propUS)]);
disp(['Cond Stim (dig)    : ' '-0007']);
disp(['UCond Stim (dig)   : ' '-0006']);
disp('----------------------------------------------');
disp(' ');

%% BEHAVIOR ANALYSIS: EXTRACT TIMESTAMPS AND STANDARD ALIGNMENTS FOR ACTIVE DATASTRUCTURE
% function details for reference:
% [tonicDataStructure] = TNC_ReadToStdRecStruct(dataStructure,unitArrayToLoad,CSid,USid,ELid)

% what data fields should be loaded?
% propUA = [3 17 18 19 20 21];    % use previously detected
% propUA = [3 4 17 18 19 20 21];    % use previously detected
% propCS = ;            % use previously detected
digUS = '-0006';
digCS = '-0007';
digEL = '-0001';
% propEL = ;            % use previously detected

% Define the smoothing to apply to the timeseries data
smthParams.rise     = 1;
smthParams.decay    = 50;
% [kernel]  = TNC_CreateCausalKernel(smthParams.rise,smthParams.decay,1);
[kernel]  = TNC_CreateGaussian(smthParams.decay.*4,smthParams.decay,smthParams.decay.*8,1);

winParams.prior     = 5e3;
winParams.after     = 5e3;
exFlag =0;

% digFlag is a 3 place vector that gives the datatype for CS, US, and EL (1==digital)
[TD_Active] = TNC_ReadToStdRecStruct(RD_Active,propUA,propCS,propUS,digEL,kernel,winParams,[0,0,1],extFlag)

TD_Active.local     = RunTonicInfo.fileNames(RunTonicInfo.activeSet).name;
TD_Active.source    = RD_Active.fileName;

%% BEHAVIOR ANALYSIS: TAG DATA STRUCTURE
% Provide some useful set of tags that can allow Ovation import later and
% also allow me to search this data again. Perhaps provide a reference to
% an analysis log that keeps command histories and other details and is
% searchable (i.e. electronic)

TD_Active.tags.reward.probDist = 'Uniform';  
TD_Active.tags.reward.probMean = 0
TD_Active.tags.reward.timeDist = 'Fixed'; 
TD_Active.tags.reward.timeMean = 2000;
TD_Active.tags.reward.timeVar  = 0;
TD_Active.tags.collector       = 'WXP';
TD_Active.tags.analyzer        = 'JTD';
TD_Active.tags.mouseID         = 'DA-MA-05';
TD_Active.tags.dateTime        = [2010,01,07];

TD_Active.tags

disp('Finished running __TAG DATA STRUCTURE__');

%% CALCULATE INTERSPIKE INTERVAL DISTRIBUTION

kMax = size(TD_Active.unit,2);

figure(10); clf; clear isiCorrMat;
subplot(1,3,1:2);

colorGradient = ([0.043, 0.518, 0.78]-[0.847, 0.161, 0.00]);

for k=1:kMax
    
    [isi] = TNC_QuantISI(TD_Active.unit(k).ts);
    TD_Active.unit(k).isi = isi;

    semilogx(TD_Active.unit(k).isi.hist.logTimes,TD_Active.unit(k).isi.hist.logCount,'Color',[0.043, 0.518, 0.78]-(colorGradient.*(k./kMax)),'LineWidth',2);
    drawnow; hold on;
    isiCorrMat(:,k) = TD_Active.unit(k).isi.hist.logCount';
end
[mapName] = TNC_CreateRBColormap(1024,'rb');
colormap(mapName);

subplot(1,3,3);
colormap(mapName);
imagesc(corr(isiCorrMat));

disp('Finished running __INTERSPIKE INTERVAL DISTRIBUTION__');


%% SHOW WAVEFORM DATA
kMax = size(TD_Active.unit,2);
numSamps = size(TD_Active.unit(1).wf,1);

figure(13); clf;

for k=1:kMax
    subplot(1,kMax,k);
    TD_Active.unit(k).WfMean = mean(TD_Active.unit(k).wf,2);
    TD_Active.unit(k).WfStd = std(TD_Active.unit(k).wf,0,2);
    plot(1:numSamps,TD_Active.unit(k).WfMean,'k','LineWidth',2);
    hold on; plot(1:numSamps,TD_Active.unit(k).WfMean+TD_Active.unit(k).WfStd,'k--',1:numSamps,TD_Active.unit(k).WfMean-TD_Active.unit(k).WfStd,'k--');
%     axis([0 numSamps -750 500]);
    axis off;
end


%% TRIAL CLASSIFIERS

sortType = '';
totalTrials = size(TD_Active.events.CS.ts,1);

switch sortType
    
    case 'interCurrent'
        tmpIntervals(1,1)           = 20e3;
        tmpIntervals(2:totalTrials) = TD_Active.events.CS.ts(2:totalTE__Drug_rials) - TD_Active.events.CS.ts(1:totalTrials-1);
        [y,sortInds] = sort(tmpIntervals,1,'ascend');

    case 'intraCurrent'
        tmpIntervals = TD_Active.events.US.ts - TD_Active.events.CS.ts;
        [y,sortInds] = sort(tmpIntervals,1,'ascend');

    case 'intraLast'
        tmpIntervals = TD_Active.events.US.ts - TD_Active.events.CS.ts;
        [y,sortInds] = sort(tmpIntervals,1,'ascend');
        sortInds = sortInds+1;
        sortInds(find(sortInds>totalTrials),1) = 1;

    case 'firstLickCS'
        for m=1:totalTrials
            firstLickCStimes(m) = find(TD_Active.events.EL.ts>TD_Active.events.CS.ts(m),1);
        end
        tmpIntervals = firstLickCStimes - TD_Active.events.CS.ts;
        [y,sortInds] = sort(tmpIntervals,1,'ascend');

    case 'firstLickUS'
        for m=1:totalTrials
            firstLickUStimes(m) = find(TD_Active.events.EL.ts>TD_Active.events.US.ts(m),1);
        end
        tmpIntervals = firstLickUStimes - TD_Active.events.US.ts;
        [y,sortInds] = sort(tmpIntervals,1,'ascend');

    otherwise
        disp('Just putting trials in standard order');
        sortInds = 1:totalTrials;

end

TD_Active.sorting.currSortInds = sortInds;
TD_Active.sorting.currSortType = sortType;

disp('Finished running __TRIAL CLASSIFIERS__');

%% CREATE DERIVED ACTIVITY VARIABLES
% Z-score integral over a time window
% peak rate
% spiking precision (i.e. height of a psth compared to gaussian background)

kMax = size(TD_Active.unit,2);

for k=1:kMax
    
    [spkCntPerTrialCS] = TNC_QuantSpksPerTrial(1,TD_Active.unit(k).responseCS.image.aligned,[winParams.prior,winParams.prior+0.4e3]);
    TD_Active.unit(k).spkCnt.spkCntPerTrialCS = spkCntPerTrialCS;

    [spkCntPerTrialCSS] = TNC_QuantSpksPerTrial(1,TD_Active.unit(k).responseCSS.image.aligned,[winParams.prior,winParams.prior+0.4e3]);
    TD_Active.unit(k).spkCnt.spkCntPerTrialCSS = spkCntPerTrialCSS;

    if TD_Active.USexists == 1 

        [spkCntPerTrialUS] = TNC_QuantSpksPerTrial(1,TD_Active.unit(k).responseUS.image.aligned,[winParams.prior,winParams.prior+0.4e3]);
        TD_Active.unit(k).spkCnt.spkCntPerTrialUS = spkCntPerTrialUS;

        [spkCntPerTrialUSS] = TNC_QuantSpksPerTrial(1,TD_Active.unit(k).responseUSS.image.aligned,[winParams.prior,winParams.prior+0.4e3]);
        TD_Active.unit(k).spkCnt.spkCntPerTrialUSS = spkCntPerTrialUSS;
    end

    % build the correlation matrix
    if k==1
        forCorrCS = spkCntPerTrialCS;
        forCorrCSS = spkCntPerTrialCSS;
        if TD_Active.USexists == 1 
            forCorrUS = spkCntPerTrialUS;     
            forCorrUSS = spkCntPerTrialUSS;     
        end
    else
        forCorrCS = [forCorrCS,spkCntPerTrialCS];     
        forCorrCSS = [forCorrCSS,spkCntPerTrialCSS];     
        if TD_Active.USexists == 1 
            forCorrUS = [forCorrUS,spkCntPerTrialUS];     
            forCorrUSS = [forCorrUSS,spkCntPerTrialUSS];     
        end
    end
    
    % take the histogram of the spike count for use in ROC analysis
    TD_Active.unit(k).spkCnt.hCntCS = hist(spkCntPerTrialCS,1:1:50);
    if TD_Active.USexists == 1 
        TD_Active.unit(k).spkCnt.hCntUS = hist(spkCntPerTrialUS,1:1:50);
    end
end

[corrCS.rho,corrCS.p] = corr(forCorrCS);
[corrCSS.rho,corrCSS.p] = corr(forCorrCSS);
if TD_Active.USexists == 1 
    [corrUS.rho,corrUS.p] = corr(forCorrUS);
    [corrUSS.rho,corrUSS.p] = corr(forCorrUSS);
end

TD_Active.derived.spkCnt.corrCS = corrCS;
TD_Active.derived.spkCnt.corrCSS = corrCSS;

if TD_Active.USexists == 1 
    TD_Active.derived.spkCnt.corrUS = corrUS;
    TD_Active.derived.spkCnt.corrUSS = corrUSS;
end

[mapName] = TNC_CreateRBColormap(1024,'rb');
figure(9); colormap(mapName);
imagesc(TD_Active.derived.spkCnt.corrCS.rho,[-1 1])

disp('Finished running __CREATE DERIVED ACTIVITY VARIABLES__');
clear forCorr* corr* spkCnt*

%% KEEP A RUNNING UPDATE OF THE MEAN PSTH FOR EACH CELL

learning = 0;

if exist('PSTHrespCS')
    totalCScells = size(PSTHrespCS.cells,2)
    totalUScells = size(PSTHrespUS.cells,2)
else
    totalCScells = 0;
    totalUScells = 0;
end

kMax = size(TD_Active.unit,2);

for k=1:kMax

        if learning
            PSTHrespCS.cells(totalCScells+k).alignedTime(1,:) = TD_Active.unit(k).responseCS.image.alignedtime;
            PSTHrespCS.cells(totalCScells+k).respCS(1,:) = TD_Active.unit(k).responseCS.image.psth;
            PSTHrespCS.cells(totalCScells+k).respCSS(1,:) = TD_Active.unit(k).responseCSS.image.psthZ;
    
            PSTHrespUS.cells(totalUScells+k).alignedTime(1,:) = TD_Active.unit(k).responseUS.image.alignedtime;
            PSTHrespUS.cells(totalUScells+k).respUS(1,:) = TD_Active.unit(k).responseUS.image.psth;
            PSTHrespUS.cells(totalUScells+k).respUSS(1,:) = TD_Active.unit(k).responseUSS.image.psthZ;
        else
            if learnAdd>0
                PSTHrespCS.cells(totalCScells-learnAdd+k).respCS(2,:) = TD_Active.unit(k).responseCS.image.psth;
                PSTHrespCS.cells(totalCScells-learnAdd+k).respCSS(2,:) = TD_Active.unit(k).responseCSS.image.psthZ;            
            end
            
        end
        
        disp(k);

end

if learning
    learnAdd = kMax
else
    learnAdd = 0
end

%% DISPLAY ALL PSTHS COLLECTED SO FAR
USplot = 1
kMax = size(PSTHrespCS.cells,2)

figure(11);
clf;
hold on;

for k=1:kMax
    if USplot
        plot(PSTHrespUS.cells(k).alignedTime(1,:),PSTHrespUS.cells(k).respUSS(1,:),'Color',[0.043, 0.518, 0.78]-(colorGradient.*(k./kMax)),'LineWidth',2);    
    else
        plot(PSTHrespCS.cells(k).alignedTime(1,:),PSTHrespCS.cells(k).respCSS(1,:),'Color',[0.043, 0.518, 0.78]-(colorGradient.*(k./kMax)),'LineWidth',2);    
    end

    corrMatUS(:,k) = PSTHrespUS.cells(k).respUSS(1,:)';
    corrMatCS(:,k) = PSTHrespCS.cells(k).respCSS(1,:)';

end

figure(12);
clf; colormap(mapName);
subplot(121);
imagesc(cov(corrMatCS),[-1 1]);
subplot(122);
imagesc(cov(corrMatUS),[-1 1]);

%% SEGMENTED TRIAL DATA
kMax = size(TD_Active.unit,2);
trialsPerSegment = round(size(TD_Active.unit(k).responseCS.image.aligned,1)./5)

% Need to decide what kind of kernel is appropriate here (causal, gaussian, ?)
[kernel] = TNC_CreateGaussian(150,10,300,1);
% [kernel] = TNC_CreateCausalKernel(1,10,1)

for k=1:kMax
    
    [segResponseCS] = TNC_QuantSegmentedPSTH(TD_Active.unit(k).responseCS.image.aligned,trialsPerSegment,kernel);
    TD_Active.unit(k).segResponseCS = segResponseCS;
    clear segResponseCS;

    if TD_Active.USexists == 1 
        [segResponseUS] = TNC_QuantSegmentedPSTH(TD_Active.unit(k).responseUS.image.aligned,trialsPerSegment,kernel);
        TD_Active.unit(k).segResponseUS = segResponseUS;
        clear segResponseUS;
    end

end

disp('Finished running __SEGMENTED TRIAL DATA__');

%% VISUALIZATION OF BASIC TRIAL STRUCTURE DATA
[mapName] = TNC_CreateRBColormap(1024,'rb');
numTrials = size(TD_Active.events.CS.ts,1); 
numUnits = size(TD_Active.unit,2);


if TD_Active.USexists == 1 
    figure(1);
else
    figure(3);
end

colormap(mapName);
axisVector = [winParams.prior-0.25e3 winParams.prior+2e3 1 numTrials];

for j = 1:numUnits
    subplot(4,numUnits,[j,j+numUnits,j+(2.*numUnits)]);   
    imagesc(TD_Active.unit(j).responseCSS.image.aligned(sortInds,:));
    title(['Unit ' num2str(j) ' aligned to CS']);
    axis(axisVector);

    subplot(4,numUnits,j+(3.*numUnits)); 
    plot(TD_Active.unit(j).responseCSS.image.psthZ,'LineWidth',2);
%     title(['PSTH of Unit ' num2str(j) ' aligned to CS']);
    axis([winParams.prior-0.25e3 winParams.prior+2e3 -5 10]);
end

if TD_Active.USexists == 1 
    figure(2);
    colormap(mapName);
    axisVector = [winParams.prior-2e3 winParams.prior+4e3 1 numTrials]; 
    for j = 1:numUnits
        subplot(4,numUnits,[j,j+numUnits,j+(2.*numUnits)]);   
        imagesc(TD_Active.unit(j).responseUSS.image.aligned(sortInds,:));
        title(['Unit ' num2str(j) ' aligned to US']);
        axis(axisVector);

        subplot(4,numUnits,j+(3.*numUnits)); 
        plot(TD_Active.unit(j).responseUSS.image.psthZ,'LineWidth',2);
        title(['PSTH of Unit ' num2str(j) ' aligned to US']);
        axis([winParams.prior-2e3 winParams.prior+4e3 -5 10]);
    end
end

%% POPULATION VECTOR PLOTS AND ANALYSIS
% my standard color scheme...
% cyanish : [0.043, 0.518, 0.78]
% deep red: [0.847, 0.161, 0.00]

figure(4); clf;
unitArray = [4,7,10];
plot3(TD_Active.unit(unitArray(1)).responseCSS.image.psthZ',TD_Active.unit(unitArray(2)).responseCSS.image.psthZ',TD_Active.unit(unitArray(3)).responseCSS.image.psthZ','Color',[0.043, 0.518, 0.78],'LineWidth',2);
hold on
if TD_Active.USexists == 1 
    plot3(TD_Active.unit(unitArray(1)).responseUSS.image.psthZ',TD_Active.unit(unitArray(2)).responseUSS.image.psthZ',TD_Active.unit(unitArray(3)).responseUSS.image.psthZ','Color',[0.847, 0.161, 0.00],'LineWidth',2);
end
patch([-2,-2,10,10],[-2,10,10,-2],[0 0 0 0],[0.8 0.8 0.8],'FaceAlpha',0.75);
grid on
title([TD_Active.local ' | units: ' num2str(unitArray)]);
xlabel(['unit ' num2str(unitArray(1))]);
ylabel(['unit ' num2str(unitArray(2))]);
zlabel(['unit ' num2str(unitArray(3))]);

figure(4); clf;
plot(TD_Active.unit(unitArray(1)).responseCSS.image.psthZ,'Color',[0.043, 0.518, 0.78],'LineWidth',2);
hold on
axisVector = [winParams.prior-0.1e3 winParams.prior+1e3 -3 12];
if TD_Active.USexists == 1 
    plot(TD_Active.unit(unitArray(1)).responseUSS.image.psthZ,'Color',[0.847, 0.161, 0.00],'LineWidth',2);
end
axis(axisVector);

%% DISPLAY OF SEGMENTED TRIAL DATA
figure(5);
[mapName] = TNC_CreateRBColormap(1024,'wblue');
colormap(mapName);
axisVector = [0 2750 1 size(TD_Active.unit(j).segResponseCS.psthS,1) 0 2.5]; 
numUnits = size(TD_Active.unit,2);

for j = 1:numUnits
    subplot(1,numUnits,j);   
    h = waterfall(TD_Active.unit(j).segResponseUS.psthS(:,winParams.prior-1.5e3:winParams.prior+0.5e3));
    title(['Unit ' num2str(j) ' aligned to US']);
    axis tight;
    grid off
    box off
end

%% COMPILE MEAN DATA ACROSS STRUCTURES AND DATAFILES

%% LOAD AND EXTRACT VIDEO DATA

%% ALIGN VIDEO DATA TO PHYSIOLOGY DATA USING RECORDED TRIGGERS

%% MOVE ACTIVE DATASTRUCTURE TO A NAME SUITABLE FOR PERSISTING
disp(['Executing:    ' 'TD_' RunTonicInfo.fileNames(RunTonicInfo.activeSet).name ' = TD_Active;']);
eval(['TD_' RunTonicInfo.fileNames(RunTonicInfo.activeSet).name ' = TD_Active;']);
disp(['Saving:    ' '~/Data_Analyzed_Matlab/TONIC_DATA/TD_' RunTonicInfo.fileNames(RunTonicInfo.activeSet).name]);
tmpFileName      = ['~/Data_Analyzed_Matlab/TONIC_DATA/TD_' RunTonicInfo.fileNames(RunTonicInfo.activeSet).name];
tmpDataStructure = ['TD_' RunTonicInfo.fileNames(RunTonicInfo.activeSet).name];
save(tmpFileName, tmpDataStructure);
ls -FGoS ~/Data_Analyzed_Matlab/TONIC_DATA

%% PERSIST THE DATA STRUCTURE IN MAT FORMAT OR THE OVATION DATABASE
