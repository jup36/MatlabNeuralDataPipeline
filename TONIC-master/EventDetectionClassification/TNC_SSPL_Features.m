% _________________________________________________________________________
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / Janelia
% 
% BUG REPORTING: josh@dudmanlab.org
% FURTHER INFORMATION: www.dudmanlab.org
% Copyright (C) 2015 by Howard Hughes Medical Institute.
% _________________________________________________________________________
% _________________________________________________________________________
% 
% DETAILS: extract spike waveforms and features suitable for clustering with the TNC_SS_GUI interface
% 
% INPUTS:
% fileName: name of the ns5 file used from processing spike data
% arrayType: name of the array - depends upon function TNC_RemapElecPos
%  *** NOTE: if the array type string is not recognized all sites are processed independently.
% chunk: how many seconds of data should be processed per chunk
% dispOn: display online progress in figures [[[SLOW]]]
% targetName: base name of target files containg the featStruct (*_ft) and sessionStruct (*_ss) outputs
% maxSegs: if 0 process all chunks, otherwise process maxSegs number of chunks
% 
% EXAMPLE:
% [sessionStruct,featStruct] = TNC_SSPL_Features('15_05_27m28003.ns5','NN_h64_dhs',90,0,'g15_05_27m28003?,0);
% 

function [sessionStruct,featStruct] = TNC_SSPL_Features(fileName,arrayType,chunk,dispOn,targetName,maxSegs)

sampPrior = 25;

%% PRE-PROCESSING
    % GET INFO ABOUT THE CURRENT RECORDING
    Ns5DATA = openNSx('report',fileName);
    numSegs = ceil(Ns5DATA.MetaTags.Duration ./ chunk);
    disp(' ');
    disp(' ');
    disp(['Data contains ' num2str(numSegs) ' x ' num2str(chunk) ' seconds long segments.']);

    sessionStruct.chunk     = chunk;
    sessionStruct.arrayType = arrayType;
    sessionStruct.fileName  = fileName;
    if maxSegs==0
        sessionStruct.numSegs   = numSegs;
    else
        numSegs = maxSegs;
        sessionStruct.numSegs   = numSegs;        
    end
    disp([num2str(numSegs) ' segments to be processed.']);
    disp(' ');
    disp(' ');

    totChar = numel(fileName);
    
    if dispOn==0
        loadBehavior=0;
    else
        d = dir([fileName(1:totChar-4) '_bh.mat']);
        if numel(d)>0
            loadBehavior = 1;
        else 
            disp('No associated behavioral data was found.');
            loadBehavior = 0;
        end
    end
    
%% EXTRACT MULTIDIMENSIONAL EVENTS

    wfs = [];
    wts = [];
    chan = [];
    inds = [];
    vals = [];

    % LOAD BEHAVIORAL DATA TO DISPLAY WITH THE SPIKING DATA
    totChar = numel(fileName);
    ppBehavFileName = [fileName(1:totChar-4) '.mat']
    
    if loadBehavior
        eval(['load ' fileName(1:totChar-4) '.mat']);
        behaviorTrue = 1;
    else
        behaviorTrue = 0;
    end
    
    for k=1:numSegs

    % LOAD ALL CHANNEL DATA
        timeStr = ['t:' num2str((k-1)*chunk) ':' num2str(k*chunk)];
        disp(' ');disp(' ');disp(' ');disp(' ');disp(' ');disp(' ');
        disp(['Loading data over the range of ' timeStr]);
        disp(' ');
        Ns5DATA = openNSx('report','read',fileName,timeStr,'sec');
        
    % SMOOTH/FILTER DATA
        avgSeg      = median(Ns5DATA.Data,1);
%         numChan     = size(Ns5DATA.Data,1); 
        numChan     = numel(find(Ns5DATA.MetaTags.ChannelID<129))

    % SOME CHECKING OF DATA FOR PROBES WHERE MULTIPLE SITES ARE USED FOR SORTING
        if ~strcmp(arrayType,'SingleSites')
            if numChan~=64
                for zz=1:numChan
                    allRecChan(zz) = str2double(Ns5DATA.MetaTags.ElecLabel(zz,5:8));
                end
                
                if min(allRecChan)>64
                    allRecChan = allRecChan-64;
                end
            else
                allRecChan = 1:numChan;
            end
        else
            if k==1 % pre-allocate memory for the structure
                sessionStruct.seg(numSegs).shank(numChan).spikes    = [];
                featStruct.seg(numSegs).shank(numChan).col          = 0;
            end
        end

        [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);
        numRows = size(electrodeMatrix,1);
        
        if numRows>1

            %________________________________________________________
            %________________________________________________________
            %________________________________________________________
            %_______CODE FOR MULTISITE ARRAYS________________________
                    
                %=============================
                % GET THE ELECTRODE MAPPINGS
                [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);
                
                for j=1:size(electrodeMatrix,2) %shanks                                       
                    
                    wfs = [];
                    wts = [];
                    chan = [];
                    inds = [];
                    vals = [];
                    
                    spkData = zeros(size(electrodeMatrix,1),numel(Ns5DATA.Data(1,:)));

                    disp('_________________________________________________________________________');
                    disp('_________________________________________________________________________');
                    disp('_________________________________________________________________________');
                    disp('_________________________________________________________________________');
                    disp(['...BEGIN FILTERING DATA for seg ' num2str(k) ' of ' num2str(numSegs) ' | shank ' num2str(j)]);
                    disp(' ');

                    for i=1:size(electrodeMatrix,1) %electrodes
                        if numel( find(allRecChan==electrodeMatrix(i,j)) ) > 0

                            disp(['Filtering data on channel ' num2str(electrodeMatrix(i,j)) '|' num2str(find(allRecChan==electrodeMatrix(i,j))) ' ...']);
                            
                            rawData = Ns5DATA.Data(find(allRecChan==electrodeMatrix(i,j)),:)-avgSeg;
%                             [lowBandData,hiBandData] = TNC_FilterData(rawData,Ns5DATA.MetaTags.SamplingFreq,0,0,[1 0]);
%                             spkData(i,:) = hiBandData.values;                                                
            
                        %=============================
                        % USE A WAVELET FILTERING METHOD FOR SPEED
                            wname = 'db4'; 
                            maxlevel = 6; % 7.5kHz low pass
                            [c,l] = wavedec(rawData, maxlevel, wname);
                            c = wthcoef('a', c, l);
                            spkData(i,:) = sgolayfilt(waverec(c, l, wname),3,9);
                                                    
                        else
                            disp(['***ALERT*** Channel ' num2str(electrodeMatrix(i,j)) ' was not recorded in this segment...']);
                            spkData(i,:) = zeros(size(Ns5DATA.Data(1,:)));
                        end
                    end
                    disp('...completed');
                    
                    disp(' ');
                    disp(['Detecting events...']);
                    
                    for i=1:size(electrodeMatrix,1) % electrodes
                    
                        %=============================
                        % DETECT EVENTS ON THIS CHANNEL
                        if numel(find(allRecChan==electrodeMatrix(i,j)))==1
%                             [events] = TNC_SSPL_EventDetect(spkData(i,:),30,4,0);
                            [events] = TNC_SSPL_EventDetect(spkData(i,:),30,3.3,1);
                        else
                            disp(['***ALERT*** Channel ' num2str(electrodeMatrix(i,j)) ' has 0 events in this segment...']);
                            events.inds = [];
                        end
                        
                        %=============================
                        % ADD TO TOTAL LIST OF EVENTS
                        inds = [inds events.inds];
                        vals = [vals spkData(i,events.inds)];
                        chan = [chan ones(1,numel(events.inds)).*i];
                    
                        %=============================
                        % DISPLAY PHYS AND BEHAVIOR
                        if dispOn==1
                            disp('Displaying data for the current channel...');
                            if i==1
                                figure(1); clf;
                            end

                            figure(1);

                            subplot(10,1,1:size(electrodeMatrix,1));
                            if rem(j,2)==0
                                plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*size(electrodeMatrix,1))+(i-1)).*1000,'Color',[0 0.67 1]); hold on; axis off; axis tight;
                            else
                                plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*size(electrodeMatrix,1))+(i-1)).*1000,'Color',[1 0    0]); hold on; axis off; axis tight;
                            end

                            if behaviorTrue
                                subplot(10,1,9);
                                plot(ContData.behavior.sLeverV((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
                                axis off; axis tight;
                                subplot(10,1,10);
                                plot(ContData.behavior.rawLick((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
                                axis off; axis tight;
                            end
                            
                            drawnow;
                        end
                                            
                    end
                
                    %=============================
                    % FIND UNIQUE EVENTS ACROSS THE SHANK (ALL EVENT MUST BE A FIXED TIME APART)
                    [allInds,origIs,uIs]    = unique(inds);
                    allVals                 = vals(origIs);
                    allChan                 = chan(origIs);
                    
                    % all indices more than ~0.67 ms apart should be independent
                    valids = find(diff(allInds)>20);
                    valids = [ 1 , valids + 1]; % deal with the N-1 size of diff, first spike is by definition well separated

                    % all indices less than 1 ms apart need to be evaluated
                    invalids = find(diff(allInds)<=20);
                    invalids = invalids + 1; % deal with the N-1 size of diff
                    
                    %=============================
                    % RECONCILE THE INVALID SAMPLES
                    % first look for sequential invalids with no intervening valid indices (these should be treated as a group).
                    % for isolated invalid indices just need to compare those to the previously valid index, but it probably won't matter because the valid index will, by the definition of invalid, capture the invalid waveform.
                    m=2; startK=1; addValids=[];
                    while m < numel(invalids)
                        
                        if invalids(m)-invalids(m-1)==1
                            % in a sequence of short gaps with no
                            % intervening valid gaps
                            
                        else
                            endK = m-1;
                            
                            % evaluate the range invalids(startK)-1 (i.e. last valid) to invalids(endK) to find the
                            values = allVals(invalids(startK)-1:invalids(endK));
                            chans  = allChan(invalids(startK)-1:invalids(endK));
                            indices= allInds(invalids(startK)-1:invalids(endK));
                            lSpkLoc = find(values==min(values)); % look for the largest spike out of the set
                            if lSpkLoc == 1 %largest spike was at the valid point
                                if abs(chans(1)-chans(2)) > 2
                                  addValids = [addValids invalids(startK)];
                                end
                            else
                                addValids = [addValids invalids(lSpkLoc+startK-2)];
                            end
                            
                            startK = m;
                        end
                        
                        m=m+1;

                    end
                
                    eventInds = unique(allInds([valids addValids])); % use unique to sort the indices
                    disp(['Total events: ' num2str(size(eventInds))]);
                    
                    %=============================
                    % EXTRACT WAVEFORMS FROM ALL TRODES FOR EACH UNIQUE EVENT
                    disp(' '); disp('Extract spike waveforms...');
                    [spikes]        = TNC_SSPL_EventExtractME(spkData,eventInds,[sampPrior,45]);
                    
                    %=============================
                    % CALCULATE FEATURES FOR SPIKE EVENTS                   
                    if k==1
                        [features] = TNC_SSPL_EventQuantME(spikes.wfs, [], round(spikes.inds./30), sampPrior, size(electrodeMatrix,1));
                        pcaStruct = features.pca;
                    else
                        [features] = TNC_SSPL_EventQuantME(spikes.wfs, pcaStruct, round(spikes.inds./30), sampPrior, size(electrodeMatrix,1));
                    end                 
                    
                    %=============================
                    % UPDATE STRUCTURE FOR SPIKE SORTING
                    disp('4) Updating the feature structure with new data...');
                        
                        featStruct.seg(k).shank(j).inds     = spikes.inds + ((k-1).*chunk.*30000);
                        featStruct.seg(k).shank(j).ts       = round(featStruct.seg(k).shank(j).inds./30);
                        featStruct.seg(k).shank(j).id       = zeros(numel(spikes.inds),1);
                        featStruct.seg(k).shank(j).params   = features.params;
                        
                        if k==1 && j==1
                            featStruct.paramNames = features.paramNames;
                            featStruct.rSp = rSp;
                            featStruct.cSp = cSp;
                            featStruct.electrodeMatrix = electrodeMatrix;
                            
                            sessionStruct.resolution = spikes.resolution;
                            sessionStruct.winL       = spikes.winL;
                            sessionStruct.winR       = spikes.winR;                        
                            
                        end

                        sessionStruct.seg(k).shank(j).wfs = spikes.wfs;
                        sessionStruct.seg(k).shank(j).inds = spikes.inds;
                    disp('_________________________________________________________________________');
                    disp(' ');
                    disp(' ');                    
                    disp(' ');
                    
                end

            %________________________________________________________
            %________________________________________________________
            %________________________________________________________
            %_______CODE FOR ISOLATED SITE ARRAYS____________________

        else
                
                for eInd=1:numChan % shanks | sites
                    
                    %=============================
                    % INITIALIZE ARRAYS
                    wfs = [];
                    wts = [];
                    chan = [];
                    inds = [];
                    vals = [];
                    pcaStruct = [];                    
                    spkData = zeros(1,numel(Ns5DATA.Data(1,:)));


                    %=============================
                    % UPDATE USER
                    disp(' ');
                    disp(' ');
                    disp('_________________________________________________________________________');
                    disp(['...BEGIN FILTERING DATA for seg ' num2str(k) ' of ' num2str(numSegs) ' | site ' num2str(eInd) ' of ' num2str(numChan)]);

                    %=============================
                    % FILTER DATA ON THIS CHANNEL
                    disp(['1) Filtering data on channel ' num2str(eInd) ' | shank:' num2str(eInd)]);

                        rawData = Ns5DATA.Data(eInd,:)-avgSeg;

                        % IDEAL ZERO PHASE FILTERING
                        % [lowBandData,hiBandData] = TNC_FilterData(rawData,Ns5DATA.MetaTags.SamplingFreq,0,0,[1 0]);
                        % spkData = hiBandData.values;
            
                        % USE A WAVELET FILTERING METHOD FOR SPEED
                        wname = 'db4'; 
                        maxlevel = 6; % 7.5kHz low pass
                        [c,l] = wavedec(rawData, maxlevel, wname);
                        c = wthcoef('a', c, l);
                        spkData = sgolayfilt(waverec(c, l, wname),3,9);
                            
                    %=============================
                    % EXTRACT WAVEFORMS FROM ALL TRODES FOR EACH UNIQUE EVENT
                    disp(['2) Detecting and extracting events...']);

%                         [events] = TNC_SSPL_EventDetect(spkData,30,4,1);
                        [events] = TNC_SSPL_EventDetect(spkData,30,4,1);

                        if numel(events.ts)>10
                            [spikes]        = TNC_SSPL_EventExtract(spkData,rawData,events.inds,[15,45],dispOn);
                                                
                    %=============================
                    % CALC SCALAR QUANT FOR SPIKES
                            disp(['3) Quantify spikes...']);

                            if k==1
                                [features]   = TNC_SSPL_EventQuantSE(spikes,[]);
                                templates(eInd).template = features.template;
                            else
                                [features]   = TNC_SSPL_EventQuantSE(spikes,templates(eInd).template);
                            end

                        else
                            
                            events.inds = [];
                            features.params = [];
                            
                        end
                        
                    %=============================
                    % UPDATE STRUCTURE FOR SPIKE SORTING
                    disp('4) Updating the feature structure with new data...');
                        
                        featStruct.seg(k).shank(eInd).inds     = spikes.inds + ((k-1).*chunk.*30000);
                        featStruct.seg(k).shank(eInd).id       = zeros(numel(spikes.inds),1);
                        featStruct.seg(k).shank(eInd).ts       = round(featStruct.seg(k).shank(eInd).inds./30)';
                        featStruct.seg(k).shank(eInd).params   = features.params;
                        
                        if k==1 && eInd==1
                            [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);
                            featStruct.paramNames = features.paramNames;
                            featStruct.rSp = rSp;
                            featStruct.cSp = cSp;
                            featStruct.electrodeMatrix = electrodeMatrix;

                            sessionStruct.resolution = spikes.resolution;
                            sessionStruct.winL       = spikes.winL;
                            sessionStruct.winR       = spikes.winR;
                        end

                        sessionStruct.seg(k).shank(eInd).wfs        = spikes.wfs;
                        sessionStruct.seg(k).shank(eInd).inds       = spikes.inds;

                    disp('_________________________________________________________________________');
                    disp(' ');
                    disp(' ');                    
                    disp(' ');
                    
                end

        end
        
    end

%% SAVE FEATURE STRUCTURE

% clear paramVals featStruct wfStruct newValues
offset                  = 0;
featStruct.chunk        = chunk;
featStruct.arrayType    = arrayType;

disp(['save ~/' targetName '_ft featStruct']); 
save([targetName '_ft.mat'],'featStruct');

%% SAVE SESSION STRUCTURE
 
disp(['save ~/' targetName '_ss sessionStruct']);
save([targetName '_ss.mat'],'sessionStruct');

%% DEPRECATED CODE

% BUILD THE WF STRUCTURE 
% 
% for j=1:size(sessionStruct.seg(1).shank,2)
%         
%     for i=1:sessionStruct.numSegs
%         
%         wfdata.seg(i).wfs(numel(sessionStruct.seg(i).shank(j).spikes.wfs)).values = zeros(8,61);
% 
%         for k=1:numel(sessionStruct.seg(i).shank(j).spikes.wfs)
%             wfdata.seg(i).wfs(k).values = int16(sessionStruct.seg(i).shank(j).spikes.wfs(k).values);
%         end
%     end
%     
%     % For getting waveform data 
%     wfdata.chunk    = chunk;
%     wfdata.shank    = j;
%     wfdata.filename = fileName;
%     wfdata.datatype = 'int16';
%     
%     disp(['save ' fileName(1:totChar-4) '_shank' num2str(j) '_wf wfdata']);
%     eval(['save ' fileName(1:totChar-4) '_shank' num2str(j) '_wf wfdata']);
% 
%     clear wfdata    
%     
% end


% 
%                         % DISPLAY PHYS AND BEHAVIOR
%                         if dispOn==1
%                             disp('Displaying data for the current channel...');
%                             if i==1
%                                 figure(1); clf;
%                             end
% 
%                             figure(1);
% 
%                             subplot(10,1,1:8);
%                             if rem(j,2)==0
%     %                             plot(spkData(i,:)+((((j-1).*8)+(i-1)).*1000),'Color',[0 0.67 1]); hold on; axis tight;
%                                 plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*8)+(i-1)).*1000,'Color',[0 0.67 1]); hold on; axis off; axis tight;
%                             else
%     %                             plot(spkData(i,:)+((((j-1).*8)+(i-1)).*1000),'Color',[1 0 0]); hold on; axis tight;
%                                 plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*8)+(i-1)).*1000,'Color',[1 0    0]); hold on; axis off; axis tight;
%                             end
% 
%                             subplot(10,1,9);
%                             plot(ContData.behavior.sLeverV((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
%                             axis off; axis tight;
%                             subplot(10,1,10);
%                             plot(ContData.behavior.rawLick((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
%                             axis off; axis tight;
%                                 
%                             drawnow;
%                         end