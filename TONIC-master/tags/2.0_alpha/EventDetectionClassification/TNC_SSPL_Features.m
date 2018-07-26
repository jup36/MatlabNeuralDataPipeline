function [sessionStruct,featStruct] = TNC_SSPL_Features(fileName,arrayType,chunk,dispOn,targetName)

%% PRE-PROCESSING
    % GET INFO ABOUT THE CURRENT RECORDING
    Ns5DATA = openNSx('report',fileName);
    numSegs = ceil(Ns5DATA.MetaTags.Duration ./ chunk);
    disp(' ');
    disp(' ');
    disp(['Data will be loaded and processed as ' num2str(numSegs) ' x ' num2str(chunk) ' seconds long segments.']);
    disp(' ');
    disp(' ');

    sessionStruct.chunk     = chunk;
    sessionStruct.arrayType = arrayType;
    sessionStruct.fileName  = fileName;
    sessionStruct.numSegs   = numSegs;

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
        numChan     = size(Ns5DATA.Data,1);        

    % SOME CHECKING OF DATA FOR PROBES WHERE MULTIPLE SITES ARE USED FOR SORTING
        if ~strcmp(arrayType,'SingleSites')
            if numChan~=64
                for zz=1:numChan
                    allRecChan(zz) = str2double(Ns5DATA.MetaTags.ElecLabel(zz,5:8));
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

        switch arrayType

            %________________________________________________________
            %________________________________________________________
            %________________________________________________________
            %_______CODE FOR BUZSAKI 64 ARRAY________________________

            case 'NN_b64'
                    
                % GET THE ELECTRODE MAPPINGS
                [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);
                
                for j=1:8 %shanks                                       
                    
                    wfs = [];
                    wts = [];
                    chan = [];
                    inds = [];
                    vals = [];
                    
                    spkData = zeros(8,numel(Ns5DATA.Data(1,:)));

                    disp('_________________________________________________________________________');
                    disp('_________________________________________________________________________');
                    disp('_________________________________________________________________________');
                    disp('_________________________________________________________________________');
                    disp(['...BEGIN FILTERING DATA for seg ' num2str(k) ' of ' num2str(numSegs) ' | shank ' num2str(j)]);
                    disp(' ');

                    for i=1:8 %electrodes
                        if numel( find(allRecChan==electrodeMatrix(i,j)) ) > 0

                            disp(['Filtering data on channel ' num2str(electrodeMatrix(i,j)) ' ...']);
                            
                            rawData = sgolayfilt(Ns5DATA.Data(electrodeMatrix(i,j),:)-avgSeg,11,21);
%                             [lowBandData,hiBandData] = TNC_FilterData(rawData,Ns5DATA.MetaTags.SamplingFreq,0,0,[1 0]);
%                             spkData(i,:) = hiBandData.values;                                                
            
                            % USE A WAVELET FILTERING METHOD FOR SPEED
                            wname = 'db4'; 
                            maxlevel = 5; % 7.5kHz low pass
                            [c,l] = wavedec(rawData, maxlevel, wname);
                            c = wthcoef('a', c, l);
                            spkData(i,:) = waverec(c, l, wname);
                        
                        else
                            disp(['***ALERT*** Channel ' num2str(electrodeMatrix(i,j)) ' was not recorded in this segment...']);
                            spkData(i,:) = zeros(size(Ns5DATA.Data(1,:)));
                        end
                    end
                    disp('...completed');
                    
                    disp(' ');
                    disp(['Detecting events...']);
                    
                    for i=1:8
                    
                        % DETECT EVENTS ON THIS CHANNEL
                        if numel(find(allRecChan==electrodeMatrix(i,j)))==1
                            [events] = TNC_EventDetect(spkData(i,:),30,4.5);
                        else
                            disp(['***ALERT*** Channel ' num2str(electrodeMatrix(i,j)) ' has 0 events in this segment...']);
                            events.inds = [];
                        end
                        
                        % ADD TO TOTAL LIST OF EVENTS
                        inds = [inds events.inds];
                        vals = [vals spkData(i,events.inds)];
                        chan = [chan ones(1,numel(events.inds)).*i];
                    
                        % DISPLAY PHYS AND BEHAVIOR
                        if dispOn==1
                            disp('Displaying data for the current channel...');
                            if i==1
                                figure(1); clf;
                            end

                            figure(1);

                            subplot(10,1,1:8);
                            if rem(j,2)==0
    %                             plot(spkData(i,:)+((((j-1).*8)+(i-1)).*1000),'Color',[0 0.67 1]); hold on; axis tight;
                                plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*8)+(i-1)).*1000,'Color',[0 0.67 1]); hold on; axis off; axis tight;
                            else
    %                             plot(spkData(i,:)+((((j-1).*8)+(i-1)).*1000),'Color',[1 0 0]); hold on; axis tight;
                                plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*8)+(i-1)).*1000,'Color',[1 0    0]); hold on; axis off; axis tight;
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
                    %disp('...done.');
                    %disp(' ');
                
                % FIND UNIQUE EVENTS ACROSS THE SHANK (ALL EVENT MUST BE 0.5 MS APART)
                    [allInds,origIs,uIs]    = unique(inds);
                    allVals                 = vals(origIs);
                    allChan                 = chan(origIs);
                    
                    % all indices more than ~0.5 ms apart should be independent
                    valids = find(diff(allInds)>20);
                    valids = [ 1 , valids + 1]; % deal with the N-1 size of diff, first spike is by definition well separated

                    % all indices less than 1 ms apart need to be evaluated
                    invalids = find(diff(allInds)<=20);
                    invalids = invalids + 1; % deal with the N-1 size of diff
                    
                % RECONCILE the invalid samples. 
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
                    
                % EXTRACT WAVEFORMS FROM ALL TRODES FOR EACH UNIQUE EVENT
                    disp(' '); disp('Extract spike waveforms...');
                    [spikes]        = TNC_EventExtractME(spkData,eventInds,[15,45]);
                    
                % CALCULATE FEATURES FOR SPIKE EVENTS                   
%                     if k==1
                        [features] = TNC_EventQuantME(spikes.wfs, [], round(spikes.inds./30), 15, 8);
%                         pcaStruct = features.pca;
%                     else
%                         [features] = TNC_EventQuantME(spikes.wfs, pcaStruct, round(spikes.inds./30), 15, 8);
%                     end                 
                    
                % UPDATE STRUCTURE FOR SPIKE SORTING
                    disp('4) Updating the feature structure with new data...');
                        
                        featStruct.seg(k).shank(j).inds     = spikes.inds + ((k-1).*chunk.*30000);
                        featStruct.seg(k).shank(j).id       = zeros(numel(spikes.inds),1);
                        featStruct.seg(k).shank(j).ts       = round(spikes.inds./30);
                        featStruct.seg(k).shank(j).params   = features.params;
                        
                        if k==1 && j==1
                            featStruct.paramNames = features.paramNames;
                            featStruct.rSp = rSp;
                            featStruct.cSp = cSp;
                            featStruct.electrodeMatrix = electrodeMatrix;
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

            otherwise
                
                for eInd=1:numChan % shanks | sites
                    
                    wfs = [];
                    wts = [];
                    chan = [];
                    inds = [];
                    vals = [];
                    pcaStruct = [];                    
                    spkData = zeros(1,numel(Ns5DATA.Data(1,:)));


                    disp(' ');
                    disp(' ');
                    disp('_________________________________________________________________________');
                    disp(['...BEGIN FILTERING DATA for seg ' num2str(k) ' of ' num2str(numSegs) ' | site ' num2str(eInd) ' of ' num2str(numChan)]);

                    [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(eInd,arrayType);
                        
                    % FILTER DATA ON THIS CHANNEL
                    disp('_________________________________________________________________________');
                    disp(['1) Filtering data on channel ' num2str(eInd) ' | shank:' num2str(col) ', site:' num2str(row)]);

                        rawData = sgolayfilt(Ns5DATA.Data(eInd,:)-avgSeg,9,21);

%                         [lowBandData,hiBandData] = TNC_FilterData(rawData,Ns5DATA.MetaTags.SamplingFreq,0,0,[1 0]);
%                         spkData = hiBandData.values;
            
                        % USE A WAVELET FILTERING METHOD FOR SPEED
                        wname = 'db4'; 
                        maxlevel = 5; % 7.5kHz low pass
                        [c,l] = wavedec(rawData, maxlevel, wname);
                        c = wthcoef('a', c, l);
                        spkData = waverec(c, l, wname);
                   

                    % EXTRACT WAVEFORMS FROM ALL TRODES FOR EACH UNIQUE EVENT
                    disp(['2) Detecting and extracting events...']);

                        [events] = TNC_EventDetect(spkData,30,4);

                        if numel(events.ts)>10
                            [spikes]        = TNC_EventExtract(spkData,rawData,events.inds,[15,45],dispOn);
                            spikes.inds     = events.inds;
                            sessionStruct.seg(k).shank(eInd).spikes = spikes;
                                                
                            % CALC SCALAR QUANT FOR SPIKES
                            disp(['3) Quantify spikes...']);

                            if k==1
                                [features]   = TNC_EventQuantSE(spikes,[]);
                                templates(eInd).template = features.template;
                            else
                                [features]   = TNC_EventQuantSE(spikes,templates(eInd).template);
                            end

                        else
                            
                            events.inds = [];
                            features.params = [];
                            
                        end
                        
                    % UPDATE STRUCTURE FOR SPIKE SORTING
                    disp('4) Updating the feature structure with new data...');
                        
                        featStruct.seg(k).shank(eInd).inds     = events.inds + ((k-1).*chunk.*30000);
                        featStruct.seg(k).shank(eInd).id       = zeros(numel(events.inds),1);
                        featStruct.seg(k).shank(eInd).ts       = round(featStruct.seg(k).shank(eInd).inds./30);
                        featStruct.seg(k).shank(eInd).params   = features.params;

                        
                        if k==1 && eInd==1
                            featStruct.paramNames = features.paramNames;
                            featStruct.rSp = rSp;
                            featStruct.cSp = cSp;
                            featStruct.electrodeMatrix = electrodeMatrix;
                        end

                        sessionStruct.seg(k).shank(eInd).wfs    = spikes.wfs;
                        sessionStruct.seg(k).shank(eInd).inds   = spikes.inds;
                        
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

disp(['save ' targetName '_ft featStruct']); 
save([targetName '_ft.mat'],'featStruct');

%% SAVE SESSION STRUCTURE
 
disp(['save ' targetName '_ss sessionStruct']);
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