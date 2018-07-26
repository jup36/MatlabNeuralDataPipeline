function [sessionStruct] = TNC_SSP_SingleSite(fileName,arrayType,chunk,dispOn)

%% PRE-PROCESSING
    % GET INFO ABOUT THE CURRENT RECORDING
    Ns5DATA = openNSx('report',fileName);
    numSegs = floor(Ns5DATA.MetaTags.Duration ./ chunk);
    disp(' ');
    disp(' ');
    disp(['Data will be loaded and processed as ' num2str(numSegs) ' x ' num2str(chunk) ' seconds long segments.']);
    disp(' ');
    disp(' ');

    sessionStruct.chunk     = chunk;
    sessionStruct.arrayType = arrayType;
    sessionStruct.fileName  = fileName;
    sessionStruct.numSegs   = numSegs;

    parfor q=1:10 % for some reason this initializes the parallel computing toolbox correctly
        q;
    end
    
%% EXTRACT MULTIDIMENSIONAL EVENTS

    wfs = [];
    wts = [];
    chan = [];
    inds = [];
    vals = [];

    % LOAD BEHAVIORAL DATA TO DISPLAY WITH THE SPIKING DATA
    totChar = numel(fileName);
    eval(['load ' fileName(1:totChar-4) '.mat']);
    
    % GET THE ELECTRODE MAPPINGS
    [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);    
    
    for k=1:numSegs

    % LOAD ALL CHANNEL DATA
        timeStr = ['t:' num2str((k-1)*chunk) ':' num2str(k*chunk)];
        disp(' ');disp(' ');disp(' ');disp(' ');disp(' ');disp(' ');
        disp(['Loading data over the range of ' timeStr]);
        disp(' ');
        Ns5DATA = openNSx('report','read',fileName,timeStr,'sec');
        
    % SMOOTH/FILTER DATA
        avgSeg      = mean(Ns5DATA.Data,1);
        numChan     = size(Ns5DATA.Data,1);
        
        switch arrayType

            %________________________________________________________
            %________________________________________________________
            %________________________________________________________
            %_______CODE FOR BUZSAKI 64 ARRAY________________________

            case 'NN_b64'

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
                    disp(['...BEGIN FILTERING DATA for seg ' num2str(k) ' | shank ' num2str(j)]);
                    disp(' ');

                    parfor i=1:8 %electrodes
                        %disp(['Filtering data on channel ' num2str(electrodeMatrix(i,j)) ' ...']);
                        rawData = sgolayfilt(Ns5DATA.Data(electrodeMatrix(i,j),:)-avgSeg,11,21);
                        [lowBandData,hiBandData] = TNC_FilterData(rawData,Ns5DATA.MetaTags.SamplingFreq,0,0,[1 0]);
                        spkData(i,:) = hiBandData.values;                                                
                    end
                    disp('...completed');
                    
                    disp(' ');
                    disp(['Detecting events...']);
                    
                    for i=1:8
                    
                        % DETECT EVENTS ON THIS CHANNEL
                        %disp('Looking for events...');
                        [events] = TNC_EventDetect(spkData(i,:),30,5);
                        
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

                            subplot(10,1,9);
                            plot(ContData.behavior.sLeverV((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
                            axis off; axis tight;
                            subplot(10,1,10);
                            plot(ContData.behavior.rawLick((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
                            axis off; axis tight;
                                
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
                    valids = find(diff(allInds)>14);
                    valids = [ 1 , valids + 1]; % deal with the N-1 size of diff, first spike is by definition well separated

                    % all indices less than 1 ms apart need to be evaluated
                    invalids = find(diff(allInds)<=14);
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
                
                    shank(j).inds = unique(allInds([valids addValids])); % use unique to sort the indices
                    
                % EXTRACT WAVEFORMS FROM ALL TRODES FOR EACH UNIQUE EVENT
                    disp('Extract spike waveforms...');
                    [spikes]        = TNC_EventExtractME(spkData,shank(j).inds,[10,30]);
                    spikes.inds     = shank(j).inds;
                    shank(j).spikes = spikes;
                    
                % CALCULATE FEATURES FOR SPIKE EVENTS                   
                    if k==1
                        wfsStruct = shank(j).spikes.wfs;
                        [features] = TNC_EventQuantME(wfsStruct, [], 'dot');
                        pcaStruct = features.pca;
                    else
                        wfsStruct = shank(j).spikes.wfs;
                        [features] = TNC_EventQuantME(wfsStruct, pcaStruct, 'dot');
                    end
                    
                    shank(j).features = features;                    
                    
                end

            %________________________________________________________
            %________________________________________________________
            %________________________________________________________
            %_______CODE FOR WIDE 64 ARRAY___________________________

            case 'NN_w64'

        
        end
        
        sessionStruct.seg(k).shank = shank;
        
    end

%% BUILD FEATURE STRUCTURE

clear paramVals featStruct wfStruct newValues
offset                  = 0;
featStruct.chunk        = chunk;
featStruct.arrayType    = arrayType;

featStruct.seg(sessionStruct.numSegs).shank(size(sessionStruct.seg(1).shank,2)).id   = 0;
wfStruct(size(sessionStruct.seg(1).shank,2)).chan(8).v = [];

for i=1:sessionStruct.numSegs
    for j=1:size(sessionStruct.seg(1).shank,2)
                
        if i==1
            paramVals(j).values = [];
        end
        
        featStruct.seg(i).shank(j).pca  = sessionStruct.seg(i).shank(j).features.scl.dot.values;
        featStruct.seg(i).shank(j).amp  = sessionStruct.seg(i).shank(j).features.scl.values;
        featStruct.seg(i).shank(j).ts   = sessionStruct.seg(i).shank(j).spikes.times  + (chunk.*(i-1).*30000);
        featStruct.seg(i).shank(j).id   = zeros(numel(sessionStruct.seg(i).shank(j).spikes.times),1);
        featStruct.seg(i).shank(j).cnt  = zeros(1,size(sessionStruct.seg(i).shank(j).features.scl.dot.values,2));
        
        newValues                       = [featStruct.seg(i).shank(j).ts + (chunk.*(i-1).*30000) , featStruct.seg(i).shank(j).pca , featStruct.seg(i).shank(j).amp.minV' , featStruct.seg(i).shank(j).amp.maxV' , featStruct.seg(i).shank(j).amp.energy' , featStruct.seg(i).shank(j).amp.peVa' ];
        paramVals(j).values             = [paramVals(j).values ; newValues];
        
        
    end
end

%% SAVE OUTPUT

totChar = numel(fileName);

disp(['save ' fileName(1:totChar-4) '_ft featStruct']); 
eval(['save ' fileName(1:totChar-4) '_ft featStruct']); 

for j=1:size(sessionStruct.seg(1).shank,2)
            
    % For use with matclust
    filedata.params     = paramVals(j).values;
    filedata.chunk      = chunk;
    filedata.shank      = j;
    filedata.filename   = fileName;
    filedata.paramnames = {'ts','pc1','pc2','pc3','pc4','pc5','pc6','pc7','pc8','minV1','minV2','minV3','minV4','minV5','minV6','minV7','minV8','maxV1','maxV2','maxV3','maxV4','maxV5','maxV6','maxV7','maxV8','energy1','energy2','energy3','energy4','energy5','energy6','energy7','energy8','peVa1','peVa2','peVa3','peVa4','peVa5','peVa6','peVa7','peVa8'}';

    disp(['save ' fileName(1:totChar-4) '_shank' num2str(j) '_mc filedata']);
    eval(['save ' fileName(1:totChar-4) '_shank' num2str(j) '_mc filedata']);

    clear fileData

end

%% BUILD THE WF STRUCTURE 

for j=1:size(sessionStruct.seg(1).shank,2)
        
    for i=1:sessionStruct.numSegs
        
        wfdata.seg(i).wfs(numel(sessionStruct.seg(i).shank(j).spikes.wfs)).values = zeros(8,41);
        
        for k=1:numel(sessionStruct.seg(i).shank(j).spikes.wfs)
            wfdata.seg(i).wfs(k).values = int16(sessionStruct.seg(i).shank(j).spikes.wfs(k).values);
        end
    end
    
    % For getting waveform data 
    wfdata.chunk    = chunk;
    wfdata.shank    = j;
    wfdata.filename = fileName;
    wfdata.datatype = 'int16';
    
    disp(['save ' fileName(1:totChar-4) '_shank' num2str(j) '_wf wfdata']);
    eval(['save ' fileName(1:totChar-4) '_shank' num2str(j) '_wf wfdata']);

    clear wfdata    
    
end

