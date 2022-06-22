function [sessionStruct] = TNC_SSpipe(fileName,arrayType,chunk,dispOn)

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
    
    
%% PRE-PROCESSING

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
    
    for k=1:1%numSegs

    % LOAD ALL CHANNEL DATA
        timeStr = ['t:' num2str((k-1)*chunk) ':' num2str(k*chunk)];
        disp(' ');
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
                                       
                    if j==1
                        figure(1); clf;
                        figure(2); clf;
                    end
                    
                    wfs = [];
                    wts = [];
                    chan = [];
                    inds = [];
                    vals = [];
                    
                    for i=1:8 %electrodes

                        tCh = electrodeMatrix(i,j);
                        disp(' '); disp(' ');
                        disp(['Filtering data on channel ' num2str(tCh) ' ...']);
                        rawData = sgolayfilt(Ns5DATA.Data(tCh,:)-avgSeg,11,21);
                        [lowBandData,hiBandData] = TNC_FilterData(rawData,Ns5DATA.MetaTags.SamplingFreq,0,0);
                        if i==1
                            spkData = zeros(8,numel(hiBandData.values));
                        end
                        spkData(i,:) = hiBandData.values;                        
                        
                    % DETECT EVENTS ON THIS CHANNEL
                        disp('Looking for events...');
                        [events] = TNC_EventDetect(spkData(i,:),30,5);
                        
                    % ADD TO TOTAL LIST OF EVENTS
                        inds = [inds events.inds];
                        vals = [vals spkData(i,events.inds)];
                        chan = [chan ones(1,numel(events.inds)).*i];
                                                
                    % DISPLAY PHYS AND BEHAVIOR
                        if dispOn==1
                            if j<=4
                                figure(1);
                            else
                                figure(2);
                            end

                            subplot(10,1,1:8);
                            if rem(j,2)==0
    %                             plot(spkData(i,:)+((((j-1).*8)+(i-1)).*1000),'Color',[0 0.67 1]); hold on; axis tight;
                                plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*8)+(i-1)).*1000,'Color',[0 0.67 1]); hold on; axis off; axis tight;
                            else
    %                             plot(spkData(i,:)+((((j-1).*8)+(i-1)).*1000),'Color',[1 0 0]); hold on; axis tight;
                                plot([events.inds ; events.inds],[spkData(i,events.inds).*0 ; spkData(i,events.inds)]+(((j-1).*8)+(i-1)).*1000,'Color',[1 0    0]); hold on; axis off; axis tight;
                            end
                            if i==1 & ( j==1 | j==4 )
                                subplot(10,1,9);
                                plot(ContData.behavior.sLeverV((120.*1000.*(k-1))+1:120.*1000.*k),'k');
                                axis off; axis tight;
                                subplot(10,1,10);
                                plot(ContData.behavior.rawLick((120.*1000.*(k-1))+1:120.*1000.*k),'k');
                                axis off; axis tight;
                            end
                            drawnow;
                        end                        
                    end
                
                % FIND UNIQUE EVENTS ACROSS THE SHANK (ALL EVENT MUST BE 0.5 MS APART)
                    [allInds,origIs,uIs]    = unique(inds);
                    allVals                 = vals(origIs);
                    allChan                 = chan(origIs);
                    
                    % all indices more than 1 ms apart should be independent
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
                    [spikes]        = TNC_EventExtractME(spkData,shank(j).inds,[15,30]);
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

%% CALCULATE FEATURES FROM WT COEFFS


for k=1:1%numSegs

    % [features] = TNC_EventQuant(aSpikes.SH_wfs,aSpikes.Sh_x,'scalar','dot');
    % aSpikes.features = features;
    
    [features] = TNC_EventQuant(aSpikes.WT_wfs,aSpikes.Sh_x,'pca','dot');
    aSpikes.features = features;
    [val,inds] = sort(aSpikes.features.pca.dot.mmStat,'descend');    
    forClust = [aSpikes.features.pca.dot.values(:,inds(1)),aSpikes.features.pca.dot.values(:,inds(2)),aSpikes.features.pca.dot.values(:,inds(3))];
    
    % distMat = linkage(aSpikes.features.pca.dot.values(:,inds(1)),'average','euclidean');
    % clusterIds = cluster(distMat,'cutoff',1.2);
    clusterIds = kmeans(forClust,10);

    figure(10); clf;
    colormap(winter);
    scatter(aSpikes.features.pca.dot.values(:,inds(1)),aSpikes.features.pca.dot.values(:,inds(2)),4,clusterIds,'.');
%     scatter3(aSpikes.features.pca.dot.values(:,inds(1)),aSpikes.features.pca.dot.values(:,inds(2)),aSpikes.features.pca.dot.values(:,inds(3)),14,clusterIds,'filled');
    axis tight; hold on;

        % calculate cluster means, errors and locations in feature space
    for m=1:10
        theseSpks = find(clusterIds==m);
        if numel(theseSpks>0)
            clstCntrs(m,1) = mean(aSpikes.features.pca.dot.values(theseSpks,inds(1)));
            clstCntrs(m,2) = mean(aSpikes.features.pca.dot.values(theseSpks,inds(2)));
            clstCntrs(m,3) = mean(aSpikes.features.pca.dot.values(theseSpks,inds(3)));
            clstErrs(m,1) = std(aSpikes.features.pca.dot.values(theseSpks,inds(1)));
            clstErrs(m,2) = std(aSpikes.features.pca.dot.values(theseSpks,inds(2)));
            clstErrs(m,3) = std(aSpikes.features.pca.dot.values(theseSpks,inds(3)));
    %         text(clstCntrs(m,1),clstCntrs(m,2),clstCntrs(m,3),['CL: ' num2str(m)]);
            plot(clstCntrs(m,1),clstCntrs(m,2),'ko','MarkerSize',6);
            plot([clstCntrs(m,1),clstCntrs(m,1)+clstErrs(m,1)],[clstCntrs(m,2),clstCntrs(m,2)],'k-');
            plot([clstCntrs(m,1),clstCntrs(m,1)],[clstCntrs(m,2),clstCntrs(m,2)+clstErrs(m,2)],'k-');
            plot([clstCntrs(m,1),clstCntrs(m,1)-clstErrs(m,1)],[clstCntrs(m,2),clstCntrs(m,2)],'k-');
            plot([clstCntrs(m,1),clstCntrs(m,1)],[clstCntrs(m,2),clstCntrs(m,2)-clstErrs(m,2)],'k-');
            text(clstCntrs(m,1),clstCntrs(m,2),['    ' num2str(m)]);
        end
    end
    
end
    
    
%% MERGE CLUSTERS AND THEN SUBTRACT CLUSTERS FROM SPKDATA

    % choose by number which clusters to merge
    [mClstIds] = TNC_MergeClusters(clusterIds, [2,4]);
    [mClstIds1] = TNC_MergeClusters(mClstIds, [3,6]);
    [mClstIds2] = TNC_MergeClusters(mClstIds1, [1,4,5,7,8,9,10]);

    figure(10); clf;
    colormap(jet)
    scatter(aSpikes.features.pca.dot.values(:,inds(1)),aSpikes.features.pca.dot.values(:,inds(2)),4,mClstIds2);
%     scatter3(aSpikes.features.pca.dot.values(:,inds(1)),aSpikes.features.pca.dot.values(:,inds(2)),aSpikes.features.pca.dot.values(:,inds(3)),14,clusterIds,'filled');
    axis tight; hold on;

    for m=1:10
        theseSpks = find(mClstIds2==m);
        if numel(theseSpks>0)
            clstCntrs(m,1) = mean(aSpikes.features.pca.dot.values(theseSpks,inds(1)));
            clstCntrs(m,2) = mean(aSpikes.features.pca.dot.values(theseSpks,inds(2)));
            clstCntrs(m,3) = mean(aSpikes.features.pca.dot.values(theseSpks,inds(3)));
            clstErrs(m,1) = std(aSpikes.features.pca.dot.values(theseSpks,inds(1)));
            clstErrs(m,2) = std(aSpikes.features.pca.dot.values(theseSpks,inds(2)));
            clstErrs(m,3) = std(aSpikes.features.pca.dot.values(theseSpks,inds(3)));
    %         text(clstCntrs(m,1),clstCntrs(m,2),clstCntrs(m,3),['CL: ' num2str(m)]);
            plot(clstCntrs(m,1),clstCntrs(m,2),'ko','MarkerSize',8,'LineWidth',2);
            [clX clY] = calculateEllipse(clstCntrs(m,1), clstCntrs(m,2), clstErrs(m,1).*2, clstErrs(m,2).*2, 0, 50);
            plot(clX,clY,'k-','LineWidth',2);            
            text(clstCntrs(m,1),clstCntrs(m,2),['    ' num2str(m)]);
        end
    end
    
    % calculate raw voltage templates
    [clstTemp] = TNC_GetTemplate(aSpikes.Sh_wfs,mClstIds2);
    figure(11); clf;
    plot(clstTemp.x,clstTemp.wf);
    
    
    % subtract from the filtered data and look for residuals
    

%% 
    if arrayType

        for j=1:numShanks
        end

    else

        for j=1:numSites
        end

    end

%% CONCATENATE
    if arrayType 
        [row,col,rSp,cSp] = TNC_RemapElecPos(electrodeNumber,arrayType);
    end
    
%% CUT CLUSTERS

%% COMPUTE TEMPLATES

    for k =1:numClusters               
    end
    
%% REMOVE TEMPLATES FROM RAW

%% DETECT WITH LOWER THRESHOLD AND CHECK RESIDUAL DATA (LOOPED)

%% WRITE OUT DATA STRUCTURE

