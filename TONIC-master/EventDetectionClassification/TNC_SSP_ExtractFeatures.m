%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function [sessionStruct,featStruct] = TNC_SSP_ExtractFeatures(dataFileName, arrayType, numSegs, varargin)

% INPUT ARGUMENTS
%
% REQUIRED:
% 1 = dataFileName = input (binary) file name
% 2 = arrayType
% 3 = number of time segments ('numSeg') 
%
% OPTIONAL (POSITION SPECIFIC):
% 4 = iNode (index of the cluster node to be used for data processing)
%
% OPTIONAL PARAMETERS (to be specified using name-value pairs)
% display: whether or not display data, e.g. 'display, 1 (default = 0)
% verbose: whether or not produce more output, e.g. 'verbose', 1 (default = 0)
%
%% PARSE INPUT ARGUMENTS LIST
    try
        p = TNC_OOP_FeatureExtractionInputParser; % call constructor
        p.parse(dataFileName, arrayType, numSegs, varargin{:});

        % Required inputs:
        dataFileName      = p.Results.dataFileName;
        options.arrayType = p.Results.arrayType;
        options.numSegs   = p.Results.numSegs;

        % Optional, position specific
        options.inode     = uint32(p.Results.iNode);

        % Optional parameters (to be specified as parameter-value pairs)
        options.targetName   = p.Results.targetName;
        if isnumeric(numSegs)  % original Matlab version
            options.iSeg     = p.Results.iSeg;
            options.iShank   = p.Results.iShank;
            options.dispOn   = p.Results.dispOn; 
            options.verbose  = p.Results.verbose;
            options.fakeData = p.Results.fakeData;
            options.rawData  = p.Results.rawData;
            options.snr      = p.Results.snr;
            options.scm      = p.Results.scm;

        else                  % compiled version
            options.iSeg     = int32(str2double(p.Results.iSeg));
            options.iShank   = int32(str2double(p.Results.iShank)); 
            options.dispOn   = int32(str2double(p.Results.dispOn));
            options.verbose  = int32(str2double(p.Results.verbose));
            options.fakeData = int32(str2double(p.Results.fakeData));
            options.rawData  = int32(str2double(p.Results.rawData));
            options.snr      =       str2double(p.Results.snr);
            options.scm      =       str2double(p.Results.scm);
        end
    catch
        output_usage_message();
        return
    end
 
    % CHECK IF INPUT DATA EXIST
    if ~(exist(dataFileName, 'file') == 2)
        disp(' ');
        disp(['Input file ' dataFileName ' is not found']);
        return;
    end

    disp(['options.verbose =' num2str(options.verbose ) ' options.scm=' num2str(options.scm)]);
    if options.verbose
        disp(' ');
        disp(['dataFileName=' dataFileName ' arrayType=' options.arrayType ...
              ' numSegs=' num2str(options.numSegs) ' inode=' num2str(options.inode) ' iSeg=' num2str(options.iSeg) ...
              ' iShank=' num2str(options.iShank) ' targetName=' options.targetName ' dispOn=' num2str(options.dispOn) ...
              ' verbose=' num2str(options.verbose) ' fakeData=' num2str(options.fakeData) ' snr=' num2str(options.snr) ...
              ' scm=' num2str(options.scm) ]);
    end

    % DEFINE CONSTANTS
    samplingRate  = 30; % samples per ms
    winL          = 16;
    winR          = 47; % winL+winR+1 must be = 64, in order to make wavelet analysis applicable
    maxSpikeDelay = 0.25;  % in ms

%% PRE-PROCESSING
    sessionStruct = [];
    featStruct    = [];

    % GET INFO ABOUT THE CURRENT RECORDING
    if ~isempty(strfind(dataFileName, '.ns'))
        try
            sigDATA = openNSx('report',dataFileName);
            chunk   = double(sigDATA.MetaTags.Duration) ./ options.numSegs;
            disp(['fileName=' dataFileName ' size(sigDATA)=' num2str(size(sigDATA)) ' Duration=' num2str(sigDATA.MetaTags.Duration)]);
            channelCount = uint32(sigDATA.MetaTags.ChannelCount);
        catch
            disp('Unable to read input data. Case 1.');
            return;
        end
    elseif ~isempty(strfind(dataFileName,'.mat'))
        samplingRate = 24;  % samples per ms
        sigDATA = extract_signal_data_from_mat_file(dataFileName,1,options);
        chunk   = double(sigDATA.MetaTags.Duration) ./ options.numSegs;
        channelCount = uint32(sigDATA.MetaTags.ChannelCount);    
    else
        disp([dataFileName ': unsupported input data file format']);
        return;
    end

    disp(' ');
    disp(['Data will be loaded and processed as ' num2str(uint32(options.numSegs))...
          ' x ' num2str(uint32(chunk)) ' seconds long segments.']);
    disp(' ');
    sampleDuration = 1/samplingRate;  % time between two adjacent scans (in ms)
    disp(['sampleDuration=' num2str(sampleDuration) ' chunk=' num2str(chunk)...
          ' samplingRate=' num2str(samplingRate)]);
    if options.verbose
        disp(' ');
        disp(['file_name=' dataFileName ' size(sigDATA)=' num2str(size(sigDATA)) ' chunk=' num2str(chunk)]);
    end
    clear sigDATA;

%% PROCESS OPTIONAL INPUT PARAMETERS 
    if findstr(options.arrayType,'NN_b64')
        if channelCount ~= 64
            disp('Incorrectly specified arrayType');
            return   
        end
        shankSize = 8;
        numShanks = channelCount / shankSize;
        numNodes = uint32(options.numSegs) .* numShanks; 
    elseif findstr(options.arrayType,'NN_b32')
        if channelCount ~= 32
            disp('Incorrectly specified arrayType');
            return
        end
        shankSize = 8;
        numShanks = channelCount / shankSize;
        numNodes = uint32(options.numSegs) .* numShanks;
    elseif strcmp(options.arrayType,'tetrode')
        shankSize    = 4;
        channelCount = 4;
        numShanks    = channelCount / shankSize;
        numNodes    = 4;
    else % shankSize == 1
        numNodes = uint32(options.numSegs) .* uint32(channelCount);  
        disp(['numSegs=' num2str(uint32(options.numSegs)) ...
              ' channelCount=' num2str(uint32(channelCount)) ...
              ' numNodes=' num2str(numNodes)]);
    end

    my_j = 0;
    my_k = 0;     
    if options.inode > 0           % running on cluster
        % my_j = the shank id which is in fact specified by user 
        %        (by specifying inode = index of a cluster node)
        my_j = mod(options.inode, numShanks);
        if my_j == 0
            my_j = numShanks;
        end
        if options.inode > numShanks
            my_k = (options.inode - my_j)/numShanks + 1;
        else
            my_k = 1;
        end
   
        if options.verbose
            disp(' ');
            disp(['numShanks=' num2str(numShanks) ' my_j=' num2str(my_j) ...
                  ' my_k=' num2str(my_k) ' inode=' num2str(options.inode) ...
                  ' numNodes=' num2str(numNodes)]);
            disp(' ');
        end
    elseif options.iShank > 0 || options.iSeg > 0
        if options.iShank > 0
            my_j = options.iShank;
        end
        if options.iSeg > 0
            my_k = options.iSeg;
        end
        if options.verbose
            disp(' ');
            disp([' my_j=' num2str(my_j) ' my_k=' num2str(my_k) ' inode=' num2str(options.inode) ]);
            disp(' ');
        end
    end

%% EXTRACT MULTIDIMENSIONAL EVENTS

% LOAD BEHAVIORAL DATA TO DISPLAY WITH THE SPIKING DATA
    totChar = numel(dataFileName);
    if options.dispOn==0
        loadBehavior=0;
    else
        d = dir([dataFileName(1:totChar-4) '.mat']);
        if numel(d)>0
            loadBehavior = 1;
        else
            disp('No associated behavioral data was found.');
            loadBehavior = 0;
        end
    end
   
% LOAD BEHAVIORAL DATA TO DISPLAY WITH THE SPIKING DATA
    totChar = numel(dataFileName);
    ppBehavFileName = [dataFileName(1:totChar-4) '.mat'];

    if loadBehavior
        eval(['load ' dataFileName(1:totChar-4) '.mat']);
        behaviorTrue = 1;
    else
        behaviorTrue = 0;
    end
 
    min_k = 1;
    max_k = options.numSegs;
    if options.inode > 0 || options.iSeg > 0
        min_k = my_k;
        max_k = my_k;
    end

    disp(['verbose=' num2str(options.verbose) ' min_k=', num2str(min_k), ' max_k=', num2str(max_k)]);
    if options.verbose
        disp (['min_k=', num2str(min_k), ' max_k=', num2str(max_k), ' numSegs=', num2str(options.numSegs)]);
    end

    for k=min_k:max_k % loop through time segments 

        % LOAD ALL CHANNEL DATA
        timeStr = strcat('t:',num2str(double(k-1)*double(chunk)),':',num2str(double(k)*double(chunk)));
        disp(['k=' num2str(k) ' chunk=' num2str(chunk) ' timeStr=' timeStr ]);
        if ~isempty(strfind(dataFileName, '.ns'))
            disp(' ');disp(' ');
            disp(['k=' num2str(k) ' chunk=' num2str(chunk) ' Loading data over the range of ' timeStr]);
            sigDATA = struct('MetaTags',[],'Data',[]);
            disp(['k=', num2str(k) ' args to openNSx: dataFileName=' dataFileName ' timeStr=' timeStr])
            try
                sigDATA = openNSx('report','read',dataFileName,timeStr,'sec');
            catch
                disp('Unable to read input data. Case 2. Please, check the file name or specify a larger # of segments');
                return
            end
        elseif ~isempty(strfind(dataFileName, '.mat'))
            disp(['file_name=' dataFileName ' k=' num2str(k) ' numSegs=' num2str(options.numSegs)]);
            try
                sigDATA = extract_signal_data_from_mat_file(dataFileName, k, options);
            catch
                disp('Unable to read input data. Case 3.');
                return   
            end
        end

        % SMOOTH/FILTER DATA
        disp(['dataFileName=' dataFileName ' timeStr=' timeStr ])
%       numChan = size(sigDATA.Data,1);
        numChan     = numel(find(sigDATA.MetaTags.ChannelID<129))
        disp(['size(sigDATA.Data)=' num2str(size(sigDATA.Data)) ' numChan=' num2str(numChan)]);

        % SOME CHECKING OF DATA FOR PROBES WHERE MULTIPLE SITES ARE USED FOR SORTING
        disp(['arrayType=' options.arrayType ' k(segment #)=' num2str(k)]);

        if k== 1 || k==my_k 
            pcaStruct = [];
            featStruct = [];
            sessionStruct = [];
        end

        %% PROCESS DATA FROM DIFFERENT TYPES OF ARRAYS    
        snrThresh  =  5.5; % the lower this value, the more events will be detected
        if options.snr > 0
            snrThresh = options.snr;
        end
        disp(['snrThresh=' num2str(snrThresh)]);
        switch options.arrayType
                          
            % Handle Buzsaki-64 array data
            case {'NN_b64' 'NN_b64_dhs' 'NN_b32' 'NN_b32_dhs'}
                featStruct.snr = snrThresh;
                [pcaStruct,featStruct,sessionStruct] = ...
                process_multichannel_data(dataFileName,sigDATA,k,...
                                          pcaStruct,featStruct,sessionStruct,...
                                          numChan,shankSize,numShanks,my_j,my_k,chunk,...
                                          samplingRate,snrThresh,winL,winR,...
                                          maxSpikeDelay,sampleDuration,options);
            case 'tetrode'
                snrThresh  =  36.0;
                if options.snr > 0
                    snrThresh = options.snr;
                end
                featStruct.snr = snrThresh;
                [pcaStruct,featStruct,sessionStruct] = ...
                process_multichannel_data(dataFileName,sigDATA,k,...
                                          pcaStruct,featStruct,sessionStruct,...
                                          numChan,shankSize,numShanks,my_j,my_k,chunk,...
                                          samplingRate,snrThresh,winL,winR,...
                                          maxSpikeDelay,sampleDuration,options);       

            % Handle probes of all other types
            otherwise
                featStruct.snr = snrThresh;
                [featStruct,sessionStruct] = ...
                process_isolated_site_data(dataFileName, sigDATA, k, min_k,...
                                           featStruct,sessionStruct,...
                                           numChan,my_j,my_k,...
                                           chunk,samplingRate,snrThresh,...
                                           winL,winR,options);        
        end % switch arrayType
        
    end % loop in k

    %% SAVE  FEATURE STRUCTURE

%   clear paramVals featStruct wfStruct newValues
    offset                  = 0;
    featStruct.chunk        = chunk;
    featStruct.arrayType    = options.arrayType;

    output_ft_data(my_j, my_k, options.targetName, featStruct);

    %% SAVE SESSION STRUCTURE
    output_ss_data(my_j, my_k, options.targetName, sessionStruct);
%   disp(['save ' targetName '_ss sessionStruct']);
%   save([targetName '_ss.mat'],'sessionStruct');

    %% SAVE OUTPUT
    % OUTPUT MC DATA
%   output_mc_data(my_j, min_k, my_k,  dataFileName, chunk, sessionStruct, paramVals);    

    % OUTPUT WF DATA
%   output_wf_data(my_j, min_k, max_k, dataFileName, chunk, sessionStruct, totChar);

    return

%-------------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: TNC_SSP_ExtractFeatures(dataFileName,arrayType,numSegs [,node] [,parameters])');
    disp('Required arguments:');
    disp('    dataFileName - name of input file');
    disp('    arrayType    - string specifying an array type used to generate data');
    disp('                   (may be NN_b64, tetrode or anything else)');
    disp('    numSegs      - number of the time segments into which to split data');
    disp('Optional argument:')
    disp('    node         - cluster node id (used only in cluster processing');
    disp('Optional parameters, specified as parameter-value pairs:');
    disp('    iSeg         - id of the time segment to be processed; <= numSegs');
    disp('                   default=-1 means processing all segments');
    disp('    iShank       - id of the shank to be processed;');
    disp('                   default=-1 means processing all shanks');
    disp('    targetName   - name prefix of an output file;');
    disp('                   default: same as prefix of input file');
    disp('    snr          - a threshold value of signal-to-noise ratio');
    disp('                   defaults: 5.0 for real data and 2.8 for fake data');
    disp('    scm          - subtract common mean, i.e. remove non-physiological ');
    disp('                   noise common for all channels; default= 0');
    disp('                   (Ludwig ea - J Neurophysiol 2009, v.101, p. 1679);');
    disp('    verbose  - weather or not to increase verbosity of output (default=0)');
    disp('    dispOn   - weather or not to display results (default=0)');
    disp('    fakeData - weather or not the input data is fake (default=0)');
    disp('    rawData  - weather or not to perform filtering of data (default=1)');
    return

%-------------------------------------------------------------------------------

function [pcaStruct,featStruct,sessionStruct] = ...
          process_multichannel_data(dataFileName, sigDATA, k, pcaStruct, ... 
                                    featStruct,sessionStruct, numChan, ...
                                    shankSize,numShanks,my_j,my_k,chunk,...
                                    samplingRate,snrThresh,winL,winR,...
                                    maxSpikeDelay,sampleDuration,options)
    sessionStruct.chunk        = chunk;
    sessionStruct.arrayType    = options.arrayType;
    sessionStruct.dataFileName = dataFileName;
    sessionStruct.numSegs      = options.numSegs;

    allRecChan = [1:numChan];
    disp(['numChan=' num2str(numChan)]);

    % scm: subtract common mean
    % remove non-physiological noise common for all channels
    % (Ludwig ea - J Neurophysiol 2009, v.101, p. 1679)
    if options.scm    % subtract common mean              
        avgSig  = mean(sigDATA.Data,1);  % mean across all channels
    else
        avgSig  = mean(sigDATA.Data,2);  % mean across one channel
    end

    % GET THE ELECTRODE MAPPINGS
    [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,options.arrayType);

    min_i = 1;
    max_i = shankSize;
    min_j = 1;
    max_j = numShanks;

    if options.inode > 0 || options.iShank > 0
        min_j = my_j;  % user-specified shank id
        max_j = my_j;
    end

    disp(['min_j=' num2str(min_j) ' max_j(shank #)=' num2str(max_j)]);
    for j=min_j:max_j % loop through shanks
        wfs = [];
        wts = [];
        chan = [];
        inds = []; % index of signal time series corresponding to spike max
        vals = []; % value of signal at the max
       
        spkData = zeros(shankSize,numel(sigDATA.Data(1,:)));
        disp('_________________________________________________________________________');
        disp('_________________________________________________________________________');
        disp('_________________________________________________________________________');
        disp('_________________________________________________________________________');
        disp(['...BEGIN FILTERING DATA for seg ' num2str(k) ' of ' num2str(options.numSegs) ' | shank ' num2str(j)]);
        disp(' ');

        disp(['Filtering data for shank ' num2str(j) ' ...']);
        disp(['size(sigDATA.Data)=' num2str(size(sigDATA.Data))]);
        if options.verbose
            disp(['min_i=' num2str(min_i) ' max_i=' num2str(max_i)]);
            disp('electrodeMatrix=');
            electrodeMatrix
        end
        for i=min_i:max_i % loop through channels in a given shank
            if numel( find(electrodeMatrix(i,j) == allRecChan) ) > 0
                if strcmp(options.arrayType, 'tetrode') && ~options.rawData % already filtered data
                    spkData(i,:) = sigDATA.Data(electrodeMatrix(i,j),:);
%                   disp(['i=' num2str(i) ' spkData(i,:)=' num2str(spkData(i,:))]); 
                else % real data or unfiltered fake data
                    disp(['   Filtering data on channel ' num2str(electrodeMatrix(i,j)) ' ...']);
                    if options.scm % subtract common mean 
                        rawData = sgolayfilt(sigDATA.Data(electrodeMatrix(i,j),:)-avgSig   ,11,21);
                    else
                        rawData = sgolayfilt(sigDATA.Data(electrodeMatrix(i,j),:)-avgSig(i),11,21);
                    end

                    if options.verbose
%                       disp(['   rawData=' num2str(rawData)]);
                        disp(['   sigDATA.MetaTags.SamplingFreq=' num2str(sigDATA.MetaTags.SamplingFreq)]);
                        disp(' ');
                    end
                    [lowBandData, hiBandData] = TNC_FilterData(rawData,sigDATA.MetaTags.SamplingFreq,0,0,[1 0]);
                    spkData(i,:) = hiBandData.values;  
%                   disp(['i=' num2str(i) ' hiBandData.values=' num2str(hiBandData.values)]); 
                end                                             
            else
                disp(['    ***ALERT*** Channel ' num2str(electrodeMatrix(i,j)) ...
                      ' was not recorded in this segment...']);
                spkData(i,:) = zeros(size(sigDATA.Data(1,:)));
            end
        end
        disp('...done');

        % REMOVE ARTIFACTS
%       len = numel(spkData(1,:));
%       K   = numel(spkData(:,1));
%       for i1=1:len
%           for j1=1:K
%               if i1 > 1 && i1 < len && abs(spkData(j1,i1)) < 0.75*abs(spkData(j1,i1-1)) &&  abs(spkData(j1,i1)) < 0.75*abs(spkData(j1,i1+1))    
%                   spkData(j1,i1) = (spkData(j1,i1-1) + spkData(j1,i1+1))/2;
%               end
%           end
%       end
        
        disp(' ');
        disp('_________________________________________________________________________');
        disp(' ');
        disp('Detecting events...');
       
        % DETECT EVENTS IN EACH CHANNEL AND CREATE JOINT ARRAY OF INDICES
        for i=min_i:max_i % loop through channels in a given shank
        
            % DETECT EVENTS ON THIS CHANNEL
%           disp(['numel(find(allRecChan==electrodeMatrix(i,j)))=' num2str(numel(find(allRecChan==electrodeMatrix(i,j))))]);
%           disp(['size(spkData)=' num2str(size(spkData))]);
            if numel(find(allRecChan==electrodeMatrix(i,j)))==1
                % Detect indices of the signal minima corresponding to spikes
                [spkEvents] = TNC_EventDetect(spkData(i,:),samplingRate,snrThresh);
            else
                disp(['***ALERT*** Channel ' num2str(electrodeMatrix(i,j)) ' has 0 events in this segment...']);
                spkEvents.inds = [];
            end
            
            % ADD TO TOTAL LIST OF EVENTS
%           disp(['Channel=' num2str(i) ' num_event_inds=' num2str(length(spkEvents.inds))]);
            inds = [inds spkEvents.inds];
            vals = [vals spkData(i,spkEvents.inds)];
            chan = [chan ones(1,numel(spkEvents.inds)).*i];
        
            % DISPLAY PHYS AND BEHAVIOR
            if options.dispOn==1
                if i==min_i
                    figure(1); clf;
                end
                fg = display_data(i, j, k, spkData, spkEvents, chunk, figure(1));       
            end
        end % loop through channels for a given shank

        if options.dispOn==1
            waitfor(figure(1));
            disp('ok');
        end

        disp('...done.');
        disp(' ');

    % GENNADY'S CODE
        Spikes = TNC_OOP_Spikes(num2cell(inds), num2cell(vals), num2cell(chan), shankSize);
        SortedSpikes = sort(Spikes);
        [ eventInds ] = find_events_across_shank(SortedSpikes, spkData, maxSpikeDelay, sampleDuration);

        disp(['After find_events_across_shank: size(eventInds)=' num2str(size(eventInds))]);
        if size(eventInds, 2) == 0
            disp('No spikes confirmed across different channels were found.');
            disp('Please, try to reduce snr or increase maxSpikeDelay');
            disp(' ');
            return; 
        end

    % EXTRACT WAVEFORMS FROM ALL TRODES FOR EACH UNIQUE EVENT
        disp(' '); 
        disp('_________________________________________________________________________');
        disp(' ');
        disp('Extracting spike waveforms...');
%       [spikes1]       = TNC_EventExtractME(spkData, eventInds1,[winL,winR]);
%       disp(['size(spikes1)=' num2str(size(spikes1))]);  
        [spikes]       = TNC_EventExtractME(spkData, eventInds,[winL,winR]);
        disp(['size(spikes)=' num2str(size(spikes))]);

    % CALCULATE FEATURES FOR SPIKE EVENTS                
        if options.verbose
            disp(['k=' num2str(k) ' my_k=' num2str(my_k)]);   
        end
        disp(['size(spikes.wfs)=' num2str(size(spikes.wfs)) ' size(spikes.wfs(1).values)=' num2str(size(spikes.wfs(1).values))]);
        if k==1 || k==my_k
            [features] = TNC_EventQuantME(spkData, spikes.wfs, [], spikes.indsM,...
                                          samplingRate, winL, shankSize);
            pcaStruct  = features.pca;
        else
            [features] = TNC_EventQuantME(spkData, spikes.wfs, pcaStruct, spikes.indsM,...
                                          samplingRate, winL, shankSize);
        end
        disp('...done.');
        disp(' ');
        
    % DETERMINE THE 'TRUE' CLASS LABELS FOR SPIKE EVENTS IN FAKE DATA
        if options.fakeData
            features.Yt = determine_true_class_labels(sigDATA.Spiketimes, ...
                                                      options);      
            disp(['Y_true=' features.Yt]);
        else
            features.Yt = [''];
        end
 
    % UPDATE STRUCTURE FOR SPIKE SORTING
        disp('_________________________________________________________________________');
        disp(' ');
        disp('Updating the feature structure with new data...');
            disp( ['j=', num2str(j) ' k=' num2str(k)]); 
            featStruct.seg(k).shank(j).inds   = spikes.inds + (double(k-1).*chunk.*double(30000));
            featStruct.seg(k).shank(j).id     = zeros(size(spikes.inds,1),1); % needed for visualization
            featStruct.seg(k).shank(j).ts     = round(spikes.inds./samplingRate);
            featStruct.seg(k).shank(j).params = features.params;
            featStruct.seg(k).shank(j).Y_true = features.Yt; 
            disp(['size(id)=' num2str(size(featStruct.seg(k).shank(j).id))]);

            % Splitting params randomly into s given number of folds in the range 2 to 50
            numData = size(features.params, 1);
            oldInds = [1:numData];
            featStruct.seg(k).shank(j).random_split(1).folds = oldInds;
            for nf=2:50
               res     = mod(numData, nf);
               if res > 0
                   % Compliment data with (nf - res) NANs to make its size multiple of nf
                   numData1= numData + (nf - res); 
               else
                   numData1= numData;
               end
               oldInds = [1:numData1];
               newInds = oldInds(randperm(numData1)); % a random permutation of data
               featStruct.seg(k).shank(j).random_split(nf).folds = ...
                         reshape(newInds, numel(newInds)/nf, nf); % each col is a 'set'
            end

            if k==1 && j==1
                featStruct.paramNames = features.paramNames;
                featStruct.rSp = rSp;
                featStruct.cSp = cSp;
                featStruct.electrodeMatrix = electrodeMatrix;
            end
            if options.fakeData > 0
                featStruct.x_neurons = sigDATA.x_neurons;
                featStruct.y_neurons = sigDATA.y_neurons;
                featStruct.z_neurons = sigDATA.z_neurons;
                featStruct.models    = sigDATA.models;
            end

            sessionStruct.seg(k).shank(j).wfs = spikes.wfs;
            sessionStruct.seg(k).shank(j).inds = spikes.inds;

        disp('_________________________________________________________________________');
        disp('...done.');
        disp(' ');
        disp(' ');

    end % loop in j 

%-------------------------------------------------------------------------------

function [featStruct,sessionStruct] = ...
         process_isolated_site_data(FileName, sigDATA, k, min_k,  ...
                                    featStruct,sessionStruct,...
                                    numChan,my_j,my_k,chunk,samplingRate,...
                                    snrThresh,winL,winR,options)             
    persistent templates;

    sessionStruct.chunk     = chunk;
    sessionStruct.arrayType = options.arrayType;
    [pathstr,name,ext]      = fileparts(FileName);
    sessionStruct.FileName  = strcat(FileName,'.',ext);
    sessionStruct.numSegs   = options.numSegs;

    if options.scm     % subtract common mean
        avgSig  = mean(sigDATA.Data,1);
    else
        avgSig  = mean(sigDATA.Data,2);
    end

    numChan = size(sigDATA.Data,1);

    min_eInd = 1;
    max_eInd = numChan;
    if options.inode > 0 || options.iShank > 0
        min_eInd = my_j;
        max_eInd = my_j;
    end

    for eInd=min_eInd:max_eInd % shanks | sites

        wfs = [];
        wts = [];
        chan = [];
        inds = [];
        vals = [];
        spkData = zeros(1,numel(sigDATA.Data(1,:)));


        disp(' ');
        disp(' ');
        disp('_________________________________________________________________________');
        disp(['...BEGIN FILTERING DATA for seg ' num2str(k) ' of ' num2str(options.numSegs) ' | site ' num2str(eInd) ' of ' num2str(numChan)]);

        [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(eInd,options.arrayType);

        % FILTER DATA ON THIS CHANNEL
        disp('_________________________________________________________________________');
        disp(['1) Filtering data on channel ' num2str(eInd) ' | shank:' num2str(col) ', site:' num2str(row)]);

            if options.scm      % subtract common mean
                rawData = sgolayfilt(sigDATA.Data(eInd,:)-avgSig      ,11,21);
            else
                rawData = sgolayfilt(sigDATA.Data(eInd,:)-avgSig(eInd),11,21);
            end

            [lowBandData,hiBandData] = TNC_FilterData(rawData,sigDATA.MetaTags.SamplingFreq,0,0,[1 0]);
            spkData0 = hiBandData.values;

        % REMOVE ARTIFACTS
        len = numel(spkData(1,:));
        K   = numel(spkData(:,1));
        for i=1:len                         
            for j=1:K
                if i > 1 && i < len && abs(spkData0(j,i)) < 1.e-8 && spkData0(j,i-1) < 0 &&  spkData0(j,i-1) < 0
                    spkData(j,i) = (spkData0(j,i-1) + spkData0(j,i+1))/2;
                else
                    spkData(j,i) = spkData0(j,i);
                end
            end
        end

        % EXTRACT WAVEFORMS FROM ALL TRODES FOR EACH UNIQUE EVENT
        disp('2) Detecting and extracting events...');

            % Detect indices of the signal minima corresponding to spikes  
            [spkEvents] = TNC_EventDetect(spkData,samplingRate,snrThresh);

            try 
                assert(numel(spkEvents.inds)>10);
                [spikes]        = TNC_EventExtract(spkData,rawData,spkEvents.inds,[winL,winR],options.dispOn);
                spikes.inds     = spkEvents.inds;

                % CALC SCALAR QUANT FOR SPIKES
                disp('3) Quantify spikes...');
                disp(['k=' num2str(k) ' my_k=' num2str(my_k) ]);
                if k==1 || k==my_k
                    [features]   = TNC_EventQuantSE(spikes,[]);
                    templates(eInd).template = features.template;
                else
                    [features]   = TNC_EventQuantSE(spikes,templates(eInd).template);
                end
            catch 
                spkEvent.inds = [];
                features.params = [];
                features.paramNames = [];
                spikes.wfs = [];
                spikes.inds = [];
            end

        % UPDATE STRUCTURE FOR SPIKE SORTING
        disp('4) Updating the feature structure with new data...');

            featStruct.seg(k).shank(eInd).inds   = spkEvents.inds + (double(k-1).*chunk.*double(30000));
            featStruct.seg(k).shank(eInd).id     = zeros(numel(spkEvents.inds),1); % needed for visualization
            featStruct.seg(k).shank(eInd).ts     = round(featStruct.seg(k).shank(eInd).inds./30);
            featStruct.seg(k).shank(eInd).params = features.params;

            % Splitting params randomly into s given number of folds in the range 2 to 50
            numData = size(features.params, 1);
            oldInds = [1:numData];
            featStruct.seg(k).shank(eInd).random_split(1).folds = oldInds;
            for nf=2:50
               res     = mod(numData, nf); % residual
               if res > 0
                   % Compliment data with (nf - res) NANs to make its size multiple of nf
                   numData1= numData + (nf - res); 
               else
                   numData1= numData;
               end
               oldInds = [1:numData1];
               newInds = oldInds(randperm(numData1)); % a random permutation of data
               featStruct.seg(k).shank(eInd).random_split(nf).folds = ...
                         reshape(newInds, numel(newInds)/nf, nf); % each col is a 'set'
            end

            if k==1 && eInd==1
                featStruct.paramNames = features.paramNames;
                featStruct.rSp = rSp;
                featStruct.cSp = cSp;
                featStruct.electrodeMatrix = electrodeMatrix;
            end

            sessionStruct.seg(k).shank(eInd).wfs  = spikes.wfs;
            sessionStruct.seg(k).shank(eInd).inds = spikes.inds;

        disp('_________________________________________________________________________');
        disp(' ');
        disp(' ');
        disp(' ');


    end % loop in eInd

%-------------------------------------------------------------------------------

function output_ss_data(my_j, my_k, targetName, sessionStruct)
    sessionStruct.program = 'TNC_SSP_ExtractFeatures.m';
    if my_j == 0 && my_k == 0
        disp(['save ' targetName '_ss sessionStruct']);
        eval(['save ' targetName '_ss sessionStruct']);
     elseif my_j >  0 && my_k == 0
        disp(['save ' targetName '_shank' num2str(my_j) '_ss sessionStruct']);
        eval(['save ' targetName '_shank' num2str(my_j) '_ss sessionStruct']);
     elseif my_j == 0 && my_k >  0
        disp(['save ' targetName '_seg' num2str(my_k) '_ss sessionStruct']);
        eval(['save ' targetName '_seg' num2str(my_k) '_ss sessionStruct']);
     elseif my_j >  0 && my_k >  0
        disp(['save ' targetName '_shank' num2str(my_j) '_seg' num2str(my_k) '_ss sessionStruct']);
        eval(['save ' targetName '_shank' num2str(my_j) '_seg' num2str(my_k) '_ss sessionStruct']);
     end

%-------------------------------------------------------------------------------

function output_ft_data(my_j, my_k, targetName, featStruct)
    featStruct.program = 'TNC_SSP_ExtractFeatures.m';
    if my_j == 0 && my_k == 0
        disp(['save ' targetName '_ft featStruct']);
        eval(['save ' targetName '_ft featStruct']);
     elseif my_j >  0 && my_k == 0
        disp(['save ' targetName '_shank' num2str(my_j) '_ft featStruct']);
        eval(['save ' targetName '_shank' num2str(my_j) '_ft featStruct']);
     elseif my_j == 0 && my_k >  0
        disp(['save ' targetName '_seg' num2str(my_k) '_ft featStruct']);
        eval(['save ' targetName '_seg' num2str(my_k) '_ft featStruct']);
     elseif my_j >  0 && my_k >  0
        disp(['save ' targetName '_shank' num2str(my_j) '_seg' num2str(my_k) '_ft featStruct']);
        eval(['save ' targetName '_shank' num2str(my_j) '_seg' num2str(my_k) '_ft featStruct']);
     end

%-------------------------------------------------------------------------------

function output_wf_data(my_j, min_k, max_k, dataFileName, chunk, sessionStruct, wfdata, totChar)
    totChar = numel(dataFileName);
    min_j = 1;
    max_j = size(sessionStruct.seg(min_k).shank,2);
    if my_j > 0
        min_j = my_j;
        max_j = my_j;
    end
    for j=min_j:max_j
        for k=min_k:max_k
            wfdata.seg(k).wfs(numel(sessionStruct.seg(k).shank(j).spikes.wfs)).values = zeros(8,41);
            for i=1:numel(sessionStruct.seg(k).shank(j).spikes.wfs)
                wfdata.seg(k).wfs(i).values = int16(sessionStruct.seg(k).shank(j).spikes.wfs(i).values);
            end
        end
    
        % For getting waveform data 
        wfdata.chunk    = chunk;
        wfdata.shank    = j;
        wfdata.filename = dataFileName;
        wfdata.datatype = 'int16';
   
        if min_k < max_k 
            disp(['save ' dataFileName(1:totChar-4) '_shank' num2str(j) '_wf wfdata']);
            eval(['save ' dataFileName(1:totChar-4) '_shank' num2str(j) '_wf wfdata']);
        else
            disp(['save ' dataFileName(1:totChar-4) '_shank' num2str(j) '_seg' num2str(min_k) '_wf wfdata']);
            eval(['save ' dataFileName(1:totChar-4) '_shank' num2str(j) '_seg' num2str(min_k) '_wf wfdata']);
        end
        clear wfdata    
    end
return

%-------------------------------------------------------------------------------

function output_mc_data(my_j, min_k, my_k, dataFileName, chunk, sessionStruct, paramVals) 
    totChar = numel(dataFileName);
    min_j = 1;
    max_j = size(sessionStruct.seg(min_k).shank,2);
    if my_j > 0
        min_j = my_j;
        max_j = my_j;
    end
    for j=min_j:max_j

        % For use with matclust
        filedata.params     = paramVals(j).values;
        filedata.chunk      = chunk;
        filedata.shank      = j;
        filedata.filename   = dataFileName;
        filedata.paramnames = {'ts','pc1','pc2','pc3','pc4','pc5','pc6','pc7','pc8','minV1','minV2','minV3','minV4','minV5','minV6','minV7','minV8','maxV1','maxV2','maxV3','maxV4','maxV5','maxV6','maxV7','maxV8','energy1','energy2','energy3','energy4','energy5','energy6','energy7','energy8','peVa1','peVa2','peVa3','peVa4','peVa5','peVa6','peVa7','peVa8'}';
        if my_k == 0
            disp(['save ' dataFileName(1:totChar-4) '_shank' num2str(j) '_mc filedata']);
            eval(['save ' dataFileName(1:totChar-4) '_shank' num2str(j) '_mc filedata']);
        else 
            disp(['save ' dataFileName(1:totChar-4) '_shank' num2str(j) '_seg' num2str(my_k) '_mc filedata']);
            eval(['save ' dataFileName(1:totChar-4) '_shank' num2str(j) '_seg' num2str(my_k) '_mc filedata']);
        end
        clear fileData

    end

return

%-------------------------------------------------------------------------------

function [figure] = display_data(i, j, k, spkData, spkEvents, chunk, figure)
    % i = channel id
    % j = shank id
    disp('Displaying data for the current channel...');
    if i==1
        figure(1); clf;
    end

    figure(1);

    subplot(10,1,1:8);
    if rem(j,2)==0
        plot(spkData(i,:)+(double(((j-1).*8)+(i-1)).*1000.),'Color',[0 0.67 1]); hold on; axis tight;
        plot([spkEvents.inds ; spkEvents.inds],[spkData(i,spkEvents.inds).*0 ; spkData(i,spkEvents.inds)]+double(((j-1).*8)+(i-1)).*1000.,'Color',[0 0.67 1]);
        hold on; axis off; axis tight;
    else
        plot(spkData(i,:)+(double(((j-1).*8)+(i-1)).*1000.),'Color',[1 0 0]); hold on; axis tight;
        plot([spkEvents.inds ; spkEvents.inds],[spkData(i,spkEvents.inds).*0 ; spkData(i,spkEvents.inds)]+double(((j-1).*8)+(i-1)).*1000.,'Color',[1 0    0]);
        hold on; axis off; axis tight;
    end

%   subplot(10,1,9);
%   plot(ContData.behavior.sLeverV((chunk.*1000.*(k-1))+1:chunk.*1000.*k),'k');
%   axis off; axis tight;
%   subplot(10,1,10);

%   axis off; axis tight;

    drawnow;
    zoom xon;

    return

%-------------------------------------------------------------------------------

function sigDATA = extract_signal_data_from_mat_file(dataFileName, iSeg, options)                    
    disp(['Loading file ' dataFileName ]);
    sd = load(dataFileName);
    len = size(sd.data, 2);
    ibeg = int32((round(len/options.numSegs)*(iSeg-1))+1);
    iend = int32((round(len/options.numSegs)* iSeg   ));
    disp(['iSeg=' num2str(iSeg) ' ibeg=' num2str(ibeg) ' iend=' num2str(iend) ' len=' num2str(len)]);
    sigDATA.Data = sd.data(:,ibeg:iend)*150;    
    Par_sim.sr = 24000; % Default sampling rate of NeuroCube software
    sigDATA.MetaTags.Duration = size(sigDATA.Data, 2)/Par_sim.sr;
    sigDATA.MetaTags.ChannelCount = size(sigDATA.Data, 1);
    sigDATA.MetaTags.SamplingFreq = Par_sim.sr;
    sigDATA.MetaTags.ElecLabel = ['elec1', 'elec2', 'elec3', 'elec4'];
    if options.fakeData
        sigDATA.Spiketimes = sd.Spiketimes 
        sigDATA.x_neurons  = sd.x_neurons;
        sigDATA.y_neurons  = sd.y_neurons;
        sigDATA.z_neurons  = sd.z_neurons;
        sigDATA.models     = sd.models;
    end
return

%-------------------------------------------------------------------------------
%
% Compute matrix of event indices (one index per channel for each event)
% (Gennady's code)
%
function [ eventInds ] = find_events_across_shank(Spikes, data, maxSpikeDelay, sampleDuration)
    N = numel(Spikes);
    numChannels = size(data, 1);
    spikeBeg = 1;
    spikeVec = [];        % array of vectors of (presumably) related spikes from different channels
    spikeVecInd = 0;      % index of vector in spikeVec
    maxSpikeDelay_scans = maxSpikeDelay/sampleDuration;
    disp(['maxSpikeDelay_scans=' num2str(maxSpikeDelay_scans)]);
    for i=1:(numel(Spikes)-1)
        if spikeBeg < 0
            spikeBeg = i;
        end
        % Ending a group of closely spaced spikes from different channels
        if spikeBeg > 0 && (Spikes(i+1).position - Spikes(i).position)*sampleDuration > maxSpikeDelay  
            spikeEnd = i;
            % Check if a spike is confirmed by a closely spaced spike from another channel
            spikeConfirmed = 0;
            channelsPresent = [];
%           disp(['spikeBeg_pos=' num2str(Spikes(spikeBeg).position) ' spikeEnd=' num2str(Spikes(spikeEnd).position) ' maxSpikeDelay=' num2str(maxSpikeDelay)]);
            for j=spikeBeg:spikeEnd
                c = Spikes(j).channel;
                if length(find(channelsPresent == c)) == 0
                    channelsPresent = [channelsPresent c];
                end
            end
            if numel(channelsPresent) > 1
                spikeConfirmed = 1;
            end
%           disp(['spikeBeg=' num2str(spikeBeg) ' spikeEnd=' num2str(spikeEnd) ...
%                 ' channelsPresent=' num2str(channelsPresent) ' spikeConfirmed=' num2str(spikeConfirmed) ...
%                 ' spikeVecInd=' num2str(spikeVecInd) ' numel(spikeVec)=' num2str(numel(spikeVec))]);
            if spikeConfirmed % closely spaced spikes are present in at least two channels
                spikeVecInd = spikeVecInd + 1;
                spikeVec(spikeVecInd).spikePosition = ...
                   round((Spikes(spikeBeg).position+...
                          Spikes(spikeEnd).position)/2);
                spikeVec(spikeVecInd).channelPos = -1*ones(1, numChannels);  % initialization
                channelsPresent = [];
                for j=spikeBeg:spikeEnd
                    c = Spikes(j).channel;
                    spikeVec(spikeVecInd).channelPos(c) = Spikes(j).position;
                end
                % Update spikeVec: set negative channelPos to approximate positive values
                for c=1:numChannels
                    if spikeVec(spikeVecInd).channelPos(c) < 0
                        ind_pos_found = find(spikeVec(spikeVecInd).channelPos > 0);
                        if length(ind_pos_found) > 0
                            spikeVec(spikeVecInd).channelPos(c) = ...
                                round(mean(spikeVec(spikeVecInd).channelPos(ind_pos_found)));
                        else
                            sbeg = max(1,             spikeVec(spikeVecInd).spikePosition - round(maxSpikeDelay/2/sampleDuration));
                            send = min(numel(Spikes), spikeVec(spikeVecInd).spikePosition + round(maxSpikeDelay/2/sampleDuration));
                            inds = find(data(c,sbeg:send) == min(data(c,sbeg:send)));
                            pos  = min(inds);
                            spikeVec(spikeVecInd).channelPos(c) = sbeg + pos - 1;
                        end
                    end
                end
            end % spike confirmed
            spikeBeg = -1;
        end
    end
    eventInds = zeros(numChannels, numel(spikeVec));
    for n=1:numel(spikeVec)
         eventInds(:,n) = spikeVec(n).channelPos(:);
    end

%-------------------------------------------------------------------------------

function Y = determine_true_class_labels(Spiketimes, options)
    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
    Y = [];
    if options.fakeData
        numSpikeTimes = 0;
        disp(['size(Spiketimes)=' num2str(size(Spiketimes))]);
        for j=1:size(Spiketimes,2)
            num_spikes = length(Spiketimes{j});
            numSpikeTimes = numSpikeTimes + num_spikes;
        end
        prev_spiketime = 0;
        for i = 1:numSpikeTimes
            best_class = 0;
            best_spike_time = Inf;            
            for j=1:size(Spiketimes,2)
                for k = 1:length(Spiketimes{j})
                    if Spiketimes{j}(k) > prev_spiketime 
                        if  best_spike_time  > Spiketimes{j}(k)
                            best_spike_time = Spiketimes{j}(k);
                            best_class = j;
                        end  
                    end
                end 
            end 
            if best_class == 0
                break;
            end
            if best_spike_time < Inf
                prev_spiketime = best_spike_time;
            end
            Y = [Y Alphabet(best_class)];
        end
        disp(['numSpikeTimes=' num2str(numSpikeTimes) ...
              ' length(Y)=' num2str(length(Y)) ...
              ' num duplicated SpikeTimes=' num2str(numSpikeTimes - length(Y))]);
    end


