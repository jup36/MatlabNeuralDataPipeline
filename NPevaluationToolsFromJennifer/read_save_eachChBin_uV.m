function read_save_eachChBin_uV(saveName, varargin )

%saveName = 'PT08_062218_eachChRawBit';

if (length(varargin) == 0)
    % probe type: NP1.0 = 1; NP2.0 SS = 21, NP2.0 = 24;
    probeType = 1;
    % If calling with no parameters, specify here
    % channels to exclude from histograms and summary statistics.
    % should exclude reference and "dead" channels
    exChan = [141,191]; % Should always exclude the reference channel (191 for NP 1.0, 127 for NP 2.0.
    % z range on probe to include in histogram and summary stats.
    % NP1.0 full z range = 0-3840 um
    % NP2.0 rull z range = 0-2880 um
    % To compare NP1.0 data to a similar tissue volume measured with NP2.0
    % set zMax = 3000.
    zMin = -inf;
    zMax = inf;    
else
    inputCell = varargin(1);
    probeType = inputCell{1};
    zMin = -inf;
    zMax = inf; 
    exChan = [141,191];
end

% set dataType
% 0 for neuropixels 1.0, 1 for NP 2.0. 
if probeType == 1
    dataType = 0;
elseif (probeType == 21) || (probeType == 24)
    dataType = 1;
else
    fprintf( "unknown probe type\n" );
    return
end

% %thresholding params
% bUseConstThresh = 1;
% constThresh = 80;    %in uV
% qqFactor = 5; % for JRClust-style thresholding; ignored if bUseConstThresh = 1

% %other params for JRClust type merging
% hCfg.evtDetectRad = 50; %in um, distance to look for duplicate peaks
% hCfg.refracIntSamp = 7; %in samples, distance in time to look for duplicat peaks

% for either data type
nchan = 385;
dataChan = 384;

maxTime = 600; %in sec; takes the first maxTime seconds in the file

% get data file from user. 
[fileName,fileDir]=uigetfile('*.bin', 'Select file' );

cd(fileDir);

% Build name for output file, which will be a matlab structure
[~,inName,~] = fileparts(fileName);
outName = [inName,'_allDataChMicroV.mat'];

% get coordinate file from user, chan, X, Y, shank index, tab delimited
% X coordiantes across shanks should reflect the distance between the
% shanks. should include all channels in teh file (e.g. include the ref
% channels, if present).
% Get site coords file from user
[coordName,coordDir]=uigetfile('*.txt', 'Select coords file' );
cID = fopen( fullfile(coordDir,coordName), 'r' );
hCfg.siteLoc = zeros(dataChan,2);
shank = zeros(dataChan,1);
for i = 1:dataChan
    tline = fgetl(cID);
    currDat = sscanf(tline, '%d%d%d%d');
    hCfg.siteLoc(i,1) = currDat(2);
    hCfg.siteLoc(i,2) = currDat(3);
    shank(i) = currDat(4);
end
fclose(cID);

switch dataType
    
    case 0
        fs = 30000;
        nBit = 1024; %2^10   %amplitudes will be histogrammed over all bits
        uVPerBit = 2.34375; % bitPerMicrovolt = 2^10 / (1.2 / gain * 1000000); ADC resolution 10 bits
     
    case 1
        fs = 30000;
        nBit = 16384; %2^14    %amplitudes will be histogrammed over all bits
        uVPerBit = 0.7629;
        
    otherwise
        fprintf( 'unknown dataType\n' );
        return;
end

edges = (0:nBit);   %nBit+1 edges -> nBit bins. Overflow goes to top/bottom bin


    fileStats = dir(fullfile(fileDir,fileName));

    fileSize = fileStats.bytes;
    
    fm = memmapfile(fullfile(fileDir,fileName),'Format','int16');
    % the size of fm.Data = fileSize/2
    % the # of samples (fileSize/2)/385, this should be time(s) multiplied by sampling rate (Hz)
    
    maxBatchBytes = 1e8;
    maxBatchSamples = floor(maxBatchBytes/(nchan*2));
    batchBytes = maxBatchSamples*nchan*2;
    
    nneighBelow = 1;

    runSec = fileSize/(nchan*2*fs);
    fprintf('Run length in seconds: %.2f\n',runSec);

    fileSamples = fileSize/(nchan*2);
    
    if runSec < maxTime
        nBatch = floor(fileSamples/maxBatchSamples) + 1;
    else
        useSamples = maxTime*fs;
        nBatch = floor(useSamples/maxBatchSamples) + 1;
    end
        
    batchSamples = maxBatchSamples;

    fprintf('nBatch, batchSamples: %d, %d\n', nBatch, batchSamples );
    
    %will be processing the nBatch-1 "full" batches (saves complexity)
    analyzedSec = (nBatch-1)*batchSamples/fs;

    sizeA = [nchan,batchSamples];

    %set up arrays to hold amplitude histograms for each channel
    %for NP1.0 data, max amplitude is 10 bits (2^10); for NP 2.0, max is 14
    %bits (2^14)
    ampHist = zeros(dataChan,nBit,'int32');   
    tempHist = zeros(1,nBit,'int32');
    
    %Note -- should we try to preallocate space for these?
    %Would need a way to estimate the number of spikes
    %Could consider stopping after some number of spikes are reached.
    allTimes = [];
    allAmps = [];
    allSites = [];

dataArrayC = []; 
    batch16bit = batchBytes/2; % batchBytes = maxBatchSamples*nchan*2;
    for i = 1:nBatch-1
         fprintf( 'Reading batch: %d of %d\n', i, nBatch-1 );
         batchOffset = int64((i-1)*batch16bit) + 1;
         dataArray = fm.Data(batchOffset:batchOffset+batch16bit-1); % data are in bit value here, transform them to uV by multiplying (1/bitPerMicrovolt)
         dataArray = reshape(dataArray,sizeA); % sizeA = [nchan,batchSamples];
         dataArrayC = [dataArrayC, dataArray]; 
    end
    
    for c = 1:dataChan
        rawDataInBit = dataArrayC(c,:);    
        save(fullfile(fileDir,strcat(saveName,sprintf('Ch#%d',c))),'rawDataInBit','uVPerBit')
        fprintf( 'Saved Ch#: %d\n', c);
    end
 % common average subtraction 
 
 
%  dataInMs = floor(length(dataArrayC)/30);  
%  reach0Cut = reach0(1:dataInMs); % reach0 corresponding to dataArrayC
%  interval = linspace(1,length(reach0Cut),length(dataArrayC));
%  reach0CutInt = interp1(1:length(reach0Cut), reach0Cut, interval); 
%  scaleFactor = 300/max(abs(reach0(:))); 
% 
%  
%  figure; plot(dataArrayC(346,:).*uVPerBit); hold on; plot(reach0CutInt.*double(scaleFactor))
%  ylim([-300, 300])
%  
 
 