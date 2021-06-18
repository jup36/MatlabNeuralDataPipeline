function blockSVD(opts)
% example code to get raw data from a given imaging experiment and perform 
% blockwise dimensionality reduction.
% Data channels are either blue (1) or violet (2) which contains intrinsic
% signals.
% opts.fPath denotes the path to the example data folder, containing raw video 
% data as .mj2 files and according timestamps for each frame. 
% nrBlocks is the total number of blocks used for the svd. Sqrt of blocks
% has to be an even number - nrBlocks is rounded down if needed.
% overlap determines the number of pixels with which individual blocks are
% overlapping to avoid edge effects.

if ~strcmpi(opts.fPath(end),filesep)
    opts.fPath = [opts.fPath filesep];
end

% these should be provided as inputs
nrBlocks = opts.nrBlocks;
overlap = opts.overlap;

%construct path to data folder and give some basic info
opts.fName = 'Frames'; %name of imaging data files.
opts.plotChans = false; %plot results for blue and violet channel when loading raw data. Will make the code break on the HPC.
opts.loadRaw = true; %flag to use raw data instead of .mj2 files.
opts.verbosity = false; %flag to supress warnings from 'splitChannel / splitVideo' code.
opts.sRate = 30; %single-channel sampling rate in Hz

% analog data lines
opts.stimLine = 6; %analog line that contains stimulus trigger.
opts.trigLine = [8 9]; %analog lines for blue and violet light triggers.

if nrBlocks ~= floor(sqrt(nrBlocks))^2
    fprintf('Chosen nrBlocks (%d) cannot be squared. Using %d instead.\n', nrBlocks, floor(sqrt(nrBlocks))^2)
    nrBlocks = floor(sqrt(nrBlocks))^2;
end

disp('==============='); tic;
disp(opts.fPath); 
fprintf('Dual-channel recording - Using %d blocks for SVD\n', nrBlocks);
disp(datestr(now));

%% check if .dat or .mj2 files are present. Use slower splitvideo code for latter case.
rawCheck = dir([opts.fPath filesep opts.fName '*dat']);
vidCheck = dir([opts.fPath filesep opts.fName '*mj2']);
analogCheck = dir([opts.fPath filesep 'Analog_*.dat']);

if size(rawCheck,1) == size(analogCheck,1)
    fileCnt = size(rawCheck,1);
    for iFiles = 1 : fileCnt %go through files and check for trialNr. this is used for transformation to .mat files later.
        temp = textscan(rawCheck(iFiles).name,'%s%f%s','Delimiter','_');
        trials(iFiles) = temp{2};
    end
elseif size(vidCheck,1) == size(analogCheck,1)
    fileCnt = size(vidCheck,1);
    for iFiles = 1 : fileCnt %go through files and check for trialNr. this is used for transformation to .mat files later.
        temp = textscan(vidCheck(iFiles).name,'%s%f%s','Delimiter','_');
        trials(iFiles) = temp{2};
    end
    opts.loadRaw = false;
    disp(['Using ' num2str(length(trials)) ' .mj2 files.']);
else
    error('Unequal number of imaging and analog data files. Aborted')
end
trials = sort(trials);

%% get reference images for motion correction
if opts.loadRaw
    [blueData,~,hemoData] = Widefield_SplitChannels(opts,trials(1));
else
    [blueData,~,hemoData] = Widefield_SplitVideo(opts,trials(1));
end
blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment
save([opts.fPath 'blueRef.mat'],'blueRef');
hemoData = single(squeeze(hemoData));
hemoRef = fft2(median(hemoData,3)); %violet reference for alignment
save([opts.fPath 'hemoRef.mat'],'hemoRef');

%% get index for individual blocks
indImg = reshape(1:numel(blueRef),size(blueRef)); %this is an 'image' with the corresponding indices
blockSize = ceil((size(blueRef) + repmat(sqrt(nrBlocks) * overlap, 1, 2))/sqrt(nrBlocks)); %size of each block
blockInd = cell(1, nrBlocks);

Cnt = 0;
colSteps = (0 : blockSize(1) - overlap : size(blueRef,1)) + 1; %steps for columns
rowSteps = (0 : blockSize(2) - overlap : size(blueRef,2)) + 1; %steps for rows
for iRows = 1 : sqrt(nrBlocks)
    for iCols = 1 : sqrt(nrBlocks)
        
        Cnt = Cnt + 1;
        % get current block and save index as vector
        colInd = colSteps(iCols) : colSteps(iCols) + blockSize(1) - 1; 
        rowInd = rowSteps(iRows) : rowSteps(iRows) + blockSize(2) - 1;
        
        colInd(colInd > size(blueRef,1)) = [];
        rowInd(rowInd > size(blueRef,2)) = [];
        
        cBlock = indImg(colInd, rowInd);
        blockInd{Cnt} = cBlock(:);
        
    end
end
save([opts.fPath 'blockInd.mat'],'blockInd');

%% perform image alignement for separate channels and collect data in mov matrix
blueAvg = zeros([size(blueData,1), size(blueData,2), fileCnt],'uint16'); %average for mean correction. Collect single session averages to verify correct channel separation.
hemoAvg = zeros([size(blueData,1), size(blueData,2), fileCnt],'uint16');
blueFrameTimes = cell(1, fileCnt);
hemoFrameTimes = cell(1, fileCnt);
if ~exist([opts.fPath 'blockData'], 'dir')
    mkdir([opts.fPath 'blockData']);
end

try
    blueRef = gpuArray(blueRef);
    hemoRef = gpuArray(hemoRef);
end

for iTrials = 1:fileCnt
    if opts.loadRaw
        [blueData,blueTimes,hemoData,hemoTimes] = Widefield_SplitChannels(opts,trials(iTrials));
    else
        [blueData,blueTimes,hemoData,hemoTimes] = Widefield_SplitVideo(opts,trials(iTrials));
    end
    
    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
    end
    
    try
        blueData = gpuArray(blueData);
        hemoData = gpuArray(hemoData);
    end
    
    %perform image alignment for both channels
    for iFrames = 1:size(blueData,3)
        [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
        blueData(:, :, iFrames) = abs(ifft2(temp));
        
        [~, temp] = Widefield_dftregistration(hemoRef, fft2(hemoData(:, :, iFrames)), 10);
        hemoData(:, :, iFrames) = abs(ifft2(temp));
    end
    blueData = gather(blueData);
    hemoData = gather(hemoData);
    
    % keep avg for each trial to check if channels were separated correctly
    blueAvg(:,:,iTrials) = mean(blueData,3);
    hemoAvg(:,:,iTrials) = mean(blueData,3);
    
    %keep timestamps for all frames
    blueFrameTimes{iTrials} = blueTimes;
    hemoFrameTimes{iTrials} = hemoTimes;

    if rem(iTrials,10) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
    
    blueData = reshape(blueData, [], size(blueData,3));
    hemoData = reshape(hemoData, [], size(hemoData,3));
    
    % save data in individual blocks. single file for each trial/block. Will delete those later.
    for iBlocks = 1:nrBlocks
        bBlock = blueData(blockInd{iBlocks}, :);
        hBlock = hemoData(blockInd{iBlocks}, :);
        save([opts.fPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'bBlock', '-v6');
        save([opts.fPath 'blockData' filesep 'hemoBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'hBlock', '-v6');
    end
end
clear blueData hemoData blueTimes hemoTiems blueRef hemoRef 
save([opts.fPath 'trials.mat'],'trials'); %save trials so order of analysis is consistent

%save frametimes for blue/hemo trials
save([opts.fPath 'blueFrameTimes.mat'],'blueFrameTimes', 'trials');
save([opts.fPath 'hemoFrameTimes.mat'],'hemoFrameTimes', 'trials');

%save averages in case you need them later
save([opts.fPath 'blueAvg.mat'],'blueAvg');
save([opts.fPath 'hemoAvg.mat'],'hemoAvg');

%take average over all trials for subsequent mean correction
blueAvg = mean(single(blueAvg),3);
hemoAvg = mean(single(hemoAvg),3);
    
%% subtract and divide each block by and compress with SVD
bU = cell(nrBlocks,1); bV = cell(nrBlocks,1);
for iBlocks = 1 : nrBlocks
    
    % rebuild current block from all trials
    Cnt = 0;
    allBlock = NaN(size(blockInd{iBlocks},1), size(cat(1,blueFrameTimes{:}),1), 2, 'single');
    for iTrials = 1:fileCnt
        load([opts.fPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'bBlock');
        load([opts.fPath 'blockData' filesep 'hemoBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'hBlock');
        
        allBlock(:, Cnt + (1:size(bBlock,2)), 1) = single(bBlock);
        allBlock(:, Cnt + (1:size(bBlock,2)), 2) = single(hBlock);
        Cnt = Cnt + size(bBlock,2);
    end
    
    % compute dF/F   
    allBlock(:,:,1) = bsxfun(@minus, allBlock(:,:,1), blueAvg(blockInd{iBlocks}));
    allBlock(:,:,1) = bsxfun(@rdivide, allBlock(:,:,1), blueAvg(blockInd{iBlocks}));
    allBlock(:,:,2) = bsxfun(@minus, allBlock(:,:,2), hemoAvg(blockInd{iBlocks}));
    allBlock(:,:,2) = bsxfun(@rdivide, allBlock(:,:,2), hemoAvg(blockInd{iBlocks}));
    
    % run SVD on current block
    [bU{iBlocks}, s, bV{iBlocks}] = svd(reshape(allBlock,size(allBlock,1),[]), 'econ');
    bV{iBlocks} = s * bV{iBlocks}'; %multiply S into V, so only U and V from here on
    bU{iBlocks} = bU{iBlocks}(:, 1:opts.blockDims); %reduce number of components
    bV{iBlocks} = bV{iBlocks}(1:opts.blockDims, :); %reduce number of components
    
    if rem(iBlocks, round(nrBlocks / 5)) == 0
        fprintf(1, 'Loading block %d out of %d\n', iBlocks, nrBlocks);
    end
end

% save blockwise SVD data from both channels
save([opts.fPath 'bV.mat'], 'bU', 'bV', 'blockInd', 'opts', '-v7.3');
toc;

