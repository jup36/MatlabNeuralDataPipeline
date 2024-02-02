function behaviorTimestampsPP(p)
%behaviorTimestamps

%addpath(genpath(''))
%filePath = '/Volumes/RAID2/parkj/NeuralData/js2.0/WR25/101718';
%p = parse_input_Js(filePath, varargin ); % parse input
%p = parse_input_Js(filePath, {} ); % use this line instead when running line-by-line
% To convert datenum to datetime use, datetime(vFronFileStartDatenum(tempVFronI),'ConvertFrom','datenum')

cd(p.Results.filePath)
% 
% check the folder whether there's rez file already
if ~isempty(dir(fullfile(p.Results.filePath,'BehVariablesPP.mat'))) % if the BehVariablesPP.mat file already exists in the filePath
    answer = questdlg('BehVariablesPP.mat already exists, Would you like to replace it?','Choice','Replace','Cancel','Cancel');
    switch answer
        case 'Replace'
        case 'Cancel'
            return
    end
end

binFile = dir(fullfile(p.Results.filePath,'*.nidq.bin')); % look for nidq.bin file
binFileI = cell2mat(cellfun(@(a) isletter(a(1)), {binFile.name}, 'un', 0));
binFile = binFile(binFileI);

if length(binFile)>1 || isempty(binFile)
    error('File could not be found or multiple nidq.bin files exist!');
end

binName = binFile.name;

% Parse the corresponding metafile
meta  = ReadMeta(binName, p.Results.filePath); % get the meta data (structure)
channels = textscan(meta.acqMnMaXaDw,'%n %n %n %n','Delimiter',',');
nSamp = floor(SampRate(meta));          % sampling rate (default: 25kHz)
totalTimeSecs = str2double(meta.fileTimeSecs); % total duration of file in seconds

if ~isempty(dir(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'))) && p.Results.reReadBin==false % if the gainCorrectRawTraces.mat file already exists in the filePath
    load(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat')) % if there are gaincorrectedrawtraces already saved, just load them
else
    % Specify the relevant behavioral channel numbers
    cmosExpCh    = channels{1}+p.Results.cmosExpCh;   % ai0
    cmosTrigCh   = channels{1}+p.Results.cmosTrigCh;  % ai1
    speakerCh    = channels{1}+p.Results.speakerCh;   % ai3
    lickCh       = channels{1}+p.Results.lickCh;      % ai4
    lick2Ch      = channels{1}+p.Results.lick2Ch;     % ai6
    bodyCamCh    = channels{1}+p.Results.bodyCamCh;   % ai5
    faceCamCh    = channels{1}+p.Results.faceCamCh;   % ai7
    photoDiodeCh = channels{1}+p.Results.photoDiodeCh; % ai14
    digitCh       = channels{1}+p.Results.digitCh;     % 17 when there are 16 analog inputs recorded

    % preallocate the behavioral data arrays
    cmosExp    = zeros(1,floor(totalTimeSecs*nSamp));
    cmosTrig  = zeros(1,floor(totalTimeSecs*nSamp));
    speaker    = zeros(1,floor(totalTimeSecs*nSamp));
    lick       = zeros(1,floor(totalTimeSecs*nSamp));
    lick2      = zeros(1,floor(totalTimeSecs*nSamp));
    bodyCam    = zeros(1,floor(totalTimeSecs*nSamp));
    faceCam    = zeros(1,floor(totalTimeSecs*nSamp));
    photoDiode = zeros(1,floor(totalTimeSecs*nSamp));
    digit      = zeros(1,floor(totalTimeSecs*nSamp));

    for i = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
        start = floor(i * nSamp);
        last  = floor((i+1) * nSamp);

        tempDataArray   = ReadBin(start, nSamp, meta, binName, p.Results.filePath); % read bin data for each second
        tempCmosExp     = tempDataArray(cmosExpCh, :);    
        tempCmosTrig    = tempDataArray(cmosTrigCh, :); 
        tempSpeaker     = tempDataArray(speakerCh, :);    
        tempLick        = tempDataArray(lickCh, :);       
        tempLick2       = tempDataArray(lick2Ch, :);     
        tempBodyCam     = tempDataArray(bodyCamCh, :);    
        tempFaceCam     = tempDataArray(faceCamCh, :);    
        tempPhotoDiode  = tempDataArray(photoDiodeCh, :); 
        tempDigit       = tempDataArray(digitCh, :); 

        cmosExp(1,start+1:last)     = tempCmosExp; % accumulated the decimated data second-by-second
        cmosTrig(1,start+1:last)    = tempCmosTrig;
        speaker(1,start+1:last)     = tempSpeaker;
        lick(1,start+1:last)        = tempLick;
        lick2(1,start+1:last)       = tempLick2;
        bodyCam(1,start+1:last)     = tempBodyCam;
        faceCam(1,start+1:last)     = tempFaceCam;
        photoDiode(1,start+1:last)  = tempPhotoDiode;
        digit(1, start+1:last)      = tempDigit; 
        fprintf('processed %d\n', i+1)
    end
    clearvars i
    
    %digitBin = dec2bin(digit);
    %water = str2num(digitBin(:, 1)); %'str2double' doesn't work here! 
    %airpuff = str2num(digitBin(:, 2)); 
    
    % note that the water strobe is connected to p0.2, and the airpuff
    % strobe is connected to p0.1. As the airpuff channel comes first
    % ordinally, it's counterintuitive that water is assigned to be '2' and
    % the airpuff is assigned to be '4'. (That said, it happens to be the case.)  
    digit_norm = digit - mode(digit); 
    digit_norm(digit_norm<0)=0; % replace negative values, a rare artifact, with 0. 
    
    if size(unique(digit_norm),2)==3 % p0.1 -> 2^1, p0.2 -> 2^2, others -> zero
      water = digit_norm==2;   % p0.2 -> 2^2
      airpuff = digit_norm==4; % p0.1 -> 2^1
    elseif size(unique(digit_norm),2)==2 % most likely, p0.2 -> 2^1, others -> zero
      water = digit_norm==2;   % p0.2 -> 2^2
      airpuff=[];
    else
      error("Digit channels have an unknown input(s) or artifact(s)!"); 
    end

%     figure; hold on; 
%     plot(water); 
%     plot(airpuff); 
% 
    % Gain correction for analog channnels 
    cmosExp     = GainCorrectNI(cmosExp, 1, meta);      
    cmosTrig    = GainCorrectNI(cmosTrig, 1, meta); 
    speaker     = GainCorrectNI(speaker, 1, meta);   
    lick        = GainCorrectNI(lick, 1, meta);   
    lick2       = GainCorrectNI(lick2, 1, meta);  
    bodyCam     = GainCorrectNI(bodyCam, 1, meta);
    faceCam     = GainCorrectNI(faceCam, 1, meta);
    photoDiode  = GainCorrectNI(photoDiode, 1, meta);
    clearvars temp*
    save('gainCorrectRawTraces', 'meta', 'cmosExp', 'cmosTrig', 'speaker', 'lick', 'lick2', ...
        'bodyCam', 'faceCam', 'photoDiode', 'digit', 'water', 'airpuff')
end


end
