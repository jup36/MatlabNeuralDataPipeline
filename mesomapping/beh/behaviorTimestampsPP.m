function behaviorTimestampsPP(p)
%behaviorTimestamps

%addpath(genpath(''))
%filePath = '/Volumes/RAID2/parkj/NeuralData/js2.0/WR25/101718';
%p = parse_input_Js(filePath, varargin ); % parse input
%p = parse_input_Js(filePath, {} ); % use this line instead when running line-by-line
% To convert datenum to datetime use, datetime(vFronFileStartDatenum(tempVFronI),'ConvertFrom','datenum')

cd(p.Results.filePath)

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
    cmosTrigCh   = channels{1}+p.Results.cmosTrigCh;  % ai0
    faceCamCh    = channels{1}+p.Results.faceCamCh;   % ai1
    speakerCh    = channels{1}+p.Results.speakerCh;   % ai2
    treadMillCh  = channels{1}+p.Results.treadMillCh; % ai4
    sideGreenCh  = channels{1}+p.Results.sideGreenCh; % ai6
    topRedCh     = channels{1}+p.Results.topRedCh;    % ai7

    digitCh      = channels{1}+p.Results.digitCh;     % 9 when there are 8 analog inputs recorded

    % preallocate the behavioral data arrays
    cmosTrig   = zeros(1,floor(totalTimeSecs*nSamp));
    faceCam    = zeros(1,floor(totalTimeSecs*nSamp));
    speaker    = zeros(1,floor(totalTimeSecs*nSamp));
    treadMill  = zeros(1,floor(totalTimeSecs*nSamp));
    sideGreen  = zeros(1,floor(totalTimeSecs*nSamp));
    topRed     = zeros(1,floor(totalTimeSecs*nSamp));

    digit      = zeros(1,floor(totalTimeSecs*nSamp));

    for i = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
        start = floor(i * nSamp);
        last  = floor((i+1) * nSamp);

        tempDataArray = ReadBin(start, nSamp, meta, binName, p.Results.filePath); % read bin data for each second
        tempCmosTrig  = tempDataArray(cmosTrigCh, :);
        tempFaceCam   = tempDataArray(faceCamCh, :);
        tempSpeaker   = tempDataArray(speakerCh, :);
        tempTreadMill = tempDataArray(treadMillCh, :);
        tempSideGreen = tempDataArray(sideGreenCh, :);
        tempTopRed    = tempDataArray(topRedCh, :);
        tempDigit       = tempDataArray(digitCh, :);

        cmosTrig(1,start+1:last)    = tempCmosTrig;
        faceCam(1,start+1:last)     = tempFaceCam;
        speaker(1,start+1:last)     = tempSpeaker;
        treadMill(1,start+1:last)   = tempTreadMill;
        sideGreen(1,start+1:last)   = tempSideGreen;
        topRed(1,start+1:last)      = tempTopRed;
        digit(1, start+1:last)      = tempDigit;
        fprintf('processed %d\n', i+1)
    end
    clearvars i

    % parse digits
    digit_norm = digit - mode(digit);
    digit_norm(digit_norm<0)=0; % replace negative values, a rare artifact, with 0.

    digitC = digit_decompose(digit_norm, nSamp);
    fprintf("Finished parsing digital channels!\n")

    water = digitC{1,1};
    airpuff = digitC{2,1};
    blueLED = digitC{3,1};
    limeLED = digitC{4,1};
    greenLED = digitC{5,1};
    redLED = digitC{6,1};
    lick = digitC{7,1};
    manualWater = digitC{8,1};

    % Gain correction for analog channnels
    cmosTrig    = GainCorrectNI(cmosTrig, 1, meta);
    faceCam     = GainCorrectNI(faceCam, 1, meta);
    speaker     = GainCorrectNI(speaker, 1, meta);
    treadMill   = GainCorrectNI(treadMill, 1, meta);
    sideGreen   = GainCorrectNI(sideGreen, 1, meta);
    topRed      = GainCorrectNI(topRed, 1, meta);

    clearvars temp*
    save(fullfile(p.Results.filePath, 'gainCorrectRawTraces'), 'meta', ...
        'cmosTrig', 'faceCam', 'speaker', 'treadMill', 'topRed', 'sideGreen', ...
        'digit_norm', 'water', 'airpuff', 'blueLED', 'limeLED', 'greenLED', 'redLED', ...
        'lick', 'manualWater');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function parses the data stream in the digit channel
    function [out_timeseries] = digit_decompose(digit_norm, nSamp)
        % Decompose a time series with sums of 2^0, 2^1, ..., 2^7 into 8 time series
        % and filter out values that last shorter than nSamp/10.
        % Inputs:
        %   digit_norm - input time series (vector)
        %   nSamp - sampling rate
        % Outputs:
        %   out_timeseries - cell array of 8 time series vectors (binary), where each
        %                    indicates the presence of each respective power of 2.

        % Define the minimum duration (in samples) that a value must last
        min_duration = round(nSamp / 200); % to cutoff short noises that last less than 5ms

        % Initialize a filtered version of digit_norm
        digit_filtered = digit_norm;

        % Find changes in the input time series
        changes = [1, find(diff(digit_norm) ~= 0) + 1, length(digit_norm) + 1];

        % Filter out any segment that lasts shorter than min_duration
        for ii = 1:(length(changes) - 1)
            segment_start = changes(ii);
            segment_end = changes(ii+1) - 1;
            if (segment_end - segment_start + 1) < min_duration
                % Set the values of this short segment to 0 (treated as noise)
                digit_filtered(segment_start:segment_end) = 0;
            end
        end

        % Define the 8 powers of 2 (2^0, 2^1, ..., 2^7)
        powers_of_2 = 2.^(0:7);

        % Initialize output timeseries as a matrix with 8 rows (one for each power of 2)
        out_timeseries = zeros(8, length(digit_filtered));

        % Iterate only at change points
        for ii = 1:(length(changes) - 1)
            segment_start = changes(ii);
            segment_end = changes(ii+1) - 1;
            value = digit_filtered(segment_start);

            if value == 0
                % If value is zero, all powers of 2 are zero
                current_decomposition = zeros(8, 1);
            elseif all(ismember(dec2bin(value), ['0', '1'])) && value >= 0 && value < 256
                % Decompose the value only if it's valid
                current_decomposition = bitand(value, powers_of_2) ~= 0;
            else
                % Treat it as noise, skip decomposition
                current_decomposition = zeros(8, 1);
            end

            % Assign the current decomposition to the timeseries for the current segment
            for t = segment_start:segment_end
                if t == segment_start
                    % For the first point in the segment, use current decomposition
                    out_timeseries(:, t) = current_decomposition;
                else
                    % For subsequent points in the segment, inherit the previous values
                    out_timeseries(:, t) = out_timeseries(:, t-1);
                end
            end
        end

        % Convert the output matrix to a cell array with 8 separate timeseries
        out_timeseries = mat2cell(out_timeseries, ones(1,8), length(digit_filtered));
    end


end
