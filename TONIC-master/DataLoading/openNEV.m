function NEV = openNEV(varargin)

% Opens an .nev file for reading, returns all file information in a NEV
% structure. Works with File Spec 2.1 & 2.2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use OUTPUT = openNEV(fname, 'read', 'report', 'noparse', 'nowarning', 'precision').
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   fname:        Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file using an open 
%                 file user interface. 
%                 DEFAULT: Will open Open File UI.
%
%   'read':       Will read the data if user passes this argument.
%                 DEFAULT: will not read data.
%
%   'report':     Will show a summary report if user passes this argument.
%                 DEFAULT: will not show report.
%
%   'noparse':    The code will not parse the experimental parameters.
%                 See below for format.
%                 DEFAULT: will parse the parameters.
%
%   'nowarning':  The code will not give a warning if there is an error in
%                 parsing.
%                 DEFAULT: will give warning message.
%
%   'nosave':     The code will not save a copy of the NEV structure as a
%                 MAT file. By default the code will save a copy in the same
%                 folder as the NEV file for easy future access.
%                 DEFAULT: will save the MAT file.
%
%   'nomat':      Will not look for a MAT file. This option will force
%                 openNEV to open a NEV file instead of any available MAT
%                 files.
%                 DEFAULT: will load the MAT file if available.
%
%   'compact':    If specified the spike data is stored in 'int16' type 
%                 instead of default double precision.
%
%   OUTPUT:       Contains the NEV structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   USAGE EXAMPLE: 
%   
%   openNEV('report','read','c:\data\sample.nev');
%
%   In the example above, the file c:\data\sample.nev will be used. A
%   report of the file contents will be shown. The digital data will be
%   parsed. The data needs to be in the proper format (refer below).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIGITAL PARAMETERS/MARKERS FORMAT:
%
%   The experimental parameters need to be in the following format for the
%   code to properly parse them:
%
%   *Label:Parameter1=value1;Parameter2=value2;Parameter3=value3;#
%
%   EXAMPLES:
%   *ExpParameter:Intensity=1.02;Duration=400;Trials=1;PageSegment=14;#
%   *Stimulation:StimCount=5;Duration=10;#
%
%   The above line is an "ExpParameter". The parameters are, "Intensity,
%   Duration, Trials, and PageSement." The values of those parameters are,
%   respectively, "1.02, 400, 1, and 14." The second example is a 
%   "Stimulation". The name of the parameters are "StimCount" and
%   "Duration" and the values are "5" and "10" respectively.
%
%   It can also read single value markers that follow the following format.
%
%   *MarkerName=Value;#
%
%   EXAMPLE:
%   *WaitSeconds=10;# OR
%   *JuiceStatus=ON;#
%
%   The above line is a "Marker". The marker value is 10 in the first 
%   and it's ON in the second example.
%
%   The label, parameter name, and values are flexible and can be anything.
%   The only required formatting is that the user needs to have a label
%   followed by a colon ':', followed by a field name 'MarkerVal', followed
%   by an equal sign '=', followed by the parameter value '10', and end
%   with a semi-colon ';'. 
%
%   NOTE:
%   Every parameter requires a pound-sign '#' at the very end. 
%   Every parameter requires a star sign '*' at the very beginning. If you
%   use my LabVIEW SendtoCerebus VI then there is no need for a '*' in the
%   beginning.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 3.5.0 - March 22, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NEVver = '3.5.0'; % optimize for memory, supporting general usage
disp(['openNEV version ' NEVver])

    %%  Validating existance of parseCommand
if exist('parseCommand.m', 'file') ~= 2
    disp('This version of openNEV requires function parseCommand.m to be placed in path.');
    return;
end

%% Defining structures
NEV = struct('MetaTags',[],'IOLabels',[], 'ElectrodesInfo', [], 'Data',[]);
NEV.MetaTags = struct('Subject', [], 'Experimenter', [], 'DateTime',[],...
    'SampleRes',[],'Comment',[],'FileTypeID',[],'Flags',[]);
NEV.Data = struct('SerialDigitalIO', [], 'Spikes', []);
NEV.Data.Spikes = struct('Timestamps',[],'Electrode',[],...
    'Unit',[],'Waveform',[]);
NEV.Data.SerialDigitalIO = struct('InputType',[],'TimeStamp',[],...
    'TimeStampSec',[],'Type',[],'Value',[]);

%% Validating input arguments
for i=1:length(varargin)
    switch varargin{i}
        case 'report'
            Report = varargin{i};
        case 'read'
            ReadData = varargin{i};
        case 'nosave'
            SaveFile = varargin{i};
        case 'nomat'
            NoMAT = varargin{i};
        case 'nowarning'
            WarningStat = varargin{i};
        case 'noparse'
            ParseData = 'noparse';
        case 'compact'
            Compact = 'compact';
        otherwise
            temp = varargin{i};
            if length(temp)>3 && strcmpi(temp(end-3),'.')
                fname = varargin{i};
                if exist(fname, 'file') ~= 2
                    disp('The file does not exist.');
                    clear variables;
                    if nargout; NEV = []; end;
                    return;
                end
            else
                if ~isnumeric(varargin{i})
                    disp(['Invalid argument ''' varargin{i} ''' .']);
                else
                    disp(['Invalid argument ''' num2str(varargin{i}) ''' .']);
                end
                clear variables;
                return;
            end
    end
end

%% Defining and validating variables
if ~exist('fname', 'var')
    [fname, fpath] = uigetfile('D:\Data\*.nev');
else
    [fpath,fname, fext] = fileparts(fname);
    fname = [fname fext];
    fpath = [fpath '/'];
end
if ~exist('Report', 'var'); Report = 'noreport'; end
if ~exist('WarningStat', 'var'); WarningStat = 'warning'; end;
if ~exist('ReadData', 'var'); ReadData = 'noread'; end
if ~exist('ParseData', 'var'); ParseData = 'parse'; end
if ~exist('SaveFile', 'var'); SaveFile = 'save'; end;
if ~exist('NoMAT', 'var'); NoMAT = 'yesmat'; end;
if ~exist('Compact', 'var'); Compact = 'non-compact'; end;
if fname==0; clear variables; if nargout; NEV = []; end; disp('No file was selected.'); return; end;
tic;
matPath = [fpath fname(1:end-4) '.mat'];

%% Check for a MAT file and load that instead of NEV
if exist(matPath, 'file') == 2 && strcmpi(NoMAT, 'yesmat')
    disp('MAT file corresponding to selected NEV file already exists. Loading MAT instead...');
    load(matPath);
    if isempty(NEV.Data.Spikes.Waveform) && strcmpi(ReadData, 'read')
        disp('The MAT file does not contain waveforms. Loading NEV instead...');
    else
        if ~nargout
            assignin('base', 'NEV', NEV);
            clear variables;
        end
        return;
    end
end

%% Reading BasicHeader information from file
FID                       = fopen([fpath fname], 'r', 'ieee-le');
BasicHeader               = fread(FID, 336, '*uint8');
NEV.MetaTags.FileTypeID   = char(BasicHeader(1:8)');
fileSpec                  = num2str((BasicHeader(9:10)'));
NEV.MetaTags.Flags        = dec2bin(typecast(BasicHeader(11:12), 'uint16'),16);
fExtendedHeader           = double(typecast(BasicHeader(13:16), 'uint32'));
PacketBytes               = double(typecast(BasicHeader(17:20), 'uint32'));
TimeRes                   = double(typecast(BasicHeader(21:24), 'uint32'));
NEV.MetaTags.SampleRes    = typecast(BasicHeader(25:28), 'uint32');
t                         = typecast(BasicHeader(29:44), 'uint16');
NEV.MetaTags.Comment      = char(BasicHeader(77:332)');
NEV.MetaTags.Filename     = fname;
ExtHeaderCount            = typecast(BasicHeader(333:336), 'uint32');
if strcmpi(NEV.MetaTags.FileTypeID, 'NEURALEV') 
    disp(['Current NEV Filespec: ' fileSpec]);
    if exist([fpath fname(1:end-8) '.sif'], 'file') == 2
        METATAGS = textread([fpath fname(1:end-8) '.sif'], '%s');
        NEV.MetaTags.Subject      = METATAGS{3}(5:end-5);
        NEV.MetaTags.Experimenter = [METATAGS{5}(8:end-8) ' ' METATAGS{6}(7:end-7)];
    end
elseif strcmpi(fileSpec, '2  1')
    disp('Current NEV Filespec: 2.1');
else
    disp('Unknown filespec. Cannot open file.')
    clear variables;
    return;
end

%% Parsing and validating FileSpec and DateTime variables
NEV.MetaTags.DateTime = [num2str(t(2)) '/'  num2str(t(4)) '/' num2str(t(1))...
    ' ' datestr(double(t(3)), 'dddd') ' ' num2str(t(5)) ':'  ...
    num2str(t(6)) ':'  num2str(t(7)) '.' num2str(t(8))] ;

%% Removing extra garbage characters from the Comment field.
NEV.MetaTags.Comment(find(NEV.MetaTags.Comment==0,1):end) = 0;

%% Recording after BasicHeader file position
fBasicHeader = ftell(FID); %#ok<NASGU>

%% Reading ExtendedHeader information
for ii=1:ExtHeaderCount
    ExtendedHeader = fread(FID, 32, '*uint8');
    PacketID = char(ExtendedHeader(1:8)');
    switch PacketID
        case 'ARRAYNME'
            NEV.ArrayInfo.ElectrodeName    = char(ExtendedHeader(9:end));
        case 'ECOMMENT'
            NEV.ArrayInfo.ArrayComment     = char(ExtendedHeader(9:end));
        case 'CCOMMENT'
            NEV.ArrayInfo.ArrayCommentCont = char(ExtendedHeader(9:end));
        case 'MAPFILE'
            NEV.ArrayInfo.MapFile          = char(ExtendedHeader(9:end));
        case 'NEUEVWAV'
            PacketID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.ElectrodeID     = PacketID;
            NEV.ElectrodesInfo{PacketID, 1}.ConnectorBank   = char(ExtendedHeader(11)+64);
            NEV.ElectrodesInfo{PacketID, 1}.ConnectorPin    = ExtendedHeader(12);
            NEV.ElectrodesInfo{PacketID, 1}.DigitalFactor   = typecast(ExtendedHeader(13:14),'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.EnergyThreshold = typecast(ExtendedHeader(15:16),'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.HighThreshold   = typecast(ExtendedHeader(17:18),'int16');
            NEV.ElectrodesInfo{PacketID, 1}.LowThreshold    = typecast(ExtendedHeader(19:20),'int16');
            NEV.ElectrodesInfo{PacketID, 1}.Units           = ExtendedHeader(21);
            NEV.ElectrodesInfo{PacketID, 1}.WaveformBytes   = ExtendedHeader(22);
        case 'NEUEVLBL'
            PacketID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.ElectrodeLabel = char(ExtendedHeader(11:26));
        case 'NEUEVFLT'
            PacketID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.HighFreqCorner = typecast(ExtendedHeader(11:14),'uint32');
            NEV.ElectrodesInfo{PacketID, 1}.HighFreqOrder  = typecast(ExtendedHeader(15:18),'uint32');
            NEV.ElectrodesInfo{PacketID, 1}.HighFilterType = typecast(ExtendedHeader(19:20),'uint16');
            NEV.ElectrodesInfo{PacketID, 1}.LowFreqCorner  = typecast(ExtendedHeader(21:24),'uint32');
            NEV.ElectrodesInfo{PacketID, 1}.LowFreqOrder   = typecast(ExtendedHeader(25:28),'uint32');
            NEV.ElectrodesInfo{PacketID, 1}.LowFilterType  = typecast(ExtendedHeader(29:30),'uint16');
        case 'DIGLABEL'
            Label                                 = char(ExtendedHeader(9:24));
            Mode                                  = ExtendedHeader(25);
            NEV.IOLabels{Mode+1, 1} = Label;
        case 'NSASEXEV' %% Not implemented in the Cerebus firmware. 
                        %% Needs to be updated once implemented into the 
                        %% firmware by Blackrock Microsystems.
            NEV.NSAS.Freq          = typecast(ExtendedHeader(9:10),'uint16');
            NEV.NSAS.DigInputConf  = char(ExtendedHeader(11));
            NEV.NSAS.AnalCh1Conf   = char(ExtendedHeader(12));
            NEV.NSAS.AnalCh1Detect = typecast(ExtendedHeader(13:14),'uint16');
            NEV.NSAS.AnalCh2Conf   = char(ExtendedHeader(15));
            NEV.NSAS.AnalCh2Detect = typecast(ExtendedHeader(16:17),'uint16');
            NEV.NSAS.AnalCh3Conf   = char(ExtendedHeader(18));
            NEV.NSAS.AnalCh3Detect = typecast(ExtendedHeader(19:20),'uint16');
            NEV.NSAS.AnalCh4Conf   = char(ExtendedHeader(21));
            NEV.NSAS.AnalCh4Detect = typecast(ExtendedHeader(22:23),'uint16');
            NEV.NSAS.AnalCh5Conf   = char(ExtendedHeader(24));
            NEV.NSAS.AnalCh5Detect = typecast(ExtendedHeader(25:26),'uint16');
        otherwise
            display(['PacketID ' PacketID ' is unknown.']);
    end
end

%% Recording after ExtendedHeader file position and calculating Data Length
%  and number of data packets
fseek(FID, 0, 'eof');
fData = ftell(FID);
DataPacketCount = (fData - fExtendedHeader)/PacketBytes;
DataLen = PacketBytes - 8; %#ok<NASGU>
UnparsedDigitalDataFlag = 0;

%% Reading data packets if 'read' is passed as an argument
fseek(FID, fExtendedHeader, 'bof');
Timestamps        = fread(FID, DataPacketCount, '*uint32', PacketBytes-4);
fseek(FID, fExtendedHeader+4, 'bof');
PacketIDs         = fread(FID, DataPacketCount, '*uint16', PacketBytes-2);
fseek(FID, fExtendedHeader+6, 'bof');
tempClassOrReason = fread(FID, DataPacketCount, '*uchar', PacketBytes-1);
fseek(FID, fExtendedHeader+8, 'bof');
tempDigiVals      = fread(FID, DataPacketCount, '*uint16', PacketBytes-2);

%% Populate the NEV structure with spike timestamps, electrode numbers
% and unit numbers
nonNeuralIndices  = find(PacketIDs == 0);
neuralIndices     = find(PacketIDs ~= 0);
nonNeuTimestamps  = Timestamps(nonNeuralIndices);
NEV.Data.Spikes.Timestamps   = Timestamps(neuralIndices);
clear Timestamps; % save memory
NEV.Data.Spikes.Electrode = PacketIDs(neuralIndices);
NEV.Data.Spikes.Unit      = tempClassOrReason(neuralIndices);
clear neuralIndices; % save memory

%% Parse read digital data. Please refer to help to learn about the proper
% formatting if the data.
if strcmp(ParseData, 'parse')
    try
        InsertionReason   = tempClassOrReason(find(PacketIDs == 0));
        DigiValues        = char(tempDigiVals(nonNeuralIndices)');
    % %   This section needs to be uncommented out for Justin
    %     if int16(DigiValues(1)) > 128
    %         DigiValues = char(DigiValues-128);
    %     end
        AsteriskIndices   = find(DigiValues == '*');
        DataBegTimestamps = nonNeuTimestamps(AsteriskIndices);
        Inputs            = {'Digital'; 'AnCh1'; 'AnCh2'; 'AnCh3'; 'AnCh4'; 'AnCh5'; 'PerSamp'; 'Serial'};
        if ~isempty(DigiValues)
            if ~int8(DigiValues(2))
                DigiValues(find(DigiValues == DigiValues(2))) = [];
            end
            if strcmp(ParseData, 'parse') && ~isempty(DigiValues)
                splitDigiValues = regexp(DigiValues(2:end), '*', 'split')';
                for idx = 1:length(splitDigiValues)
                    try
                        if isempty(find(splitDigiValues{idx} == ':', 1))
                            splitDigiValues{idx}(find(splitDigiValues{idx} == '#')) = [];
                            NEV.Data.SerialDigitalIO(idx).Value = splitDigiValues{idx};
                            NEV.Data.SerialDigitalIO(idx).Type = 'Marker';
                        else
                            [tempParsedCommand error] = parseCommand(splitDigiValues{idx});
                            if ~error
                                pcFields = fields(tempParsedCommand);
                                for fidx = 1:length(pcFields)
                                    NEV.Data.SerialDigitalIO(idx).(pcFields{fidx}) = tempParsedCommand.(pcFields{fidx});
                                end
                            else
                               NEV.Data.SerialDigitalIO(idx).Value = splitDigiValues{idx};
                               NEV.Data.SerialDigitalIO(idx).Type = 'UnparsedData';
                               UnparsedDigitalDataFlag = 1;
                            end
                        end
                    catch
                        disp(['Error parsing: ' splitDigiValues{idx}]);
                        disp('Please refer to the help for more information on how to properly format the digital data for parsing.');
                    end
                end            
                % Populate the NEV structure with timestamps and inputtypes for the
                % digital data
                if ~isempty(DataBegTimestamps)
                    c = num2cell(DataBegTimestamps); [NEV.Data.SerialDigitalIO(1:length(NEV.Data.SerialDigitalIO)).TimeStamp] = deal(c{1:end});
                    c = num2cell(DataBegTimestamps/NEV.MetaTags.SampleRes); [NEV.Data.SerialDigitalIO.TimeStampSec] = deal(c{1:end});
                    c = {Inputs{InsertionReason(AsteriskIndices)}}; [NEV.Data.SerialDigitalIO.InputType] = deal(c{1:end});
                end
            elseif ~isempty(DigiValues)
                NEV.Data.SerialDigitalIO.TimeStamp = nonNeuTimestamps;
                NEV.Data.SerialDigitalIO.UnparsedData = DigiValues;
            end
        else
            disp('No digital data to read.');
        end
    catch
        disp('An error occured during reading digital data. This is due to a problem with formatting digital data.');
        disp('Refer to help ''help openNEV'' for more information on how to properly format the digital data.');
    end
else
    clear tempClassOrReason;
    clear PacketIDs;
    NEV.Data.SerialDigitalIO.UnparsedData = tempDigiVals(nonNeuralIndices);
    clear tempDigiVals;
    NEV.Data.SerialDigitalIO.TimeStamp    = nonNeuTimestamps;
    NEV.Data.SerialDigitalIO.TimeStampSec = double(nonNeuTimestamps) / double(NEV.MetaTags.SampleRes);
end

%% Reads the waveforms if 'read' is passed to the function
if strcmp(ReadData, 'read')
    fseek(FID, fExtendedHeader + 8, 'bof'); % go to where spikes are located
    if (strcmp(Compact, 'compact'))
        SpikeWaveformType='int16';
    else
        SpikeWaveformType='double';
    end
    NEV.Data.Spikes.Waveform  = fread(FID, [(PacketBytes-8)/2 DataPacketCount], ...
                                [num2str((PacketBytes-8)/2) '*int16=>' SpikeWaveformType], 8)';
    NEV.Data.Spikes.Waveform(nonNeuralIndices, :) = [];
    clear nonNeuTimestamps;
end

%% Calculating the length of the data
fseek(FID, -PacketBytes, 'eof');
DataDuration = fread(FID, 1, 'uint32=>double');
%% Show a report if 'report' is passed as an argument
if strcmp(Report, 'report')
    disp( '*** FILE INFO **************************');
    disp(['File Path           = ' fpath(1:end-1)]);
    disp(['File Name           = ' fname]);
    disp(['Filespec            = ' fileSpec]);
    disp(['Data Duration (min) = ' num2str(round(DataDuration/NEV.MetaTags.SampleRes/60))]);
    disp(['Packet Counts       = ' num2str(DataPacketCount)]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');    
    disp(['File Type ID        = '         NEV.MetaTags.FileTypeID]);
    disp(['Sample Resolution   = ' num2str(NEV.MetaTags.SampleRes)]);
    disp(['Date and Time       = '         NEV.MetaTags.DateTime]);
    disp(['Comment             = '         NEV.MetaTags.Comment(1:64)   ]);
    disp(['                      '         NEV.MetaTags.Comment(65:128) ]);
    disp(['                      '         NEV.MetaTags.Comment(129:192)]);
    disp(['                      '         NEV.MetaTags.Comment(193:256)]);
end


%% Display how fast the function was executed.
if strcmp(Report, 'report')
    disp(['The load time was for NEV file was ' num2str(toc, '%0.1f') ' seconds.']);
end

%% Closing and clearing memory
if strcmp(ParseData, 'parse')
    if UnparsedDigitalDataFlag && strcmp(WarningStat, 'warning')
        fprintf(2, 'WARNING: The NEV file contains unparsed digital data.\n'); % where the file is opened?
        pause;
    end
end

%% Saving the NEV structure as a MAT file for easy access
if strcmp(SaveFile, 'save')
    if exist(matPath, 'file') == 2
        disp(['File ' matPath ' already exists.']);
        overWrite = input('Would you like to overwrite (Y/N)? ', 's');
        if strcmpi(overWrite, 'y')
            disp('Saving MAT file. This may take a few seconds...');
            save(matPath, 'NEV');
        else
            disp('File was not overwritten.');
        end
    else
        disp('Saving MAT file. This may take a few seconds...');
        save(matPath, 'NEV');
    end
end

%% Closing and clearing memory
if ~nargout
    assignin('base', 'NEV', NEV);
    fclose(FID);
    clear variables;
else
    fclose(FID);
end
end
