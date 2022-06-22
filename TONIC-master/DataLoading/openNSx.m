function NSx = openNSx(varargin)

%%Opens an NSx file for reading, returns all file information in a NSx
% structure. Works with File Spec 2.1 and 2.2. 
% Use OUTPUT = openNSx(fname, 'read', 'report', 'electrodes', 'duration', 'mode', 'precision').
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   fname:        Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file. 
%                 DEFAULT: Will open Open File UI.
%
%   'read':       Will read the data in addition to the header information
%                 if user passes this argument.
%                 DEFAULT: will only read the header information.
%
%   'report':     Will show a summary report if user passes this argument.
%                 DEFAULT: will not show report.
%
%   'electrodes': User can specify which electrodes need to be read. The
%                 number of electrodes can be greater than or equal to 1
%                 and less than or equal to 128. The electrodes can be
%                 selected either by specifying a range (e.g. 20:45) or by
%                 indicating individual electrodes (e.g. 3,6,7,90) or both.
%                 This field needs to be followed by the prefix 'e:'. See
%                 example for more details.
%                 DEFAULT: will read all existing channels.
%
%   'duration':   User can specify the beginning and end of the data
%                 segment to be read. If the start time is greater than the
%                 length of data the program will exit with an error
%                 message. If the end time is greater than the length of
%                 data the end packet will be selected for end of data. The
%                 user can specify the start and end values by comma 
%                 (e.g. [20,50]) or by a colon (e.g. [20:50]). To use this
%                 argument the user must specify the [electrodes] or the
%                 interval will be used for [electrodes] automatically.
%                 This field needs to be followed by the prefix 't:'. See
%                 example for more details.
%                 DEFAULT: will read the entire file.
%
%   'mode':       The user can specify the mode of duration in [duration],
%                 such as 'sec', 'min', 'hour', or 'sample'. If 'sec' is
%                 specified the numbers in [duration] will correspond to
%                 the number of seconds. The same is true for 'min', 'hour'
%                 and 'sample'.
%                 DEFAULT: will be set to 'sample'.
%
%   'precision':  The data storage class can be any format known by
%                 MATLAB such as 'double', 'int16', 'int'.
%                 This field needs to be followed by the prefix 'p:'. See
%                 example for more details.
%                 DEFAULT: will be set to 'double'
%   
%   OUTPUT:       Contains the NSx structure.
%
%   Example: 
%   
%   openNSx('report','read','c:\data\sample.ns5', 'e:15:30', 't:3:10', 'min', 'p:int16');
%
%   In the example above, the file c:\data\sample.ns5 will be used. A
%   report of the file contents will be shown. The data will be read from
%   electrodes 15 through 50 in the 3-10 minute time interval. 
%   The data is saved in 'int16' type.
%   If any of the arguments above are omitted the default values will be used.
%
%   Kian Torab
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 2.1.0 - March 22, 2010

NSxver = '2.1.0';
disp(['openNSx version ' NSxver])

%% Defining the NSx data structure and sub-branches.
NSx             = struct('MetaTags',[],'Data',[]);
NSx.MetaTags    = struct('FileTypeID',[],'ChannelCount',[],'SamplingFreq',[], ...
                         'ChannelID',[],'Version',[],'ElecLabel',[],'CreateDateTime',[], ...
                         'Resolution',[], 'Comments',[], 'VoltsPerDigit', []);

%% Validating the input arguments. Exit with error message if error occurs.
for i=1:length(varargin)
    inputArgument = varargin{i};
    if strcmpi(inputArgument, 'report')
        Report = inputArgument;
    elseif strcmpi(inputArgument, 'read')
        ReadData = inputArgument;
    elseif strncmp(varargin{i}, 't:', 2)
        colonIndex = find(inputArgument(3:end) == ':');
        StartPacket = str2num(inputArgument(3:colonIndex+1));
        EndPacket = str2num(inputArgument(colonIndex+3:end));    
        if min(varargin{i})<1 || max(varargin{i})>128
            display('The electrode number cannot be less than 1 or greater than 128.');
            clear variables;
            if nargout; NSx = []; end
            return;
        end
    elseif strncmp(varargin{i}, 'e:', 2)
        Elec = str2num(inputArgument(3:end)); %#ok<ST2NM>        
    elseif strncmp(varargin{i}, 'p:', 2)
        Precision = inputArgument(3:end); % precision for storage
    elseif strfind(' hour min sec sample ', [' ' inputArgument ' ']) ~= 0
        TimeScale = inputArgument;
    else
        temp = inputArgument;
        if length(temp)>3 && strcmpi(temp(end-3),'.')
            fname = inputArgument;
            if exist(fname, 'file') ~= 2
                display('The file does not exist.');
                clear variables;
                if nargout; NSx = []; end
                return;
            end
        else
            display(['Invalid argument ''' inputArgument ''' .']);
            clear variables;
            if nargout; NSx = []; end
            return;
        end
    end
end

%% Popup the Open File UI. Also, process the file name, path, and extension
%  for later use, and validate the entry.
if ~exist('fname', 'var')
    [fname, path] = uigetfile('D:\Data\*.ns*');
    if fname == 0
        clear variables;
        if nargout; NSx = []; end
        return;
    end
    fext = fname(end-3:end);
else
    [path, fname, fext] = fileparts(fname);
    fname = [fname fext];
end
if fname==0; return; end;

tic;

%% Give all input arguments a default value. All input arguments are
%  optional.
if ~exist('Report', 'var');      Report = 'noreport'; end
if ~exist('ReadData', 'var');    ReadData = 'noread'; end
if ~exist('StartPacket', 'var'); StartPacket = 0;     end
if ~exist('TimeScale', 'var');   TimeScale = 1;       end
if ~exist('Precision', 'var');   Precision = 'double';       end

%% Reading Basic Header from file into NSx structure.
% we use fullfile instead of [path '\' fname] to support nix platforms
% 
% General format of calling fopen: 
%     fileID = fopen(filename [,permission [,machinefmt,encodingIn]])
% where
%     permission = 'r' (default) | 'w' | 'a' | 'r+' | 'w+' | 'a+' | 'A' | 'W' | ...
%     machinefmt = Order for reading or writing bytes or bits'n' (default) | 'b' | 'l' | 's' | 'a' | ...
%                  'n' or 'native'           Your system byte ordering (default)
%                  'b' or 'ieee-be'          Big-endian ordering
%                  'l' or 'ieee-le'          Little-endian ordering
%                  's' or 'ieee-be.l64'      Big-endian ordering, 64-bit long data type
%                  'a' or 'ieee-le.l64'      Little-endian ordering, 64-bit long data type
%                  (for bionary files, default = 'b' or 'l')
%     encodingIn = character encoding
%                  'Big5'         'ISO-8859-1'   'windows-932'
%                  'GB2312'       'ISO-8859-2'   'windows-936'
%                  'EUC-KR'       'ISO-8859-3'   'windows-949'
%                  'EUC-JP'       'ISO-8859-4'   'windows-950'
%                  'GBK'          'ISO-8859-9'   'windows-1250'
%                  'KSC_5601'     'ISO-8859-13'  'windows-1251'
%                  'Macintosh'    'ISO-8859-15'  'windows-1252'
%                  'Shift_JIS'    'windows-1253' 'US-ASCII'
%                  'windows-1254' 'UTF-8'        'windows-1257'
%
machinefmt = {'native', 'ieee-be', 'ieee-le', 'ieee-be.l64', 'ieee-le.l64'};
encodingIn = {'Macintosh', ...
              'Big5'     ,   'ISO-8859-1'  ,'windows-932',...
              'GB2312'   ,   'ISO-8859-2'  ,'windows-936',...
              'EUC-KR'   ,   'ISO-8859-3'  ,'windows-949',...
              'EUC-JP'   ,   'ISO-8859-4'  ,'windows-950',...
              'GBK'      ,   'ISO-8859-9'  ,'windows-1250',...
              'KSC_5601' ,   'ISO-8859-13' ,'windows-1251',...
                             'ISO-8859-15' ,'windows-1252',...
              'Shift_JIS',   'windows-1253','US-ASCII',...
              'windows-1254','UTF-8'       ,'windows-1257' };
for j=1:numel(encodingIn)
    for i = 1:numel(machinefmt)
        FID = fopen(fullfile(path,fname), 'r', machinefmt{i}, encodingIn{j}); 
        if int16(FID) > 0
            disp(['encodingIn=' encodingIn{j} ' machinefmt='  machinefmt{i} ]);
            my_machinefmt = machinefmt{i};
            my_encodingIn = encodingIn{j};
            break;
        end
    end
    if int16(FID) > 0
        break;
    end
end
% General format of calling fread:
%     FileTypeID = fread(fileID [,sizeA [,precision [,skip [,machinefmt]]]])
% where
%     sizeA      = dimensions of output array: Inf (default) | integer two-element vector
%     precision  = class and size of values to read: 'uint8=>double' (default) | string
%     skip       = number of bites to skip: default = 0
%     machinefmt
NSx.MetaTags.FileTypeID = fread(FID, [1,8], '*char');

if strcmp(NSx.MetaTags.FileTypeID, 'NEURALSG')                              % 2.1
    NSx.MetaTags.Version = '2.1';
    NSx.MetaTags.Label        = fread(FID, [1,16]  , '*char'          );
    NSx.MetaTags.SamplingFreq = 30000/fread(FID, 1 , 'uint32=>double' );
    ChannelCount                 = fread(FID, 1       , 'uint32=>double' );
    NSx.MetaTags.ChannelID    = fread(FID, [ChannelCount 1], '*uint32');
    NSx.MetaTags.ChannelCount = ChannelCount;
    fHeader = ftell(FID);
elseif strcmp(NSx.MetaTags.FileTypeID, 'NEURALCD')                          % 2.2
    Major        = num2str(fread(FID, 1  , 'uint8=>double'   ));
    Minor        = num2str(fread(FID, 1  , 'uint8=>double'   ));
    NSx.MetaTags.Version = [Major '.' Minor];
    fHeader      = fread(FID, 1  , 'uint32=>double'          );
    NSx.MetaTags.Label        = fread(FID, [1,16]  , '*char'          );
    NSx.MetaTags.Comments        = fread(FID, [1,256]  , '*char'          );
    NSx.MetaTags.SamplingFreq = 30000/fread(FID, 1 , 'uint32=>double' );
    NSx.MetaTags.Resolution = fread(FID, 1 , 'uint32=>double' );
    NSx.MetaTags.CreateDateTime = fread(FID, [1,8] , 'uint16=>double' );
    ChannelCount                 = fread(FID, 1       , 'uint32=>double' );
    NSx.MetaTags.ChannelCount = ChannelCount;
    NSx.MetaTags.ChannelID = zeros(ChannelCount, 1);
    NSx.MetaTags.ElecLabel = char(zeros(ChannelCount, 16));         % Electrode label
    NSx.MetaTags.VoltsPerDigit = ones(ChannelCount, 1); % volet per digit
    % now read external header
    for ii = 1:ChannelCount
        CC = fread(FID, [1,2], '*char' );
%       disp(['CC=' CC ]);
        if ~strcmp(CC, 'CC') 
            disp('Wrong extension header');
            fclose(FID);
            clear variables;
            return; 
        end
        NSx.MetaTags.ChannelID(ii) = fread(FID, 1 , 'uint16=>double' );
        NSx.MetaTags.ElecLabel(ii,:) = fread(FID, [1,16]  , '*char' );
        dummy = fread(FID, 4, 'uint8=>double'   ); % dummy
        maxDig = fread(FID, 1, 'int16=>double'   ); % minimum digital value
        dummy = fread(FID, 2, 'uint8=>double'   ); % dummy
        maxAna = fread(FID, 1, 'int16=>double'   ); % minimum analog value
        unit = fread(FID, [1,16]  , '*char' );
        if (strcmpi(unit(1:2), 'uV')) % micro volts
            NSx.MetaTags.VoltsPerDigit(ii) = maxAna / maxDig * 1e-6;
        elseif (strcmpi(unit(1:2), 'mV')) % milli volts
            NSx.MetaTags.VoltsPerDigit(ii) = maxAna / maxDig * 1e-3;
        else
            disp('Wrong unit');
            fclose(FID);
            clear variables;
            return; 
        end
        % We do not care about the rest now
        dummy = fread(FID, [20,1], 'uint8=>double'   ); % dummy
    end
    clear dummy;
    clear unit;
    clear maxDig;
    clear maxAna;
    if fHeader ~= ftell(FID)
        display('Header file corrupted!');
        fHeader = ftell(FID);
    end
    fHeader = fHeader + 9; % to account for the data header
else
    display('This version of openNSx can only read File Specs 2.1 or 2.2');
    display(['The selected file label is ' NSx.MetaTags.FileTypeID '.']);
    fclose(FID);
    clear variables;
    if nargout; NSx = []; end;
    return; 
end;
% find out number of data points
fseek(FID, 0, 'eof');
fData = ftell(FID);
fseek(FID, fHeader, 'bof');

%% Adjusts StartPacket and EndPacket based on what time setting (sec, min,
%  hour, or packets) the user has indicated in the input argument.
switch TimeScale
    case 'sec'
        StartPacket = round(StartPacket * NSx.MetaTags.SamplingFreq);
        EndPacket = round(EndPacket * NSx.MetaTags.SamplingFreq);
    case 'min'
        StartPacket = StartPacket * NSx.MetaTags.SamplingFreq * 60;
        EndPacket = EndPacket * NSx.MetaTags.SamplingFreq * 60;
    case 'hour'
        StartPacket = StartPacket * NSx.MetaTags.SamplingFreq * 3600;
        EndPacket = EndPacket * NSx.MetaTags.SamplingFreq * 3600;
end

%% Validate StartPacket and EndPacket to make sure they do not exceed the
%  length of packets in the file. If EndPacket is over then the last packet
%  will be set for EndPacket. If StartPacket is over then will exit with an
%  error message.
NumofPackets = (fData-fHeader)/(2*ChannelCount);
if exist('EndPacket', 'var') && (EndPacket >= NumofPackets)
    disp('The time interval specified is longer than the data duration.');
    if StartPacket >= NumofPackets
        disp('The starting packet is greater than the total data duration.');
        clear variables;
        if nargout; NSx = []; end
        return;
    end
    disp('The time interval specified is longer than the data duration.');
    disp('Last data point will be used instead.');
%     disp('Press enter to continue...');
%     pause;
    EndPacket = NumofPackets - 1;
elseif ~exist('EndPacket', 'var')
    EndPacket = NumofPackets - 1;
    disp('The time interval specified is longer than the data duration. Using last packet');
end
DataLength = EndPacket - StartPacket + 1;
clear TimeScale

%% Displaying a report of basic file information and the Basic Header.
if strcmp(Report, 'report')
    disp( '*** FILE INFO **************************');
    disp(['File Path          = '  path]);
    disp(['File Name          = '  fname   ]);
    disp(['File Version       = '  NSx.MetaTags.Version   ]);
    disp(['Duration (seconds) = '  num2str(NumofPackets/NSx.MetaTags.SamplingFreq)]);
    disp(['Total Data Points  = '  num2str(NumofPackets)                   ]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');
    disp(['File Type ID       = '          NSx.MetaTags.FileTypeID      ]);
    disp(['Label              = '          NSx.MetaTags.Label           ]);
    disp(['Sample Resolution  = '  num2str(NSx.MetaTags.SamplingFreq)         ]);
    disp(['Electrodes Read    = '  num2str(NSx.MetaTags.ChannelCount)   ]);
end
NSx.MetaTags.Duration = double(round(NumofPackets/NSx.MetaTags.SamplingFreq));

%%
if ~exist('Elec', 'var'); 
    Elec = 1:ChannelCount; 
end
ReadElec = max(Elec)-min(Elec)+1;
if (ReadElec <= ChannelCount)
    if strcmp(ReadData, 'read')
        fseek(FID, StartPacket * 2 * ChannelCount + fHeader, 'bof');
        fseek(FID, (min(Elec)-1) * 2, 'cof');
        NSx.Data = fread(FID, [ReadElec DataLength-1], [num2str(ReadElec) '*int16=>' Precision], (ChannelCount-ReadElec) * 2);    
    end
end
%% If user does not specify an output argument it will automatically create a structure.
outputName = ['NS' fext(4)];
if (nargout == 0),
    assignin('caller', outputName, NSx);
end

if strcmp(Report, 'report')
    disp(['The load time for ' outputName ' file was ' num2str(toc, '%0.1f') ' seconds.']);
end
fclose(FID);

end
