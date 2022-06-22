function [outStruct] = TNC_PackBlackrockSession(fileNameBase,recParams)
%% FUNCTION DETAILS: Load and pack a given session of Blackrock Data into a single standard structure organization
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: dudmanlab.org/projects.html
% _________________________________________________________________________
% SOME EXMAPLE DATA...
% 
% F01_ACQ_004_NEVDATA.Data
% 
%     SerialDigitalIO: [1x1 struct]
%              Spikes: [1x1 struct]
% 
% F01_ACQ_004_NEVDATA.Data.Spikes
% 
%     Timestamps: [4283704x1 uint32]
%      Electrode: [4283704x1 uint16]
%           Unit: [4283704x1 uint8]
%       Waveform: [4283704x30 double]
%     
% F01_ACQ_004_NS4DATA = 
% 
%     MetaTags: [1x1 struct]
%         Data: [37x30272809 double]
%     fileName: 'F01-100923-trace-P70-004.ns4'
% 

%% PARAMETERS
% recParams is a structure containing some essential information about the recording 
% including the number of electrode arrays, event channel mappings, etc. -- detailed below

% this will be rig specific. Settings for Weixing's BlackRock system
recParams.rewChan;
recParams.lickChan = 138;
recParams.eventChan = 137;
recParams.arrays;

session.fileNameBase    = fileNameBase;
session.description     = []; % would like someway to input metadata here... perhaps load from a file? xml? hmm..
session.cranitomy.date  = [year,month,day];
session.cranitomy.side  = 'left';
session.cranitomy.num   = 1;
session.recording.logic = 1;
session.recording.numArray  = numel(recParams.arrays);

for j = 1:session.recording.numArray
    session.array(j).elecType   = [];
end

paramsC.tapers   = [3 5];
paramsC.pad      = 0;
paramsC.Fs       = lowBandData.sampleRate;
paramsC.fpass    = [0 100];
paramsC.err      = 0;
session.paramsC = paramsC;

movingwin        = [0.5 0.05];
session.movingwin = movingwin;

%% NS4
aC = 1; eC = 1;
fileNameStr         = [fileNameBase '.ns4'];
% [Ns4DATA]           = TNC_LoadData(0, 1, fileNameStr);

disp(' ');
disp('_____________ ns4 loading _____________');
disp(' ');
disp('Check header...');
Ns4DATA             = openNSx('read',fileNameStr, 't:1:200');
session.creation    = Ns4DATA.MetaTags.CreateDateTime;

% first electrode

% how many arrays
% arrays = session.recording.numArray;

% Get the number of the lowest recorded electrode
if strcmp(Ns4DATA.MetaTags.ElecLabel(1,6),' ')
    howLong = 5;
else
    howLong = 6;
end
firstElecName   = Ns4DATA.MetaTags.ElecLabel(1,5:howLong);
firstElecNum    = str2num(firstElecName);

disp('Begin repacking...');

for i=1:size(Ns4DATA.Data,1)

    if strcmp(Ns4DATA.MetaTags.ElecLabel(i,6),' ')
        howLong = 5;
    else
        howLong = 6;
    end
    brNum   = str2num(Ns4DATA.MetaTags.ElecLabel(i,5:howLong));

    electrode = ['e:' num2str(i)];

    disp(['Loading ... from ' Ns4DATA.MetaTags.ElecLabel(i,:)]);    
%     dataSpks = openNSx('read',fileNameStr, electrode); 
    dataSpks = openNSx('read',fileNameStr, electrode,'t:1:18000000'); %load first 10 minutes for debuggingx
    disp('...loaded');    

    if strcmp(Ns4DATA.MetaTags.ElecLabel(i,1:4),'elec')

        session.elec(eC).elecName               = Ns4DATA.MetaTags.ElecLabel(i,1:howLong);
        
        session.elec(eC).brNum                  = brNum;
        session.elec(eC).arNum                  = brNum-firstElecNum-1;
%         session.elec(eC).contData.tenk(1,:)     = dataSpks.Data(1,:);
        disp('...decimating');
        session.elec(eC).contData.onek(1,:)     = decimate(dataSpks.Data(1,:),10);

        % Find position on the electrode array
        session.elec(eC).array                  = 1;
        eC = eC+1;

    else
        
        session.ainp(aC).ainpName               = Ns4DATA.MetaTags.ElecLabel(i,1:howLong);
        session.ainp(aC).contData.tenk(1,:)     = dataSpks.Data(1,:);
        aC = aC+1;
        
    end

    disp(['...packed ' Ns4DATA.MetaTags.ElecLabel(i,:)]);
    disp(' ');
    
end

disp(' ');
disp('_____________ ns4 completed _____________');
disp(' ');

%% NS5
fileNameStr = [fileNameBase '.ns5'];
% test if the file exists
NS5true = 0 % default setting

if NS5true
    [Ns5DATA] = TNC_LoadData(0, 1, fileNameStr);

    % how many arrays
    arrays = session.recording.numArray;

    % Get the number of the lowest recorded electrode
    howLong         = numel(Ns5DATA.MetaTags.ElecLabel(1,:));
    firstElecName   = Ns5DATA.MetaTags.ElecLabel(1,5:howLong);
    firstElecNum    = str2num(firstElecName);

    for i=1:numChan

        disp(['Loading data from electrode ' num2str(i) '...']);
        electrode = ['e:' num2str(i)];
        dataSpks = openNSx('report','read',fileNameStr, electrode);

        % Filter the data to select for spikes
        disp('Filtering data...');
        [lowBandData,hiBandData] = TNC_FilterData(dataSpks.Data(1,:),dataSpks.MetaTags.SamplingFreq,1,0);
        data = rmlinesmovingwinc(lowBandData.values,movingwin,10,paramsC,0.05,'n',60);
        lowBandData.values=data';

    %     % Find threshold crossings with a specific SNR
    %     disp('Finding threshold crossings...');
    %     [events] = TNC_EventDetect(hiBandData.values,hiBandData.sampleRate,snrThresh);
    % 
    %     % Store the delta function sampled at 1kHz and the lfp amplitude
    %     tmpTs = events.inds .* events.sampleRate;
    %     delta = zeros(1,round(size(dataSpks.Data(1,:),2).*(1000./events.sampleRate)));
    %     delta(1,round(events.inds.*(1000./events.sampleRate))) = 1;
    %     outputStructure.channel(i).delta = delta;

        % Apply chronux multitaper method to extract power spectrum in time
        disp('Calculating the spectrogram of the lowBand data...');
        movingwin       = [1 0.01];
        outputStructure.channel(i).lfp.params      = paramsC;
        outputStructure.channel(i).lfp.movingwin   = movingwin;
        [S,t,f] = mtspecgramc(lowBandData.values,movingwin,paramsC);

        session.array(j).elec(i).contData.t = t;
        session.array(j).elec(i).contData.f = f;
        session.array(j).elec(i).contData.S = S;
    end
end

%% NEV
disp(' ');
disp('_____________ nev loading _____________');
disp(' ');
disp('Loading...');

fileNameStr = [fileNameBase '.nev'];
[NevDATA] = TNC_LoadData(0, 1, fileNameStr);

nevRes = NevDATA.MetaTags.SampleRes;

% pass in the first electrode name
firstTrodeNum = str2num(firstElecName);

% walk through the session.elec(i).elecName space and then grab all ts for individual units
% totalTS = numel(NevDATA.Data.Spikes.Electrode);
totalTS = size(find(NevDATA.Data.Spikes.Timestamps<3*18000000),1); % for loader debugging

session.event(1).true = 1;
session.event(1).unit(1).true = 1; 

disp('Sorting events...');


% walk through all spikes and sort into containers
for j=1:totalTS
    
    if rem(j,round(totalTS./10))==0
        disp(['Fraction complete: ' num2str(j./totalTS) ' ']);
    end
    
    currStamp   = NevDATA.Data.Spikes.Timestamps(j);
    eN          = NevDATA.Data.Spikes.Electrode(j);
    uN          = NevDATA.Data.Spikes.Unit(j) + 1;
    
    if eN < 129
        session.elec(eN - firstTrodeNum + 1).true = 1;
        session.elec(eN - firstTrodeNum + 1).unit(uN).true = 1; 
        if isfield(session.elec(eN - firstTrodeNum + 1).unit(uN),'ts')
            session.elec(eN - firstTrodeNum + 1).unit(uN).ts = [session.elec(eN - firstTrodeNum + 1).unit(uN).ts,currStamp];
        else
            session.elec(eN - firstTrodeNum + 1).unit(uN).ts = currStamp;
        end
    else
        session.elec(eN).true = 1;
        session.elec(eN).unit(uN).true = 1; 
        if isfield(session.elec(eN).unit(uN),'ts')
            session.elec(eN).unit(uN).ts = [session.elec(eN).unit(uN).ts,currStamp];
        else
            session.elec(eN).unit(uN).ts = currStamp;
        end
    end
end

disp(' ');
disp('_____________ nev completed _____________');
disp(' ');

%% End structure

% session
% 
% session = 
% 
%     creation: [2011 9 5 23 22 50 51 895]
%         elec: [1x138 struct]
%         ainp: [1x7 struct]
%        event: [1x1 struct]
% 
% session.elec
% 
% ans = 
% 
% 1x138 struct array with fields:
%     elecName
%     brNum
%     arNum
%     contData
%     array
%     true
%     unit
% 
% session.elec(3).unit
% 
% ans = 
% 
% 1x2 struct array with fields:
%     true
%     ts
% 
% session.elec(3).unit(2)
% 
% ans = 
% 
%     true: 1
%       ts: [1x7297 uint32]
% 
