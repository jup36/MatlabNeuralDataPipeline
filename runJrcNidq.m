function runJrcNidq(filePath,binFileList,probeFileName)
%RUNJRC Executes JRCLUST program
%   RUNJRC(BINFILELIST) runs JRCLUST program. Makes parameter file
%   (*.prm) using 'jrc makeprm'. Then, it executes spike sorting using 'jrc
%   spikesort'. Finally, it starts manual clustering using  'jrc manual'.

%   input variable 'binFileList' should be cell list that includes exact
%   file names.

%   Dohoung Kim
%   Howard Hughes Medical Institute
%   Janelia Research Campus
%   19700 Helix Drive
%   Ashburn, Virginia 20147
%   kimd11@janelia.hhmi.org

% default data directory
cd(filePath)
DATA_PATH = 'E:\SGL_SSD_Data';
PROBE_PATH = 'Z:\parkj\Probe'; 

%% 1-1. Check whether there's a bin file in the filePath. If not, pop up the window to select manually
if nargin < 1 || isempty(binFileList) || ~iscell(binFileList)
    binList = dir(fullfile(filePath,'*.nidq.bin'));
    
    if isempty(binList)
        dataPath = uigetdir(DATA_PATH);
        if ~ischar(dataPath); return; end
        binList = dir(fullfile(dataPath,'*.nidq.bin'));
    else
        dataPath = filePath;
    end

    nBin = length(binList);
    binFile = {};
    for iBin = 1:nBin
        if binList(iBin).bytes > 10^10
            binFile = [binFile; {fullfile(dataPath, binList(iBin).name)}];
        end
    end
else
    nBin = length(binFileList);
    binFile = {};
    for iBin = 1:nBin
        if exist(binFileList{iBin}, 'file')
            binFile = [binFile; binFileList{iBin}];
        end
    end
end

%% 2. Run JRC to make prm file and do spike sorting
nBin = length(binFile);
[prmFile, spkwavFile] = deal(cell(nBin, 1));
for iBin = 1:nBin
    binFileName = binFile{iBin};
    %option = readOption(binFileName);
    prmFile{iBin} = replace(binFileName, '.bin', ['_', replace(probeFileName, '.prb', '') ,'.prm']);
    spkwavFile{iBin} = replace(prmFile{iBin}, '.prm', '_spkwav.jrc');

    % make prm
    if exist(prmFile{iBin}, 'file') ~= 2
        if exist(probeFileName, 'file') ~= 2
            error('Probe file not found!! Probe is required to be in the current folder to make the .prm file!!')
        end
        
        jrc('makeprm', binFileName, probeFileName);
    end
    
    if exist(spkwavFile{iBin}, 'file') ~= 2
        jrc('spikesort', prmFile{iBin});
    end
end

%% 3. After automated spike sorting, do manual spike sorting
iFile = listdlg('PromptString', 'Select a file for manual clustering', ...
    'SelectionMode', 'single', ...
    'ListSize', [400, 200], ...
    'ListString', prmFile);

if ~isempty(iFile)
    jrc('manual', prmFile{iFile});
end

function option = readOption(binFileName)
% Parse ini file into cell entries C{1}{i} = C{2}{i}
metaFile = replace(binFileName, '.bin', '.meta');
fid = fopen(metaFile, 'r');
C = textscan(fid, '%[^=] = %[^\r\n]');
fclose(fid);

% New empty struct
meta = struct();

% Convert each cell entry into a struct entry
for i = 1:length(C{1})
    tag = C{1}{i};
    if tag(1) == '~'
        % remake tag excluding first character
        tag = sprintf('%s', tag(2:end));
    end
    meta.(tag) = C{2}{i};
end
option = meta.imProbeOpt(1);


end

end
