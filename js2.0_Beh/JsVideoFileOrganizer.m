function [S] = jsVideoFileOrganizer(filePath)
%This function inspects the trial-by-trial front and side videos, and
% assign them to corresponding trials. The output is a structure named 'jsTime1k_KV'.  
% This function has been updated in september/2019 to add a field named
% 'vUseFrameIdx', which is a logical important to determine which frames
% from the video comprise the hand trajectory of the current trial. This is
% especially critical when using a video file for multiple trials (e.g. S:\Junchol_Data\JS2.0\WR40_081419\jsTime1k_Kinematics_VideoFiles.mat\jsTime1k_KV(17).vUseFrameIdx). 

cd(filePath)
if exist('jsTime1k_K','var')==1
    S=jsTime1k_K;
else
    pathJsTime1k_K = dir('**/*_kinematics.mat');
    S=load(fullfile(pathJsTime1k_K.folder,pathJsTime1k_K.name),'jsTime1k_K'); % just building up on the outcome of jsKinematicsAnalysis.m to create one all-inclusive file
    S = S.('jsTime1k_K');
end
S = rmfield(S,{'baseJsTrajmm','baseSmJsVel','basePeriodicAbsVelSum'}); % rmfield remove fields from a structure array

pathBehVar = dir('**/*BehVariablesJs.mat');
load(fullfile(pathBehVar.folder,pathBehVar.name), 'evtIdx1k', 'p')
trStartIdx = evtIdx1k.trStartIdx; 
%load(fullfile(filePath,'evtIndices.mat'), 'trStartIdx', 'trEndIdx')

behFilePath = dir(fullfile(filePath,'201*')); % dir where the trial-by-trial behavioral csv files are saved
tbytCsvList = dir(fullfile(behFilePath.folder,behFilePath.name,'trial_*'));    % trial-by-trial files
allTrialCsv = dir(fullfile(behFilePath.folder,behFilePath.name,'trials.csv')); % all trial file
if length(allTrialCsv)==1
    trialsFileName = fullfile(allTrialCsv.folder,allTrialCsv.name);
    trialsCsv = readtable(trialsFileName);
else
    error('More than one trials.csv file detected!')
end

[~,tbytCsvdateSort] = sort(datenum({tbytCsvList(:).date}, 'dd-mmm-yyyy hh:MM:ss'), 1, 'ascend'); % sorted fileList

%% organize the trial-by-trial csv and avi files
% spot the trial-by-trial video files
vFiles = dir('**/*.avi'); % list all the video files
if isempty(vFiles)
    error('No Video files were found!!')
end
vFronFiles = vFiles(cellfun(@(c)contains(c,'cam0'), {vFiles(:).name})); % front cam files
[vFronFiles(:).fileCalled] = deal(0); % accumulate the number of times the video file called before
[vFronFiles(:).framesUsed] = deal(0); % the number of frames previously assigned (to a trigger pulse) before
vSideFiles = vFiles(cellfun(@(c)contains(c,'cam1'), {vFiles(:).name})); % side cam files
[vSideFiles(:).fileCalled] = deal(0); % accumulate the number of times the video file called before
[vSideFiles(:).framesUsed] = deal(0); % the number of frames previously assigned before

vFileDateExp = '(20\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)'; % the expression of the date in the video file names

for v = 1:length(vFronFiles)
    vFronFileStartDatenum(v,1) = datenum(regexp(vFronFiles(v).name, vFileDateExp, 'match'),'yyyy_mm_dd_HH_MM_SS'); % front video file start date read from the video fileName
end

for v = 1:length(vSideFiles)
    vSideFileStartDatenum(v,1) = datenum(regexp(vSideFiles(v).name, vFileDateExp, 'match'),'yyyy_mm_dd_HH_MM_SS'); % front video file start date read from the video fileName
end

[S(:).csvFile] = deal(NaN);
[S(:).fVideo] = deal(NaN);
[S(:).sVideo] = deal(NaN);
[S(:).framePulseId] = deal(NaN);
[S(:).vFrameTime] = deal(NaN);
[S(:).origVideoFramesCnt] = deal(NaN);
[S(:).usedFrameCnt] = deal(NaN);
[S(:).framePulses] = deal(NaN);
[S(:).fVideoInfo] = deal(NaN);
[S(:).sVideoInfo] = deal(NaN);
[S(:).framePulses] = deal(NaN);
[S(:).fVideoStart] = deal(NaN);
[S(:).sVideoStart] = deal(NaN);
[S(:).fVideoCmplt] = deal(NaN);
[S(:).sVideoCmplt] = deal(NaN);
[S(:).tbytCsvFileCmplt] = deal(NaN);
[S(:).vUseFrameIdx] = deal(NaN); % logical indicating which video frames match current trial's trigger pulses

if ~isempty(vFronFiles)&&~isempty(vSideFiles)
    if abs(length(trStartIdx)-length(tbytCsvList))<=1 % the length of the tbytCsvList is supposed to be smaller than that of the trStartIdx (or S) as the last trial's csv file is not generated
        for t = 1:length(S)
            % find the front and side video files of the current trial
            tempTrCamPulseTRise = find(evtIdx1k.camTrigRiseIdx<S(t).trStart,1,'last'); % identify the pulse before the trStart the relevant train of pulses should encompass the trStart and trEnd
            tempTrCamPulseTFall = find(evtIdx1k.camTrigFallIdx>S(t).trEnd,1,'first');  % identify the pulse after the trEnd
            if ~isempty(tempTrCamPulseTRise)&&~isempty(tempTrCamPulseTFall)&& evtIdx1k.camPulseTrainIdx(tempTrCamPulseTRise)==evtIdx1k.camPulseTrainIdx(tempTrCamPulseTFall) % this ensures that the train pulses encompass the current trStart and trEnd
                
                tempPulseTrainId = evtIdx1k.camPulseTrainIdx(tempTrCamPulseTRise); % pulse train id
                tempPulseTrainLength = sum(evtIdx1k.camPulseTrainIdx==evtIdx1k.camPulseTrainIdx(tempTrCamPulseTRise)); % pulse train length (the # of pulses)
                tempPulseTStartIdx = find(evtIdx1k.camPulseTrainIdx==tempPulseTrainId,1,'first'); % locate the 1st cam trig pulse of this trial
                tempPulseTStopIdx = find(evtIdx1k.camPulseTrainIdx==tempPulseTrainId,1,'last');   % locate the last cam trig pulse of this trial
                
                if t<=length(S)-2 % for all trials except the last two
                    S(t).csvFile = fullfile(p.Results.filePath, behFilePath.name, tbytCsvList(tbytCsvdateSort(t)).name); % trial-by-trial csv file
                    
                    % identify the relevant front video file
                    tempCsvVFronTD = abs(datetime(tbytCsvList(tbytCsvdateSort(t)).date)-datetime({vFronFiles(:).date})); % absolute CSV and front Video file completion time difference
                    [tempMinCsvVFronTD,~] = min(tempCsvVFronTD); % the front video file index that potentially matches the current trial
                    tempMinCsvVFronTDI = find(tempCsvVFronTD==tempMinCsvVFronTD);
                    if length(tempMinCsvVFronTDI)==1
                        tempVFronI = tempMinCsvVFronTDI;
                    else % in case there're multiple files with the minimum time difference, choose the one with later fileStart
                        tempVFronI = tempMinCsvVFronTDI(find(vFronFileStartDatenum(tempMinCsvVFronTDI) <= tbytCsvList(tbytCsvdateSort(t)).datenum,1,'last'));
                    end
                    
                    % identify the relevant side video file
                    tempCsvVSideTD = abs(datetime(tbytCsvList(tbytCsvdateSort(t)).date)-datetime({vSideFiles(:).date})); % absolute CSV and side Video file completion time difference
                    [tempMinCsvVSideTD,~] = min(tempCsvVSideTD); % the side video file index that potentially matches the current trial
                    tempMinCsvVSideTDI = find(tempCsvVSideTD==tempMinCsvVSideTD);
                    if length(tempMinCsvVSideTDI)==1
                        tempVSideI = tempMinCsvVSideTDI;
                    else % in case there're multiple files with the minimum time difference, choose the one with later fileStart
                        tempVSideI = tempMinCsvVSideTDI(find(vSideFileStartDatenum(tempMinCsvVSideTDI) <= tbytCsvList(tbytCsvdateSort(t)).datenum,1,'last'));
                    end
                    
                    % read Video with mmread
                    tempVF = mmread(fullfile(vFronFiles(tempVFronI).folder, vFronFiles(tempVFronI).name),1); % get the trial-by-trial front cam video info, just read the first frame only for speed
                    tempVF.path = fullfile(vFronFiles(tempVFronI).folder, vFronFiles(tempVFronI).name);
                    tempVF.name = vFronFiles(tempVFronI).name;
                    
                    tempVS = mmread(fullfile(vSideFiles(tempVSideI).folder, vSideFiles(tempVSideI).name),1); % get the trial-by-trial side cam video info, just read the first frame only for speed
                    tempVS.path = fullfile(vSideFiles(tempVSideI).folder, vSideFiles(tempVSideI).name);
                    tempVS.name = vSideFiles(tempVSideI).name;
                    
                    if vFronFiles(tempVFronI).fileCalled==1 && abs(tempVF.totalDuration - vFronFiles(tempVFronI).framesUsed) < 4 % if video file used before, and most of the frames were assigned in previous calls
                        tempVFronI = tempVFronI + 1; % go try the next video file
                        tempVF = mmread(fullfile(vFronFiles(tempVFronI).folder, vFronFiles(tempVFronI).name),1); % reload the trial-by-trial front cam video info
                        tempVF.path = fullfile(vFronFiles(tempVFronI).folder, vFronFiles(tempVFronI).name);
                        tempVF.name = vFronFiles(tempVFronI).name;
                    end
                    
                    if vSideFiles(tempVSideI).fileCalled==1 && abs(tempVS.totalDuration - vFronFiles(tempVSideI).framesUsed) < 4 % if video file used before, and most of the frames were assigned in previous calls
                        tempVSideI = tempVSideI + 1; % go try the next video file
                        tempVS = mmread(fullfile(vSideFiles(tempVSideI).folder, vSideFiles(tempVSideI).name),1); % reload the trial-by-trial side cam video info
                        tempVS.path = fullfile(vSideFiles(tempVSideI).folder, vSideFiles(tempVSideI).name);
                        tempVS.name = vSideFiles(tempVSideI).name;
                    end
                    
                    if tempVF.totalDuration>1 && tempVS.totalDuration>1
                        [tempFrameDur,tempFrameDurI] = min([tempVF.totalDuration, tempPulseTrainLength]);
                        if tempFrameDurI==1
                        elseif tempFrameDurI==2
                            tempFrameDur=tempFrameDur-2;
                        end
                        
                        if vFronFiles(tempVFronI).fileCalled==0 && vSideFiles(tempVSideI).fileCalled==0 % if not used before take the pulse start as the frame reference
                            tempFrameTime = evtIdx1k.camTrigFallIdx(tempPulseTStartIdx+1:tempPulseTStartIdx+tempFrameDur); %evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempVF.totalDuration+1:tempPulseTStopIdx); % get the frame times
                            tempFrameIdx = false(tempVF.totalDuration,1);
                            tempFrameIdx(1:tempFrameDur)=true;
                        elseif vFronFiles(tempVFronI).fileCalled>=1 && vSideFiles(tempVSideI).fileCalled>=1 % if used before take the pulse stop as the frame reference
                            tempFrameTime = evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempFrameDur+1:tempPulseTStopIdx); %evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempVF.totalDuration+1:tempPulseTStopIdx); % get the frame times
                            tempFrameIdxPrevEnd = find(S(t-1).vUseFrameIdx==true,1,'last');
                            tempFrameIdx(tempFrameIdxPrevEnd+3:min(tempFrameIdxPrevEnd+2+tempFrameDur,tempVF.totalDuration))=true;
                            tempFrameIdx(1:tempFrameIdxPrevEnd)=false;
                        end
                        
                        % ensure that the vFile start time is before the tbytCsv file completion, also ensure that the vFile completion time is before the next tbytCsv file completion
                        %                     if vFronFileStartDatenum(tempVFronI)<=tbytCsvList(tbytCsvdateSort(t)).datenum && tbytCsvList(tbytCsvdateSort(t+1)).datenum>vFronFiles(tempVFronI).datenum ...
                        %                             && vSideFileStartDatenum(tempVSideI)<=tbytCsvList(tbytCsvdateSort(t)).datenum && tbytCsvList(tbytCsvdateSort(t+1)).datenum>vSideFiles(tempVSideI).datenum
                        
                        if vFronFileStartDatenum(tempVFronI)<=tbytCsvList(tbytCsvdateSort(t)).datenum && vSideFileStartDatenum(tempVSideI)<=tbytCsvList(tbytCsvdateSort(t)).datenum
                            if tempFrameTime(1)<S(t).trJsReady && S(t).trEnd<tempFrameTime(end) && abs(tempVF.totalDuration-tempVS.totalDuration)<=1 && vFronFiles(tempVFronI).fileCalled<3
                                % mark the video file usage
                                vFronFiles(tempVFronI).fileCalled = vFronFiles(tempVFronI).fileCalled + 1;
                                vSideFiles(tempVSideI).fileCalled = vSideFiles(tempVSideI).fileCalled + 1;
                                
                                % mark the number of frames already assigned
                                vFronFiles(tempVFronI).framesUsed = tempFrameDur;
                                vSideFiles(tempVSideI).framesUsed = tempFrameDur;
                                
                                S(t).fVideo = tempVF.path; % front video path
                                S(t).sVideo = tempVS.path; % side video path
                                S(t).vFrameTime = tempFrameTime; % mark the frame time points
                                S(t).framePulseId = tempPulseTrainId;
                                S(t).origVideoFramesCnt = tempVF.totalDuration;
                                S(t).fVideoInfo = tempVF;
                                S(t).sVideoInfo = tempVS;
                                S(t).framePulses = tempPulseTrainLength;
                                S(t).usedFrameCnt = tempFrameDur;
                                S(t).fVideoStart = datetime(vFronFileStartDatenum(tempVFronI),'ConvertFrom','datenum');
                                S(t).sVideoStart = datetime(vSideFileStartDatenum(tempVSideI),'ConvertFrom','datenum');
                                S(t).fVideoCmplt = datetime(vFronFiles(tempVFronI).datenum,'ConvertFrom','datenum');
                                S(t).sVideoCmplt = datetime(vSideFiles(tempVSideI).datenum,'ConvertFrom','datenum');
                                S(t).tbytCsvFileCmplt = datetime(tbytCsvList(tbytCsvdateSort(t)).datenum,'ConvertFrom','datenum');
                                S(t).vFronFileCalled = vFronFiles(tempVFronI).fileCalled; % # of file called
                                S(t).vSideFileCalled = vSideFiles(tempVSideI).fileCalled; % # of file called
                                S(t).vUseFrameIdx = tempFrameIdx;
                            end
                        end
                    end
                else % for the last two trials
                    if t==length(S)-1 % for the second last trial
                        % identify the relevant front video file
                        tempCsvVFronTD = abs(datetime(tbytCsvList(tbytCsvdateSort(t)).date)-datetime({vFronFiles(:).date})); % absolute CSV and front Video file completion time difference
                        [tempMinCsvVFronTD,~] = min(tempCsvVFronTD); % the front video file index that potentially matches the current trial
                        tempMinCsvVFronTDI = find(tempCsvVFronTD==tempMinCsvVFronTD);
                        if length(tempMinCsvVFronTDI)==1
                            tempVFronI = tempMinCsvVFronTDI;
                        else % in case there're multiple files with the minimum time difference
                            tempVFronI = tempMinCsvVFronTDI(find(vSideFileStartDatenum(tempMinCsvVFronTDI) <= tbytCsvList(tbytCsvdateSort(t)).datenum,1,'last'));
                        end
                        % identify the relevant side video file
                        tempCsvVSideTD = abs(datetime(tbytCsvList(tbytCsvdateSort(t)).date)-datetime({vSideFiles(:).date})); % absolute CSV and side Video file completion time difference
                        [tempMinCsvVSideTD,~] = min(tempCsvVSideTD); % the side video file index that potentially matches the current trial
                        tempMinCsvVSideTDI = find(tempCsvVSideTD==tempMinCsvVSideTD);
                        if length(tempMinCsvVSideTDI)==1
                            tempVSideI = tempMinCsvVSideTDI;
                        else % in case there're multiple files with the minimum time difference
                            tempVSideI = tempMinCsvVSideTDI(find(vSideFileStartDatenum(tempMinCsvVSideTDI) <= tbytCsvList(tbytCsvdateSort(t)).datenum,1,'last'));
                        end
                        
                        S(t).csvFile = fullfile(p.Results.filePath, behFilePath.name, tbytCsvList(tbytCsvdateSort(t)).name); % trial-by-trial csv file
                        S(t).tbytCsvFileCmplt = datetime(tbytCsvList(tbytCsvdateSort(t)).datenum,'ConvertFrom','datenum');
                        %tempVFronI = find(vFronFileStartDatenum <tbytCsvList(tbytCsvdateSort(t)).datenum,1,'last'); % the front video file index that potentially matches the current trial
                        %tempVSideI = find(vSideFileStartDatenum <tbytCsvList(tbytCsvdateSort(t)).datenum,1,'last'); % the side video file index that potentially matches the currnet trial
                    elseif t==length(S) % for the last trial
                        S(t).csvFile = NaN; % trial-by-trial csv file not saved for the final trial
                        tempVFronI = length(vFronFiles);
                        tempVSideI = length(vSideFiles);
                    end
                    tempVF = mmread(fullfile(vFronFiles(tempVFronI).folder, vFronFiles(tempVFronI).name),1); % get the trial-by-trial front cam video info
                    tempVF.path = fullfile(vFronFiles(tempVFronI).folder, vFronFiles(tempVFronI).name);
                    tempVF.name = vFronFiles(tempVFronI).name;
                    
                    tempVS = mmread(fullfile(vSideFiles(tempVSideI).folder, vSideFiles(tempVSideI).name),1); % get the trial-by-trial side cam video info
                    tempVS.path = fullfile(vSideFiles(tempVSideI).folder, vSideFiles(tempVSideI).name);
                    tempVS.name = vSideFiles(tempVSideI).name;
                    
                    if tempVF.totalDuration>1 && tempVS.totalDuration>1
                        [tempFrameDur,tempFrameDurI] = min([tempVF.totalDuration, tempPulseTrainLength]);
                        if tempFrameDurI==1
                        elseif tempFrameDurI==2
                            tempFrameDur=tempFrameDur-2;
                        end
                        %tempFrameTime = evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempVF.totalDuration+1:tempPulseTStopIdx); % get the frame times
                        if vFronFiles(tempVFronI).fileCalled==0 && vSideFiles(tempVSideI).fileCalled==0 % if not used before take the pulse start as the frame reference
                            tempFrameTime = evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempVF.totalDuration+1:tempPulseTStopIdx); %evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempVF.totalDuration+1:tempPulseTStopIdx); % get the frame times
                            tempFrameIdx = false(tempVF.totalDuration,1);
                            tempFrameIdx(1:tempFrameDur)=true;
                        elseif vFronFiles(tempVFronI).fileCalled>=1 && vSideFiles(tempVSideI).fileCalled>=1 % if used before take the pulse stop as the frame reference
                            tempFrameTime = evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempFrameDur+1:tempPulseTStopIdx); %evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempVF.totalDuration+1:tempPulseTStopIdx); % get the frame times
                            tempFrameIdxPrevEnd = find(S(t-1).vUseFrameIdx==true,1,'last');
                            tempFrameIdx(tempFrameIdxPrevEnd+3:min(tempFrameIdxPrevEnd+2+tempFrameDur,tempVF.totalDuration))=true;
                            tempFrameIdx(1:tempFrameIdxPrevEnd)=false;
                        end
                        
                        if vFronFiles(tempVFronI).fileCalled<3 && vSideFiles(tempVSideI).fileCalled<3 && tempFrameTime(1)<S(t).trJsReady && S(t).trEnd<tempFrameTime(end) && tempVF.totalDuration==tempVS.totalDuration && abs(tempPulseTrainLength-tempVF.totalDuration)<4 % the video frame counts happen to be consistently fewer than the # of frame pulses by 2
                            % mark the video file usage
                            vFronFiles(tempVFronI).fileCalled = vFronFiles(tempVFronI).fileCalled + 1;
                            vSideFiles(tempVSideI).fileCalled = vSideFiles(tempVSideI).fileCalled + 1;
                            
                            % mark the number of frames already assigned
                            vFronFiles(tempVFronI).framesUsed = tempFrameDur;
                            vSideFiles(tempVSideI).framesUsed = tempFrameDur;
                            
                            S(t).fVideo = tempVF.path; % front video path
                            S(t).sVideo = tempVS.path; % side video path
                            S(t).vFrameTime = tempFrameTime; % mark the frame time points
                            S(t).framePulseId = tempPulseTrainId;
                            S(t).origVideoFramesCnt = tempVF.totalDuration;
                            S(t).fVideoInfo = tempVF;
                            S(t).sVideoInfo = tempVS;
                            S(t).framePulses = tempPulseTrainLength;
                            S(t).usedFrameCnt = tempFrameDur;
                            S(t).fVideoStart = datetime(vFronFileStartDatenum(tempVFronI),'ConvertFrom','datenum');
                            S(t).sVideoStart = datetime(vSideFileStartDatenum(tempVSideI),'ConvertFrom','datenum');
                            S(t).fVideoCmplt = datetime(vFronFiles(tempVFronI).datenum,'ConvertFrom','datenum');
                            S(t).sVideoCmplt = datetime(vSideFiles(tempVSideI).datenum,'ConvertFrom','datenum');
                            S(t).vFronFileCalled = vFronFiles(tempVFronI).fileCalled; % the # of file called
                            S(t).vSideFileCalled = vSideFiles(tempVSideI).fileCalled; % the # of file called
                            S(t).vUseFrameIdx = tempFrameIdx;
                        end
                    end
                end
            end
            fprintf('organized csv and video files for trial #%d\n', t);
            clearvars tempVF tempVS temp*
        end
    else
        error('Some trial-by-trial csv files might be missing!')
    end
else
    warning('No video files were found!')
end

% rename and save the structure
jsTime1k_KV = S;
save(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles'),'jsTime1k_KV');

