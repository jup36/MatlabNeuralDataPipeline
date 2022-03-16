function [S] = JsVideoFileOrganizer_vUseFrameIdx(filePath)
%This function inspects the trial-by-trial front and side videos, and
% assign them to corresponding trials. The output is a structure named 'jsTime1k_KV'.
% This function has been updated in september/2019 to add a field named
% 'vUseFrameIdx', which is a logical important to determine which frames
% from the video comprise the hand trajectory of the current trial. This is
% especially critical when using a video file for multiple trials (e.g. S:\Junchol_Data\JS2.0\WR40_081419\jsTime1k_Kinematics_VideoFiles.mat\jsTime1k_KV(17).vUseFrameIdx).

cd(filePath)

pathJsTime1k_KV = dir('**/*_kinematics_VideoFiles.mat');
S=load(fullfile(pathJsTime1k_KV(1).folder,pathJsTime1k_KV(1).name),'jsTime1k_KV'); % just building up on the outcome of jsKinematicsAnalysis.m to create one all-inclusive file
S = S.('jsTime1k_KV');

pathBehVar = dir('**/*BehVariablesJs.mat');
load(fullfile(pathBehVar(1).folder,pathBehVar(1).name), 'evtIdx1k', 'p')
trStartIdx = evtIdx1k.trStartIdx;

behFilePath = dir(fullfile(filePath,'20*-*')); % dir where the trial-by-trial behavioral csv files are saved
tbytCsvList = dir(fullfile(behFilePath(1).folder,behFilePath(1).name,'trial_*'));    % trial-by-trial files
allTrialCsv = dir(fullfile(behFilePath(1).folder,behFilePath(1).name,'trials.csv')); % all trial file
if length(allTrialCsv)==1
    trialsFileName = fullfile(allTrialCsv.folder,allTrialCsv.name);
    trialsCsv = readtable(trialsFileName);
else
    error('More than one trials.csv file detected!')
end

%[~,tbytCsvdateSort] = sort(datenum({tbytCsvList(:).date}, 'dd-mmm-yyyy hh:MM:ss'), 1, 'ascend'); % sorted fileList

%%
[S(:).vUseFrameIdx] = deal(NaN); % logical indicating which video frames match current trial's trigger pulses

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
            
            if isstruct(S(t).fVideoInfo) && isstruct(S(t).sVideoInfo)
                if S(t).fVideoInfo.totalDuration>1 && S(t).sVideoInfo.totalDuration>1
                    [tempFrameDur,tempFrameDurI] = min([S(t).fVideoInfo.totalDuration, tempPulseTrainLength]);
                    if tempFrameDurI==1
                    elseif tempFrameDurI==2
                        tempFrameDur=tempFrameDur-2;
                    end
                    
                    if S(t).vFronFileCalled==1 && S(t).vSideFileCalled==1 % if not used before take the pulse start as the frame reference
                        tempFrameTime = evtIdx1k.camTrigFallIdx(tempPulseTStartIdx+1:tempPulseTStartIdx+tempFrameDur); %evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempVF.totalDuration+1:tempPulseTStopIdx); % get the frame times
                        tempFrameIdx = false(S(t).fVideoInfo.totalDuration,1);
                        tempFrameIdx(1:tempFrameDur)=true;
                    elseif S(t).vFronFileCalled>=2 && S(t).vSideFileCalled>=2 % if used before take the pulse stop as the frame reference
                        tempFrameTime = evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempFrameDur+1:tempPulseTStopIdx); %evtIdx1k.camTrigFallIdx(tempPulseTStopIdx-tempVF.totalDuration+1:tempPulseTStopIdx); % get the frame times
                        tempFrameIdxPrevEnd = find(S(t-1).vUseFrameIdx==true,1,'last');
                        tempFrameIdx = false(S(t).fVideoInfo.totalDuration,1);
                        tempFrameIdx(tempFrameIdxPrevEnd+3:min(tempFrameIdxPrevEnd+2+tempFrameDur,S(t).fVideoInfo.totalDuration))=true;
                        tempFrameIdx(1:tempFrameIdxPrevEnd)=false;
                    end
                    
                    %if vFronFileStartDatenum(tempVFronI)<=tbytCsvList(tbytCsvdateSort(t)).datenum && vSideFileStartDatenum(tempVSideI)<=tbytCsvList(tbytCsvdateSort(t)).datenum
                    if tempFrameTime(1)<S(t).trJsReady && S(t).trEnd<tempFrameTime(end) && abs(S(t).fVideoInfo.totalDuration-S(t).sVideoInfo.totalDuration)<=1 && S(t).vFronFileCalled<4
                        S(t).vUseFrameIdx = tempFrameIdx;
                    end
                    %end
                end
            end
        end
        fprintf('organized csv and video files for trial #%d\n', t);
        %clearvars tempVF tempVS temp*
    end
else
    error('Some trial-by-trial csv files might be missing!')
end


% rename and save the structure
jsTime1k_KV = S;
save(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles'),'jsTime1k_KV');

