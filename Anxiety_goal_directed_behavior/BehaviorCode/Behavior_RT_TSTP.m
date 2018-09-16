%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MUST RUN THIS SCRIPT TO PUT THE BEHAVIORAL DATA THE TSTP.EVENT CELLS  %
%   IN THE THIRD COLUMN OF THE TSTP.EVENT                                 %      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to get the following behavioral indices.
% 1) islatency: latency from cue onset to instrumental poke. 
% 2) ftlatency: latency from food delivery to food retrieval.
% 3) isimmobile: absolute and percent time spent in immobility between cue
% onset and instrumental pokes.
% 5) immobileITI: 1st col: length of prior ITI, 2nd col: abs. immobile time
% during the prior ITI, 3rd col: per. immobile time during the prior ITI. 
% The output STRUCTURE 'tstpbeh' contains: 
%   1.name (filename) 2.cue tstp 3.ispk tstp 4.pellet tstp 5.ftpk tstp 6.ispk latency 
%   7.immobile ispklatency 8.ftpklatency 9.normalized ispklatency 10.normalized immobile ispklatency 
%   11.normalized ftpklatency   

% B1: P(O|A) = 1, P(S|A) = 0
% B2: P(O|A) = 1, P(S|A) = 0.06
% B3: P(O|A) = 1, P(S|A) = 0.1

% Ascending order session: B1 - B2 - B3
% Descending order session: B3 - B2 - B1 

clear all; close all; clc

%addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

%% Extract data from txt files 
tstpDirectory = '/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/MAT';
cd(tstpDirectory)
load('SUA_AGB_tstp_022315.mat'); clearvars tstpDirectory
FileDirectory = '/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/Behavior/TXT';      % directory for the txt files, this folder must only contain the valid txt files to be analyzed, all the other files must be stored in TXT_Whole  
cd(FileDirectory)
SaveDirectory = '/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/MAT';       % directory to save the outcome of this script
SaveName = 'SUA_AGB_tstpbeh_022315';
FileList = dir('*GS*');

GLMcell = cell(size(FileList,1),1);         % GLMcell (Session by Rat) is a cell containing the BEHmat of all the animals 

for i = 1:size(FileList,1)          % numbrat
    cd(FileDirectory)
    FileName = FileList(i).name;
    [~,~,~,~,~,~,bv.subject,bv.time,~,bv.state,bv.event,bv.motion,~,bv.instpoke,bv.FTpoke,bv.motionOff,bv.instpokeOff,~,bv.FTpokeOff,~,~,~,~,~,~,~,~] = textread(FileName,'%s %s %d %d %d %d %s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d', 'headerlines',1);
    % bv.subject = animal ID, bv.time = time in unit (20 ms), bv.state = state (changes), bv.motion = thermodetector, bv.instpoke = instrumental pokes (A3), bv.FTpoke = food trough pokes
    %% Calculate the latency from each cue onset to instrumental poke
    [B1.ISstatebegin, B1.ISstateshift, B1.islatency, B1.isimmobile, B1.isimmobilecount, B1.immobileITI, B1.isimmobilecountITI] = statedetect( 2, 3, 9, bv.state, bv.time, 20/1000, bv.motion );           % state 2 = BLOCK1 CUE, state 3 = BLOCK1 ISPOKE, For ispokes, use state 2(cue) as index state to avoid complications due to reminder cues. 20/1000 is the time constant.
    [B2.ISstatebegin, B2.ISstateshift, B2.islatency, B2.isimmobile, B2.isimmobilecount, B2.immobileITI, B2.isimmobilecountITI] = statedetect( 12, 13, 19, bv.state, bv.time, 20/1000, bv.motion );        % state 12= BLOCK2 CUE, state 13= BLOCK2 ISPOKE
    [B3.ISstatebegin, B3.ISstateshift, B3.islatency, B3.isimmobile, B3.isimmobilecount, B3.immobileITI, B3.isimmobilecountITI] = statedetect( 22, 23, 29, bv.state, bv.time, 20/1000, bv.motion );        % state 22= BLOCK3 CUE, state 23= BLOCK3 ISPOKE
    [B4.ISstatebegin, B4.ISstateshift, B4.islatency, B4.isimmobile, B4.isimmobilecount, B4.immobileITI, B4.isimmobilecountITI] = statedetect( 32, 33, 39, bv.state, bv.time, 20/1000, bv.motion );        % state 22= BLOCK3 CUE, state 23= BLOCK3 ISPOKE
      
    %% Calculate the latency from each pellet drop to pellet take
    [~, ~, B1.ftlatency, ~, ~, ~, ~] = statedetect( 6, 6, 9, bv.state, bv.time, 20/1000, bv.motion );           % bv.state 6 = TAKE PELLET1
    [~, ~, B2.ftlatency, ~, ~, ~, ~] = statedetect( 16, 16, 19, bv.state, bv.time, 20/1000, bv.motion );        % bv.state 16 = TAKE PELLET2
    [~, ~, B3.ftlatency, ~, ~, ~, ~] = statedetect( 26, 26, 29, bv.state, bv.time, 20/1000, bv.motion );        % bv.state 26 = TAKE PELLET3
    [~, ~, B4.ftlatency, ~, ~, ~, ~] = statedetect( 36, 36, 39, bv.state, bv.time, 20/1000, bv.motion );        % bv.state 36 = TAKE PELLET4
       
    %% construct the BEHmat
    B1.trialnumb = length(B1.ISstatebegin);         % the # of trials for B1
    B2.trialnumb = length(B2.ISstatebegin);         % the # of trials for B2
    B3.trialnumb = length(B3.ISstatebegin);         % the # of trials for B3
    B4.trialnumb = unique(~isnan(B4.ISstatebegin)*length(B4.ISstatebegin));         % the # of trials for B4
    
    B1.BEHmat = nan(B1.trialnumb,6);        % preallocate the BEHmat for B1, the number of columns must correspond the number of variables in GLM
    B1.BEHmat(:,1) = B1.ISstateshift;       % the state shift info is to be used for sorting the BEHmat
    B1.BEHmat(:,2) = 1;                     % put the number of blocks
    B1.BEHmat(:,3:5) = [ B1.islatency, B1.isimmobile(:,1), B1.ftlatency ];      % put islatency, isimmobile (absolute time), ftlatency to the BEHmat
    B1.BEHmat(:,6) = B1.ISstatebegin;
    
    B2.BEHmat = nan(B2.trialnumb,6);        % preallocate the BEHmat for B2, the number of columns must correspond the number of variables in GLM
    B2.BEHmat(:,1) = B2.ISstateshift;       % the state shift info is to be used for sorting the BEHmat
    B2.BEHmat(:,2) = 2;                     % put the number of blocks
    B2.BEHmat(:,3:5) = [ B2.islatency, B2.isimmobile(:,1), B2.ftlatency ];      % put islatency, isimmobile (absolute time), ftlatency to the BEHmat
    B2.BEHmat(:,6) = B2.ISstatebegin;
    
    B3.BEHmat = nan(B3.trialnumb,6);        % preallocate the BEHmat for B3, the number of columns must correspond the number of variables in GLM
    B3.BEHmat(:,1) = B3.ISstateshift;       % the state shift info is to be used for sorting the BEHmat
    B3.BEHmat(:,2) = 3;                     % put the number of blocks
    B3.BEHmat(:,3:5) = [ B3.islatency, B3.isimmobile(:,1), B3.ftlatency ];      % put islatency, isimmobile (absolute time), ftlatency to the BEHmat
    B3.BEHmat(:,6) = B3.ISstatebegin;
    
    if sum(isnan(B4.ISstatebegin)) == 1       % in case that B4 is invalid
        B4.BEHmat = NaN;        % preallocate the BEHmat for B4, the number of columns must correspond the number of variables in GLM
        tempGLMmat = [ B1.BEHmat; B2.BEHmat; B3.BEHmat ];        % concatenate the BEHmat across 4 blocks
    else      % if B4 is valid
        B4.BEHmat = nan(B4.trialnumb,5);        % preallocate the BEHmat for B4, the number of columns must correspond the number of variables in GLM
        B4.BEHmat(:,1) = B4.ISstateshift;       % the state shift info is to be used for sorting the BEHmat
        B4.BEHmat(:,2) = 0;                     % put the number of blocks
        B4.BEHmat(:,3:5) = [ B4.islatency, B4.isimmobile(:,1), B4.ftlatency ];      % put islatency, isimmobile (absolute time), ftlatency to the BEHmat
        B4.BEHmat(:,6) = B4.ISstatebegin;
        tempGLMmat = [ B4.BEHmat; B1.BEHmat; B2.BEHmat; B3.BEHmat ];        % concatenate the BEHmat across 4 blocks
    end
    
    tempGLMmat = sortrows(tempGLMmat,1);            % sort the concatenated BEHmat based on the number of states
      
%% match the trial numbers between new and old files, due to the trial number reduction from 180 to 150    
    BEHmat = nan(150,6);    % preallocate the BEHmat, the number of trials must be fixed at 150 trials, the number of columns (variables) may vary
    if size(tempGLMmat,1) == 150    % if the number of trials = 150: no change occurs
        BEHmat = tempGLMmat;
    elseif size(tempGLMmat,1) == 180    % if the number of trials = 180: reorganize trials
        BEHmat(1:50,:) = tempGLMmat(1:50,:);       % take B1 trials: first 50 trials: B0(30) + B1(20)
        BEHmat(51:100,:) = tempGLMmat(81:130,:);     % take B2 trials
        BEHmat(101:150,:) = tempGLMmat(131:180,:);   % take B3 trials
    else
        warning('Check the number of trials')
    end
        
%% match the trial numbers, mismatch happened due to missing trials in the recording file, and put the matched information to the tstp cell for the future analyses (decoding) 
    temptstpispk = tstp.ispk{i,1};
    if length(temptstpispk) < 150 && length(BEHmat) == 150    % in case, there's a missing trial in recorded file
        gstime = bv.time(BEHmat(:,1))*(20/1000);     % graphic state time, length = 150
        plxtime = temptstpispk;    % plexon time, length < 150 (due to a missing trial
        
        numbmismatch = length(BEHmat) - length(temptstpispk);
        % get the trial-by-trial time difference in the plexon time
        for u = 1:length(plxtime)
            if u < length(plxtime)
                plxtime(u,2) = plxtime(u+1,1)-plxtime(u,1);
            else
                plxtime(u,2) = nan;
            end
        end
        
        % get the trial-by-trial time difference in the graphic state time
        for o = 1:length(gstime)
            if o < length(gstime)
                gstime(o,2) = gstime(o+1,1)-gstime(o,1);
            else
                gstime(o,2) = nan;
            end
        end
        
        tempmissingref = nan(150,1);
        for j = 51:150      % only deal with the B2 or B3 trials, as there's no missing trials in B1
            tempmissingref(j,1) = abs(nansum(gstime(j:end,2)-plxtime(j-numbmismatch:end,2)));       % use the 'block match' to find the missing point
        end 
        
        missingpt = min(find(tempmissingref < 1));      % the point where the unmissed point resumes (one (or two) points before this missing point were omitted in the recorded file)
        
        crtGLMmat = nan(150-numbmismatch,size(BEHmat,2));
        crtGLMmat(1:missingpt-numbmismatch-1,:) = BEHmat(1:missingpt-numbmismatch-1,:);
        crtGLMmat(missingpt-numbmismatch:end,:) = BEHmat(missingpt:end,:);
        crtGLMmat = crtGLMmat(:,3:5);       % select variables: ISletency, immobile ISlatency, FT latency
        %tstp.ispk{i,3} = crtGLMmat; tstp.pellet{i,3} = crtGLMmat; tstp.iti{i,3} = crtGLMmat;     % store the BEHmat into the GLMcell in the 3rd column of the tstp cell
    else    % no change occurs
        crtGLMmat = BEHmat(:,3:5);       % select variables: ISletency, immobile ISlatency, FT latency
        %tstp.ispk{i,3} = crtGLMmat; tstp.pellet{i,3} = crtGLMmat; tstp.iti{i,3} = crtGLMmat;      % store the BEHmat into the GLMcell in the 3rd column of the tstp cell
    end
clearvars u o j missingpt

 %% Create a structure containing behavioral data
    tstpbeh(i).name = tstp.ispk{i,2};       % label first!
    tstpbeh(i).cue = tstp.cue{i,1};         % get the cue timestamps 
    tstpbeh(i).ispk = tstp.ispk{i,1};       % get the ispk timestamps
    tstpbeh(i).pellet = tstp.pellet{i,1};   % get the cue timestamps
    tstpbeh(i).ftpk = tstp.iti{i,1};         % get the cue timestamps    
    tstpbeh(i).islate = crtGLMmat(:,1);     % raw instrumental poke latency
    tstpbeh(i).imbislate = crtGLMmat(:,2);  % raw immobile instrumental poke latency
    tstpbeh(i).ftlate = crtGLMmat(:,3);     % raw food trough poke latency   
    
%% Normalize the islatency and ftlatency to the entire session mean and std
    [tempmean,tempstd,~] = meanstdsem(crtGLMmat(1:50,1:3));
    tempnormcrtGLMmat = nan(size(crtGLMmat,1),size(crtGLMmat,2));
    for jj = 1:length(crtGLMmat)
        tempnormcrtGLMmat(jj,:) = (crtGLMmat(jj,1:3) - tempmean)./tempstd;
    end
    tstpbeh(i).normislate = tempnormcrtGLMmat(:,1);
    tstpbeh(i).normimbislate = tempnormcrtGLMmat(:,2);
    tstpbeh(i).normftlate = tempnormcrtGLMmat(:,3);
    
    clearvars jj temp*
    
%% Correct the mismatch of the numbers cue and ispk timestamps ONLY FOR 'AD02'  
    switch strcmp(tstpbeh(i).name,'AD02_Day4')
        case 1      % in case of 'AD02_Day4' the numbers of cue and ispk timestamps must be matched
            mismatchpt = 51;    % start off screening from the B2 trials
            while tstpbeh(i).cue(mismatchpt) > tstpbeh(i).ispk(mismatchpt-1)    % mismatch is defined by the point where the cue_timestamp(t) < ispk_timestamp(t-1)
                  mismatchpt = mismatchpt + 1;  
            end
            % correct cue timestamps
            mismatchcuetstp = tstpbeh(i).cue;   % the old mismatched cue timestamps
            tstpbeh(i).cue = [];       % empty the cue field to correct it
            tstpbeh(i).cue = cat(1,mismatchcuetstp(1:mismatchpt-1,1),mismatchcuetstp(mismatchpt+1:end,1));
            % correct ftpk timestamps
            mismatchftpktstp = tstpbeh(i).ftpk;   % the old mismatched ftpk timestamps
            tstpbeh(i).ftpk = [];       % empty the ftpk field to correct it
            tstpbeh(i).ftpk = cat(1,mismatchftpktstp(1:mismatchpt-1,1),mismatchftpktstp(mismatchpt+1:end,1));
        case 0      % all the other files must bypass this
    end
    clearvars mismatchpt mismatchcuetstp mismatchftpktstp   
    
%%  Correct the mismatch between the numbers of ftpk and ispk timestamps
    if length(tstpbeh(i).ftpk) > length(tstpbeh(i).ispk)
        extftpktstp = tstpbeh(i).ftpk;      % extended ftpk
        tstpbeh(i).ftpk = [];
        tstpbeh(i).ftpk = extftpktstp(1:length(tstpbeh(i).ispk),1);
    end
    
end
clearvars bv i ii tempGLMmat 

%% Save GLMcell
cd(SaveDirectory)
save(SaveName,'tstpbeh')       % update the tstp with the behavioral data added in the 3rd column of the tstp cell

