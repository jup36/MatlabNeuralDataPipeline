% This script is to get the following behavioral indices.
% 1) islatency: latency from cue onset to instrumental poke. 
% 2) ftlatency: latency from food delivery to food retrieval.
% 3) isimmobile: absolute and percent time spent in immobility between cue
% onset and instrumental pokes.
% 4) isimmobilecount: # of immobile times counted if longer than 0.5
% sec. 
% 5) immobileITI: 1st col: length of prior ITI, 2nd col: abs. immobile time
% during the prior ITI, 3rd col: per. immobile time during the prior ITI. 
% 6) isimmobilecountITI: # of immobile times during the prior ITI counted
% if longer than 0.5 sec. 
% 7) shockdistmat: calculates distance from prev(1st col) and next(2nd col) shock trials. 
%    To generate shockdistmat 'distanceshock.m' must be run, and the cell
%    containing shock distance information e.g. 'shockcellbatch3_1' must be
%    ready. 

% B1: P(O|A) = 1, P(S|A) = 0
% B2: P(O|A) = 1, P(S|A) = 0.06
% B3: P(O|A) = 1, P(S|A) = 0.1

% Ascending order session: B1 - B2 - B3
% Descending order session: B3 - B2 - B1 

clear all; close all; clc

addpath('/Users/parkj/Desktop/OldPC/Matlab/');
addpath('/Users/parkj/Desktop/OldPC/Matlab/Functions/'); 
addpath('/Users/parkj/Desktop/OldPC/Matlab/Functions/Anxiety_goal-directed_behavior/Behavior_package/');

%% Extract data from txt files 
FileDirectory = '/Users/parkj/Google Drive/Chowdhury/GS_DATA (TGC)/TGC47-52_GS12_noinjectionDay5_022616/';  % directory for the txt files, this folder must only contain the valid txt files to be analyzed, all the other files must be stored in TXT_Whole  
cd(FileDirectory)
SaveDirectory = '/Users/parkj/Google Drive/Chowdhury/';  % directory to save the outcome of this script
SaveGLMcellname = 'AGB_TGCcell_061817';
FileList = dir('*GS*');

% experiment variables. Modify this!
numbses = 1;       % # of sessions
numbrat = 6;      % # of animals (files)
numbtrial = 50;    % # of trials for B1, B2, B3

%FileList = reshape(FileList, numbses, numbrat);         % reshape the FileList to arrange sessions of each animal in each column (e.g. JP01 in the first column)

% get shock trials 
%cd('C:\Documents and Settings\jup36\Desktop\Data\Anxiety_goal_directed_behavior\RecBatch_data\MAT')       % temporarily change directory just to load 'Batch3-1_shocktrials'. 
%load('Recbatch_shocktrials','shockcellrecbatch')        % load shocktrials cell: this cell contains shock trial information of B2 & B3 of every session of every animal in the current batch. This cell must be created and saved in the directory specified in the above line before running this script. 

GLMcell = cell(numbrat,numbses);         % GLMcell (Session by Rat) is a cell containing the GLMmat of all the animals 

for i = 1:size(FileList,1)          % numbrat
    for ii = 1:size(FileList,2)     % numbses
        cd(FileDirectory)
        FileName = FileList(i,ii).name;
        fid   = fopen(FileName);
        SaveName = [FileName(1:9) '.mat'];
        
        textdata = textscan(fid,'%s %s %d %d %d %d %s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d', 'headerlines',1);
        bv.subject = textdata{1,7};
        bv.time = textdata{1,8};
        bv.state = textdata{1,10};
        bv.event = textdata{1,11};
        bv.instpoke = textdata{1,12};
        bv.motion = textdata{1,14};
        bv.FTpoke = textdata{1,15};
        bv.instpokeOff = textdata{1,16};
        bv.motionOff = textdata{1,18};
        bv.FTpokeOff = textdata{1,19};
          
        % bv.subject = animal ID, bv.time = time in unit (20 ms), bv.state = state (changes), bv.motion = thermodetector, bv.instpoke = instrumental pokes (A3), bv.FTpoke = food trough pokes
        %% Calculate the latency from each cue onset to instrumental poke
        [B1.ISstatebegin, B1.ISstateshift, B1.islatency, B1.isimmobile, B1.isimmobilecount, B1.immobileITI, B1.isimmobilecountITI] = statedetect( 2, 3, 9, bv.state, bv.time, 20/1000, bv.motion );           % state 2 = BLOCK1 CUE, state 3 = BLOCK1 ISPOKE, For ispokes, use state 2(cue) as index state to avoid complications due to reminder cues. 20/1000 is the time constant.      
        [B2.ISstatebegin, B2.ISstateshift, B2.islatency, B2.isimmobile, B2.isimmobilecount, B2.immobileITI, B2.isimmobilecountITI] = statedetect( 12, 13, 19, bv.state, bv.time, 20/1000, bv.motion );        % state 12= BLOCK2 CUE, state 13= BLOCK2 ISPOKE
        [B3.ISstatebegin, B3.ISstateshift, B3.islatency, B3.isimmobile, B3.isimmobilecount, B3.immobileITI, B3.isimmobilecountITI] = statedetect( 22, 23, 29, bv.state, bv.time, 20/1000, bv.motion );        % state 22= BLOCK3 CUE, state 23= BLOCK3 ISPOKE
        [B4.ISstatebegin, B4.ISstateshift, B4.islatency, B4.isimmobile, B4.isimmobilecount, B4.immobileITI, B4.isimmobilecountITI] = statedetect( 32, 33, 39, bv.state, bv.time, 20/1000, bv.motion );        % state 22= BLOCK3 CUE, state 23= BLOCK3 ISPOKE       
                
        [B1.avg_ispoke,~,B1.sem_ispoke] = meanstdsem(B1.islatency);
        [B2.avg_ispoke,~,B2.sem_ispoke] = meanstdsem(B2.islatency);
        [B3.avg_ispoke,~,B3.sem_ispoke] = meanstdsem(B3.islatency);
        [B4.avg_ispoke,~,B4.sem_ispoke] = meanstdsem(B4.islatency);
               
        %% z-score normalization of the instrumental RTs using the mean and std of the RTs in the B4.
        B1.normislatency = (B1.islatency-B1.avg_ispoke)./nanstd(B1.islatency,0,1);
        B2.normislatency = (B2.islatency-B1.avg_ispoke)./nanstd(B1.islatency,0,1);
        B3.normislatency = (B3.islatency-B1.avg_ispoke)./nanstd(B1.islatency,0,1);
        B4.normislatency = (B4.islatency-B1.avg_ispoke)./nanstd(B1.islatency,0,1);
        
        %% Calculate the latency from each pellet drop to pellet take
        [~, ~, B1.ftlatency, ~, ~, ~, ~] = statedetect( 6, 6, 9, bv.state, bv.time, 20/1000, bv.motion );           % bv.state 6 = TAKE PELLET1
        [~, ~, B2.ftlatency, ~, ~, ~, ~] = statedetect( 16, 16, 19, bv.state, bv.time, 20/1000, bv.motion );        % bv.state 16 = TAKE PELLET2
        [~, ~, B3.ftlatency, ~, ~, ~, ~] = statedetect( 26, 26, 29, bv.state, bv.time, 20/1000, bv.motion );        % bv.state 26 = TAKE PELLET3
        [~, ~, B4.ftlatency, ~, ~, ~, ~] = statedetect( 36, 36, 39, bv.state, bv.time, 20/1000, bv.motion );        % bv.state 36 = TAKE PELLET4
        
        [B1.avg_ftpoke,~,B1.sem_ftpoke] = meanstdsem(B1.ftlatency);
        [B2.avg_ftpoke,~,B2.sem_ftpoke] = meanstdsem(B2.ftlatency);
        [B3.avg_ftpoke,~,B3.sem_ftpoke] = meanstdsem(B3.ftlatency);
        [B4.avg_ftpoke,~,B4.sem_ftpoke] = meanstdsem(B4.ftlatency);
                               
        %% z-score normalization of the FT RTs using the mean and std of the RTs in the B4.
        B1.normftlatency = (B1.ftlatency-B4.avg_ftpoke)./std(B4.ftlatency,0,1);
        B2.normftlatency = (B2.ftlatency-B4.avg_ftpoke)./std(B4.ftlatency,0,1);
        B3.normftlatency = (B3.ftlatency-B4.avg_ftpoke)./std(B4.ftlatency,0,1);
        B4.normftlatency = (B4.ftlatency-B4.avg_ftpoke)./std(B4.ftlatency,0,1);
        
        %% get the behavioral indices.
        contrastIS_B3B1 = (B3.avg_ispoke-B1.avg_ispoke)/(B3.avg_ispoke+B1.avg_ispoke);
        contrastIS_B2B1 = (B2.avg_ispoke-B1.avg_ispoke)/(B2.avg_ispoke+B1.avg_ispoke);
        contrastFT_B3B1 = (B3.avg_ftpoke-B1.avg_ftpoke)/(B3.avg_ftpoke+B1.avg_ftpoke);
        contrastFT_B2B1 = (B2.avg_ftpoke-B1.avg_ftpoke)/(B2.avg_ftpoke+B1.avg_ftpoke);
        
        %% construct the GLMmat 
        B1.trialnumb = length(B1.ISstatebegin);         % the # of trials for B1
        B2.trialnumb = length(B2.ISstatebegin);         % the # of trials for B2
        B3.trialnumb = length(B3.ISstatebegin);         % the # of trials for B3
        B4.trialnumb = unique(~isnan(B4.ISstatebegin)*length(B4.ISstatebegin));         % the # of trials for B4
        
        B1.GLMmat = nan(B1.trialnumb,5);        % preallocate the GLMmat for B1, the number of columns must correspond the number of variables in GLM 
        B1.GLMmat(:,1) = B1.ISstateshift;       % the state number info is to be used for sorting the GLMmat 
        B1.GLMmat(:,2) = 1;                     % put the number of blocks 
        B1.GLMmat(:,3:5) = [ B1.islatency, B1.isimmobile(:,1), B1.ftlatency ];      % put islatency, isimmobile (absolute time), ftlatency to the GLMmat
                
        B2.GLMmat = nan(B2.trialnumb,5);        % preallocate the GLMmat for B2, the number of columns must correspond the number of variables in GLM 
        B2.GLMmat(:,1) = B2.ISstateshift;       % the state number info is to be used for sorting the GLMmat 
        B2.GLMmat(:,2) = 2;                     % put the number of blocks 
        B2.GLMmat(:,3:5) = [ B2.islatency, B2.isimmobile(:,1), B2.ftlatency ];      % put islatency, isimmobile (absolute time), ftlatency to the GLMmat
        
        B3.GLMmat = nan(B3.trialnumb,5);        % preallocate the GLMmat for B3, the number of columns must correspond the number of variables in GLM 
        B3.GLMmat(:,1) = B3.ISstateshift;       % the state number info is to be used for sorting the GLMmat 
        B3.GLMmat(:,2) = 3;                     % put the number of blocks 
        B3.GLMmat(:,3:5) = [ B3.islatency, B3.isimmobile(:,1), B3.ftlatency ];      % put islatency, isimmobile (absolute time), ftlatency to the GLMmat
        
        if sum(isnan(B4.ISstatebegin)) == 1       % in case that B4 is invalid
            B4.GLMmat = NaN;        % preallocate the GLMmat for B4, the number of columns must correspond the number of variables in GLM
            tempGLMmat = [ B1.GLMmat; B2.GLMmat; B3.GLMmat ];        % concatenate the GLMmat across 4 blocks
        else      % if B4 is valid
            B4.GLMmat = nan(B4.trialnumb,5);        % preallocate the GLMmat for B4, the number of columns must correspond the number of variables in GLM
            B4.GLMmat(:,1) = B4.ISstateshift;       % the state number info is to be used for sorting the GLMmat
            B4.GLMmat(:,2) = 0;                     % put the number of blocks
            B4.GLMmat(:,3:5) = [ B4.islatency, B4.isimmobile(:,1), B4.ftlatency ];      % put islatency, isimmobile (absolute time), ftlatency to the GLMmat
            tempGLMmat = [ B4.GLMmat; B1.GLMmat; B2.GLMmat; B3.GLMmat ];        % concatenate the GLMmat across 4 blocks
        end
        
        tempGLMmat = sortrows(tempGLMmat,1);            % sort the concatenated GLMmat based on the number of states                     
             
        % replace the number of states with the number of trials in the 1st column  
        %GLMmat(1:size(GLMmat,1),1) = [1:size(GLMmat,1)]';                     
                       
        GLMmat = nan(150,5);    % preallocate the GLMmat, the number of trials must be fixed at 150 trials, the number of columns (variables) may vary
        if size(tempGLMmat,1) == 150    % if the number of trials = 150: no change occurs
            GLMmat = tempGLMmat;
        GLMcell{i,ii} = GLMmat;    
        elseif size(tempGLMmat,1) == 180    % if the number of trials = 180: reorganize trials
            GLMmat(1:50,:) = tempGLMmat(1:50,:);       % take B1 trials: first 50 trials: B0(30) + B1(20)
            GLMmat(51:100,:) = tempGLMmat(81:130,:);     % take B2 trials
            GLMmat(101:150,:) = tempGLMmat(131:180,:);   % take B3 trials
        GLMcell{i,ii} = GLMmat;    
        else
        warning('Check the number of trials')  
        GLMcell{i,ii} = GLMmat;     %  store the GLMmat into the GLMcell
        
        end       
    end
end
clearvars i ii tempGLMmat GLMmat


%% Save GLMcell
save(SaveGLMcellname,'GLMcell')

%% Just to creat a structure femaleBeh
for i = 1:length(GLMcell)
    
    tempFileName        = [FileList(i).name(1:5) '_' FileList(i).name(23:26)]; 
    femaleBeh(i).name   = tempFileName;
    femaleBeh(i).RT     = GLMcell{i,1}(:,3);
    femaleBeh(i).foodRT = GLMcell{i,1}(:,5);
    
end

