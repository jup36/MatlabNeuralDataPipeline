function [stimR] = stimEffectReach( filePath, varargin )
%'stimEffectReach' examines the effect of laser stim on the reach behavior. 
%  It generates/stores a structure 'stimR' which contains relevant quantities such as max velocity, position.  
%  

cd(filePath)

fileList = dir('BehVariables*'); % get the behavioral data BehVariables*.mat 

%% Preprocessing of the behavioral data
stimR.reachPos = cell(length(fileList),1);       % reach position
stimR.laserReachPos = cell(length(fileList),1); % reach position detected near laser stim
stimR.reachPos1Dcm   = cell(length(fileList),1);
stimR.laserReachPos1Dcm = cell(length(fileList),1);
stimR.stimReachIdx   = cell(length(fileList),1);

for f = 1:length(fileList) % increment files
    load(fileList(f).name,'pos1'); % moving-averaged reach position (amplitude)
    load(fileList(f).name,'positionData'); % xPos; yPos
    load(fileList(f).name,'ts');   % timestamps
    
    tempPos1 = pos1; clearvars pos1 
    tempPositionData = positionData; clearvars positionData
    tempTS   = ts; clearvars ts
    
    stimR.reachPos{f,1} = tempPos1;
    tempDcmPos1 = zeros(size(tempPos1,1), size(tempPos1,2)/10); % preallocate tempDcmStmPos1 (decimated by a factor of 10)
    for t = 1:size(tempPos1,1)
        tempDcmPos1(t,:) = decimate(tempPos1(t,:),10); % decimated Pos1 (moving-averaged reach position) 
    end
    clearvars t
    stimR.reachPos1Dcm{f,1} = tempDcmPos1;    % decimated non-stim reach amplitudes
    
    % get stimPos1 of laserReaches
    [ tempLaserReachStart, ~, tempLaserPos1 ] = getPseudoReachTimesForStimTrials( tempPositionData, tempTS.stmLaser ); % detect pseudo-reach trajectories around the time of laser stimulation deliveries. 
    [ tempPseudoLaserReachStart, ~, tempPseudoLaserPos1 ] = getPseudoReachTimesForStimTrials( tempPositionData, tempTS.pseudoLaser ); %  detect the pseudo-reach trajectories around the time of pseudo Laser TTL pulses (without actual laser delivery)
    
    stimR.laserReachPos{f,1} = tempLaserPos1; % laser reach amplitudes 
    stimR.pseudoLaserReachPos{f,1} = tempPseudoLaserPos1; % pseudo laser reach amplitudes
    
    tempDcmStmPos1 = zeros(size(tempLaserPos1,1), size(tempLaserPos1,2)/10); % preallocate tempDcmStmPos1 (decimated by a factor of 10)
    for t = 1:size(tempLaserPos1,1)       % increment laserReaches
        tempDcmStmPos1(t,:) = decimate(tempLaserPos1(t,:),10);
    end
    clearvars t
    stimR.laserReachPos1Dcm{f,1} = tempDcmStmPos1; % decimated stim pseudo-reach amplitudes
    
    % get the tempStmReachIdx by matching ts.stmLaser and ts.stmLaserReach 
    tempStmReachIdx = zeros(length(tempLaserReachStart),1); % index for stmReaches (laserReaches actually satisfying reach criterion)
    for r = 1:length(tempTS.stmLaserReach) % increment stmReach Trials
        tempStmReachIdx(find(tempTS.stmLaser==tempTS.stmLaserReach(r),1)) = true; % get the tempStmReachIdx
    end
    clearvars r
    
    stimR.stimReachIdx{f,1} = tempStmReachIdx;
    
end

%% pca on the decimted pos1 trajectories 
stimR.reachPos1Dcm = cell2mat(stimR.reachPos1Dcm);   % pca on decimated reach amp trajectories
[stimR.posCoeff, stimR.posScore] = pca(stimR.reachPos1Dcm(:,1:100)); % plot(stimR.posCoeff(:,2))
stimR.stmPos1Score = stimR.laserReachPos1Dcm{1,1}(:,1:100)*stimR.posCoeff; % project stmPos1 data to the principal components of the no-stim pos1 trajectories

%figure; hold on; plot(stimR.posScore(:,1),stimR.posScore(:,2),'.','MarkerSize',20); plot(stimR.stmPos1Score(:,1),stimR.stmPos1Score(:,2),'r.','MarkerSize',20); hold off;
figure; hold on; plot3(stimR.posScore(:,1),stimR.posScore(:,2),stimR.posScore(:,3),'.','MarkerSize',20); plot3(stimR.stmPos1Score(:,1),stimR.stmPos1Score(:,2),stimR.stmPos1Score(:,3),'r.','MarkerSize',20); hold off;
grid on

%% Further quantification of reaches
stimR.velReaches    = diff([zeros(size(stimR.reachPos1Dcm,1),1)  stimR.reachPos1Dcm(:,1:end)],1,2); % velocity of reaches on decimated position data
stimR.maxVelReaches = max(stimR.velReaches,[],2); % max velocity of full reaches

stimR.velStimReaches = diff([zeros(size(stimR.laserReachPos1Dcm{1,1}(:,1:end),1),1)  stimR.laserReachPos1Dcm{1,1}(:,1:end)],1,2); % velocity of full reaches
stimR.maxVelStimReaches = max(stimR.velStimReaches,[],2); % max velocity of stim reaches

stimR.maxPosReaches = max(stimR.reachPos1Dcm(:,1:end),[],2); % max position of full reaches
stimR.maxPosStimReaches = max(stimR.laserReachPos1Dcm{1,1}(:,1:end),[],2); % max velocity of stim reaches

%% Plot with Gramm
stimR.expLabel = [zeros(length(stimR.maxPosReaches),1)+1; zeros(length(stimR.maxPosStimReaches),1)+2]; % zeros(length(maxPosConStimReaches),1)+3]; % experiment label
stimR.maxPosArr = [ stimR.maxPosReaches; stimR.maxPosStimReaches ]; % concatenated max amp (pos1) data
stimR.maxVelArr = [ stimR.maxVelReaches; stimR.maxVelStimReaches ]; % concatenated max vel (diff(pos1)) data

color = [zeros(length(stimR.maxPosReaches),1); stimR.stimReachIdx{1}]; % to mark the stim reach trials with a different color

clear g
g(1,1)=gramm('x',stimR.expLabel,'y',stimR.maxPosArr,'color',color);
g(1,2)=copy(g(1));

g(1,1).geom_point();
g(1,1).set_names('x','label','y','maxPos');
g(1,1).set_title('No groups');

g(1,2).geom_jitter('width',0.4,'height',0);
g(1,2).set_title('geom_jitter()');

gf = copy(g);

figure('Position',[100 100 800 400]);
g.draw();

gf.coord_flip();
figure('Position',[100 100 800 400]);
gf.draw();

%% save 
if ismac
    delimitIdx = strfind(filePath,'/');
elseif ispc
    delimitIdx = strfind(filePath,'\');
else
    disp('Platform not supported')
end

saveName = strcat(filePath(delimitIdx(end)+1:end),'_stimR');
save(saveName,'stimR')

end