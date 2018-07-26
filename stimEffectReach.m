function [stimR] = stimEffectReach( filePath, experimentType, varargin )
%'stimEffectReach' examines the effect of laser stim on the reach behavior (kinematics). 
%  It generates/stores a structure 'stimR' which contains relevant quantities such as max velocity, position.  

cd(filePath)

fileList = dir('BehVariables*'); % get the behavioral data BehVariables*.mat 

%% Preprocessing of the behavioral data
stimR.reachPos = cell(length(fileList),1);      % reach position
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
    
    % get stimPos1 of laser Reaches and pseudo-laser reaches
    [ tempLaserReachStart, ~, tempLaserPos1 ] = getPseudoReachTimesForStimTrials( tempPositionData, tempTS.stmLaser ); % detect pseudo-reach trajectories around the time of laser stimulation deliveries. 
    [ ~, ~, tempPseudoLaserPos1 ] = getPseudoReachTimesForStimTrials( tempPositionData, tempTS.pseudoLaser ); %  detect the pseudo-reach trajectories around the time of pseudo Laser TTL pulses (without actual laser delivery)
    
    stimR.laserReachPos{f,1} = tempLaserPos1; % laser reach amplitudes 
    stimR.pseudoLaserReachPos{f,1} = tempPseudoLaserPos1; % pseudo laser reach amplitudes
    
    tempDcmStmPos1 = zeros(size(tempLaserPos1,1), round(size(tempLaserPos1,2)/10)); % preallocate tempDcmStmPos1 (decimated by a factor of 10)
    for t = 1:size(tempLaserPos1,1)       % increment laserReaches
        tempDcmStmPos1(t,:) = decimate(tempLaserPos1(t,:),10); % decimate reach position data of stim trials
    end
    clearvars t
    stimR.laserReachPos1Dcm{f,1} = tempDcmStmPos1; % decimated stim reach amplitudes
    
    tempDcmPseudoStmPos1 = zeros(size(tempPseudoLaserPos1,1), round(size(tempPseudoLaserPos1,2)/10)); % preallocate tempDcmPseudoStmPos1 (decimated by a factor of 10)
    for t = 1:size(tempPseudoLaserPos1,1)
        tempDcmPseudoStmPos1(t,:) = decimate(tempPseudoLaserPos1(t,:),10); % decimate pseudo reach position data of the pseudo-stim trials 
    end
    clearvars t
    stimR.pseudoLaserReachPos1Dcm{f,1} = tempDcmPseudoStmPos1; % decimated pseudo-stim reach amplitudes

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
[stimR.posPCcoeff, stimR.posPCscore] = pca(stimR.reachPos1Dcm(:,1:100)); % plot(stimR.posCoeff(:,2))
stimR.stmPos1PCscore = stimR.laserReachPos1Dcm{1,1}(:,1:100)*stimR.posPCcoeff; % project stmPos1 data to the principal components of the no-stim pos1 trajectories
stimR.pseudoStmPos1PCscore = stimR.pseudoLaserReachPos1Dcm{1,1}(:,1:100)*stimR.posPCcoeff; % project stmPos1 data to the principal components of the no-stim pos1 trajectories

%figure; hold on; plot(stimR.posPCscore(:,1),stimR.posPCscore(:,2),'.','MarkerSize',20); plot(stimR.stmPos1PCscore(:,1),stimR.stmPos1PCscore(:,2),'r.','MarkerSize',20); hold off;
figure; hold on; 
plot3(stimR.posPCscore(:,1),stimR.posPCscore(:,2),stimR.posPCscore(:,3),'.','MarkerSize',20); 
plot3(stimR.stmPos1PCscore(:,1),stimR.stmPos1PCscore(:,2),stimR.stmPos1PCscore(:,3),'r.','MarkerSize',20); 
plot3(stimR.pseudoStmPos1PCscore(:,1), stimR.pseudoStmPos1PCscore(:,2), stimR.pseudoStmPos1PCscore(:,3),'g.','MarkerSize',20); 
hold off;
grid on

%% Further quantification of reaches
stimR.velReaches    = diff([zeros(size(stimR.reachPos1Dcm,1),1)  stimR.reachPos1Dcm(:,1:end)],1,2); % velocity of reaches on decimated position data
stimR.maxVelReaches = max(stimR.velReaches,[],2); % max velocity of full reaches

stimR.velStimReaches = diff([zeros(size(stimR.laserReachPos1Dcm{1,1}(:,1:end),1),1)  stimR.laserReachPos1Dcm{1,1}(:,1:end)],1,2); % velocity of stim-reaches
stimR.maxVelStimReaches = max(stimR.velStimReaches,[],2); % max velocity of stim reaches

stimR.velPseudoStimReaches = diff([zeros(size(stimR.pseudoLaserReachPos1Dcm{1,1}(:,1:end),1),1)  stimR.pseudoLaserReachPos1Dcm{1,1}(:,1:end)],1,2); % velocity of pseudo-stim reaches
stimR.maxVelPseudoStimReaches = max(stimR.velPseudoStimReaches,[],2); % max velocity of pseudo-stim reaches

stimR.maxPosReaches = max(stimR.reachPos1Dcm(:,1:end),[],2); % max position of full reaches
stimR.maxPosStimReaches = max(stimR.laserReachPos1Dcm{1,1}(:,1:end),[],2); % max velocity of stim reaches
stimR.maxPosPseudoStimReaches = max(stimR.pseudoLaserReachPos1Dcm{1,1}(:,1:end),[],2); % max velocity of pseudo-stim reaches

%% save stimR
if ismac
    delimitIdx = strfind(filePath,'/');
elseif ispc
    delimitIdx = strfind(filePath,'\');
else
    disp('Platform not supported')
end

saveName = strcat(filePath(delimitIdx(end)+1:end),'_stimR');
save(saveName,'stimR')

%% Plot with Gramm
stimR.expLabel = [zeros(length(stimR.maxPosReaches),1)+1; zeros(length(stimR.maxPosPseudoStimReaches),1)+2; zeros(length(stimR.maxPosStimReaches),1)+3]; % zeros(length(maxPosConStimReaches),1)+3]; % experiment label
stimR.maxPosArr = [ stimR.maxPosReaches; stimR.maxPosPseudoStimReaches; stimR.maxPosStimReaches ]; % concatenated max amp (pos1) data
stimR.maxVelArr = [ stimR.maxVelReaches; stimR.maxVelPseudoStimReaches; stimR.maxVelStimReaches ]; % concatenated max vel (diff(pos1)) data

color = [zeros(length(stimR.maxPosReaches),1); zeros(length(stimR.maxPosPseudoStimReaches),1)+1; stimR.stimReachIdx{1}+2]; % to mark the stim reach trials with a different color

clear g
g(1,1)=gramm('x',stimR.expLabel,'y',stimR.maxPosArr,'color',color);
g(1,2)=gramm('x',stimR.expLabel,'y',stimR.maxVelArr,'color',color);

g(1,1).geom_jitter('width',0.4,'height',0);
g(1,1).set_names('x','label','y','maxPos');
g(1,1).set_title('Max Pos');
 
g(1,2).geom_jitter('width',0.4,'height',0);
g(1,2).set_names('x','label','y','maxVel');
g(1,2).set_title('Max Vel');

fig1 = figure('Position',round([100 100 800 400]./1.5));
g.draw();

print(fig1, strcat(saveName, experimentType), '-dpdf'); % print figure as pdf

%gf.coord_flip();
%figure('Position',[100 100 800 400]);
%gf.draw();


end