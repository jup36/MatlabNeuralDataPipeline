function [ reachStart, reachMW, pos1, xpos1, ypos1 ] = getPseudoReachTimesForStimTrials( positionData, stmLaser )
%This function detects pseudo-reach trajectories around the time of laser
% stimulation deliveries.

TbetweenR = 2000; % minimum time between reaches, in ms
maxDur=1500;      % max reach duration

xPos = positionData(1,:); % Xposition
yPos = positionData(2,:); % Yposition

reachSG = sgolayfilt(sqrt((yPos-median(yPos)).^2 + (xPos-median(xPos)).^2),3,201); % smoothing with sgolay filter
reachMW = moveavg(sqrt((yPos-median(yPos)).^2 + (xPos-median(xPos)).^2),10);       % smoothing with moving average filter
XposMW  = moveavg(xPos,10); % smoothing with moving average filter
YposMW  = moveavg(yPos,10); % smoothing with moving average filter

vel=[0, diff(reachSG)]; % velocity - differentiation of the position trajectory

sortPos=sort(reachSG);  % sort position data
%sortVel=sort(vel);      % sort velocity data
sortPos(sortPos<median(sortPos))=median(sortPos); % replace less-than median position values with the median position value
cutoff  = sortPos(round(length(sortPos)*.90));    % set the 90% cutoff for the position data
%cutoffd = sortVel(round(length(sortVel)*.98));    % set the 98% cutoff for the velocity data (velocity cutoff is required, since reaches are expected to occur at high velocity)

%% detect reach start
reachStart=[]; % the goal here is to find the near-zero point that is immediately prior to each high-velocity point
for i  = 1:length(stmLaser)
    if stmLaser(i)+1000 < length(reachSG)
        currHighAmpPoints = find(reachSG(1,stmLaser(i):stmLaser(i)+1000)>cutoff) + stmLaser(i);  % current valid sample (high-velocity data point)
        
        if ~isempty(currHighAmpPoints) % in case there's no point above the position cutoff
            
            currHighAmpPointsIdx = zeros(length(currHighAmpPoints),1);
            for ii = 1:length(currHighAmpPoints)
                currHighAmpPointsIdx(ii,1) = sum(reachSG(1,currHighAmpPoints(ii):currHighAmpPoints(ii)+200)>cutoff)>20;
            end
            
            if ~isempty(find(currHighAmpPointsIdx==1,1))
                currSamp = currHighAmpPoints(find(currHighAmpPointsIdx==1,1));
            else
                [~,currSamp] = max(reachSG(1,stmLaser(i):stmLaser(i)+1000),[],2);
                currSamp = currSamp + stmLaser(i);
            end
            
        else
            [~,currSamp] = max(reachSG(1,stmLaser(i):stmLaser(i)+1000),[],2);
            currSamp = currSamp + stmLaser(i);
        end
        
        nearZero= find(reachSG(1:currSamp)<cutoff*.4,1,'last'); % find the most recent near zero-point
        if nearZero>10
            if currSamp-nearZero>1000  % sometimes the joystick doesn't reset to zero and induces 20 second long reaches. this prevents that by setting the second (additional) cutoff
                sortPos2=sort((reachSG(nearZero:currSamp)));
                cutoff2 = sortPos2(round(length(sortPos2)*.90)); % reset the cutoff2 to be a point immediatly before the reach peak (farthest point)
                nearZero=find(reachSG(1:currSamp)<cutoff2,1,'last');
            end
            realStartMinus10=find(vel(1:nearZero)<0,10,'last'); % find 10 elements whose velocity is negative prior to the current near zero point
        else % if there's no near zero point
            realStartMinus10=currSamp; % take the current valid sample (high-velocity data point) as a reach start point
        end
        reachStart = [reachStart, realStartMinus10(1)];
    end
    
    % this will eliminate pseudo starts, you might not want to do eliminate them
    pos1=[];    % reach position data aligned to reachStart
    xpos1=[];   % x position data aligned to reachStart normalized to the x position at the reach start
    ypos1=[];   % y position data aligned to reachStart normalized to the y position at the reach start
    for i = reachStart
        if i > 1000
            if i+1499 > length(reachMW) % in case the reach trajectory end exceeds the length of the time series (this must happen only for the last reach)
                pos1 =[pos1; nan(1,size(pos1,2))]; % position data corresponding to each reach from 200-ms before to 1500-ms after each reach start
                xpos1=[xpos1; nan(1,size(xpos1,2))];
                ypos1=[ypos1; nan(1,size(ypos1,2))];
                
            else % in case the reach trajectory is within the length of the time series
                pos1  = [pos1; reachMW(i-200:i+1499)-reachMW(i)]; % position data corresponding to each reach from 200-ms before to 1500-ms after each reach start
                xpos1 = [xpos1; XposMW(i-200:i+1499)-XposMW(i)];  % x position data aligned to reachStart normalized to the x position at the reach start
                ypos1 = [ypos1; YposMW(i-200:i+1499)-YposMW(i)];  % y position data aligned to reachStart normalized to the y position at the reach start
            end
        end
    end
    % figure; imagesc(pos1); colorbar
    % title('start times, aligned at 200ms')
    
end
end

