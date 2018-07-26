% Create the filename list from the current directory

fileDirectory = '';
type = 2;

disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(['________________Running RetrieveAllBehaviorData... Script_________________________'])
disp(' ');
disp(['Directory: ' fileDirectory]);
disp(['Timestamp: ' datestr(now)]);
disp(' ');

% create a list of all the files
allFiles = dir(sprintf('%s*_p.csv',fileDirectory));

for k = 1:size(allFiles,1)

    switch type
        
        case 1

            fileNameBase = allFiles(k).name(1:strfind(allFiles(k).name,'_')-1);
            year = fileNameBase(strfind(fileNameBase,'201'):strfind(fileNameBase,'201')+3);
            month = fileNameBase(strfind(fileNameBase,'201')+4:strfind(fileNameBase,'201')+4);
            if month==1
                month = fileNameBase(strfind(fileNameBase,'201')+4:strfind(fileNameBase,'201')+5);
                day = fileNameBase(strfind(fileNameBase,'201')+6:strfind(fileNameBase,'201')+7);
            else
                month = ['0' fileNameBase(strfind(fileNameBase,'201')+4:strfind(fileNameBase,'201')+4)];
                day = fileNameBase(strfind(fileNameBase,'201')+5:strfind(fileNameBase,'201')+6);
            end

        case 2
            
    end
            
    disp(' ');
    disp(['File ' num2str(k) ' base name: ' fileNameBase]);
    disp(['Date created: ' year '/' month '/' day]);    
    disp(' ');

    % Script to load csv behavior data
    clear disp* trial* 

    skipInitTrials = 0;

    fileNameStr = [fileDirectory fileNameBase '_p.csv'];
    trialParams = dlmread(fileNameStr,',',1+skipInitTrials,0);

    fileNameStr = [fileDirectory fileNameBase '.csv'];
    trialTimes = dlmread(fileNameStr,',',2+skipInitTrials,0);

    fileNameStr = [fileDirectory fileNameBase 'tX.csv'];
    tmpX = dlmread(fileNameStr,',',1+skipInitTrials,0);
    trajectories.dispX = tmpX(:,1:100);

    fileNameStr = [fileDirectory fileNameBase 'tY.csv'];
    tmpY = dlmread(fileNameStr,',',1+skipInitTrials,0);
    trajectories.dispY = tmpY(:,1:100);

    % Analyze the properties of behavioral response times
    trialData.trialParams = trialParams;
    trialData.trialTimes = trialTimes;

    indArray = 2:size(trialTimes,1);

    trialData.startToEvent  = trialTimes(indArray,3) - trialTimes(indArray,2);
    trialData.endToStart    = trialTimes(indArray,3) - trialTimes(indArray-1,5);

    % Smooth the trajectories to get rid of the 8-bit discretization ugliness

    for i=1:size(trajectories.dispX,1)
        trajectories.dispXnorm(i,:) = trajectories.dispX(i,:) - mean(trajectories.dispX(i,1:5),2);
        trajectories.dispYnorm(i,:) = trajectories.dispY(i,:) - mean(trajectories.dispY(i,1:5),2); 
        trajectories.dispXsg(i,:)   = sgolayfilt(trajectories.dispXnorm(i,:),3,11);
        trajectories.dispYsg(i,:)   = sgolayfilt(trajectories.dispYnorm(i,:),3,11);
    end

    % Create vector sums of the trajectories and report angles, polar histogram, and lengths
    figure(5); clf;

    for i=1:size(trajectories.dispX,1)

        % get peak and average velocity
        j = 2:size(trajectories.dispX,2);
        velocity(j-1) = sqrt(trajectories.dispXsg(i,j).^2 + trajectories.dispYsg(i,j).^2);
        trajectories.velocity(i,:) = velocity;

        tmpXe = cumsum(trajectories.dispXsg(i,1:50),2);
        tmpYe = cumsum(trajectories.dispYsg(i,1:50),2);
        tmpXe=tmpXe/50;
        tmpYe=tmpYe/50;

        tmpXl = cumsum(trajectories.dispXsg(i,96:100),2);
        tmpYl = cumsum(trajectories.dispYsg(i,96:100),2);
        tmpXl=tmpXl/5;
        tmpYl=tmpYl/5;

        % should normalize to the length and also store the length
        currVectorXe = tmpXe(size(tmpXe,2));
        currVectorYe = tmpYe(size(tmpXe,2));

        % should normalize to the length and also store the length
        currVectorXl = tmpXl(size(tmpXl,2));
        currVectorYl = tmpYl(size(tmpXl,2));

        if (currVectorXe.^2+currVectorYe.^2)==0
            trajectories.lengthE(i,:) = 0;
        else
            trajectories.lengthE(i,:) = sqrt(currVectorXe.^2+currVectorYe.^2);    
        end
        trajectories.vectorXe(i,:) = currVectorXe./trajectories.lengthE(i,:);
        trajectories.vectorYe(i,:) = currVectorYe./trajectories.lengthE(i,:);

        if (currVectorXl.^2+currVectorYl.^2)==0
            trajectories.lengthL(i,:) = 0;
        else
            trajectories.lengthL(i,:) = sqrt(currVectorXl.^2+currVectorYl.^2);
        end
        trajectories.vectorXl(i,:) = currVectorXl./trajectories.lengthL(i,:);
        trajectories.vectorYl(i,:) = currVectorYl./trajectories.lengthL(i,:);

    end

    [trajectories.thetaE,trajectories.rhoE] = cart2pol(trajectories.vectorXe,trajectories.vectorYe);
    needZeros = find(isnan(trajectories.thetaE)); % CATCH THE CASES WITH 0 LENGTH VECTORS
    disp(['Number of zero vectors: ' num2str(length(needZeros))]);
    trajectories.thetaE(needZeros) = 0;
    trajectories.polarHistE = hist(trajectories.thetaE,-pi:0.1:pi);

    [trajectories.thetaL,trajectories.rhoL] = cart2pol(trajectories.vectorXl,trajectories.vectorYl);
    needZeros = find(isnan(trajectories.thetaL));
    trajectories.thetaL(needZeros) = 0;
    trajectories.polarHistL = hist(trajectories.thetaL,-pi:0.1:pi);

    trajectories.vacillate = trajectories.thetaL-trajectories.thetaE;

    figure(3); clf;
    subplot(141);
    compass(trajectories.vectorXe,trajectories.vectorYe,'r'); hold on;
    compass(trajectories.vectorXl,trajectories.vectorYl,'b'); hold on;
    subplot(142);
    rose(trajectories.thetaL); hold on;
    rose(trajectories.thetaE);
    subplot(1,4,3:4);
    plot(trajectories.lengthL,'b.-'); hold on;
    plot(trajectories.lengthE,'r.-');
    axis([0 200 0 50]);

    figure(5); clf;
    subplot(211);
    plot(trajectories.velMean,'k.-');
    axis([0 200 0 20]);
    subplot(212);
    plot(trajectories.velMax,'k.-');
    axis([0 200 0 50]);

    % Combine the session data into a single structure

    if k>1
        disp(['Current number of sessions: ' num2str(size(behavData.session,2))]);
        z = size(behavData.session,2)+1;
    else
        disp(['Current number of sessions: 0']);
        z = 1;
    end
    
    behavData.session(z).trajectories = trajectories;
    behavData.session(z).trialData = trialData;
    behavData.session(z).dateStamp.year = year;
    behavData.session(z).dateStamp.month = month;
    behavData.session(z).dateStamp.day = day;
    
    % want some simple number (like days since January 1, 2008) that allows sequencing of behavior data for every mouse
    behavData.session(z).dateStamp.nDate = datenum([day '-' month '-' year]); % returns n = 730625
    
end

% Write the data structure to a specific name
disp(' ');
disp(' ');
eval([fileNameBase(1,1:3) '= behavData;']);
clear behavData;

disp(['Output data in the structure: ' fileNameBase(1,1:3)]);
disp(' ');
disp(['________________Completed Script_______________________________________________________'])
disp(' ');
disp(' ');
disp(' ');
disp(' ');


%     % Plot the inidivudal trajectories
%     figure(2); clf;
% 
%     for i=1:size(trajectories.dispX,1)
%         figure(2); hold off;
%         plot(trajectories.dispXsg(i,:),trajectories.dispYsg(i,:),'k.-',[0 trajectories.vectorXe(i).*10], [0 trajectories.vectorYe(i).*10],'r',[0 trajectories.vectorXl(i).*10], [0 trajectories.vectorYl(i).*10],'b');
%         axis([-40 40 -40 40]);
%         drawnow;
%         pause(0.1);
%     end
