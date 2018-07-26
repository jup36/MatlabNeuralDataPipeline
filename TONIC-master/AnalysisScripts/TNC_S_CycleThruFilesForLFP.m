%% RUN THROUGH ONE DIRECTORY OF FILES TO EXTRACT LFP, REACHES

paramsC.tapers      = [5 9];
paramsC.pad         = 0;
paramsC.Fs          = 500; % reflects the fact that the data was decimated down to 1 kHz sampling frequency
paramsC.fpass       = [2 115];
paramsC.err         = 0;
movingwin           = [0.4 0.01];

elecType = 'short';
justMiddle = 1;     % just use the middle electrodes of the array for LFP analysis
onlyOut = 1;        % select for only outward movements

disp(' ');
disp(' ');
disp(['Initialized the parameters for the analysis of LFP data ...']);
disp(' ');
disp(' ');

fileList = dir('*.ns5');
numFiles = size(fileList,1);

for q = 1:numFiles

    filenamestr = fileList(q).name;
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Beginning file: ' filenamestr(1,1:length(filenamestr)-3)]);
    disp(' ');disp(' ');  
    
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % LOAD EVENT DATA
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    
    filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'nev']
    dataEvents = openNEV(filenamestrE,'read','nosave','nomat');

    rewIndsTmp = find(dataEvents.Data.Spikes.Electrode==142); % reward
        % eliminate any pulses that were detected twice
        rewTsTmp = round(dataEvents.Data.Spikes.Timestamps(rewIndsTmp)./30);
        validRewInds = find(diff(rewTsTmp)>10);
        % write valid data to structure
        ContData.behavior.rewardInds = rewTsTmp(validRewInds);

    thrIndsTmp = find(dataEvents.Data.Spikes.Electrode==143); % threshold crossing
        % eliminate any pulses that were detected twice
        thrTsTmp = round(dataEvents.Data.Spikes.Timestamps(thrIndsTmp)./30);
        validThrInds = find(diff(thrTsTmp)>10);
        % write valid data to structure
        ContData.behavior.threshInds = thrTsTmp(validThrInds);

    lickInds = find(dataEvents.Data.Spikes.Electrode==139); % licking
        ContData.behavior.lickInds = round(dataEvents.Data.Spikes.Timestamps(lickInds)./30);
    
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % EXTRACT THE LFP DATA
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------

    switch elecType

        case 'short'
            if justMiddle==1
                chanList = [   20 12 5 61 58 50   ];            
            else
                chanList = [28 20 12 5 61 58 50 42];
            end
            ContData.chanList   = chanList;

        case 'long'

            if justMiddle==1
                chanList = [   24 16 8 33 64 54   ];     
            else
                chanList = [27 24 16 8 33 64 54 46];
            end
            ContData.chanList   = chanList;

    end

    for j=1:numel(chanList)
        ContData.shank(j).repChannel = chanList(j);
    end

    % cycle through all channels and filter/calculate power within each band
    for j=1:numel(chanList)

        tstart = tic;

        i = chanList(j);

        electrode = ['e:' num2str(i)];
        Ns5DATA = openNSx('report','read',filenamestr, electrode);    

        rawData = decimate(Ns5DATA.Data(1,:),60);    

        % Filter the data to select for spikes
        disp(' ');
        disp(' ');
        disp(['Filtering data on channel ' num2str(i) ' ...']);

        decimateCheck = numel(Ns5DATA.Data(1,:)) / numel(rawData);
        disp(['Number of samples in lowBand: ' num2str(numel(rawData)) ' ... from ' num2str(numel(Ns5DATA.Data(1,:))) ' original samples. Downsampled by ' num2str(decimateCheck)]);

        % design the bandpass filter to use for field potentials
        dBT = fdesign.bandpass('n,fc1,fc2', 300, 18, 32, paramsC.Fs);
        HdBT = design(dBT);   
        beta = filtfilt(HdBT.Numerator,1,rawData); % zero-phase filtering

        dGH = fdesign.bandpass('n,fc1,fc2', 300, 65, 115, paramsC.Fs);
        HdGH = design(dGH);  
        gammaHi = filtfilt(HdGH.Numerator,1,rawData); % zero-phase filtering
        
        dGL = fdesign.bandpass('n,fc1,fc2', 300, 40, 55, paramsC.Fs);
        HdGL = design(dGL);   
        gammaLo = filtfilt(HdGL.Numerator,1,rawData); % zero-phase filtering

        ContData.shank(j).contData.beta     = beta;
        ContData.shank(j).contData.gammaLo  = gammaLo;
        ContData.shank(j).contData.gammaHi  = gammaHi;
        
        % Get estimates of instantaneous power by using a median filter of squared amplitudes
        % uses a window of 80 points (3-4 cycles of beta frequency), or 0.25 sec at 500 Hz sampling rate
        ContData.shank(j).contData.betaP     = medfilt1(beta.^2,80);
        ContData.shank(j).contData.gammaLoP  = medfilt1(gammaLo.^2,40);
        ContData.shank(j).contData.gammaHiP  = medfilt1(gammaHi.^2,25);

        % Apply chronux multitaper method to extract power spectrum in time
        disp('Calculating the spectrogram of the lowBand data...');
% 
%         ContData.lfp.params      = paramsC;
%         ContData.lfp.movingwin   = movingwin;
%         [S,t,f] = mtspecgramc(lowBandData.values,movingwin,paramsC);
% 
%         dLo = fdesign.bandpass('n,fc1,fc2', 300, 2, 115, paramsC.Fs);
%         HdLo = design(dLo);   
%         lowBandDataTMP = filtfilt(HdLo.Numerator,1,rawData); % zero-phase filtering
%         [S,f] = mtspectrumc( lowBandDataTMP, paramsC );
%         
%         ContData.shank(j).contData.t = t;
%         ContData.shank(j).contData.f = f;
%         ContData.shank(j).contData.S = S;
% 
%         ContData.shank(j).contData.beta = betaRaw;
%         
%         ContData.bands(1).name = '2to5';
%         ContData.bands(2).name = '5to12';
%         ContData.bands(3).name = '14to34';
%         ContData.bands(4).name = '40to55';
%         ContData.bands(5).name = '65to115';
% 
%         ContData.bands(1).inds = find(f<2,1,'last')  : find(f>5,1,'first');
%         ContData.bands(2).inds = find(f<5,1,'last')  : find(f>12,1,'first');
%         ContData.bands(3).inds = find(f<14,1,'last') : find(f>34,1,'first');
%         ContData.bands(4).inds = find(f<40,1,'last') : find(f>55,1,'first');
%         ContData.bands(5).inds = find(f<65,1,'last') : find(f>115,1,'first');
% 
% 
%         for k=1:5
%             ContData.shank(j).bands(k).values = mean(ContData.shank(j).contData.S(:, ContData.bands(k).inds)');
%         end

%         ContData.shank(j).bands(1).values = beta.^2;
%         ContData.shank(j).bands(2).values = gammaLo.^2;
%         ContData.shank(j).bands(3).values = gammaHi.^2;

        telapsed = toc(tstart);
        disp(['Filtered and extracted power spec from one channel in ' num2str(telapsed) ' seconds.']);

%         figure(1);
%         subplot(numel(chanList),1,j);
%         semilogy(f,S);

%         hold off;
%         plot(f,mean(S)); hold on;
%         peak = max(mean(S));
%         plot([2 2],[0 peak],'r--');
%         plot([5 5],[0 peak],'r--');
%         plot([14 14],[0 peak],'r--');
%         plot([34 34],[0 peak],'r--');
%         plot([40 40],[0 peak],'r--');
%         plot([55 55],[0 peak],'r--');
%         plot([65 65],[0 peak],'r--');
%         plot([130 130],[0 peak],'r--');
%         xlabel('Frequency (Hz)');
%         ylabel('Mean');
%         title('Frequency bands selected for analysis');
        
        figure(2);
        subplot(numel(chanList),1,j); hold off;
        plot(ContData.shank(j).contData.beta(10000:70000)); hold on;
        plot(ContData.shank(j).contData.gammaHi(10000:70000),'k');
        drawnow;

    end

%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % LOAD CONTINUOUS LEVER DATA
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------

    filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'ns4']
    Ns4DATA = openNSx('report','read',filenamestrE);
    clear leverData sLeverData tmpLeverData 

    leverData(1,:) = decimate(Ns4DATA.Data(1,:),10);
    leverData(2,:) = decimate(Ns4DATA.Data(2,:),10);

    sLeverData(1,:) = sgolayfilt(leverData(1,:),9,101);
    sLeverData(2,:) = sgolayfilt(leverData(2,:),9,101);

    ContData.behavior.sLeverData = sLeverData;

    if strcmp(filenamestr(1),'M')
        figure(55); clf;
    else
        figure(5); clf;
    end
    plotArray = 1:250000;
    plot3(sLeverData(1,plotArray),sLeverData(2,plotArray),plotArray,'k');

    rawLick = decimate(Ns4DATA.Data(3,:),10);
    difLick = diff(rawLick);
    evLick=zeros(1,numel(rawLick));

    ContData.behavior.evLick    = evLick;
    ContData.behavior.rawLick   = rawLick;

%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % EXTRACT VELOCITY DATA FROM CONT LEVER DATA
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------

    numSamples = size(ContData.behavior.sLeverData,2);
    tmpLeverData(1,:) = sgolayfilt(ContData.behavior.sLeverData(1,:),3,201);
    tmpLeverData(2,:) = sgolayfilt(ContData.behavior.sLeverData(2,:),3,201);

% Above is a more accurate version of a ~10Hz cutoff filter (see below)
% that is standard in the field.
%     h=fdesign.lowpass('Fp,Fst,Ap,Ast',0.01,0.015,1,60);
%     d=design(h,'equiripple'); %Lowpass FIR filter
%     tmpLeverDataFX=filtfilt(d.Numerator,1,ContData.behavior.sLeverData(1,:)); %zero-phase filtering
%     tmpLeverDataFY=filtfilt(d.Numerator,1,ContData.behavior.sLeverData(2,:)); %zero-phase filtering

    sLeverV = zeros(1,numSamples);

    disp(' ');disp(' ');disp('Extracting velocity...');

    dX = diff(tmpLeverData(1,:));
    dY = diff(tmpLeverData(2,:));

    sLeverV = sqrt( dX.^2 + dY.^2 );

    disp(' ');disp(' Complete. ');disp(' ');

    ContData.behavior.sLeverV = sLeverV;
    ContData.behavior.sLeverVm = sgolayfilt(ContData.behavior.sLeverV,3,501);
    clear sLeverV;

%     figure(2); clf;
%     subplot(211);
%     plot(ContData.behavior.sLeverData(1,:),'k'); hold on;
%     subplot(212);
%     plot(ContData.behavior.sLeverV,'k'); hold on;
%     plot(ContData.behavior.sLeverVm,'r');

%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % FIND MOVEMENTS
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    method = 'vel';
    numC = 5;
    clear progSt* reach

    currX       = ContData.behavior.sLeverData(1,:);
    currY       = ContData.behavior.sLeverData(2,:);
    currV       = ContData.behavior.sLeverVm; 
    

    switch method

        case 'vel'

            pre     = 10;
            post    = 10;       
            minSpace = 250;
            count = 1;
            
            % threshold the velocities
            allValidSamps   = find(currV(1:numel(currV)-(minSpace.*2))>2);

%             progStartTMP(count,1)  = currStamps(allValidSamps(1));
            progStartTMP(count,2)  = currX(allValidSamps(1));
            progStartTMP(count,3)  = currY(allValidSamps(1));
            progStartTMP(count,4)  = allValidSamps(1);

            for j=2:numel(allValidSamps)
                
                if allValidSamps(j)>allValidSamps(j-1)+minSpace
                    
                    postN = find(currV(allValidSamps(j-1):allValidSamps(j))<1,1,'first');
                    if numel(postN)<1
                        figure(3); clf;
                        plot(progStartTMP(count,4):allValidSamps(j),currV(progStartTMP(count,4):allValidSamps(j)),'k');
                        hold on; plot(progStartTMP(count,4),currV(progStartTMP(count,4)),'bo',allValidSamps(j-1),currV(allValidSamps(j-1)),'ro');
                        pause(0.1);
                        
                        disp('Cannot find stop');                        
                        postN=post;
                    end
%                     progStopTMP(count,1)   = currStamps(allValidSamps(j-1)+post);
                    progStopTMP(count,2)   = currX(allValidSamps(j-1)+postN);
                    progStopTMP(count,3)   = currY(allValidSamps(j-1)+postN);
                    progStopTMP(count,4)   = allValidSamps(j-1)+postN;

                    count                  = count+1;

                    preN = find(currV(allValidSamps(j-1) : allValidSamps(j))<1,1,'last');
                    if numel(preN)<1
                        disp('Cannot find start');
%                     progStartTMP(count,1)  = currStamps(allValidSamps(j)-pre);
                        progStartTMP(count,2)  = currX(allValidSamps(j)-pre);
                        progStartTMP(count,3)  = currY(allValidSamps(j)-pre);
                        progStartTMP(count,4)  = allValidSamps(j)-pre;                    
                    else
%                     progStartTMP(count,1)  = currStamps(allValidSamps(j-1)+preN);
                        progStartTMP(count,2)  = currX(allValidSamps(j-1)+preN);
                        progStartTMP(count,3)  = currY(allValidSamps(j-1)+preN);
                        progStartTMP(count,4)  = allValidSamps(j-1)+preN;
                    end
                end

                if j==numel(allValidSamps)
%                     post = find(currV(allValidSamps(j):allValidSamps(j)+minSpace)<0.5,1,'first');
%                     progStopTMP(count,1)   = currStamps(allValidSamps(j)+post);
                    progStopTMP(count,2)   = currX(allValidSamps(j)+post);
                    progStopTMP(count,3)   = currY(allValidSamps(j)+post);
                    progStopTMP(count,4)   = allValidSamps(j)+post;
                end

            end

            count = 1;
            for k = 1:size(progStartTMP,1)
                
                if k==1
                    reach.init = 1;
                end
                
                % reaches must be at least 50 ms long
                if progStopTMP(k,4)-progStartTMP(k,4)>=90 & progStartTMP(k,4)>minSpace

                    trajAngle   = atan2(progStopTMP(k,3)-progStartTMP(k,3),progStopTMP(k,2)-progStartTMP(k,2));

                    if (pdist2([progStopTMP(k,2),progStopTMP(k,3)],[mean(currX),mean(currY)]) > pdist2([progStartTMP(k,2),progStartTMP(k,3)],[mean(currX),mean(currY)]))
                        reach.out(count) = 1;
                    else
                        reach.out(count) = 0;
                    end
                    velTraj = ContData.behavior.sLeverV(progStartTMP(k,4) : progStopTMP(k,4));
                    xVals = ContData.behavior.sLeverData(1,progStartTMP(k,4) : progStopTMP(k,4));
                    yVals = ContData.behavior.sLeverData(2,progStartTMP(k,4) : progStopTMP(k,4));
                    
                    reach.start(count,:)  = progStartTMP(k,:);
                    reach.stop(count,:)   = progStopTMP(k,:);
                    reach.angle(count,1)  = trajAngle;
                    reach.dist(count,:)   = trapz(velTraj);
                    
                    tmp = findpeaks(velTraj);
                    reach.numpks(count,1) = numel(tmp.loc);
                    reach.dur(count,1)    = progStopTMP(k,4) - progStartTMP(k,4);
                    reach.vel(count,1)   = max(velTraj);
                    reach.vel(count,2)   = trapz(velTraj) ./ reach.dur(count,1);
                    reach.vel(count,3)   = var(velTraj);

                    reach.acc(count,1)   = max(diff(velTraj));
                    reach.acc(count,2)   = mean(diff(velTraj));
                    reach.acc(count,3)   = max(diff(velTraj(1:90))); % max in first 90 ms of movement
                    
                    reach.tort(count,1)  = reach.dist(count,:) ./ pdist2([progStopTMP(k,2),progStopTMP(k,3)],[progStartTMP(k,2),progStartTMP(k,3)]);
                    
                    count                 = count+1;

                end            
            end


    end

    level = ones(1,size(progStartTMP,1));
    figure(3); clf;
    plot(ContData.behavior.sLeverV,'k'); hold on;
    plot(progStartTMP(:,4),level,'r^');
    plot(progStopTMP(:,4),level,'bo');

%     [vals,inds] = sort(reach.dist+(reach.maxV(:,3).*1e3),'ascend');        
%     [vals,inds] = sort(reach.maxV(:,3),'ascend');    
    [vals,inds] = sort(reach.dur,'ascend');    

    numReaches = size(reach.start,1)
%     if strcmp(filenamestr(1),'M')
%         figure(55); clf;
%     else
        figure(5); clf;
        figure(7); clf;
%     end

    dims = floor(sqrt(numReaches));
    scaler = 3000;
    
    for p = 1:1:dims.^2
        
        l = inds(p);
        
    %   offset = reach.start(l,4).*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
    %   plot3(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4)), currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4)) , offset , 'r','LineWidth',1);
        
        offset = l.*ones(1,numel(reach.start(l,4):reach.stop(l,4)));
        [x,y] = ind2sub(dims,p);
        
        figure(5);
        plot((x.*scaler),(y.*scaler),'k.','MarkerSize',8);
        hold on;
        plot(currX(reach.start(l,4):reach.stop(l,4))-currX(reach.start(l,4))+(x.*scaler),currY(reach.start(l,4):reach.stop(l,4))-currY(reach.start(l,4))+(y.*scaler) , 'k','LineWidth',1,'Color',[p/numReaches 0.67-(p/numReaches)*0.67 1-(p/numReaches)]);
        if i==dims.^2
            axis tight; axis off;
        end
                
        figure(7);
        startInd = round(reach.start(l,4)./2);
        stopInd = round(reach.stop(l,4)./2);
        
        meanField = sqrt(ContData.shank(2).contData.betaP(startInd-(minSpace/2):startInd+600));
        for d=3:5
            meanField = meanField + sqrt(ContData.shank(d).contData.betaP(startInd-(minSpace/2):startInd+600));
        end
        meanField = meanField ./ 4;
            
        plot((x.*scaler)+([startInd-(minSpace/2):startInd+600]-startInd).*2 , meanField + (y.*scaler./20),'k','LineWidth',1,'Color',[p/numReaches 0.67-(p/numReaches)*0.67 1-(p/numReaches)]);
        hold on;
%         plot((x.*scaler)+([startInd-(minSpace/2):startInd+600]-startInd).*2 , ContData.shank(4).contData.beta(startInd-(minSpace/2):startInd+600) + (y.*scaler./20),'k--','LineWidth',1,'Color',[p/numReaches 0.67-(p/numReaches)*0.67 1-(p/numReaches)]);
        plot((x.*scaler)+[reach.start(l,4)-minSpace:reach.start(l,4)+1200]-reach.start(l,4) , ContData.behavior.sLeverV(reach.start(l,4)-minSpace:reach.start(l,4)+1200).*5 + (y.*scaler./20),'k','LineWidth',2);
        plot((x.*scaler),(y.*scaler./20),'k.','MarkerSize',8);
        if i==dims.^2
            axis tight; axis off;
        end
    end
    
    figure(6); clf;
    subplot(311);
    p = polyfit(reach.dist,reach.vel(:,1),1);
    plot(reach.dist,reach.vel(:,1),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Peak Velocity');
    
    subplot(312);
    p = polyfit(reach.dist,reach.dur,2);
    plot(reach.dist,reach.dur,'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Movement Duration (ms)');
    
    subplot(313);
    p = polyfit(reach.dist,reach.acc(:,3),1);
    plot(reach.dist,reach.acc(:,3),'k.',1:max(reach.dist),polyval(p,1:max(reach.dist)),'r-');
    xlabel('Total distance');
    ylabel('Initial acceleration');

%     plot3([0 0],[0 0],[-5 numReaches+5],'k-','LineWidth',2);
%     grid on; box on;
%     xlabel('X (a.u.)');
%     ylabel('Y (a.u.)');
%     zlabel('Time (samples)');
%     title(['Number of reaches: ' num2str(numReaches)]);
%     view([-45 10]);
%     axis tight;

    figure(8); clf;
    rose(reach.angle);
    
    reach.numReaches = numReaches;
    ContData.behavior.reach = reach;
    
%     if numReaches>0
%         [winReach] = TNC_ReachVigorWindow(reach.vel(:,1),reach.numReaches,9);
%         ContData.behaviora.winReach = winReach;
%     end
    
    drawnow;
    
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    % COMPLETED ALL ANALYSIS. SAVE AND START OVER
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
    disp(' ');    
    disp('%-------------------------------------------------------------------');
    disp(['Completed file: ' filenamestr(1,1:length(filenamestr)-3)]);

    ContData
    
    % save the data from this file
    disp(['save ' filenamestr(1,1:length(filenamestr)-3) 'mat ContData']);
    eval(['save ' filenamestr(1,1:length(filenamestr)-3) 'mat ContData']);
    
    disp('%-------------------------------------------------------------------');
    disp(' ');
    
    % clear memory and start again
    clear ContData prog* sLev* lowB* Ns* dat* reach
    
end