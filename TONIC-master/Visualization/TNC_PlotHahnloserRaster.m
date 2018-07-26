function [] = TNC_PlotHahnloserRaster(unitStruct,window,alignTs,freeTs,sortScalar,numTrials)

%% Prep the figure
    disp(['Found ' num2str(numTrials) ' valid trials']);

    h1=figure(40); clf; 
    set(h1,'Color',[1 1 1]);
    h2=figure(41); clf; 
    set(h2,'Color',[1 1 1]);

%% Sort based on 'sortScalar'

    [vals,trialInds] = sort(sortScalar);
    
%% Plot the behavioral data

    disp('Plotting routine for behavior in time...');
    figure(40);
    subplot(5,1,1); hold on;
    for k = 1:numTrials
        thisTrialRow    = trialInds(k);
        theseTimeStamps = PopData.session(i).behavior.raster.trial(thisTrialRow).ts;
        validStamps = 1:numel(PopData.session(i).behavior.raster.trial(thisTrialRow).ts);
        plot(theseTimeStamps(validStamps),ones(1,numel(validStamps)).*k,'k.','MarkerSize',1);%,'.','Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
    end

%% Plot all trials for each cell

    numUnits = size(PopData.session(i).unit,2);

    figure(40); hold on;

    for m=1:numUnits

        delta = zeros(1,7300);

        thisUnit = m;

        for k = 1:numTrials

            thisTrialRow    = ((m-1).*numTrials) + (5.*m) + k;
            theseTimeStamps = PopData.session(i).unit(thisUnit).raster.trial().ts;

           if numel(theseTimeStamps)>0
               tmp = find(theseTimeStamps<6300);

               delta(ceil(theseTimeStamps(tmp)+1000)) = delta(ceil(theseTimeStamps(tmp)+1000)) + 1;

               if rem(m,2)
                    subplot(5,1,2:5); hold on;
                    plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[m./numUnits 0.67-(0.67.*(m./numUnits)) 1-(m./numUnits)]);
               else
                    subplot(5,1,2:5);  hold on;
                    plot(theseTimeStamps,ones(1,numel(theseTimeStamps)).*thisTrialRow,'.','MarkerSize',1,'Color',[1-(m./numUnits) 0.67.*(m./numUnits) m./numUnits]);            
               end

               plot(zeros(1,numTrials),((m-1).*numTrials) + (5.*m) + [1:numTrials],'b');
               plot(PopData.session(i).behavior.flCSon(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');
               plot(PopData.session(i).behavior.USon(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'r');
               plot(PopData.session(i).behavior.flCSoff(trialInds(1:numTrials)),((m-1).*numTrials) + (5.*m) + [1:numTrials],'k');
           end

        end

        PopData.session(i).unit(thisUnit).allPhase.timePSTH.delta = delta ./ numTrials;

    end

    subplot(5,1,2:5);
    set(gca,'TickDir','out','TickLength',[0.005 0]); box off;
    axis([-1000 8500 -1 ((numUnits).*numTrials) + (5.*(numUnits+1))]);
    xlabel('Time (ms)');
    ylabel('Sorted Trials per Unit');
    drawnow;
    