function [] = TNC_POP_AlignedRaster(PopData,sessNum,alignToTimes,trialArray,unitArray,alignWindow,markTimes)

%% EXAMPLE SESSION

%     [PopData] = TNC_ConvertTSDtoPopData('',1)
%     TNC_POP_AllCrossCorr(PopData,1,[],15,2)
%     unitArray=[]
%     sessNum=1
%     trialArray = 1:200;
%     alignToTimes = ContData.behavior.reach.start(:,4);
%     load m38-behaviorData_bh
%     alignToTimes = ContData.behavior.reach.start(:,4);
%     alignWindow = [600,1200]
%     [vals,inds] = sort(ContData.behavior.reach.tort,'ascend');
%     trialArray = inds;

%% STANDARD ANALYSIS PARAMETERS

    % Define the smoothing to apply to the timeseries data
    currParams.smthParams.rise     = 1;
    currParams.smthParams.decay    = 50;
    currParams.filter.causal       = 0;

    [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
    if currParams.filter.causal
        currParams.filter.kernel(1:currParams.smthParams.decay.*15) = 0;
    end

    PopData.currParams = currParams;
    
    [mapName] = TNC_CreateRBColormap(1024,'bo');
    [eMap] = TNC_CreateRBColormap(8,'bo');

    disp(' ');
    disp(' ');
    disp('________________________________');
    disp(' ');
    disp('Initialized analysis parameters.');
    disp('________________________________');
    disp(' ');
    disp(' ');

%% COMPUTE the AlignRasters for alignToTimes and over the window

       figure(sessNum); clf;
       
    if numel(unitArray)>1
        noUnitList = 1;
        numUnits = numel(unitArray);
    else
        noUnitList = 0;
        numUnits = numel(PopData.session(sessNum).unit);
        unitArray = 1:numUnits;
    end
    
    for index=1:numUnits
        
%         currUnit = unitArray(index)
        currUnit = index;
        
        numStamps   = length(PopData.session(sessNum).unit(currUnit).ts);
        delta       = zeros(1,ceil(PopData.session(sessNum).unit(currUnit).ts(numStamps)));
        delta(1,ceil(PopData.session(sessNum).unit(currUnit).ts)) = 1;
        tmpSmooth   = conv(delta,currParams.filter.kernel,'same');

        % Only analyze events for which the unit was well isolated
        validEventsForThisCell = find(alignToTimes>PopData.session(sessNum).unit(currUnit).ts(1) & alignToTimes<PopData.session(sessNum).unit(currUnit).ts(numStamps));
        stimTimes   = double(alignToTimes(validEventsForThisCell));

        % for reference
        % [response] = TNC_AlignRasters(delta,spkStamps,stableTrials,alignStamps,window,rasterFlag,boxcar)
        [alignedData] = TNC_AlignRasters(delta , double(PopData.session(sessNum).unit(currUnit).ts) , -1 , stimTimes , alignWindow, 1, 1);
        PopData.session(sessNum).unit(currUnit).response.raster      = alignedData.raster;
        PopData.session(sessNum).unit(currUnit).response.boxcar      = alignedData.image.boxcar;

        [alignedDataS] = TNC_AlignRasters(tmpSmooth , double(PopData.session(sessNum).unit(currUnit).ts) , -1 , stimTimes , alignWindow, 0, 1);
        PopData.session(sessNum).unit(currUnit).response.psthAVG    = alignedDataS.image.psthAVG;
        PopData.session(sessNum).unit(currUnit).response.psthSEM    = alignedDataS.image.psthSEM;
        PopData.session(sessNum).unit(currUnit).response.psthZ      = alignedDataS.image.psthZ;
        PopData.session(sessNum).unit(currUnit).response.psthZe     = alignedDataS.image.psthZe;
        PopData.session(sessNum).unit(currUnit).response.trialwise  = alignedDataS.image.aligned;
        
        halfMax = max(PopData.session(sessNum).unit(currUnit).response.psthZ(alignWindow(1):alignWindow(1)+alignWindow(2))) ./ 2;
%         try
            if halfMax>0.33
                PopData.session(sessNum).uLatencies(index)  = find(PopData.session(sessNum).unit(currUnit).response.psthZ(alignWindow(1):alignWindow(1)+alignWindow(2)) > halfMax,1,'first');
            else
                PopData.session(sessNum).uLatencies(index)  = 0;                
            end
%         catch ME
%             PopData.session(sessNum).uLatencies(index)  = 0;
%         end
        
        try
            PopData.session(sessNum).duration(index)    = find(PopData.session(sessNum).unit(currUnit).response.psthZ(PopData.session(sessNum).uLatencies(index)+alignWindow(1):alignWindow(2))>halfMax,1,'last')+PopData.session(sessNum).uLatencies(index);
        catch ME
            PopData.session(sessNum).duration(index)    = 0;
        end

        PopData.session(sessNum).respAmps(index) = trapz(PopData.session(sessNum).unit(currUnit).response.psthZ(alignWindow(1):alignWindow(1)+600));
        
        figure(3);
        subplot(3,ceil(numUnits./3),index);     
%         subplot(8,1,PopData.session(sessNum).unit(currUnit).sh);     
%         patch([-alignWindow(1):alignWindow(2) ; -alignWindow(1):alignWindow(2)] , [ PopData.session(sessNum).unit(currUnit).response.psthZ + PopData.session(sessNum).unit(currUnit).response.psthZe ; PopData.session(sessNum).unit(currUnit).response.psthZ - PopData.session(sessNum).unit(currUnit).response.psthZe ], [index./numUnits 0.5 1-(index./numUnits)]); hold on;
        shadedErrorBar([-alignWindow(1):alignWindow(2)],PopData.session(sessNum).unit(currUnit).response.psthZ,PopData.session(sessNum).unit(currUnit).response.psthZe,{'k'}); hold on;

%         plot(-alignWindow(1):alignWindow(2), PopData.session(sessNum).unit(currUnit).response.psthZ ./ max(PopData.session(sessNum).unit(currUnit).response.psthZ), 'LineWidth', 1, 'Color', [index./numUnits 0.5 1-(index./numUnits)]); hold on;
%         plot([PopData.session(sessNum).uLatencies(index), PopData.session(sessNum).duration(index) + PopData.session(sessNum).uLatencies(index)],[0.5 0.5], 'o', 'LineWidth', 2, 'Color', [index./numUnits 0.5 1-(index./numUnits)]);

%         if index==numUnits
            plot([0 0],[-1 2],'k--');
            plot([-alignWindow(1) alignWindow(2)],[0 0],'k--');
%         end
%         colormap(mapName); 
        axis([-alignWindow(1) alignWindow(2) -1 2]);

    end
    
%% ORDER CELLS BY RESPONSE LATENCY
    
    [vals,trialArray] = sort(markTimes,'ascend');

    figure(600);
    subplot(311);
    uLatenciesDist.x = 0:1:120;
    uLatenciesDist.y = hist(PopData.session(sessNum).uLatencies,uLatenciesDist.x);
    bar(uLatenciesDist.x,uLatenciesDist.y);    
    
    figure(600);
    subplot(312);
    uDurationDist.x = 0:1:120;
    uDurationDist.y = hist(PopData.session(sessNum).duration,uDurationDist.x);
    bar(uDurationDist.x,uDurationDist.y);
    
    subplot(313);
    plot(PopData.session(sessNum).uLatencies,PopData.session(sessNum).respAmps,'k.');
    
    [values,inds] = sort(PopData.session(sessNum).uLatencies + PopData.session(sessNum).respAmps,'ascend');    
unitArray=inds;
%     unitArray = daUnits(inds);
    
%% DISPLAY RASTERS SORTED BY trialArray and unitArray

%     [vals,unitArray] = sort(PopData.session(sessNum).respAmps,'ascend');
figOffset = 0;
yspc=1;
    figure(5+figOffset); clf;
    figure(6+figOffset); clf;
    kludge = 0;
    numTrials = numel(trialArray);
    numUnits = numel(unitArray);

    disp(' ');
    disp(' ');
    disp('________________________________');
    disp(' ');
    disp('Generating plot.');
    disp('________________________________');
    disp(' ');
    disp(' ');

    [colorMap] = TNC_CreateRBColormap(1024,'bo');
    
% ['rgb(166,206,227)','rgb(31,120,180)','rgb(178,223,138)','rgb(51,160,44)','rgb(251,154,153)','rgb(227,26,28)','rgb(253,191,111)','rgb(255,127,0)','rgb(202,178,214)','rgb(106,61,154)']

shankColors = [ 166,206,227;
                 31,120,180;
                178,223,138;
                 51,160, 44;
                251,154,153;
                227, 26, 28;
                253,191,111;
                255,127,  0;
                202,178,214;
                106, 61,154;
                255,255,153;
                177, 89, 40] ./ 256;
            
    
    for index = 1:numUnits

        currUnit = unitArray(index);
        figure(5+figOffset);
        if rem(index,2)==0
            thisColor = colorMap(round(1024-(index./numUnits.*512)),:);            
        else
            thisColor = colorMap(round(513-(index./numUnits.*512)),:);            
        end
        
        allRows     = [];
        allStamps   = [];
        markStamps  = [];
        markRows    = [];
                                
        for jindex=1:numTrials

            currTrial = trialArray(jindex);
            
            figure(5+figOffset);
            rowOffset   = ((index-1).*numTrials.*yspc) + 10 + jindex;
            theseStamps = PopData.session(sessNum).unit(currUnit).response.raster.trial(currTrial).ts;
            theseRows   = ones(1,numel(theseStamps)) .* rowOffset;            
            theseMarks  = markTimes(currTrial);
            theseMarkRows= ones(1,numel(theseMarks)) .* rowOffset; 
            
            
            allStamps   = [allStamps theseStamps];
            allRows     = [allRows theseRows];
            markStamps  = [markStamps theseMarks];
            markRows    = [markRows theseMarkRows];
    
        end
        
        
        if kludge
            for p=1:numel(allStamps)
                text(allStamps(p),allRows(p),'|','Color',thisColor,'FontName','Times','Fontsize',4); 
                hold on;
            end
        else
                plot(allStamps,allRows,'.','Color',thisColor,'MarkerSize',6); hold on;
        end
        
%         plot(markStamps,markRows,'o','Color',[0.1 0.1 0.1],'MarkerSize',6);
        
%         plot(-alignWindow(1):alignWindow(2), (PopData.session(sessNum).unit(currUnit).response.psthZ .* (numTrials ./ 4)) + ((index-1).*numTrials) + 10 + (numTrials./4),'-','Color',[0 0 0],'LineWidth',1); hold on;
%         text(-alignWindow(1)-500, ((index-1).*numTrials) + 10 + (numTrials./4), ['e' num2str(PopData.session(sessNum).unit(currUnit).el) 'u' num2str(PopData.session(sessNum).unit(currUnit).ui)],'FontName','Arial','FontSize',12,'FontWeight','bold'); hold on;


        figure(6+figOffset);
%         patch([-alignWindow(1):alignWindow(2) ; -alignWindow(1):alignWindow(2)] , [ PopData.session(sessNum).unit(currUnit).response.psthZ + PopData.session(sessNum).unit(currUnit).response.psthZe ; PopData.session(sessNum).unit(currUnit).response.psthZ - PopData.session(sessNum).unit(currUnit).response.psthZe ] + (index), [1 0 0]); hold on;
%         plot(-alignWindow(1):alignWindow(2) , PopData.session(sessNum).unit(currUnit).response.psthZ + (index), 'Color' , shankColors(PopData.session(sessNum).unit(currUnit).sh,:) , 'LineWidth' , 3); hold on;
        
        shadedErrorBar([-alignWindow(1):alignWindow(2)],PopData.session(sessNum).unit(currUnit).response.psthZ + (index*yspc),PopData.session(sessNum).unit(currUnit).response.psthZe,{'k' 'color' eMap(PopData.session(sessNum).unit(currUnit).sh,:)}); hold on;
        
        text(-100,(index*yspc)+0.5,num2str(PopData.session(sessNum).unit(currUnit).sh));
        
        axis([-alignWindow(1) alignWindow(2) 0 (numUnits+2)*yspc]);
        set(gcf,'Color','white');
        set(gca,'TickDir','out','TickLength',[0.005 0]);
        set(gca,'YTickLabel',''); box off;
        if index==1
            plot([0 0] , [0 (numUnits+2)*yspc] , 'k-' ); 
            plot([-alignWindow(1) alignWindow(2)] , [index*yspc index*yspc] , 'k--' );             
        else
            plot([-alignWindow(1) alignWindow(2)] , [index*yspc index*yspc] , 'k--' );             
        end

    end
    
    
    figure(5+figOffset);
    plot([0 0] , [0 rowOffset] , 'k--' ); 
    axis tight;     
    set(gcf,'Color','white');
    set(gca,'TickDir','out','TickLength',[0.005 0]);
    set(gca,'YTickLabel','');
    box off;
    axis([-alignWindow(1) alignWindow(2) 0 numTrials.*numUnits.*yspc]);
%     plot2svg('allDaraster.svg');
xlabel('Time (ms)');
    
%% DISPLAY rasters averaged over data subsets
    
    firstThird = trialArray(1:round(numel(trialArray)./3));
    secondThird = trialArray(ceil(numel(trialArray)./3):round(numel(trialArray)./3)+round(numel(trialArray)./3));
    lastThird = trialArray(ceil(numel(trialArray)./3)+round(numel(trialArray)./3):numel(trialArray));
        
    figure(50+figOffset); clf;
    
    for index=1:numUnits
        thirds(1).avg = mean(PopData.session(sessNum).unit(index).response.trialwise(firstThird,:));
        thirds(1).err = std(PopData.session(sessNum).unit(index).response.trialwise(firstThird,:),[],1) ./ sqrt(numel(firstThird)-1);

        thirds(2).avg = mean(PopData.session(sessNum).unit(index).response.trialwise(secondThird,:));
        thirds(2).err = std(PopData.session(sessNum).unit(index).response.trialwise(secondThird,:),[],1) ./ sqrt(numel(secondThird)-1);

        thirds(3).avg = mean(PopData.session(sessNum).unit(index).response.trialwise(lastThird,:));
        thirds(3).err = std(PopData.session(sessNum).unit(index).response.trialwise(lastThird,:),[],1) ./ sqrt(numel(lastThird)-1);

        subplot(3,ceil(numUnits./3),index); 
%         shadedErrorBar([-alignWindow(1):alignWindow(2)],thirds(1).avg,thirds(1).err,{'k'}); hold on;
%         shadedErrorBar([-alignWindow(1):alignWindow(2)],thirds(3).avg,thirds(3).err,{'r'}); hold on;
        plot([-alignWindow(1):alignWindow(2)],thirds(1).avg,'k'); hold on;
        plot([-alignWindow(1):alignWindow(2)],thirds(2).avg,'b'); hold on;
        plot([-alignWindow(1):alignWindow(2)],thirds(3).avg,'r'); hold on;
        axis([-alignWindow(1) alignWindow(2) -0.0025 0.1]);
    end
    
%% PROJECT onto mapped dimensions to look for differences in the population state according to reward rate deviations

    dispRange = 1000:3000;
    filtWidth = 21;
    
    [mappedA, mapping_t] = compute_mapping(analysis.approach.left.pAvg','PCA',3);

    for m = 1:numel([-alignWindow(1):alignWindow(2)])
        
        for index=1:numUnits
            thirds(1).popVector(m,index) = mean(PopData.session(sessNum).unit(index).response.trialwise(firstThird,m));
            thirds(2).popVector(m,index) = mean(PopData.session(sessNum).unit(index).response.trialwise(secondThird,m));
            thirds(3).popVector(m,index) = mean(PopData.session(sessNum).unit(index).response.trialwise(lastThird,m));        
        end
        
        for n=1:3
            thirds(1).proj(n).values(m) = dot(thirds(1).popVector(m,:)',mapping_t.M(:,n));
            thirds(2).proj(n).values(m) = dot(thirds(2).popVector(m,:)',mapping_t.M(:,n));
            thirds(3).proj(n).values(m) = dot(thirds(3).popVector(m,:)',mapping_t.M(:,n));
        end    
        
    end
    
    figure(51+figOffset); clf; subplot(4,1,1:3);
    plot3(sgolayfilt(thirds(1).proj(1).values(dispRange), 3, filtWidth),sgolayfilt(thirds(1).proj(2).values(dispRange), 3, filtWidth),sgolayfilt(thirds(1).proj(3).values(dispRange), 3, filtWidth),'k-', ...
          sgolayfilt(thirds(2).proj(1).values(dispRange), 3, filtWidth),sgolayfilt(thirds(2).proj(2).values(dispRange), 3, filtWidth),sgolayfilt(thirds(2).proj(3).values(dispRange), 3, filtWidth),'b-', ...
          sgolayfilt(thirds(3).proj(1).values(dispRange), 3, filtWidth),sgolayfilt(thirds(3).proj(2).values(dispRange), 3, filtWidth),sgolayfilt(thirds(3).proj(3).values(dispRange), 3, filtWidth),'r-','LineWidth',1); hold on;
          
    plot3(thirds(1).proj(1).values(2000),thirds(1).proj(2).values(2000),thirds(1).proj(3).values(2000),'ks', ...
          thirds(2).proj(1).values(2000),thirds(2).proj(2).values(2000),thirds(2).proj(3).values(2000),'bs',...
          thirds(3).proj(1).values(2000),thirds(3).proj(2).values(2000),thirds(3).proj(3).values(2000),'rs',...
          'MarkerSize',15, 'LineWidth', 3);
      xlabel('PC 1');      ylabel('PC 2');      zlabel('PC 3');

    figure(51+figOffset); subplot(4,1,4);
        xsq = (thirds(1).proj(1).values-thirds(3).proj(1).values).^2;
        ysq = (thirds(1).proj(2).values-thirds(3).proj(2).values).^2;
        zsq = (thirds(1).proj(3).values-thirds(3).proj(3).values).^2;

        xsq2 = (thirds(2).proj(1).values-thirds(3).proj(1).values).^2;
        ysq2 = (thirds(2).proj(2).values-thirds(3).proj(2).values).^2;
        zsq2 = (thirds(2).proj(3).values-thirds(3).proj(3).values).^2;
        
        plot([-alignWindow(1):alignWindow(2)],sgolayfilt(sqrt(xsq+ysq+zsq)   , 3, filtWidth), 'k', ...
             [-alignWindow(1):alignWindow(2)],sgolayfilt(sqrt(xsq2+ysq2+zsq2), 3, filtWidth), 'b' );
         
             xlabel('Distance');      ylabel('Time from reach (ms)');
    
    
%% TEMP SANDBOX
clear firstReach*
% Find the first reach after reward
    
    for j = 2:numel(ContData.behavior.rewardInds)-5
        
        currReward = ContData.behavior.rewardInds(j);        
        nextReachInd = find(ContData.behavior.reach.start(:,4) > currReward+1000 , 1, 'first');
        firstReachInds(j)   = ContData.behavior.reach.start(nextReachInd,4);
        firstReachLats(j)   = ContData.behavior.reach.start(nextReachInd,4) - currReward;
        firstReachVel(j)    = ContData.behavior.reach.vel(nextReachInd,1);
        firstReachDur(j)    = ContData.behavior.reach.dist(nextReachInd);
        
        timeSinceReward(j)  = currReward - ContData.behavior.rewardInds(j-1);
    
    end
    