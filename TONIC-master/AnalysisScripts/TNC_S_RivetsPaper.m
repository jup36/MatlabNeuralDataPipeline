% TNC_S_RivetsPaper

%% Display example units overlaid on filtered data with simultaneous video

ExemplarUnit    = 1;
BestProj        = [6 8];
ShankToLoad     = 1;
SegToDisplay    = 3.4;
figNum          = 1;
chunk           = 120;
sortBaseName    = 'm38-20110324-buz278b_'
behDataName     = 'm38-behaviorData_bh.mat'
arrayType       = 'NN_b64'
fileName        = 'm38-20110324-buz278b-005.ns5'

[row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);

%% LOAD and filter the continuous data
        
    timeStr = ['t:' num2str((SegToDisplay-1)*chunk) ':' num2str(SegToDisplay*chunk)];
    disp(' ');disp(' ');disp(' ');disp(' ');disp(' ');disp(' ');
    disp(['Loading data over the range of ' timeStr]);
    disp(' ');
    Ns5DATA = openNSx('report','read',fileName,timeStr,'sec');
    
    spkData = zeros(8,numel(Ns5DATA.Data(1,:)));
        
    % SMOOTH/FILTER DATA
    avgSeg      = mean(Ns5DATA.Data,1);
    numChan     = size(Ns5DATA.Data,1);

    for i=1:8 %electrodes

        disp(['Filtering data on channel ' num2str(electrodeMatrix(i,ShankToLoad)) ' ...']);
        rawData = sgolayfilt(Ns5DATA.Data(electrodeMatrix(i,ShankToLoad),:)-avgSeg,7,21);
        
        dHi = fdesign.bandpass('n,fc1,fc2', 1000, 500, 7000, 30000);
        HdHi = design(dHi);        
        spkData(i,:) = filtfilt(HdHi.Numerator,1,rawData);                                              
        
    end
    disp('...completed');
    
%% LOAD the relevant behavioral data for display

    behavData = load(behDataName);

    dataRange(1)    = (SegToDisplay-1)*chunk*1000; 
    dataRange(2)    = SegToDisplay*chunk*1000;
    
    reachDisp       = behavData.ContData.behavior.sLeverData(:,dataRange(1):dataRange(2));
    lickDisp        = behavData.ContData.behavior.rawLick(dataRange(1):dataRange(2));
    
    captRewards     = behavData.ContData.behavior.rewardInds(find(behavData.ContData.behavior.rewardInds>dataRange(1) & behavData.ContData.behavior.rewardInds<dataRange(2))) - dataRange(1);
    captThresh      = behavData.ContData.behavior.threshInds(find(behavData.ContData.behavior.threshInds>dataRange(1) & behavData.ContData.behavior.threshInds<dataRange(2))) - dataRange(1);
%     lickThresh      = behavData.ContData.behavior.lickInds(find(behavData.ContData.behavior.lickInds>dataRange(1) & behavData.ContData.behavior.lickInds<dataRange(2))) - dataRange(1);
    
%% LOAD the sorting data for this shank

    sortData    = load([sortBaseName 'shank' num2str(ShankToLoad) '_tsd.mat']);
    featDataId  = load([sortBaseName 'ft_tns.mat']);
    featData    = load([sortBaseName 'ft.mat']);

    for k=1:numel(sortData.shank.unit)

        validInds = find(sortData.shank.unit(k).ts>dataRange(1) & sortData.shank.unit(k).ts<dataRange(2));
        unit.sub(k).ts      = sortData.shank.unit(k).ts(validInds) - dataRange(1);
        unit.sub(k).inds    = sortData.shank.unit(k).inds(validInds) - (dataRange(1).*30);

    end

%% GET parameters of well sorted example unit

    featDataId  = load([sortBaseName 'ft_tns.mat']);
    featData    = load([sortBaseName 'ft.mat']);
    
    xCol = 23;
    yCol = 1;
    zCol = 25;
    
    for m = 1:numel(featDataId.idList.seg)
        
        validIds = find(featDataId.idList.seg(m).shank(ShankToLoad).id==ExemplarUnit);
        
        if m==1
            exemplarPlot.x  = featData.featStruct.seg(m).shank(ShankToLoad).params(validIds,xCol)';
            noisePlot.x     = featData.featStruct.seg(m).shank(ShankToLoad).params(:,xCol)';
            exemplarPlot.y  = featData.featStruct.seg(m).shank(ShankToLoad).params(validIds,yCol)';
            noisePlot.y     = featData.featStruct.seg(m).shank(ShankToLoad).params(:,yCol)';
            exemplarPlot.z  = featData.featStruct.seg(m).shank(ShankToLoad).params(validIds,zCol)';
            noisePlot.z     = featData.featStruct.seg(m).shank(ShankToLoad).params(:,zCol)';
        else
            exemplarPlot.x  = [exemplarPlot.x   featData.featStruct.seg(m).shank(ShankToLoad).params(validIds,xCol)'];
            noisePlot.x     = [noisePlot.x      featData.featStruct.seg(m).shank(ShankToLoad).params(:,xCol)'];
            exemplarPlot.y  = [exemplarPlot.y   featData.featStruct.seg(m).shank(ShankToLoad).params(validIds,yCol)'];
            noisePlot.y     = [noisePlot.y      featData.featStruct.seg(m).shank(ShankToLoad).params(:,yCol)'];
            exemplarPlot.z  = [exemplarPlot.z   featData.featStruct.seg(m).shank(ShankToLoad).params(validIds,zCol)'];
            noisePlot.z     = [noisePlot.z      featData.featStruct.seg(m).shank(ShankToLoad).params(:,zCol)'];            
        end
        
    end    
    
%% LOAD the waveform data

    ssData = load([sortBaseName 'ss.mat']);

%% GENERATE THE EXAMPLE RECORDING FIGURE

    h = figure(figNum);
    clf;
    set(h,'Name',sortBaseName(1:numel(sortBaseName)-1),'Color',[1 1 1]);
    yspc = 1000;
    win = [1 2] .* 30;
    colorSpec               = [ 92,200,200;
                               134,223,  4;
                                   204,  3, 95;
                                    3, 95,204; 
                                    95,204,  3;
                               247,254,  0;
                               220,  0, 85;
                               255,116,  0;
                               113,  9,170;
                                 0,153,153;
                               230, 60,137;
                                 2,146,146;
                               186,239,110;                               
                               249,175,114 ];
                           
    colorSpec = colorSpec ./255;    
    
    % for reference: subplot('position',[left bottom width height])
    
    % Continuous recording data with some sorted spikes overlaid
    for i=1:8
        subplot('position',[0.3 0.5 0.6 0.45]);
        plot(spkData(i,:)+(yspc.*(i-1)),'Color',[0.7 0.7 0.7]); hold on;
    end
    
    subplot('position',[0.3 0.5 0.6 0.45]);
    for k=1:3
        unit.sub(k).wf = zeros(8,numel([-win(1):win(2)]));

        for p=1:numel(unit.sub(k).ts)
            for i=1:8

                currTS = unit.sub(k).inds(p); 
                
                plot(currTS-win(1):currTS+win(2),spkData(i,currTS-win(1):currTS+win(2))+(yspc.*(i-1)),'Color',colorSpec(k+2,:)); hold on;
                
                unit.sub(k).wf(i,:) = unit.sub(k).wf(i,:) + spkData(i,currTS-win(1):currTS+win(2));
                
%                 subplot('position',[0.925 0.7 0.05 0.25]);
%                 plot(-win(1):win(2),spkData(i,currTS-win(1):currTS+win(2))+(yspc.*(i-1)),'Color',colorSpec(k+2,:)); hold on;
            end
        end
        unit.sub(k).wf = unit.sub(k).wf ./ numel(unit.sub(k).ts);
    end
    subplot('position',[0.3 0.5 0.6 0.45]);
    title(['m38-s20110324-shank' num2str(ShankToLoad)]);
    axis tight; 
    axis off;
    ylim = get(gca,'Ylim');

    % Example unit data
    for k=1:3
        subplot('position',[0.91+(0.022.*(k-1)) 0.5 0.02 0.45]);
        for i=1:8
            plot(-win(1):win(2),unit.sub(k).wf(i,:)+(yspc.*(i-1)),'Color',colorSpec(k+2,:),'LineWidth',2); hold on;
        end
        title(['u0' num2str(k)]);
        axis([-win(1) win(2) ylim]);
        axis off;
    end                    
    
    
    % Lever data and event data
    subplot('position',[0.3 0.45 0.6 0.045]);
    plot(reachDisp(2,:),'Color',[0.3 0.3 0.3]); hold on;
    plot(captThresh,ones(1,numel(captThresh)).*14000,'+','Color',colorSpec(4,:),'LineWidth',1,'MarkerSize',8);
    axis tight; 
    axis off;
    
    % Lick data and reward delivery data
    subplot('position',[0.3 0.4 0.6 0.045]);
    plot(lickDisp,'Color',[0.3 0.3 0.3]); hold on;
    plot(captRewards,ones(1,numel(captRewards)).*20000,'+','Color',colorSpec(3,:),'LineWidth',1,'MarkerSize',8);
    axis tight; 
    axis off;
    
    % Sorting quality example for unit 01
    subplot('position',[0.92 0.4 0.055 0.095]);
    plot3(noisePlot.x,-noisePlot.z,noisePlot.y./1e3,'.','MarkerSize',8,'Color',[0.7 0.7 0.7]); hold on;
    plot3(exemplarPlot.x,-exemplarPlot.z,exemplarPlot.y./1e3,'.','MarkerSize',8,'Color',colorSpec(ExemplarUnit+2,:));
    view([-150 10]); axis([-yspc 0 0 yspc 0 200])
    axis tight;    
