function [newClustIds] = TNC_RmEvents(featureData,origClustIds,dimsToPlot)

    totalClusts = max(origClustIds);

    % display the original clustered output
    figure(201); clf;
    set(gcf, 'color', [0 0 0]);
    featureDims = size(featureData,2);
    set(gca, 'color', [0 0 0]);
    scatter(featureData(:,dimsToPlot(1)),featureData(:,dimsToPlot(2)),2,origClustIds,'filled'); hold on;
    scatter(origClustCenters(:,dimsToPlot(1)),origClustCenters(:,dimsToPlot(2)),84,'w');
    axis tight; axis off;
    
    [newX,newY] = ginput;
    
    figure(201); hold on;
    newX = [newX ; newX(1)];
    newY = [newY ; newY(1)];
    plot(newX,newY,'w');
    axis tight; 
    
    inTheHull = inpolygon(featureData(:,dimsToPlot(1)),featureData(:,dimsToPlot(2)),newX,newY);
    toChange = find(inTheHull==1);

    newClustIds = origClustIds;
    newClustIds(toChange) = -1;

    % display the updated clustered output
    figure(201); clf;
    set(gcf, 'color', [0 0 0]);
    featureDims = size(featureData,2);
    set(gca, 'color', [0 0 0]);
    scatter(featureData(:,dimsToPlot(1)),featureData(:,dimsToPlot(2)),2,newClustIds,'filled'); hold on;
    scatter(newClustCenters(:,dimsToPlot(1)),newClustCenters(:,dimsToPlot(2)),84,'w');
    axis tight; axis off;
    