function [newClustIds] = TNC_SS_AddCluster(handles,newId)

    origClustIds = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id;

    axes(handles.axes1);

    totalClusts = max(origClustIds);

    [newX,newY] = ginput;

    hold on; plot(newX,newY,'w');

    inTheHull = inpolygon(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,handles.xPlotNum),...
                          handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,handles.yPlotNum),...
                          newX,newY);

    toChange = find(inTheHull==1);
    
    newClustIds = origClustIds;
    
    if newId==0
        newId = totalClusts+1;
    end        
    
    newClustIds(toChange) = newId;
