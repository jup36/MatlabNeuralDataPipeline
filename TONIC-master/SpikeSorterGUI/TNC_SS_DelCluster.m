function [newClustIds] = TNC_SS_DelCluster(handles)

    origClustIds = handles.featureData.seg(handles.segList(handles.currAxes)).shank(handles.shankNum).id;

    eval(['axes(handles.axes' num2str(handles.currAxes) ')']);

    [newX,newY] = ginput;

    hold on; plot(newX,newY,'w');

    inTheHull = inpolygon(handles.featureData.seg(handles.segList(handles.currAxes)).shank(handles.shankNum).params(:,handles.xPlotNum),handles.featureData.seg(handles.segList(handles.currAxes)).shank(handles.shankNum).params(:,handles.yPlotNum),newX,newY);

    toChange = find(inTheHull==1);
    
    newClustIds = origClustIds;
    newClustIds(toChange) = 0;