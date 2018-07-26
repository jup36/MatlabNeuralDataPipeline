function [newIds] = TNC_SS_CropCluster(handles,inFlag)

% if cropId==0 then auto crop all clusters
% TO BE IMPLEMENTED


    [handles]   = TNC_SS_UpdateClusterCenters(handles);
    [eX,eY]     = TNC_SS_UpdateClusterBoundaries(handles,handles.segList(1),handles.cropId,handles.boundMethod);
    
    axes(handles.axes1);
    hold on; plot(eX,eY,'r','LineWidth',2); 

    newIds      = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id;    
    inTheHull   = inpolygon(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,handles.xPlotNum),handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,handles.yPlotNum),eX,eY);
    
    toChangeIn  = find(inTheHull==1);
    toChangeOut = find(inTheHull==0 & newIds==handles.cropId);
    
    if inFlag==1
        newIds(toChangeIn)  = handles.cropId;    
    end
    newIds(toChangeOut)     = 0;

