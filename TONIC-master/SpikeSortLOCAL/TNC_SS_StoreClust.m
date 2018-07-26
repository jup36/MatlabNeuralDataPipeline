function [handles] = TNC_SS_StoreClust(leaveSeg,enterSeg,handles)

clustNum = handles.clustToProp;

if clustNum==0

    clustNums = unique(handles.featureData.seg(leaveSeg).shank(handles.shankNum).id);
    for p=1:numel(clustNums)
        if clustNums(p) > 0
            handles.featureData.seg(enterSeg).shank(handles.shankNum).cnt(clustNums(p),:) = handles.featureData.seg(leaveSeg).shank(handles.shankNum).cnt(clustNums(p),:);
            handles.featureData.seg(enterSeg).shank(handles.shankNum).std(clustNums(p),:) = handles.featureData.seg(leaveSeg).shank(handles.shankNum).std(clustNums(p),:);

            Calculate the ellipse that contains the events
                cntX    = handles.featureData.seg(enterSeg).shank(handles.shankNum).cnt(clustNums(p),handles.xPlotNum);
                stdX    = handles.featureData.seg(enterSeg).shank(handles.shankNum).std(clustNums(p),handles.xPlotNum);
                cntY    = handles.featureData.seg(enterSeg).shank(handles.shankNum).cnt(clustNums(p),handles.yPlotNum);
                stdY    = handles.featureData.seg(enterSeg).shank(handles.shankNum).std(clustNums(p),handles.yPlotNum);

                [eX,eY] = TNC_SS_UpdateClusterBoundaries(handles,leaveSeg,clustNums(p),handles.boundMethod);

            Find all of the events in the new segment that fall in that cluster

                inTheHull   = inpolygon(handles.featureData.seg(enterSeg).shank(handles.shankNum).params(:,handles.xPlotNum),handles.featureData.seg(enterSeg).shank(handles.shankNum).params(:,handles.yPlotNum),eX,eY);
                toChange    = find(inTheHull==1);

                handles.featureData.seg(enterSeg).shank(handles.shankNum).id(toChange) = clustNums(p);
        end        
    end    
    handles.featureData.seg(leaveSeg).shank(handles.shankNum)
    handles.featureData.seg(enterSeg).shank(handles.shankNum)
    
else

    handles.featureData.seg(handles.segList(enterSeg)).shank(handles.shankNum).cnt(clustNum,:) = handles.featureData.seg(handles.segList(leaveSeg)).shank(handles.shankNum).cnt(clustNum,:);
    handles.featureData.seg(handles.segList(enterSeg)).shank(handles.shankNum).std(clustNum,:) = handles.featureData.seg(handles.segList(leaveSeg)).shank(handles.shankNum).std(clustNum,:);

end