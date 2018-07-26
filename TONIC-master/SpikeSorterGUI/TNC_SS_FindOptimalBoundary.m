function [newIds] = TNC_SS_FindOptimalBoundary(handles)

%% Get the best estimate of current cluster center
    [handles]   = TNC_SS_UpdateClusterCenters(handles);

%% Find cluster boundaries at a range of std values
    % use x,y coords for centers
    cntX    = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).cnt(handles.pickId,handles.xPlotNum);
    cntY    = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).cnt(handles.pickId,handles.yPlotNum);
    
    thisClustInds = find(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id==handles.pickId);
    x       = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(thisClustInds,handles.xPlotNum);
    y       = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(thisClustInds,handles.yPlotNum);
    
    % calc eig vectors
    [u,s,v] = svd([x,y]);
    dim1 = v(:,1);
    dim2 = v(:,2);

    % project into rotated dimensions
    yPrime = [x,y] * dim1;
    xPrime = [x,y] * dim2;

    % recalculate cluster std
    stdX = std(xPrime,[],1);
    stdY = std(yPrime,[],1);

    % calculate angle
    angle = atan(dim1(1)./dim1(2)) .* (180./pi);

    % generate ellipse
    [eX1,eY1] = TNC_SS_CalcClusterEllipse(cntX,stdX,cntY,stdY,handles.confBound,angle,0);
    [eX2,eY2] = TNC_SS_CalcClusterEllipse(cntX,stdX,cntY,stdY,handles.confBound,angle,1);
    [eX3,eY3] = TNC_SS_CalcClusterEllipse(cntX,stdX,cntY,stdY,handles.confBound,angle,2);

%% Look for extra points biasing towards large eccentricities
    
    axes(handles.axes1);
    hold on; plot(eX1,eY1,'r',eX2,eY2,'r',eX3,eY3,'r','LineWidth',2); 
% pause();
%     
    newIds = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id;
%     
    xVals = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,handles.xPlotNum);
    yVals = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,handles.yPlotNum);
%     
%     inTheHull1   = inpolygon(xVals,yVals,eX1,eY1);
%     inTheHull2   = inpolygon(xVals,yVals,eX2,eY2);
%     inTheHull3   = inpolygon(xVals,yVals,eX3,eY3);
%     
%     level0 = find(inTheHull1==1);
%     level1 = find(inTheHull1==0 & inTheHull2==1);
%     level2 = find(inTheHull1==0 & inTheHull2==0 & inTheHull3==1);
%     
%     shift1 = (mean([xVals(level1) yVals(level1)],1) - mean([xVals(level0) yVals(level0)],1)) ./ (numel(level0) ./ numel(level1));
%     shift2 = (mean([xVals(level2) yVals(level2)],1) - mean([xVals(level0) yVals(level0)],1)) ./ (numel(level2) ./ (numel(level0) + numel(level1)));
%         

%% Use shifted center to recalculate
%     newCenter = [cntX cntY] + shift1 + shift2;
%     [eX,eY] = TNC_SS_CalcClusterEllipse(newCenter(1),stdX,newCenter(2),stdY,3,angle,0);
%         
    inTheHull               = inpolygon(xVals,yVals,eX3,eY3);
    
    toChangeIn              = find(inTheHull==1);
    
    newIds(toChangeIn)      = handles.pickId;    
        