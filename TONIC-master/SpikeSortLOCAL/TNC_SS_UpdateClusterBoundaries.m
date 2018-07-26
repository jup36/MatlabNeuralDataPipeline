function [eX,eY] = TNC_SS_UpdateClusterBoundaries(handles,currSeg,currClust,method)

    thisClustInds = find(handles.featureData.seg(currSeg).shank(handles.shankNum).id==currClust);
    x = handles.featureData.seg(currSeg).shank(handles.shankNum).params(thisClustInds,handles.xPlotNum);
    y = handles.featureData.seg(currSeg).shank(handles.shankNum).params(thisClustInds,handles.yPlotNum);


    switch method

        case 'ellipse'

            % use x,y coords for centers
            cntX    = handles.featureData.seg(currSeg).shank(handles.shankNum).cnt(currClust,handles.xPlotNum);
            cntY    = handles.featureData.seg(currSeg).shank(handles.shankNum).cnt(currClust,handles.yPlotNum);

%             stdX    = handles.featureData.seg(currSeg).shank(handles.shankNum).std(currClust,handles.xPlotNum);
%             stdY    = handles.featureData.seg(currSeg).shank(handles.shankNum).std(currClust,handles.yPlotNum);

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
            [eX,eY] = TNC_SS_CalcClusterEllipse(cntX,stdX,cntY,stdY,handles.confBound,angle,handles.scaler);

        case 'hull'

            dt  = DelaunayTri(x,y);
            k   = convexHull(dt);
            eX  = dt.X(k,1);
            eY  = dt.X(k,2);
    end
    