function [newIds] = TNC_SS_GrowClusterBounds(handles)

    newIds = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id;
    axes(handles.axes1);
    totalClusts = numel(find(unique(newIds)>0));
    [newX,newY] = ginput(1);
   
%% Get rough boundaries from the chosen center point
    
    xVals           = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,handles.xPlotNum);
    yVals           = handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,handles.yPlotNum);

    distance        = sqrt( (xVals-newX).^2 + (yVals-newY).^2 );  
    
    distIncr        = round(max(distance)) ./ 100;
    distSpacing     = distIncr:distIncr:round(max(distance));
    
    distH           = hist(distance,distSpacing);

    smthHist = sgolayfilt(distH,3,21);
    
    inflection1 = find( diff(smthHist) < 0 , 1, 'first');
    inflection2 = find( diff(smthHist(inflection1:numel(smthHist)-1)) > 0 , 1, 'first') + inflection1 - 1;

    [eX,eY] = TNC_SS_CalcClusterEllipse(newX,distSpacing(inflection2),newY,distSpacing(inflection2),0.8,0,0);
    
    validPnts = inpolygon(xVals,yVals,eX,eY);
    inTheHull = find(validPnts==1);
    
    newIds(inTheHull) = handles.pickId;

%% For testing purposes

% newX = -3100;
% newY = 1650;
% 
% % newX = -50;
% % newY = 3300;
% 
%     xVals           = featStruct.seg(3).shank(1).pca(:,5);
%     yVals           = featStruct.seg(3).shank(1).pca(:,7);
% 
%     figure(1); semilogy(distSpacing,sgolayfilt(distH,3,21),'k',distSpacing(inflection2),smthHist(inflection2),'ro');    
%     figure(2); plot(xVals,yVals,'k.',newX,newY,'ro',eX,eY,'r--',xVals(inTheHull),yVals(inTheHull),'b.');

%% DEPRECATED VERSION    
%     distanceX = (xVals-newX);   
%     distanceY = (yVals-newY);  
%     
% 
%     
%     distIncrTmp1 = round( (max(distanceX)-min(distanceX)) ./ 250);
%     distIncrTmp2 = round( (max(distanceY)-min(distanceY)) ./ 250);
%     if distIncrTmp1>distIncrTmp2
%         distIncr = distIncrTmp2;
%     else
%         distIncr = distIncrTmp1;        
%     end
%     distSpacing = -1e4:distIncr:1e4;
%     
%     distXh = hist(distanceX,distSpacing);
%     distYh = hist(distanceY,distSpacing);
%     
%     midPoint = round(numel(distSpacing)./2);
%     peakX = max(distXh(midPoint-5:midPoint+5));
%     peakY = max(distYh(midPoint-5:midPoint+5));
%     
%     end0x   = find(distXh(midPoint:numel(distSpacing)) < peakX./3, 1, 'first');
%     start0x = find(distXh(1:midPoint) < peakX./3, 1, 'last');
% 
%     end0y   = find(distYh(midPoint:numel(distSpacing)) < peakY./3, 1, 'first');
%     start0y = find(distYh(1:midPoint) < peakY./3, 1, 'last');
%     
%     boundX = min([end0x start0x end0y start0y]);
%     boundY = min([end0x start0x end0y start0y]);
% % 
% %     pd = fitdist(distYh(midPoint-boundY:midPoint+boundY)','Normal');
% %     stdY = pd.std;
% 
% %     figure(1)
% %     subplot(221); plot(distSpacing,distXh,'k',distSpacing(midPoint-boundX:midPoint+boundX),distXh(midPoint-boundX:midPoint+boundX),'r');
% %     subplot(223); plot(distSpacing,distYh,'k',distSpacing(midPoint-boundY:midPoint+boundY),distYh(midPoint-boundY:midPoint+boundY),'r');
% %     
%     validPnts = find( distanceX>distSpacing(midPoint-boundX) & distanceX<distSpacing(midPoint+boundX) & distanceY>distSpacing(midPoint-boundY) & distanceY<distSpacing(midPoint+boundY) );
% %     subplot(2,2,[2,4]); plot(xVals,yVals,'k.',xVals(validPnts),yVals(validPnts),'ro');
% 
%     newIds(validPnts) = handles.pickId;
