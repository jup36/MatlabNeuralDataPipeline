function [] = TNC_SS_WriteData(handles)

    % expected structure for session data in other code:
    % session(I).unit(J).ts

    % going to supply shank.unit(J).ts which will be converted after the fact
    % into expected sessions structure by an external script

    fileNameBase    = handles.fileName(1:numel(handles.fileName)-7);
    numSegs         = numel(handles.featureData.seg);

    clustList = [];
    for i=1:numSegs
        nonZ        = find(handles.featureData.seg(i).shank(handles.shankNum).id>0);
        clustList   = [clustList unique(handles.featureData.seg(i).shank(handles.shankNum).id(nonZ))'];
    end

    allClusts = unique(clustList)
    numClusts = numel(allClusts);

    for j=1:numClusts
        
        disp(['Collecting timestamps for cluster ' num2str(allClusts(j))]);
        
        tsTmp = []; indTmp = [];        
        
        for i=1:numSegs
            
            inds    = find(handles.featureData.seg(i).shank(handles.shankNum).id==allClusts(j));            
%             disp([num2str(numel(inds)) ' of ' num2str(numel(handles.featureData.seg(i).shank(handles.shankNum).inds)) ' events in seg ' num2str(i)]);
            if size(handles.featureData.seg(i).shank(handles.shankNum).ts,1) < size(handles.featureData.seg(i).shank(handles.shankNum).ts,2)
                tsTmp   = [tsTmp handles.featureData.seg(i).shank(handles.shankNum).ts(inds)+ ((i-1).*handles.featureData.chunk.*1000)];
                indTmp  = [indTmp handles.featureData.seg(i).shank(handles.shankNum).inds(inds)];
            else
                tsTmp   = [tsTmp handles.featureData.seg(i).shank(handles.shankNum).ts(inds)'+ ((i-1).*handles.featureData.chunk.*1000)];
                indTmp  = [indTmp handles.featureData.seg(i).shank(handles.shankNum).inds(inds)'];
            end
        end
        
        shank.unit(j).ts    = tsTmp;
        shank.unit(j).inds  = indTmp;
        
    end

    shank.shankNum = handles.shankNum;
    
    disp(['save ' fileNameBase '_shank' num2str(handles.shankNum) '_tsd shank']); 
    eval(['save ' fileNameBase '_shank' num2str(handles.shankNum) '_tsd shank']); 