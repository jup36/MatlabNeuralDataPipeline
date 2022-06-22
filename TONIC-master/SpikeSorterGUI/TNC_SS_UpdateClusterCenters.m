function [handles] = TNC_SS_UpdateClusterCenters(handles)

temp_segList = handles.segList;
for i=1:handles.numGraphs
    if  temp_segList(i) < 1
        temp_segList(i) = 1;
    elseif temp_segList(i) > handles.numSegs
           temp_segList(i) = handles.numSegs;
    end 
end

for i=1:numel(temp_segList)
    
    clustNums = unique(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id);
    
    if max(clustNums) == 0
        
        handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).cnt=zeros(40,size(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params,2));
        handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).std=zeros(40,size(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params,2));
    
    else        
        for p=2:numel(clustNums)
            
            if clustNums(p)>0

                thisClustInds = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==clustNums(p));

                handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).cnt(clustNums(p),:) = mean(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(thisClustInds,:),1);
                handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).std(clustNums(p),:) = std(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(thisClustInds,:),[],1);
                                
            end
        end
    end    
end

handles.featureData.seg(temp_segList(i)).shank(handles.shankNum);
