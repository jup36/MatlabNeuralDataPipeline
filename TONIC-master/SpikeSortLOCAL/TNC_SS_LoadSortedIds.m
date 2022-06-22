function [handles] = TNC_SS_LoadSortedIds(handles,fileToLoad)

S = load(fileToLoad);

for j=1:numel(handles.featureData.seg(1).shank)
    for i=1:numel(handles.featureData.seg)
        handles.featureData.seg(i).shank(j).id = S.idList.seg(i).shank(j).id;
    end
end
