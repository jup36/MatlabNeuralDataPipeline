function [] = TNC_SS_SaveSortedIds(handles)

    for j=1:numel(handles.featureData.seg(1).shank)
        for i=1:numel(handles.featureData.seg)
            idList.seg(i).shank(j).id = handles.featureData.seg(i).shank(j).id;
        end
    end

    newFileName = [handles.fileName(1:numel(handles.fileName)-4) '_tns.mat'];
    save(newFileName,'idList');
    