function [ handles ] = TNC_SS_StoreAutoClustIds(handles)

    numSegs = numel(handles.featureData.seg);
    temp_segList = handles.segList;
    Alphabet = get_alphabet();

    for i=1:handles.numGraphs
        if  temp_segList(i) > 0        && ...
            temp_segList(i) <= numSegs && ...
            numel(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id)>0
            Y = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).Y;
            for j=1:numel(Alphabet)
                inds = find(Y == Alphabet(j));
                if numel(inds) > 0
                    handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id(inds) = j;    
                end
            end 
        end 
    end

% -------------------------------------------------------------------------------

function Alphabet = get_alphabet()
    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
    return

