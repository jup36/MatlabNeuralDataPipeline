function forSort = TNC_SS_CreateSortStruct(featStruct)

if isfield(featStruct,'paramNames')
    forSort = featStruct;
else

    for i=1:numel(featStruct.seg)
        for j=1:numel(featStruct.seg(i).shank)
            featStruct.seg(i).shank(j).params = [featStruct.seg(i).shank(j).ts, featStruct.seg(i).shank(j).pca, featStruct.seg(i).shank(j).amp.energy'];
            featStruct.seg(i).shank(j).cnt = zeros(40,size(featStruct.seg(i).shank(j).params,2));
            featStruct.seg(i).shank(j).std = zeros(40,size(featStruct.seg(i).shank(j).params,2));
        end
    end
    
    forSort = featStruct;

    forSort.paramNames{1} = 'ts';
    forSort.paramNames{2} = 'cPC1';
    forSort.paramNames{3} = 'cPC2';
    forSort.paramNames{4} = 'cPC3';
    forSort.paramNames{5} = 'cPC4';
    forSort.paramNames{6} = 'cPC5';
    forSort.paramNames{7} = 'cPC6';
    forSort.paramNames{8} = 'cPC7';
    forSort.paramNames{9} = 'cPC8';
    forSort.paramNames{10} = 'eng1';
    forSort.paramNames{11} = 'eng2';
    forSort.paramNames{12} = 'eng3';
    forSort.paramNames{13} = 'eng4';
    forSort.paramNames{14} = 'eng5';
    forSort.paramNames{15} = 'eng6';
    forSort.paramNames{16} = 'eng7';
    forSort.paramNames{17} = 'eng8';
    
end
