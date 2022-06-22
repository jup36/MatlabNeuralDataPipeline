function forSort = TNC_SS_CreateSortStruct(featStruct, gmStruct)

    if isfield(featStruct,'paramNames')
        forSort = featStruct;
        for i=1:numel(featStruct.seg)
            for j=1:numel(featStruct.seg(i).shank)
%               disp(['seg=' num2str(i) ' shank=' num2str(j) ...
%                     ' size(params)=' num2str(size(featStruct.seg(i).shank(j).params)) ]);
                if numel(gmStruct) > 0 && ...
                   numel(gmStruct.shank) >= j && ...
                   numel(gmStruct.shank(j).seg) >= i && ...
                   numel(gmStruct.shank(j).seg(i).Y) > 0
                    forSort.seg(i).shank(j).Y = gmStruct.shank(j).seg(i).Y;
                    forSort.seg(i).shank(j).Q = gmStruct.shank(j).seg(i).Q;
                else  
                    npar = size(featStruct.seg(i).shank(j).params, 1);
                    forSort.seg(i).shank(j).Y = char('A' *ones(1, npar));
                    forSort.seg(i).shank(j).Q =          zeros(1, npar);
                end
            end
        end
        if numel(gmStruct) == 0
            disp(' ');
            disp('...No color used to mark clusters because *_gm.mat file is not available');
        end;
    else

        for i=1:numel(featStruct.seg)
            for j=1:numel(featStruct.seg(i).shank)
                featStruct.seg(i).shank(j).params = ...
                   [featStruct.seg(i).shank(j).ts, ...
                    featStruct.seg(i).shank(j).pca, ...
                    featStruct.seg(i).shank(j).amp.energy'];
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
