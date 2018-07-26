function TNC_HPC_ShowGMSolution(ft_file, ind_seg, ind_shank, method, M)
    % Example of usage:
    % TNC_GM_SortSpikes('WT4RD1LSTR002_2012_02_21_ft.mat', 'ms', '2', '1', '1')         
    % Second argument (method) may be one of: 'ps', 'gs', 'ms', 'ga', 'sa'
    % Third argument is # of Gaussian components
    % Forth and fifth arguments are the segment # and shank #
    % Sixth argument indicates whether or not to perform continuation from solution
    %       obtained for lower M (=0 - no continuation; > 0 - continuation from a given #) 
    % Seventh argument is node ID, used only in the name of output file
    %
    if nargin ~= 5
        disp('Usage: ShowGMSolution(ft_file, ind_seg, ind_shank, method, M)');
        return
    end  
    prefixName = ft_file(1:(length(ft_file)-7));
    gm_file = strcat(prefixName,'_seg',num2str(ind_seg), ...
                         '_shank',num2str(ind_shank),'_mc',num2str(M),...
                         '_sol',method, '_gm');
    
    disp(['gm_file=' gm_file]);
    disp(' ');
    ind_node = 1;
    best_fmin = 0;
    best_fmin2 = 0;
    xmin = [];
    xmin2= [];
    iseg   = int32(str2double(ind_seg));
    ishank = int32(str2double(ind_shank));
    while ind_node < 100
        input_gmFile = strcat(gm_file, '_', num2str(ind_node));
        if exist(strcat(input_gmFile,'.mat'), 'file') == 2
            gmData = load(input_gmFile); 
            fmin = gmData.gmStruct.seg(iseg).shank(ishank).fmin;
%           disp(['fmin=' num2str(fmin) ' best_fmin=' num2str(best_fmin)]);
            if  best_fmin > fmin
                best_fmin2 = best_fmin;
                best_fmin = fmin;
                xmin2 = xmin;
                xmin = gmData.gmStruct.seg(iseg).shank(ishank).xmin;
            elseif best_fmin <= fmin && best_fmin2 > fmin
                best_fmin2 = fmin;
                xmin2 = gmData.gmStruct.seg(iseg).shank(ishank).xmin;
            end
            disp([' input_gmFile=' input_gmFile ' fmin=' num2str(fmin) ' best_fmin=' num2str(best_fmin)]);
        end
        ind_node = ind_node + 1;
    end
     
    disp(['method=' method, ' M=' num2str(M) ' ind_seg=' ind_seg ' ind_shank=' ind_shank]);
    disp(['best_fmin=' num2str(best_fmin) ]);
    disp('xmin=');
    for i=1:(length(xmin)/6)
        disp(['           ' num2str(xmin((1:6)+(i-1)*6))]);
    end
    disp(['best_fmin2=' num2str(best_fmin2) ]);
    disp('xmin2=');
    for i=1:(length(xmin2)/6)
        disp(['           ' num2str(xmin2((1:6)+(i-1)*6))]);
    end
