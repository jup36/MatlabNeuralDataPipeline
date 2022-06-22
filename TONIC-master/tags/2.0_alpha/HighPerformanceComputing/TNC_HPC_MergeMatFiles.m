function [] = TNC_HPC_MergeMatFiles(varargin)
    num_input_files = nargin - 1;
    disp(['num_input_files=' num2str(num_input_files)]);
    output_name = varargin{nargin};
    len = length(output_name);
    eval(['load ' varargin{1}]);
    for i=1:num_input_files
        disp(['...loading '   varargin{i}]);
        ind_seg   = get_index(varargin{i}, '_seg');  
        ind_shank = get_index(varargin{i}, '_shank');
        disp(['ind_seg=' num2str(ind_seg) ' ind_shank=' num2str(ind_shank)]);
        mf = load(varargin{i});
        if length(findstr('_ft.mat', output_name)) > 0 
           % Report snr value in the name of output *_ft.mat file           
            if i == 1
%               name_split = strsplit(output_name, '_ft.mat'); % getting cell array
%               output_name = strcat(name_split{1}, '_snr', ...
%                                    num2str(mf.featStruct.snr), '_ft.mat');
                featStruct.snr = mf.featStruct.snr;
            end
            featStruct.seg(ind_seg).shank(ind_shank).inds   = ...
                mf.featStruct.seg(ind_seg).shank(ind_shank).inds;
            featStruct.seg(ind_seg).shank(ind_shank).ts     = ...
                mf.featStruct.seg(ind_seg).shank(ind_shank).ts;
            featStruct.seg(ind_seg).shank(ind_shank).params = ...
                mf.featStruct.seg(ind_seg).shank(ind_shank).params;
            numFolds = numel(mf.featStruct.seg(ind_seg).shank(ind_shank).random_split);
            for f=1:numFolds
                featStruct.seg(ind_seg).shank(ind_shank).random_split(f).folds = ...
                    mf.featStruct.seg(ind_seg).shank(ind_shank).random_split(f).folds;
            end
        elseif length(findstr('_ss.mat', output_name)) > 0
            sessionStruct.seg(ind_seg).shank(ind_shank).wfs  = ...
                mf.sessionStruct.seg(ind_seg).shank(ind_shank).wfs;
            sessionStruct.seg(ind_seg).shank(ind_shank).inds = ...
                mf.sessionStruct.seg(ind_seg).shank(ind_shank).inds;        
        elseif length(findstr('_gm',  output_name)) > 0  && ...
               length(findstr('.mat', output_name)) > 0
            if i == 1
                gmStruct.model      = mf.gmStruct.model;
                gmStruct.method     = mf.gmStruct.method;
                gmStruct.snr        = mf.gmStruct.snr;
                disp(['snr=' num2str(gmStruct.snr)]);
                % Update the output name by inserting the snr value
                if length(findstr('_gm.mat',  output_name)) > 0 
                    name_split = strsplit(output_name, '_gm.mat'); % getting cell array
                    output_name = strcat(name_split{1}, '_snr', ...
                                         num2str(mf.gmStruct.snr), '_gm.mat'); 
                end
            end
            if ind_seg > 0
                gmStruct.shank(ind_shank).seg(ind_seg).LL   = ...
                   mf.gmStruct.shank(ind_shank).seg(ind_seg).LL   
                gmStruct.shank(ind_shank).seg(ind_seg).xpar = ...
                   mf.gmStruct.shank(ind_shank).seg(ind_seg).xpar;
                gmStruct.shank(ind_shank).seg(ind_seg).PI = ...
                   mf.gmStruct.shank(ind_shank).seg(ind_seg).PI;
                gmStruct.shank(ind_shank).seg(ind_seg).MU = ...
                   mf.gmStruct.shank(ind_shank).seg(ind_seg).MU;
                gmStruct.shank(ind_shank).seg(ind_seg).SIGMA = ...
                   mf.gmStruct.shank(ind_shank).seg(ind_seg).SIGMA;
                if gmStruct.model == 2
                    gmStruct.shank(ind_shank).seg(ind_seg).EPS = ...
                        mf.gmStruct.shank(ind_shank).seg(ind_seg).EPS; 
                end
            else
                gmStruct.shank(ind_shank) = mf.gmStruct.shank(ind_shank);
            end
        end
    end
    if     length(findstr('_ft.mat', output_name)) > 0
        featStruct.program = 'TNC_HPC_MergeMatFiles.m';
        eval(['save ' output_name ' featStruct']);
    elseif length(findstr('_ss.mat', output_name)) > 0
        sessionStruct.program = 'TNC_HPC_MergeMatFiles.m';
        eval(['save ' output_name ' sessionStruct']); 
    elseif length(findstr('_gm',  output_name)) > 0 && ...
           length(findstr('.mat', output_name)) > 0
        gmStruct.program        = 'TNC_HPC_MergeMatFiles.m';
        eval(['save ' output_name ' gmStruct']);
    else 
        for i=1:num_input_files
            disp(['...loading ' varargin{i}])
            mf(i) = load(varargin{i});
        end
        eval(['save ' output_name ' mf']);
    end

% ------------------------------------------------------------------------------

function [ ind ] = get_index(mystring, pattern)
    if length(findstr(pattern, mystring)) > 0
        split1 = [ strsplit(mystring, pattern) ];
        disp(['split1=' split1]);
        split2 = strsplit(char(split1(2)), '_');
        ind    = int32(str2double(split2(1)));
    else
        ind = -1;
    end
