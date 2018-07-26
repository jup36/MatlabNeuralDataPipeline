%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function [] = TNC_HPC_MergeMatFiles(varargin)
    num_input_files = nargin - 2;
    output_name = varargin{nargin-1};
    debug       = int32(str2double(varargin{nargin}));
    disp(['num_input_files=' num2str(num_input_files) ' debug=' num2str(debug) ' output_name=' output_name]);
    len = length(output_name);
    curr_shank = -1;
    best_mc_shank = [];
    eval(['load ' varargin{1}]);
    for i=1:num_input_files
        disp(['...loading '   varargin{i}]);
        ind_seg   = get_index(varargin{i}, '_seg');  
        ind_shank = get_index(varargin{i}, '_shank');
        disp(['ind_seg=' num2str(ind_seg) ' ind_shank=' num2str(ind_shank)]);
        try
            mf = load(varargin{i});
        catch
            continue;
        end
        if length(findstr('_ft.mat', output_name)) > 0 
           % Report thr value in the name of output *_ft.mat file           
            if i == 1
                try
                    featStruct.thr = mf.featStruct.thr;
                catch
                    featStruct.thr = 0;
                end
                if debug
                    name_split = strsplit(output_name, '_ft.mat'); % getting cell array
                    output_name = strcat(name_split{1}, '_thr', ...
                                         num2str(featStruct.thr), '_ft.mat');
                end
            end
            try
                featStruct.seg(ind_seg).shank(ind_shank).inds   = ...
                    mf.featStruct.seg(ind_seg).shank(ind_shank).inds;
                featStruct.seg(ind_seg).shank(ind_shank).id   = ...
                    mf.featStruct.seg(ind_seg).shank(ind_shank).id;
                featStruct.seg(ind_seg).shank(ind_shank).ts     = ...
                    mf.featStruct.seg(ind_seg).shank(ind_shank).ts;
                featStruct.seg(ind_seg).shank(ind_shank).params = ...
                    mf.featStruct.seg(ind_seg).shank(ind_shank).params;
                featStruct.seg(ind_seg).shank(ind_shank).Y_true = ...
                    mf.featStruct.seg(ind_seg).shank(ind_shank).Y_true;
                numFolds = numel(mf.featStruct.seg(ind_seg).shank(ind_shank).random_split);
                for f=1:numFolds
                    featStruct.seg(ind_seg).shank(ind_shank).random_split(f).folds = ...
                        mf.featStruct.seg(ind_seg).shank(ind_shank).random_split(f).folds;
                end
            catch
                disp(['Failed to extract information from file ' varargin{i} ...
                      '. Continue processing ...']);
                featStruct.seg(ind_seg).shank(ind_shank).inds = [];
                featStruct.seg(ind_seg).shank(ind_shank).ts   = [];
                featStruct.seg(ind_seg).shank(ind_shank).params = [];
                continue;
            end
        elseif length(findstr('_ss.mat', output_name)) > 0
            try 
                sessionStruct.seg(ind_seg).shank(ind_shank).wfs  = ...
                    mf.sessionStruct.seg(ind_seg).shank(ind_shank).wfs;
                sessionStruct.seg(ind_seg).shank(ind_shank).inds = ...
                    mf.sessionStruct.seg(ind_seg).shank(ind_shank).inds;        
            catch
                disp(['Failed to extract information from file ' varargin{i} ...
                      '. Continue processing ...']);
                sessionStruct.seg(ind_seg).shank(ind_shank).wfs = [];
                sessionStruct.seg(ind_seg).shank(ind_shank).inds = [];
                continue;
            end
        elseif length(findstr('_gm',  output_name)) > 0  && ...
               length(findstr('.mat', output_name)) > 0
            if i == 1
                gmStruct.model      = mf.gmStruct.model;
                gmStruct.method     = mf.gmStruct.method;
                try
                    gmStruct.thr    = mf.gmStruct.thr;
                catch
                    gmStruct.thr    = 0;
                end
                gmStruct.numFolds   = mf.gmStruct.numFolds;
                disp(['thr=' num2str(gmStruct.thr)]);
                % Update the output name by inserting the thr value
                if length(findstr('_gm.mat',  output_name)) > 0 
                    if debug
                        name_split = strsplit(output_name, '_gm.mat'); % getting cell array
                        output_name = strcat(name_split{1}, '_thr', ...
                                             num2str(gmStruct.thr), '_gm.mat'); 
                        disp(['debug=' num2str(debug) ' output_name=' output_name ]);
                    end
                end
            end
            % Handle the best_mc for entire shank
            disp(['ind_shank=' num2str(ind_shank) ' curr_shank=' num2str(curr_shank) ...
                  ' length(best_mc_shank)=' num2str(length(best_mc_shank))]);
            if ind_shank ~= curr_shank 
                if curr_shank > 0 && length(best_mc_shank) > 0
                    gmStruct.shank(curr_shank).best_mc  = round(mean(best_mc_shank));
                    disp(['curr_shank=' num2str(curr_shank) ' seg=1']);
                    gmStruct.shank(curr_shank).mc_range = gmStruct.shank(curr_shank).seg(1).mc_range;
                end
                curr_shank = ind_shank;
                best_mc_shank = [ mf.gmStruct.shank(ind_shank).seg(ind_seg).best_mc ];
            elseif ind_seg > 0
                best_mc_shank = [best_mc_shank mf.gmStruct.shank(ind_shank).seg(ind_seg).best_mc];
            end
            if ind_seg > 0
                try 
                    gmStruct.shank(ind_shank).seg(ind_seg).best_mc = ...
                        mf.gmStruct.shank(ind_shank).seg(ind_seg).best_mc;
                    gmStruct.shank(ind_shank).seg(ind_seg).mc_range = ...
                        mf.gmStruct.shank(ind_shank).seg(ind_seg).mc_range;
                    gmStruct.shank(ind_shank).seg(ind_seg).cvl = ...
                        mf.gmStruct.shank(ind_shank).seg(ind_seg).cvl;
                catch
                    disp('best_mc, mc_range and/or cvl is not available');
                end
                try 
                    gmStruct.shank(ind_shank).seg(ind_seg).Y_true = ...
                        mf.gmStruct.shank(ind_shank).seg(ind_seg).Y_true;
                    gmStruct.shank(ind_shank).seg(ind_seg).Y    = ...
                        mf.gmStruct.shank(ind_shank).seg(ind_seg).Y;
                    gmStruct.shank(ind_shank).seg(ind_seg).Q    = ...
                        mf.gmStruct.shank(ind_shank).seg(ind_seg).Q;
                    gmStruct.shank(ind_shank).seg(ind_seg).classError = ...
                        mf.gmStruct.shank(ind_shank).seg(ind_seg).classError;
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
                catch
                    disp(['Failed to extract information from file ' varargin{i} ...
                      '. Continue processing ...']);
                    gmStruct.shank(ind_shank).seg(ind_seg).Y_true = [];         
                    gmStruct.shank(ind_shank).seg(ind_seg).Y = [];
                    gmStruct.shank(ind_shank).seg(ind_seg).Q = [];
                    gmStruct.shank(ind_shank).seg(ind_seg).classError = [];
                    gmStruct.shank(ind_shank).seg(ind_seg).LL = [];
                    gmStruct.shank(ind_shank).seg(ind_seg).xpar = [];
                    gmStruct.shank(ind_shank).seg(ind_seg).PI = [];
                    gmStruct.shank(ind_shank).seg(ind_seg).MU = [];
                    gmStruct.shank(ind_shank).seg(ind_seg).SIGMA = [];
                    gmStruct.shank(ind_shank).seg(ind_seg).EPS = [];  
                    continue;        
                end
            else
                gmStruct.shank(ind_shank) = mf.gmStruct.shank(ind_shank);
            end
        end
    end

    % Finalize handling and output the structures
    if     length(findstr('_ft.mat', output_name)) > 0
        featStruct.program = 'TNC_HPC_MergeMatFiles.m';
        eval(['save ' output_name ' featStruct']);
    elseif length(findstr('_ss.mat', output_name)) > 0
        sessionStruct.program = 'TNC_HPC_MergeMatFiles.m';
        eval(['save ' output_name ' sessionStruct']); 
    elseif length(findstr('_gm',  output_name)) > 0 && ...
           length(findstr('.mat', output_name)) > 0
        gmStruct.program        = 'TNC_HPC_MergeMatFiles.m';
        % Handle the last shank
        disp(['curr_shank=' num2str(curr_shank)]);
        gmStruct.shank(curr_shank).best_mc = round(mean(best_mc_shank));
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
