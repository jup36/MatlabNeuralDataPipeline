function [] = TNC_HPC_MergeMatFiles(varargin)
    if length(findstr('_kk_', varargin{1})) > 0
        my_matrix = zeros(nargin, 2);
        disp('M_true best_kk_M');
    else
        my_matrix = zeros(nargin, 4);
        disp('M_true best_M best2_M best3_M');
    end
    sum_sq = 0.;
    for i=1:nargin
        if exist(varargin{i}, 'file') == 2
%           disp(['...loading '   varargin{i}]);
            M_true = get_index(varargin{i}, 'close_');  
            mf = load(varargin{i});
            best_M_predicted = mf.mmStruct.best_mc;
            sum_sq = sum_sq + (M_true - best_M_predicted)^2;
            if length(findstr('_kk_', varargin{1})) == 0
                best2_M_predicted = mf.mmStruct.best2_mc;
                best3_M_predicted = mf.mmStruct.best3_mc;
                my_matrix(i,:) = [M_true best_M_predicted best2_M_predicted best3_M_predicted];
            else
                my_matrix(i,:) = [M_true best_M_predicted]; 
            end
        end
    end
    sorted_matrix = sortrows(my_matrix, 1);
    sigma = sqrt(sum_sq / nargin);
    sorted_matrix
    disp(' ');
    disp(['sigma=' num2str(sigma)]);

% ------------------------------------------------------------------------------

function [ ind ] = get_index(mystring, pattern)
    split1 = [ strsplit(mystring, pattern) ];
    split2 = strsplit(char(split1(1)), '_');
    ind    = str2double(split2(int32(length(split2))));

