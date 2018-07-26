%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function [] = TNC_HPC_ExtractBestSolution(varargin)
%   This application accepts as input files produced by TNC_GM_SortSpikes.m
    num_input_files = nargin - 1;
    output_name = varargin{nargin};
    disp(['num_input_files=' num2str(num_input_files) ' output_name=' output_name ' class(output_name)=' class(output_name)]);
    disp(' ');
    iSeg      = int32(str2double(get_value(char(output_name), '_seg')));
    iShank    = int32(str2double(get_value(char(output_name), '_shank'))); 
    M         = str2double(get_value(char(output_name), '_mc'));
    best_LLf  = -Inf;
    best2_LLf = -Inf;
    best3_LLf = -Inf;
    best_LL   = -Inf;
    xpar      = [];
    PI        = [];
    MU        = [];
    SIGMA     = [];
    SIGMA_INV = [];
    EPS       = [];
    Y         = [];
    Q         = [];
    Y_true    = [''];
    classError= NaN;
    LL        = 0.;
    method     = ''
    model      = [];
    thr       = 0.;
    gmStruct.program        = 'TNC_HPC_ExtractBestSolution.m';
    gmStruct.num_components = M;
    gmStruct.numFolds       = 0;
    for i=1:num_input_files
        file_exists = (exist(varargin{i}, 'file') == 2);
        disp(['   file exists=' num2str(file_exists) ' file name=' varargin{i}]);
        if file_exists
            disp(['...loading '   varargin{i}]);
            gmData   = load(char(varargin{i}));
            disp(['    file=' char(varargin{i}) ]);
            % Choose the replica with the highest LLf
            LLf      = gmData.gmStruct.shank(iShank).seg(iSeg).LLf  
            if  best_LLf  < LLf  | (best_LLf == LLf & best_LL < LL)
                best3_LLf = best2_LLf;
                best2_LLf = best_LLf;
                best_LLf  = LLf ;
                xpar      = gmData.gmStruct.shank(iShank).seg(iSeg).xpar;
                PI        = gmData.gmStruct.shank(iShank).seg(iSeg).PI;
                MU        = gmData.gmStruct.shank(iShank).seg(iSeg).MU;
                SIGMA     = gmData.gmStruct.shank(iShank).seg(iSeg).SIGMA;
                SIGMA_INV = gmData.gmStruct.shank(iShank).seg(iSeg).SIGMA_INV;
                Y_true    = gmData.gmStruct.shank(iShank).seg(iSeg).Y_true;
                Y         = gmData.gmStruct.shank(iShank).seg(iSeg).Y;
                Q         = gmData.gmStruct.shank(iShank).seg(iSeg).Q;
                classError= gmData.gmStruct.shank(iShank).seg(iSeg).classError;
                best_LL   = gmData.gmStruct.shank(iShank).seg(iSeg).LL;
                method    = gmData.gmStruct.method;
                model     = gmData.gmStruct.model;
                try
                    thr   = gmData.gmStruct.thr;
                catch
                    thr   = 0;
                end
                if gmData.gmStruct.numFolds > 0
                    gmStruct.numFolds = gmData.gmStruct.numFolds;
                end
            elseif best_LLf >= LLf & best2_LLf < LLf  
                best3_LLf  = best2_LLf;
                best2_LLf  = LLf;
            elseif best_LLf >= LLf & best2_LLf >= LLf & best3_LLf < LLf   
                best2_LLf  = LLf;
            end
            disp(['    LLf =' num2str(LLf ) ' best_LLf=' num2str(best_LLf)]);
        end
        disp(['best LLf=' num2str(best_LLf)]);
    end
    gmStruct.program        = 'TNC_HPC_ExtractBestSolution.m';
    gmStruct.num_components = M;
%   disp(['iSeg=' num2str(iSeg) ' iShank=' num2str(iShank)]);
    gmStruct.shank(iShank).seg(iSeg).LLf       = best_LLf;   
    gmStruct.shank(iShank).seg(iSeg).LLf2      = best2_LLf;
    gmStruct.shank(iShank).seg(iSeg).LLf3      = best3_LLf;

    gmStruct.shank(iShank).seg(iSeg).xpar      = xpar;
    gmStruct.shank(iShank).seg(iSeg).PI        = PI;
    gmStruct.shank(iShank).seg(iSeg).MU        = MU;   
    gmStruct.shank(iShank).seg(iSeg).SIGMA     = SIGMA;
    gmStruct.shank(iShank).seg(iSeg).SIGMA_INV = SIGMA_INV;
    gmStruct.shank(iShank).seg(iSeg).Y_true    = Y_true;
    gmStruct.shank(iShank).seg(iSeg).Y         = Y;
    gmStruct.shank(iShank).seg(iSeg).Q         = Q;
    gmStruct.shank(iShank).seg(iSeg).classError= classError;
    gmStruct.shank(iShank).seg(iSeg).LL       = LL;
    gmStruct.method         = method;
    gmStruct.model          = model;
    gmStruct.thr            = thr;

    eval(['save ' output_name ' gmStruct']);
    disp([' shank=' num2str(iShank) ' seg=' num2str(iSeg) ...
          ' M=' num2str(M) ' best_LLf=' num2str(best_LLf)]);
    

% ------------------------------------------------------------------------------

function [ val ] = get_value(mystring, pattern)
    split1 = strsplit(mystring, pattern);
%   disp(['mystring=' mystring, ' split1=' split1 ]);
    split2 = strsplit(char(split1(2)), '_');
    val    = split2(1);
