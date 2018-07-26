%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function TNC_GM_SortSpikes(ftFile,iSeg,iShank,M,varargin)      
    % ftFile, iSeg, iShank and M = required positional arguments
    % iRepl                      = optional, position-specific argument
    % lowerM, verbose            = optional parameters
    %
    % Example of usage:
    % TNC_GM_SortSpikes('WT4RD1LSTR002_2012_02_21_ft.mat','2','3',3,'ms',16,'verbose',1)
    % Sixth argument indicates whether or not to perform continuation from solution
    %       obtained for lower M (=0 - no continuation; > 0 - continuation from a given #) 
    % Seventh argument is node ID, used only in the name of output file
    %  
    tic;

    try  
        p = TNC_OOP_SpikeSortingInputParser;     % call constructor
        p.parse(ftFile, iSeg, iShank, M, varargin{:});

        % Required inputs:
        options.ftFile =       p.Results.ftFile;
        options.iSeg   = int32(p.Results.iSeg);
        options.iShank = int32(p.Results.iShank);
        options.M      = int32(p.Results.M);

        % Optional, but position specific inputs
        options.iRepl  = int32(p.Results.iRepl);

        % Optional parameters
        options.method       = p.Results.method;
        options.guess_method = p.Results.guess_method;
        compiled_code = 0;
        if isnumeric(iSeg)     % original Matlab code     
            options.verbose  = p.Results.verbose;
            options.debug    = p.Results.debug;
            options.dispOn   = p.Results.dispOn;
            options.model    = p.Results.model;
            options.power    = p.Results.power;
            options.fakeData = p.Results.fakeData;
            options.maxIter  = p.Results.maxIter;
            options.numFolds = p.Results.numFolds;
            options.iFold    = p.Results.iFold;
            options.parNames = p.Results.parNames;
        else 
            compiled_code    = 1;    % compiled code    
            options.verbose  = int32(str2double(p.Results.verbose));
            options.debug    = int32(str2double(p.Results.debug));
            options.dispOn   = int32(str2double(p.Results.dispOn));
            options.model    = int32(str2double(p.Results.model));
            options.power    =       str2double(p.Results.power);
            options.fakeData = int32(str2double(p.Results.fakeData));
            options.maxIter  = int32(str2double(p.Results.maxIter));
            options.numFolds = int32(str2double(p.Results.numFolds));
            options.iFold    = int32(str2double(p.Results.iFold));
            options.parNames =                  p.Results.parNames;
 
        end
    catch
        output_usage_message();
        return
    end

    if check_inputs(options) ~= 1
        return;
    end

    if length(strfind(options.ftFile, '_ft')) == 0
        disp(['Input file is not of ft type']);
        return
    end

    if options.verbose
        disp(['ftFile=' options.ftFile ' iSeg=' num2str(options.iSeg) ' iShank=' num2str(options.iShank) ...
              ' M=' num2str(options.M) ' iRepl=' num2str(options.iRepl) ' method=' options.method ...
              ' guess_method=' options.guess_method ... 
              ' model=' num2str(options.model) ' dispOn=' num2str(options.dispOn) ...
              ' fakeData=' num2str(options.fakeData) ' debug=' num2str(options.debug)]);
    end

    % DATA_f = data from specified fold; DATA_nf = complimentary data
%   try
        [K,ftData,DATA_nf,DATA_f,Y_true,options,status] = extract_data(options);
%   catch
%       disp(['Failed to read data for segment ' num2str(options.iSeg) ' shank ' num2str(options.iShank) ]);
%       return;
%   end

    if status ~= 0
        return;
    end

    if options.fakeData
        try
            options.coord.x_neurons = ftData.featStruct.x_neurons;
            options.coord.y_neurons = ftData.featStruct.y_neurons;
            options.coord.z_neurons = ftData.featStruct.z_neurons;
        catch
            disp('Neuron coordinates are not found in the input file');
            disp('Consider re-running feature extraction with fakeData option on');
            return
        end
    end
   
    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'; 
    prefixName = options.ftFile(1:(length(options.ftFile)-7));
    options.prefixName = prefixName;
    disp('Start ...');

    if strcmp(options.method, 'kk')
        options.maxIter = 1;
    end
    LL      = -Inf;
    LLf     = -Inf;
    Iter    = 1;
    success = 0;
    disp(['options.M=' num2str(options.M) ' options.verbose='  num2str(options.verbose)]);
    while ~success && Iter < 10                 
        Iter = Iter + 1;
        try
            [init_guess, MU, SIGMA, DET] = TNC_GM_InitialGuess(DATA_nf, Iter, options);
            if options.verbose
                disp('After Initial guess:');
                MU
                SIGMA
            end
            [LL, xpar, LLf, Y, Q, PI, MU, SIGMA, SIGMA_INV, status] = ...
                get_opt_solution(K, init_guess, MU, SIGMA, DATA_nf, ...
                                 DATA_f, Y_true, Alphabet, options);
            disp(['status=' num2str(status) ' LL=' num2str(LL) ' LLf=' num2str(LLf)]);
            classError = NaN;
            assert(imag(LL)  == 0);
            assert(imag(LLf) == 0);
            if options.fakeData
                disp(['Y_true=' num2str(Y_true)]);
                disp(['Y     =' num2str(Y     )]);
                counts = get_label_counts(Y);
                disp(['Y_unique=', num2str(unique(Y)) ' counts=' num2str(counts)]);
                [map, classError, Y_best] = compare_class_labels(length(PI), Y_true, Y, ...
                                                                 Alphabet, options.verbose);
                disp(' ');
                disp(['Y_true=' num2str(Y_true)]);
                disp(['Y_best=' num2str(Y_best)]);
                classError
            end
            success = 1; 
            break;    
        catch
            if ~strcmp(options.method, 'kk')
                disp(['status=' num2str(status) ' replica=' num2str(options.iRepl) ' LL=' num2str(LL) ' LLf=' num2str(LLf)]);
                disp(['get_opt_solution FAILED in attempt# ' num2str(Iter) ' ...']);
            end
        end
    end
    if ~strcmp(options.method, 'kk') && success    
        disp(' ');
        disp(['# Attempts made=' num2str(Iter)]);
        disp(['Output to gm_data file: LL=' num2str(LL) ' LLf='  num2str(LLf) ' prefix=' prefixName ' iFold=' num2str(options.iFold)]);
        output_gm_data(xpar, LL, LLf, Y_true, Y, Q, classError, PI, MU, SIGMA, ...
                       SIGMA_INV, prefixName, ftData.featStruct, options);
    end

    toc;

    % Make sure the code will terminate
    if compiled_code
        exit;
    end

    return

% ---------------------------------------------------------------------

function counts = get_label_counts(Y)
    Y_unique = unique(Y);
    counts = zeros(1, length(Y_unique));
    for i=1:length(Y_unique)
        counts(i) = length(find(Y == Y_unique(i)));
    end

% ---------------------------------------------------------------------

function output_usage_message()
    disp('Usage: TNC_GM_SortSpikes(ftFile, iSeg, iShank, M [,iRepl] [,parameters])');
    disp('Required arguments:');
    disp('    ftFile       - name of an input file');
    disp('    iSeg         - index of the time segment to be processed');
    disp('    iShank       - index of the shank to be processed');
    disp('    M            - number of Gaussian components');
    disp('Optional argument:')
    disp('    iRepl        - index of solution replica (used only in cluster processing)');
    disp('Optional parameters, specified as parameter-value pairs:');
    disp('    debug        - weather or not to keep intermediate data files ');
    disp('                   and shell scripts after cluster processing');
    disp('    dispOn       - weather or not to display covariance matrix (default=0)');
    disp('    fakeData     - weather or not the input data is fake (default=0)');
    disp('    guess_method - method used to generate initial guess for EM algorithm');
    disp('                   default = km (k-means); alternative = rc (random cell)');
    disp('    iFold        - index of the fold to be tested in cross-validation');
    disp('                   (>= 0; <= numFolds, default = 1');
    disp('    maxIter      - max number of iterations to be performed');
    disp('                   (default = Inf)');
    disp('    method       - computational method to be used');
    disp('                   default=em (for Expectation Maximization)');
    disp('                   other available options: ms (for MultiStart)');
    disp('                   and kk (for KlustaKwick)');
    disp('    numFolds     - total # of folds to be used in cross-validation');
    disp('                   (default = 1)');
    disp('    parNames     - a string of comma-separated names of parameters/features ;');
    disp('                   (default = min; allowed values = cpc,eng,min,pc1,pc2,pc3)');
    disp('    power          when doing triangulation procedure, the external');
    disp('                   potential is assumed to drop as this power');
    disp('                   of distance (default = 1)');
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return

% ---------------------------------------------------------------------

function status = check_inputs(options)
    status = 1;
    if strcmp(options.guess_method, 'rc') && options.iRepl <= 0
        disp(' ');
        disp('With random cell initial guess, positive node # must be specified');
        status = 0;
    end
    if options.numFolds < 0 || options.numFolds > 50
        disp(' ');
        disp('numFolds must be between 0 and 50');
        status = 0;
    end
    if options.iFold < 0 || options.iFold > options.numFolds
        disp(' ');
        disp('iFold should be non-negative and should not acceed numFolds'); 
        status = 0;
    end
    parNames = char(strsplit(options.parNames, ',')');
    if size(parNames, 1) > 6
        disp(' ');
        disp('currently supported number of parameters does not exceed 6');
        status = 0;
    end
    

% ---------------------------------------------------------------------

function num_features = get_num_features(paramNames)
    paramNames_short = paramNames;
    disp(['numel(paramNames)=' num2str(numel(paramNames))]);
    for i=1:numel(paramNames)
        my_str = paramNames{i};
        len = length(my_str)-1;
        disp(['i=' num2str(i) ' my_str=' my_str ' len=' num2str(len)]);
        paramNames_short{i} = my_str(1:len);
    end
    disp(['paramNames_short='])
    paramNames_short
    num_features = numel(unique(paramNames_short)) -1;

% ---------------------------------------------------------------------

function [K,ftData,DATA_nf,DATA_f,Y_true,options,status] = extract_data(options)
    status = 0;
    ftData = load(char(options.ftFile));
    try
        numSegs = length(ftData.featStruct.seg);
    catch
        disp('No detected spikes found in the input file.');
        disp('May need to decrease thr');
        status = -1;
        K = 0;
        ftData =[];
        DATA_nf=[];
        DATA_f =[];
        Y_true =[''];
        return
    end 
    K = 0;
    DATA_f  = []; % data from the specified fold
    DATA_nf = []; % complimentary data (not from the fold)
    Y_true  = [''];
    if options.iSeg > numSegs
        disp(' ');
        disp('Optional parameter iSeg exceeds numSegs'); 
        status = 1;
        return
    end

    shankSize = length(ftData.featStruct.seg(options.iSeg).shank);
    if options.iShank > shankSize
        disp(' ');
        disp('Optional parameter iShank exceeds shank size');
        status = 1;
        return
    end 

    % Determine the indices of the columns of parameters
    % that will be used in spike sorting
    npar = size(ftData.featStruct.seg(options.iSeg).shank(options.iShank).params,2);
    num_features = get_num_features(ftData.featStruct.paramNames);
    disp(['num_features=' num2str(num_features)]);
    KMAX = (npar-1)/num_features;
    options.iE = (1:KMAX);
    parNames = char(strsplit(options.parNames, ',')');
    keys  = {'cpc', 'eng', 'min', 'minT', 'pc1', 'pc2', 'pc3'};
    inds  = {1,     2,     3,      4,     5,     6,      7};
    params_hash = containers.Map(keys, inds);
    ipar = [];
    for i=1:size(parNames,1)
        try
            ipar = [ipar params_hash(parNames(i,:))];
        catch
            disp(' ');
            disp(['NOTE: ignoring  invalid parameter name ' parNames(i,:)]);
            disp(' ');
        end
    end
    cols = [];
    for i=1:size(parNames,1)
        add_cols = 1+KMAX*(ipar(i)-1)+options.iE;
        cols     = [cols add_cols];
    end
    disp(['ipar=' num2str(ipar)]);

    if options.verbose
        disp(['iSeg=' num2str(options.iSeg) ' iShank=' num2str(options.iShank)]);
        disp(['size(ftData.featStruct.seg(iSeg).shank(iShank).params)=' ...
               num2str(size(ftData.featStruct.seg(options.iSeg).shank(options.iShank).params))]);
        disp(['size(params)=' num2str(size(ftData.featStruct.seg(options.iSeg).shank(options.iShank).params)) ...
              ' KMAX=' num2str(KMAX)]);
    end
    sg = options.iSeg;
    sh = options.iShank;
    DATA0 = ftData.featStruct.seg(sg).shank(sh).params;
    N     = size(DATA0, 1);
    N
    KMAX
    DATA  = zeros(N, KMAX); % KMAX columns = normalized 'distances' between spikes

    % Compute DATA as 'distances' between spikes
    if length(ipar) > 1
        for j=1:length(ipar)
            ind = 1 + KMAX*(ipar(j)-1)+(1:KMAX);
            ind
            min_data = min(min(DATA0(:,ind)));
            std_data = std(std(DATA0(:,ind)));
            min_data
            std_data
            for i=1:KMAX
                c = 1 + KMAX*(ipar(j)-1)+i;
                DATA(:,i) = DATA(:,i) + ((DATA0(:,c) - min_data)/std_data).^2;
            end
        end
        for i=1:KMAX
            DATA(:,i) = sqrt(DATA(:,i));
        end
    else
        min_data = min(min(DATA0(:,cols)))
        DATA     = abs(DATA0(:, cols)-min_data);
    end
%   DATA

    if options.numFolds == 1 || options.iFold == 0
        DATA_nf = DATA;
        % Extract the 'true' class labels
        if options.fakeData
            Y_true = ftData.featStruct.seg(options.iSeg).shank(options.iShank).Y_true;
        else
            Y_true = [''];
        end
    else
        rows   = [];
        rows_f = [];
        nf   = options.numFolds;
        rs   = ftData.featStruct.seg(sg).shank(sh).random_split(nf);
        for i=1:nf
            if i ~= options.iFold
                rows   = [rows   rs.folds(:,i)'];
            else
                rows_f = [rows_f rs.folds(:,i)'];
            end
        end
        num_data   = size(DATA, 1);
        res        = mod(num_data, nf);
        if res > 0
            % Compliment data with (nf - res) NANs to make its size multiple of nf
            num_data1= num_data + (nf - res);
        else
            num_data1= num_data;
        end
        DATA_compl = [DATA; NaN*ones(num_data1 - num_data, length(cols))];
        DATA_nf    =  DATA_compl(rows  ,:);
        DATA_f     =  DATA_compl(rows_f,:);
        % Eliminate the rows with NaNs
        inds_NaN   = find(isnan(DATA_nf(:, 1)));
        DATA_nf(inds_NaN, :) = [];
        inds_NaN   = find(isnan(DATA_f(:, 1)));
        DATA_f(inds_NaN, :) = [];
        disp(['num_data=' num2str(num_data) ' num_data1=', num2str(num_data1) ...
              ' nf=' num2str(nf) ' res=' num2str(res) ' size(rows)=' ...
              num2str(size(rows)) ' size(rows_f)=' num2str(size(rows_f))]);
    end
    disp(['num DATA_nf rows=' num2str(size(DATA_nf, 1))]);
%   disp(['size(params0=' num2str(size(ftData.featStruct.seg(options.iSeg).shank(options.iShank).params))]);
    disp(['options.fakeData=' num2str(options.fakeData)]);
    K = size(DATA_nf, 2);

% ---------------------------------------------------------------------

function [LL, xpar, LLf, Y, Q, PI, MU, SIGMA, SIGMA_INV, num_iter] = ...
          get_opt_solution(K, initial_guess, MU, SIGMA, DATA_nf, DATA_f, ...
                           Y_true, Alphabet, my_options)
    Y = [];
    Q = [];

    % Variables: MU, SIGMA^2, X, Y, Z for each cell 
    num_iter = 0; % # iterations performed

    % Upper and lower bounds
    [lb, ub] = get_bound_constraints(my_options.M, K, my_options.model);

    % Equality constraints
    [Aeq, beq] = get_equality_constraints(my_options.M, K, my_options.model);
 
    % Inequality constraints
    if my_options.M > 1
        [Aineq, bineq] = get_inequality_constraints(my_options.M, K, my_options.model);
    else
        Aineq = [];
        bineq = [];
    end

    % Initial guess 
    [x0] = initial_guess;

    if strcmp(my_options.method, 'em') % Simulated Annealing (only bound constraints!)
        % Only bound constraints are allowed
        disp('Applying Expectation-Maximization algorithm ...');
        my_options.TolFun = 0.1;
        my_options.TolX   = 1.e-8;  
        [xpar, LL, Y, Q, PI, MU, SIGMA, SIGMA_INV, num_iter] = ...
            ExpectationMaximization(x0, MU, SIGMA, DATA_nf, ...
            Y_true, Alphabet, my_options);     
          
        disp(['# EM iterations=' num2str(num_iter) ' LL=' num2str(LL)]);
        if my_options.iFold == 0 || my_options.numFolds == 1
            LLf = 0;
        elseif LL > -Inf
            M = length(PI);
            SIGMA_DET = zeros(1,M);
            for m=1:M
                SIGMA1   = SIGMA(:,(K*(m-1)+1):(K*m));
                SIGMA_DET(m) = det(SIGMA1);
            end
%           disp(['SIGMA_DET=' num2str(SIGMA_DET)]);
            LLf = TNC_GM_GaussMixLogLikelihood(PI, MU, SIGMA, SIGMA_DET, SIGMA_INV, DATA_f, 2);
        else
            LLf = -Inf;
        end
        disp(['my_options.iFold=' num2str(my_options.iFold) ' size(DATA_nf)=' num2str(size(DATA_nf,1)) ...
              ' size(DATA_f)=' num2str(size(DATA_f,1)) ' LL=' num2str(LL) ' LLf=' num2str(LLf)]);
        disp(' ');
    elseif strcmp(my_options.method, 'kk') % Simulated Annealing (only bound constraints!)
        % Only bound constraints are allowed
        disp('Applying KlustaKwik algorithm ...');
        disp(['size(DATA_nf)=' num2str(size(DATA_nf))]);
        PI = [];
        MU = [];
        SIGMA = [];
        SIGMA_INV = [];
        GAMMA = [];
        LLf = 0;
        if size(DATA_nf, 2) == 4
            [xpar, LL, M, Y] = KlustaKwik(DATA_nf, my_options);
            disp(' ');
            disp(['num_clusters=' num2str(M)]);
        elseif size(DATA_nf, 2) == 8
            DATA1 = DATA_nf(:,[1 3 5 7]);
            [xpar, LL, M1, Y] = KlustaKwik(DATA1, my_options);
            DATA2 = DATA_nf(:,[2 4 6 8]);
            [xpar, LL, M2, Y] = KlustaKwik(DATA2, my_options);
            disp(' ');
            disp(['num_clusters1=' num2str(M1) ' num_clusters2=' num2str(M2)]);
        end
        num_iter = 1;
    else
        disp([ 'Incorrectly specified name of global opt method: ' method ]);
        return;
    end

%-------------------------------------------------------------------------------

function  [lb, ub] = get_bound_constraints(M, K, model)

    if model ~= 2
        lb = ones(1, 6*M);
        ub = ones(1, 6*M);
    else
        lb = ones(1, 6*M+K);
        ub = ones(1, 6*M+K);
    end
    lb_one = [0 -2.0  1.e-6  -50 -50   0]; % PI, MU, SIGMA, X, Y, Z
    ub_one = [1  0   10       50 250  50]; % PI, MU, SIGMA, X, Y, Z
    for m=1:M
        lb((1+(m-1)*6):(m*6)) = lb_one;
        ub((1+(m-1)*6):(m*6)) = ub_one;
    end                              
    if model == 2       
        min_eps =  1.e-9;
        max_eps = 10;
        lb(end-K+1:end) = min_eps*ones(1,K);
        ub(end-K+1:end) = max_eps*ones(1,K);
    end    

%-------------------------------------------------------------------------------

function [Aeq, beq] = get_equality_constraints(M, K, model)
    if model == 2       
        Aeq = zeros(1,M*6+K);
    else
        Aeq = zeros(1,M*6);
    end
    Aeq(1,1+(0:(M-1))*6) = 1; % sum of PIs = 1, or A(1,1)*X(1)+A(1,7)*x(7)+...= 1
    beq = 1;

%-------------------------------------------------------------------------------

function [Aineq, bineq] = get_inequality_constraints(M, K, model)
    if model == 2       
        Aineq = zeros(M-1,M*6+K); 
    else
        Aineq = zeros(M-1,M*6);
    end
    for m=1:(M-1)
        Aineq(m,1+(m-1)*6) =  1;  
        Aineq(m,1+ m   *6) = -1; % MU_m <= MU_m+1 
    end
    bineq = -zeros(1, M-1);

%-------------------------------------------------------------------------------

function output_gm_data(xpar, LL, LLf, Y_true, Y, Q, classError, PI, MU, SIGMA, ...
                        SIGMA_INV, targetName, featStruct, options)
    K = length(MU);
    M = length(PI);
    iSeg                    = options.iSeg;
    iShank                  = options.iShank;
    gmStruct.program        = 'TNC_GM_SortSpikes.m';
    gmStruct.model          = options.model;
    gmStruct.paramNames     = {'PI', 'MU' 'SIGMA' 'XC' 'YC' 'ZC'};
    gmStruct.num_components = M;
    try
        gmStruct.thr        = featStruct.thr;
    catch
        gmStruct.thr        = 0;
    end
    gmStruct.method         = options.method;
    gmStruct.iFold          = options.iFold;
    gmStruct.numFolds       = options.numFolds;
    gmStruct.shank(iShank).seg(iSeg).Y_true= Y_true;
    gmStruct.shank(iShank).seg(iSeg).Y     = Y;
    gmStruct.shank(iShank).seg(iSeg).Q     = Q;
    gmStruct.shank(iShank).seg(iSeg).classError = classError;
    gmStruct.shank(iShank).seg(iSeg).PI    = PI;
    gmStruct.shank(iShank).seg(iSeg).MU    = MU;
    gmStruct.shank(iShank).seg(iSeg).SIGMA = SIGMA;
    gmStruct.shank(iShank).seg(iSeg).SIGMA_INV = SIGMA_INV;
    gmStruct.shank(iShank).seg(iSeg).xpar  = xpar;
    gmStruct.shank(iShank).seg(iSeg).LL    = LL;
    gmStruct.shank(iShank).seg(iSeg).LLf   = LLf;
    if options.model == 2       
        eps_beg = int32(length(xpar)-K+1);
        eps_end = int32(length(xpar));
        disp(['iSeg=' num2str(options.iSeg) ' iShank=' num2str(options.iShank) ...
              ' eps_beg=' num2str(eps_beg) ' eps_end=' num2str(eps_end) ]);
        gmStruct.shank(iShank).seg(iSeg).EPS = xpar(eps_beg:eps_end);
    end
    if options.fakeData
        gmStruct.coord.x_neurons = featStruct.x_neurons;
        gmStruct.coord.y_neurons = featStruct.y_neurons;
        gmStruct.coord.z_neurons = featStruct.z_neurons;
    end
    disp(['targetName=' targetName]);
    if options.iRepl > 0 % run on cluster, not locally
        output_name = strcat(targetName,'_shank', num2str(options.iShank),...
                             '_seg',num2str(options.iSeg), '_mc',num2str(M),...
                             '_gm_', num2str(options.iFold), ...
                             '_',    num2str(options.iRepl));
    elseif options.iFold > 0 && options.numFolds > 1
        output_name = strcat(targetName,'_shank',num2str(options.iShank),...
                             '_seg',num2str(options.iSeg), '_mc',num2str(M),...
                             '_gm_', num2str(options.iFold));
    else
        output_name = strcat(targetName, '_gm');
    end
    disp(['save ' output_name ' gmStruct']);
    eval(['save ' output_name ' gmStruct']);

%-------------------------------------------------------------------------------

function [ SIGMA_MATR ] = get_covariance_matrix(K, SIGMA0, ALPHA1, EPS)
    SIGMA_MATR = zeros(K, K);
    for k=1:K
        for j=1:K
            SIGMA_MATR(k,j) = SIGMA0*ALPHA1(k)*ALPHA1(j);
            if k == j
                SIGMA_MATR(k,j) = SIGMA_MATR(k,j) + max(EPS(k), 1.e-9);
            end
        end
    end

%-------------------------------------------------------------------------------

function [PI, SIGMA, SIGMA_DET, SIGMA_INV, EPS, status] = ...
         preprocessing(K, M, x0, MU, SIGMA, options)

    status = 0;

    % Do pre-processing
    if options.model == 1
        SIGMA_DET = zeros(1, M);
        PI        = zeros(1, M);
        SIGMA_INV = zeros(K, K*M);
        EPS       = zeros(1,K) + 1;
        PI(1:M)   = x0(1+(0:(M-1))*6);

        for m = 1:M
            SIGMA1 = SIGMA(:,((m-1)*K+1):(m*K));
            SIGMA_DET(m) = det(SIGMA1);
            SIGMA(:,((m-1)*K+1):(m*K)) = SIGMA1;
            SIGMA_INV(1:K,((m-1)*K+1):(m*K)) = inv(SIGMA1);
        end
        disp(['In preprocessing: SIGMA_DET=' num2str(SIGMA_DET)]);
    else
        [PI, ALPHA, ~, ~, SIGMA_DET, SIGMA_INV, EPS] = get_apparent_variables(x0, K, options);
    end

    if options.verbose >= 2
        disp(['After pre-processing:']);
        MU
        SIGMA
        if options.model == 2
            ALPHA
        end
        disp('size(SIGMA_DET)=');
        size(SIGMA_DET)
        SIGMA_DET
        SIGMA_INV
    end

%-------------------------------------------------------------------------------

function [best_Y,best_Q,best_PI,best_MU,best_SIGMA,best_SIGMA_INV,best_xpar,best_LL,em_iter] = ...
              perform_iterations(x0, LL0, PI, MU, SIGMA, SIGMA_INV, EPS, DATA, Y_true, Alphabet, options)    
    K         = size(MU, 1);
    M         = size(PI, 1);
    N         = size(DATA, 1);
    best_Y    = zeros(1, N);
    best_Q    = zeros(1, N);
    best_GAMMA= zeros(N, M);
    prev_xpar = Inf*ones(1, int32(length(x0)));
    prev_LL   = -Inf;
    xpar      = x0;
    LL        = -Inf;
    best_LL   = LL0;
    best_PI   = PI;
    em_iter = 0;
    stop = 0;
    disp('Start EM iterations ...');
    while ~stop && em_iter <= options.maxIter

        em_iter = em_iter + 1;
%       disp(['Before Expectation_step (em_iter=' num2str(em_iter) ')']);
%       PI
%       MU
        [Y, Q, GAMMA] = Expectation_step(PI, MU, SIGMA, SIGMA_INV, DATA, Alphabet, options);
%       disp(['After Expectation step (em_iter=' num2str(em_iter) ')']);
        counts = get_label_counts(Y);
        if options.verbose
            disp(['Y_unique=', num2str(unique(Y)) ' counts=' num2str(counts)]);
        end
%       disp(['Before Maximization_step (em_iter=' num2str(em_iter) ')']);
%       PI
%       MU
%       SIGMA
%       SIGMA_INV
        [PI,MU,SIGMA,SIGMA_INV,DET] = Maximization_step(GAMMA, DATA, options);
        [PI,MU,SIGMA,SIGMA_INV,DET] = sort_components(PI,MU,SIGMA,SIGMA_INV,DET);
%       disp(['After Maximization_step (em_iter=' num2str(em_iter) ')']);
%       PI
%       MU
%       SIGMA
%       DET
%       SIGMA_INV
        try 
            [MU0,SIGMA0,EPS,stop] ...
                      = TNC_GM_RefineModelParameters(MU,SIGMA,options, stop);
        catch
            return;
        end
%       disp(['After Refinement (em_iter=' num2str(em_iter) ')']);
%       SIGMA
%       DET      
%       SIGMA_INV  
        if stop
            break;
        end 
        if options.verbose >= 2
            MU0
            SIGMA0
            EPS
            SIGMA
            DET
            SIGMA_INV
        end
        prev_xpar = xpar;
        prev_LL   = LL;
        xpar      = get_intrinsic_variables(PI,MU0,SIGMA0,EPS,MU,options);

        % 1 = use mvnpdf; 2 = use an efficient implementation
        MVN_VERSION = 2;
        LL        = TNC_GM_GaussMixLogLikelihood(PI, MU, SIGMA, DET, SIGMA_INV, DATA, MVN_VERSION);
        if best_LL   < LL && imag(LL) == 0 
           best_LL   = LL;
           best_xpar = xpar;
           best_PI   = PI;
           best_MU   = MU;
           best_SIGMA     = SIGMA;
           best_SIGMA_INV = SIGMA_INV;
           best_Y    = Y;
           best_Q    = Q;
           best_GAMMA= GAMMA;
        end

        disp(['Iteration# ' num2str(em_iter) ' M= ' num2str(options.M) ...
              ' method= ' options.method ' model=' num2str(options.model) ...
              ' abs(LL - prev_LL)=' num2str(abs(LL - prev_LL))]);
        disp(['counts=' num2str(counts)]);

        if isnan(LL)
            break;
        end

        if options.verbose
            disp(['prev_LL=' num2str(prev_LL) ' LL=' num2str(LL) ...
                  ' TolFun=' num2str(options.TolFun)]);
            disp(' ' );
        end

        if abs(LL - prev_LL) <=options.TolFun
            stop = 1;
        end

        if options.verbose >= 2
            MU_MATR
        end

        disp(['LL=' num2str(LL) ' best LL=' num2str(best_LL) ]);
        disp(['PI=' num2str(PI) ' DET=' num2str(DET)]);
        disp(' ');
        if options.verbose >= 2
            disp('     xpar= ');
            for i=1:int32((length(xpar)-(options.model-1)*K)/6)
                disp(['           ' num2str(xpar(int32(1:6)+(i-1)*6))]);
            end
            if options.model == 2
                disp(['           ' num2str(xpar((length(xpar)-K+1):length(xpar)))]);
            end
            disp('best xpar= ');
            for i=1:((length(best_xpar)-(options.model-1)*K)/6)
                disp(['           ' num2str(best_xpar(int32(1:6)+(i-1)*6))]);
            end
            if options.model == 2
                disp(['           ' num2str(best_xpar((length(best_xpar)-K+1):length(best_xpar)))]);
            end
            disp(' ');
        end

        if options.fakeData
            estimate_coord_accuracy(K, M, xpar, options);
            disp(' ');
        end
    end

    if options.dispOn > 0
        disp(['size(SIGMA)=' num2str(size(SIGMA))]);
        M = size(MU, 2);
        m = 1;
        SIGMA1 = SIGMA(:,(K*(m-1)+1):(K*m));
        colormap('hot');
        imagesc(SIGMA);
        colorbar;
    end

    disp(' ');
    if imag(LL) == 0 && abs(LL - prev_LL) <= options.TolFun
        disp(['Iterations converged with change in LL= ' num2str(abs(LL - prev_LL)) ...
              ' < TolFun= ' num2str(options.TolFun) ' imag(LL)=' num2str(imag(LL))]);
    end

%-------------------------------------------------------------------------------

function [xpar, LL, Y, Q, PI, MU, SIGMA, SIGMA_INV, num_iter] = ...
    ExpectationMaximization(x0, MU, SIGMA, DATA, Y_true, Alphabet, options)      

    num_iter = 0;
    K = int32(size(MU,1));           
    M = int32(options.M);            
    N = size(DATA, 1); 
    xpar = [];
    Y = zeros(1, N);
    Q = zeros(1, N);

    if options.verbose
        disp('After initial guess:');
        MU
        SIGMA    
    end

    % Do pre-processing 
    disp('Preprocessing...');
    [PI, SIGMA, DET, SIGMA_INV, EPS, status] = preprocessing(K, M, x0, MU, SIGMA, options);
    disp(['Preprocessing status=' num2str(status)]);
%   if status > 0
%       LL = -Inf; 
%       return;
%   end

    LL0 = -Inf;
    if options.verbose
        disp('After pre-processing:');
        DET
        MU
        disp(' ' );
    end

    % Do iterations
    % NOTE: don't use an additional condition 
    %       || norm(xpar - prev_xpar) > options.TolX  !!!
    disp('Performing iterations...');
    [Y,Q,PI,MU,SIGMA,SIGMA_INV,xpar,LL,num_iter] = ...
          perform_iterations(x0, LL0, PI, MU, SIGMA, SIGMA_INV, EPS, DATA, ...
                             Y_true, Alphabet, options);
    [PI,MU,SIGMA,SIGMA_INV,DET] = sort_components(PI,MU,SIGMA,SIGMA_INV,DET);

    disp(' ');
    disp(['Final LL=' num2str(LL) ' (M=' num2str(M) ', node=' num2str(options.iRepl) ' num_iter=' num2str(num_iter) ')']);
    disp('     xpar= ');
    for i=1:int32((length(xpar)-(options.model-1)*K)/6)
        disp(['           ' num2str(xpar(int32(1:6)+(i-1)*6))]);
    end
%   if options.verbose
        disp(['      PI  =' num2str(PI)]);
        disp( '      MU  =' );
        MU
        disp( '      SIGMA=');
        SIGMA
%   end
    disp(' ');

%-------------------------------------------------------------------------------

function [PI, ALPHA, MU_MATR, SIGMA, SIGMA_DET, SIGMA_INV, EPS] = ...
         get_apparent_variables(X, K, options)

    LX   = int32(length(X));
    M    = int32((length(X)-(options.model-1)*K)/6);

    SIGMA_DET = zeros(1, M);
    PI        = zeros(1, M);
    MU0       = zeros(1, M);
    SIGMA0    = zeros(1, M);
    XC        = zeros(1, M);
    YC        = zeros(1, M);
    ZC        = zeros(1, M);

    PI(1:M)     = X(1+(0:(M-1))*6);
    MU0(1:M)    = X(2+(0:(M-1))*6);
    SIGMA0(1:M) = X(3+(0:(M-1))*6);
    XC(1:M)     = X(4+(0:(M-1))*6); % cell X-coordinate
    YC(1:M)     = X(5+(0:(M-1))*6); % cell Y-coordinate
    ZC(1:M)     = X(6+(0:(M-1))*6); % cell Z-coordinate
    EPS = zeros(1,K) + 1;
    if options.model == 2       
        EPS= X(int32(6*M) + int32(1:K));
    end

    MU_MATR   = zeros(K, M);
    SIGMA     = zeros(K, K*M);
    SIGMA_INV = zeros(K, K*M);
 
    iE = [1:K];
    if K == 8
        [XE, YE, ZE] = TNC_GM_ElectrodeCoordinates('Buzsaki64', iE);
    else
        [XE, YE, ZE] = TNC_GM_ElectrodeCoordinates('tetrode', iE);
    end
    ALPHA_VECT = []; % vector of length M*K
    for m=1:M
        dist = [];
        for k=1:K
            ce_dist = sqrt((XC(m)-XE(k))^2+(YC(m)-YE(k))^2+(ZC(m)-ZE(k))^2);
            dist = [dist ce_dist];
        end
        mean_dist = mean(dist);         
%       disp(['mean_dist=' num2str(mean_dist) ' dist=' num2str(dist)]);
        ALPHA_ADD = [];
        for k=1:K
            my_alpha = mean_dist./dist(k);
            ALPHA_ADD = [ALPHA_ADD my_alpha];
%           disp(['    k=' num2str(k) ' my_alpha=' num2str(my_alpha) ' ALPHA_ADD=' num2str(ALPHA_ADD)]);
        end
        ALPHA_VECT = [ALPHA_VECT ALPHA_ADD];
    end
    ALPHA = reshape(ALPHA_VECT, K,M);

    % Matrix of means
    for m = 1:M
        MU_MATR(:,m) = MU0(m).*ALPHA(:,m);
    end

    % Covariance matrixes for all components = M matrices of size K x K
    for m = 1:M
        SIGMA(:,((m-1)*K+1):(m*K)) = get_covariance_matrix(K,SIGMA0(m),ALPHA(:,m),EPS);
%       disp(['m=' num2str(m) ' det(SIGMA)=' det(SIGMA(:,((m-1)*K+1):(m*K)))]);
        SIGMA_DET(m) = det(SIGMA(:,((m-1)*K+1):(m*K)));
        SIGMA_INV(1:K,((m-1)*K+1):(m*K)) = inv(SIGMA(1:K,((m-1)*K+1):(m*K)));
    end

%-------------------------------------------------------------------------------

% Compute "responsibilities" of different Gaussian components
%                            for each data point
% GAMMA is N x M matrix of "responsibilities"
%                       of different components for a given data point
% NOTE: use formula (9.13) from:
%       Christopher M. Bishop, "Pattern Recognition and Machine Learning"
%       Springer, 2006
%

function [Y, Q, GAMMA ] = Expectation_step(PI, MU, SIGMA, SIGMA_INV, DATA, Alphabet, options)
    N = (size(DATA, 1));
    M = (length(PI));
    K = (size(MU,   1));

    Y     = char('-'*ones(1, N));
    Q     = zeros(1, N);
    GAMMA = zeros(N, M);
    SIGMA_DET = zeros(1, M); 
    arg   = zeros(1, M);
    for m=1:M
        SIGMA_DET(m) = det(SIGMA(:,(K*(m-1)+1):(K*m)));
    end
    if options.verbose >= 2
        disp('In Expectation_step:');
        SIGMA
        SIGMA_DET
        SIGMA_INV
    end

    % Compute 'responsibilities' of different clusters for given datapoint
    % Determine class labels Y 
    for n=1:N
        Xn = DATA(n,:);
        Num = zeros(1,M);
        for m=1:M
            MU1        = MU(:,m)';
            SIGMA1_INV = SIGMA_INV(:,(K*(m-1)+1):(K*m));
            DET1       = SIGMA_DET(m);
            Num(m) = PI(m)*TNC_GM_MultiVarGauss(Xn, MU1, SIGMA1_INV, DET1); 
            SIGMA1 = SIGMA(:,(K*(m-1)+1):(K*m));
%           SIGMA1
%           det(SIGMA1)
%           Num(m) = PI(m)* mvnpdf(Xn,MU1,SIGMA1);
%           arg(m) = (Xn-MU1)*SIGMA1_INV*(Xn-MU1)';
%           disp(['n=' num2str(n) ' m=' num2str(m) ' MU1=' num2str(MU1) ' Xn=' num2str(Xn) ' arg=' num2str(arg) ' Num=' num2str(Num(m))]);
%           SIGMA1 = SIGMA(:,(K*(m-1)+1):(K*m));
%           Num(m) = PI(m)*mvnpdf(Xn, MU1, SIGMA1);
        end
        Denom = sum(Num);
%       disp(['Xn=' num2str(Xn/1.e5) ' Num=' num2str(Num)]);
        if ~isnan(Denom)
            for m=1:M 
                if ~isnan(Num(m)/Denom)
                    [val ind] = max(Num);
                    GAMMA(n,m) = Num(m)/Denom;
                end
            end
        end
        des_GAMMA = sort(GAMMA(n,:), 'descend');
        sum_GAMMA = sum(GAMMA(n,:));
        best_m    = min(find(GAMMA(n,:) == des_GAMMA(1)));
        Y(n)      = Alphabet(best_m);
        if M == 1
            Q(n)  = 100;
        else
            Q(n)  = round(100*((des_GAMMA(1)-des_GAMMA(2))/sum_GAMMA));
        end
%       disp(['n=' num2str(n) ' GAMMA=' num2str(GAMMA(n,:)) ' Y=' num2str(Y(n)) ' Q=' num2str(Q(n))]);
%       disp(['M=' num2str(M) ' Xn=' num2str(Xn) ' arg=' num2str(arg)]);
    end

%-------------------------------------------------------------------------------

% Compute next approximation to parameters of Gaussian Mixture MOdel
% using as input the "responsibilities" of different Gaussian components
% for each data point
%
% PI is M-vector of mixing coefficients
% MU is K x M matrix of mean values
% SIGMA is K x K*M matrix of covariances
%
% NOTE: use formulas (9.17), (9.18), (9.19) and (9.22) from:
%       Christopher M. Bishop, "Pattern Recognition and Machine Learning"
%       Springer, 2006
%

function [PI, MU, SIGMA, SIGMA_INV, DET]  = Maximization_step(GAMMA, DATA, options)
    N = (size(DATA,  1));
    M = (size(GAMMA, 2));
    K = (size(DATA,  2));

    Nm    = zeros(1, M);
    PI    = zeros(1, M);
    MU    = zeros(K, M);

    % Compute mean values
    for m=1:M
        for n=1:N
                Xn      = DATA(n,:);
                Nm(m)   = Nm(m)   + GAMMA(n,m);  
                if isnan(GAMMA(n,m))
                    disp(['GAMMA=NaN for n=' num2str(n) ' m=' num2str(m)]);
                else
                    for k=1:K
                        MU(k,m) = MU(k,m) + GAMMA(n,m)*Xn(k);
                    end
                end
        end
        MU(:,m) = MU(:,m)./Nm(m);
    end

    for m=1:M
        PI(m) = Nm(m) / sum(Nm);
    end

    % Compute covariance matrix
    if options.model == 1
        [SIGMA,SIGMA_INV,DET]=UnconstrainedMaximization(K,M,N,Nm,MU,GAMMA,DATA);
    elseif options.model == 2
        [SIGMA,SIGMA_INV,DET]=  ConstrainedMaximization(K,M,N,Nm,MU,GAMMA,DATA);
    end

%-------------------------------------------------------------------------------

% Sort Gaussian components in the reverse order of PI

function[new_PI,new_MU,new_SIGMA,new_SIGMA_INV,new_DET] = sort_components(...
             PI,    MU,    SIGMA,    SIGMA_INV,    DET)
    % Sort PI and MU in reverse order of PI
    K = size(MU, 1);
    M = length(PI);
    [new_PI, I]   = sort(PI, 'descend');
    new_DET       = DET(I);
    new_MU        = zeros(K, M);
    new_SIGMA     = zeros(K, K*M);
    new_SIGMA_INV = zeros(K, K*M);
    for m=1:length(I)
        new_MU(:,m) = MU(:, I(m));
        new_SIGMA(    :,(K*(m-1)+1):(K*m)) = SIGMA(    :,(K*(I(m)-1)+1):(K*I(m)));
        new_SIGMA_INV(:,(K*(m-1)+1):(K*m)) = SIGMA_INV(:,(K*(I(m)-1)+1):(K*I(m)));
    end

%-------------------------------------------------------------------------------

function [SIGMA,SIGMA_INV,DET]=UnconstrainedMaximization(K,M,N,Nm,MU,GAMMA,DATA)
    SIGMA     = zeros(K, K*M);
    SIGMA_INV = zeros(K, K*M);
    DET       = zeros(1, M);
    for m=1:M
        for n=1:N
            Xn  = DATA(n,:);
            for k1=1:K
                for k2=k1:K
                    SIGMA(k1,K*(m-1)+k2) = SIGMA(k1,K*(m-1)+k2) + GAMMA(n,m) ...
                                          *(Xn(k1)-MU(k1,m))*(Xn(k2)-MU(k2,m));
                    if k2 ~= k1
                        SIGMA(k2,K*(m-1)+k1) = SIGMA(k1,K*(m-1)+k2);
                    end
                end
            end
        end

        % Regularized covariance matrix
        SIGMA(:,(K*(m-1)+1):(K*m)) = SIGMA(:,(K*(m-1)+1):(K*m))./Nm(m);
        SIGMA_INV(:,(K*(m-1)+1):(K*m)) = inv(SIGMA(:,(K*(m-1)+1):(K*m)));
        DET(m) = det(SIGMA(:,(K*(m-1)+1):(K*m)));
    end

%-------------------------------------------------------------------------------

function [SIGMA,SIGMA_INV,DET] = ConstrainedMaximization(K,M,N,Nm,MU,GAMMA,DATA)
    SIGMA     = zeros(K, K*M);
    SIGMA_INV = zeros(K, K*M);
    DET       = zeros(1, M);
    ALPHA     = zeros(K, M);
    for m=1:M
        MU0 = mean(MU(:,m));
        ALPHA(:, m) = MU(:, m)./MU0;
        Sm = 0.;
        Gm = 0.;
        Dm = 0.;
        Am = 0.;
        for k=1:K
            Sm = Sm + ALPHA(k,m)^2;
        end
        for n=1:N
            Gm = Gm + GAMMA(n, m);
            Xn  = DATA(n,:);
            for k1=1:K
                Dm = Dm + GAMMA(n,m)*(Xn(k1)-MU(k1,m))^2;
                for k2=1:K
                    Am = Am + GAMMA(n,m)*ALPHA(k1,m)*ALPHA(k2,m)...
                                   *(Xn(k1)-MU(k1,m))...
                                               *(Xn(k2)-MU(k2,m));
                end
            end
        end
        
        rho_m     = 0.;
        kapa_m    = 0.;
        epsilon_m = 0;
        sigma_m   = 0;
%       disp(['m=' num2str(m) ' Am=' num2str(Am) ' Dm=' num2str(Dm) ' Gm=' num2str(Gm) ' Sm=' num2str(Sm)]);
        num = (Am*K - Dm*Sm);
        den = (Dm*Sm - Am);
%       disp(['num=' num2str(num) ' den=' num2str(den)]);
        if Am <= 0
            kapa_m = K/Dm;
        else
            kapa_m0 = K*Gm/Dm;
            kapa_m1 = Gm*Sm/Am;
            if kapa_m0 < kapa_m1
                kapa_m = kapa_m0;
            else
                rho_m  = (Am*K - Dm*Sm)/(Dm*Sm - Am)/Sm;
                kapa_m = Gm * (K + Sm*rho_m)/Dm;
            end
        end            
        epsilon_m = 1/kapa_m;
        sigma_m   = rho_m * epsilon_m;
 
        % Compute SIGMA
        for k1=1:K
            for k2=k1:K
                SIGMA(k1,K*(m-1)+k2)     = sigma_m*ALPHA(k1,m)*ALPHA(k2,m);
                SIGMA_INV(k1,K*(m-1)+k2) = -kapa_m*ALPHA(k1,m)*ALPHA(k2,m)...
                                           *rho_m/(1.+rho_m*Sm);
                if k2 ~= k1
                    SIGMA(k2,K*(m-1)+k1)     = SIGMA(k1,K*(m-1)+k2);
                    SIGMA_INV(k2,K*(m-1)+k1) = SIGMA_INV(k1,K*(m-1)+k2);
                else
                    SIGMA(k1,K*(m-1)+k2) = SIGMA(k1,K*(m-1)+k2) + epsilon_m;
                    SIGMA(k1,K*(m-1)+k2) = SIGMA(k1,K*(m-1)+k2) + kapa_m;
                end
            end
        end     
        DET(m) = (1. + rho_m*Sm)*epsilon_m^K;      
    end

%-------------------------------------------------------------------------------

function X = get_intrinsic_variables(PI, MU0, SIGMA0, EPS, MU, options)

    M  = (size(MU0, 2));
    K  = (size(EPS, 2));
    LX = int32(6*M + (options.model - 1)*K);
    X  = zeros(1, LX);
    
    XC = zeros(1, M);
    YC = zeros(1, M);
    ZC = zeros(1, M);

    % Compute vector MUO
    if K == 8
        [XE, YE, ZE] = TNC_GM_ElectrodeCoordinates('Buzsaki64', [1:K]);
    else
        [XE, YE, ZE] = TNC_GM_ElectrodeCoordinates('tetrode', [1:K]);
    end
    for m=1:M
        if K == 4
            tetrode = reshape([XE, YE, ZE], K, 3);
            [XC(m), YC(m), ZC(m)] = TNC_GM_Triangulation(MU(:,m)', tetrode', options);
            if options.verbose
                disp(['m=' num2str(m) ' MU=' num2str(MU(:,m)') ' XC(m), YC(m), ZC(m)=']);
                [XC(m), YC(m), ZC(m)] 
            end
        elseif K == 8
            tetrode1 = reshape([XE([1 3 5 7]), YE([1 3 5 7]), ZE([1 3 5 7])], 4, 3);
            tetrode2 = reshape([XE([2 4 6 8]), YE([2 4 6 8]), ZE([2 4 6 8])], 4, 3);
            [XC1, YC1, ZC1] = TNC_GM_Triangulation(MU([1 3 5 7],m)', tetrode1', options);
            [XC2, YC2, ZC2] = TNC_GM_Triangulation(MU([1 3 5 7],m)', tetrode1', options);
            XC(m) = (XC1 + XC2)/2;
            YC(m) = (YC1 + YC2)/2;
            ZC(m) = (ZC1 + ZC2)/2;
        else
            disp('TNC_GM_SortSpikes: K must be a multiple of 4');
            return
        end
        if options.verbose
            disp(['m=' num2str(m) ' MU=' num2str(MU(:,m)') ...
                  ' power=' num2str(options.power) ...
                  ' X,Y,Z=' num2str([XC(m), YC(m), ZC(m)])]);
        end
    end
    X(1+(0:(M-1))*6) = PI;
    X(2+(0:(M-1))*6) = MU0;
    X(3+(0:(M-1))*6) = SIGMA0;
    X(4+(0:(M-1))*6) = XC;
    X(5+(0:(M-1))*6) = YC;
    X(6+(0:(M-1))*6) = ZC;
    if options.model == 2       
        X(6*M + (1:K)) = EPS(1:K);
    end

%-------------------------------------------------------------------------------

function estimate_coord_accuracy(K, M, X, options)
    xc_true = options.coord.x_neurons;
    yc_true = options.coord.y_neurons;
    zc_true = options.coord.z_neurons;
    for m=1:M
        xc = X(4+(m-1)*6);
        yc = X(5+(m-1)*6);
        zc = X(6+(m-1)*6);
        best_dist = inf;
        best_accuracy = inf;
        best_xc  = inf;
        best_yc  = inf;
        best_zc  = inf;
        for i=1:length(options.coord.x_neurons)
            dist = sqrt((xc-xc_true(i))^2 + ...
                        (yc-yc_true(i))^2 + ...
                        (zc-zc_true(i))^2);      
            if (best_dist > dist)
                best_dist = dist;
                best_xc = xc_true(i);
                best_yc = yc_true(i);
                best_zc = zc_true(i); 
            end
        end 
        coord_error = [abs(xc-best_xc) abs(yc-best_yc) abs(zc-best_zc)];
%       disp(['m=' num2str(m) ' true xc,yc,zc=' num2str([best_xc best_yc best_zc]) ...
%             ' coord. error=' num2str(coord_error)...
%             ' % error=' num2str(100*max(coord_error)/max([abs(best_xc) abs(best_yc) abs(best_zc)]))]);
    end    

%-------------------------------------------------------------------------------

function [xpar, LL, M, Y] = KlustaKwik(DATA, options)
    K = size(DATA, 2);
    tonic_data = getenv('TONIC_DATA');
%   prefixPath = [ tonic_data '/' my_options.prefixName ];
    prefixPath = fullfile(tonic_data, options.prefixName);
    disp(['prefixPath=' prefixPath ]);
    input_tfile = create_input_textfile(DATA, options);
    disp(['input_file=' input_tfile ]); 
    output_file = run_kk(input_tfile);       
    Y = read_output_file(output_file);
    M = length(unique(Y));
    output_kk_gm_data(M, options);
    disp(' ');
    xpar = [];
    LL   = 0;

%-------------------------------------------------------------------------------

function output_kk_gm_data(M, options)
    gmStruct.program = 'TNC_GM_SortSpikes.m';
    gmStruct.best_mc = M;
    gmStruct.method  = 'KlustaKwik';
    output_name = strcat(options.prefixName,'_kk_gm.', num2str(options.iRepl),'.mat');
    disp(['outputName=' output_name]);
    disp(['save ' output_name ' gmStruct']);
    eval(['save ' output_name ' gmStruct']);

%-------------------------------------------------------------------------------

function fetFilePrefix = create_input_textfile(DATA, options)
    K = size(DATA, 2);
    N = size(DATA, 1);
    totChar = numel(options.ftFile);
    if 0
        fileName =  dir([options.ftFile(1:totChar-4) '.mat'])
        disp(['fileName.name=' fileName.name]);
        prefixName = fileName.name(1:(length(fileName.name)-4))
    else
        fileName = options.ftFile;
        disp(['fileName=' fileName ]);
        prefixName = fileName(1:(length(fileName)-4))
    end
    tonic_data = getenv('TONIC_DATA')
    fetFilePrefix = fullfile(tonic_data, prefixName);
    fetFileName = fullfile(tonic_data, strcat(prefixName,'.fet.1' ));
    disp(['fetFileName=' fetFileName ' fetFilePrefix=' fetFilePrefix]);
    disp(' ');
    try
        fid = fopen(char(fetFileName), 'w');
        fprintf(fid, '%d\n', K);
        for n=1:N
            for k=1:K
                fprintf(fid, '%f ', DATA(n, k));
            end
            fprintf(fid, '\n');
        end
        fclose(fid);
    catch
        disp(['Cannot open file ' fetFileName ]);
    end

%-------------------------------------------------------------------------------

function output_file = run_kk(prefixPath)          
    tonic_data = getenv('TONIC_DATA');
    executable_path = [ tonic_data '/KlustaKwik' ];
    command = [ executable_path ' ' prefixPath ' 1'];
    disp(['kk unix_command=' command ]);
    status = unix(command);
    output_file = [ prefixPath '.clu.1' ];
    disp(['output_file=' output_file ]); 

%-------------------------------------------------------------------------------

function Y = read_output_file(output_file)
    Y = [];
    fid = fopen(output_file,'r'); %# open csv file for reading
    while ~feof(fid)
        line = fgets(fid); %# read line by line
        y = sscanf(line,'%d'); %# sscanf can read only numeric data :(
        Y = [Y y];
    end
    fclose(fid);

%-------------------------------------------------------------------------------

function [map,classError,Y_best] = compare_class_labels(K, Y_true, Y_comp, ...
                                                        Alphabet, verbose)
    map = ones(1, K); % best map for color(Y_true) => color(Y_comp)

    % Compute a map of labels by maximizing the map score
    for k=1:K
        best_Error = Inf;
        best_i = 0;
        temp_Y_true = Y_true;
        temp_Y_true(Y_true ~= Alphabet(k)) = '0'; 
        for i=1:K
            temp_Y_comp = Y_comp;
            temp_Y_comp(Y_comp == Alphabet(i)) = Alphabet(k);
            temp_Y_comp(Y_comp ~= Alphabet(i)) = '1';
            Error = TNC_GM_AlignClassLabels(temp_Y_comp, temp_Y_true, 0);
            if   best_Error > Error.tot
                 best_Error = Error.tot;
                 best_i = i; 
            end
        end   
        map(best_i) = k; 
    end
    disp(['map=' num2str(map)]);

    Y_best = Y_comp;
    for i=1:K
        Y_best(Y_comp == Alphabet(i)) =  Alphabet(map(i));    
    end
    classError = TNC_GM_AlignClassLabels(Y_best, Y_true, 1);

