%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function TNC_MM_SortSpikes(ftFile,iSeg,iShank,M,varargin)      
    % ftFile, iSeg, iShank and M = required positional arguments
    % iRepl                      = optional, position-specific argument
    % lowerM, verbose            = optional parameters
    %
    % Example of usage:
    % TNC_MM_SortSpikes('WT4RD1LSTR002_2012_02_21_ft.mat','2','3',3,'ms',16,'verbose',1)
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
        disp('Input file is not of ft type');
        return
    end

    if options.verbose
        disp(['ftFile=' options.ftFile ' iSeg=' num2str(options.iSeg) ' iShank=' num2str(options.iShank) ...
              ' M=' num2str(options.M) ' iRepl=' num2str(options.iRepl) ' method=' options.method ...
              ' guess_method=' options.guess_method ' dispOn=' num2str(options.dispOn) ...
              ' fakeData=' num2str(options.fakeData) ' debug=' num2str(options.debug)]);
    end

    % Read input data and split it into 2 subsets:
    % DATA_f = data from the specified fold; and DATA_nf = complimentary data
    try
        [K,ftData,TS,DATA_nf,DATA_f,Y_true,options,status] = extract_data(options);
    catch
        disp(['Failed to read data for segment ' num2str(options.iSeg) ' shank ' num2str(options.iShank) ]);
        return;
    end

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

    if strcmp(options.method, 'KK')
        options.maxIter = 1;
    end
    LL      = -Inf;
    LLf     = -Inf;
    Iter    = 1;
    success = 0;
    disp(['options.M=' num2str(options.M) ' options.verbose='  num2str(options.verbose)]);
    while ~success && Iter < 2                 
        Iter = Iter + 1;
%       try
            % Initial guessi/clustering that will be used 
            % in the training / parameter estimation procedure
            [PI, MU, SIGMA] = TNC_MM_InitialGuess(DATA_nf, Iter, options);
            if options.verbose
                disp(' ');
                disp('Initial guess:');
                PI
                MU
                SIGMA
                M = size(MU,2);
                DET_SIGMA = zeros(1,M);
                for m=1:M
                    SIGMA1   = SIGMA(:,(K*(m-1)+1):(K*m));
                    DET_SIGMA(m) = det(SIGMA1);
                end
                disp(['DET_SIGMA=' num2str(DET_SIGMA)]);
                disp(' ');
            end
            % Traing the parameters PI, MU and SIGMA on DATA_nf
            % (get the likelihood LL, or LL_nf, as a by-product of the training)
            % and apply the training results to compute the likelihood LLf on DATA_f
            [LL, LLf, Y, Q, TSM, PI, MU, SIGMA, UC, status] = ...
                compute_likelihoods(PI, MU, SIGMA, TS, DATA_nf, ...
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
%       catch
%           if ~strcmp(options.method, 'kk')
%               disp(['status=' num2str(status) ' replica=' num2str(options.iRepl) ...
%                     ' LL=' num2str(LL) ' LLf=' num2str(LLf)]);
%               disp(['compute_likelihoods FAILED in attempt# ' num2str(Iter) ' ...']);
%           end
%       end
    end
    if ~strcmp(options.method, 'kk') && success    
        disp(' ');
        disp(['# Attempts made=' num2str(Iter)]);
        disp(['Output to mm_data file: LL=' num2str(LL) ' LLf='  num2str(LLf) ...
              ' prefix=' prefixName ' iFold=' num2str(options.iFold)]);
        output_mm_data(LL, LLf, K, Y_true, Y, Q, TSM, classError, PI, MU, SIGMA, ...
                       UC, prefixName, ftData.featStruct, options);
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
    disp('Usage: TNC_MM_SortSpikes(ftFile, iSeg, iShank, M [,iRepl] [,parameters])');
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
    disp('                   default = kmmeans; alternative = random_cell');
    disp('    iFold        - index of the fold to be tested in cross-validation');
    disp('                   (>= 0; <= numFolds, default = 1');
    disp('    maxIter      - max number of iterations to be performed');
    disp('                   (default = Inf)');
    disp('    method       - computational method to be used');
    disp('                   default=GEM (for Expectation Maximization)');
    disp('                   other available options: ms (for MultiStart)');
    disp('                   and KK (for KlustaKwick)');
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
    if strcmp(options.guess_method, 'random_cell') && options.iRepl <= 0
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
%       disp(['i=' num2str(i) ' my_str=' my_str ' len=' num2str(len)]);
        paramNames_short{i} = my_str(1:len);
    end
    disp(['paramNames_short='])
    paramNames_short
    num_features = numel(unique(paramNames_short)) -1;

% ---------------------------------------------------------------------

function [K,ftData,TS,DATA_nf,DATA_f,Y_true,options,status] = extract_data(options)
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
        TS = [];
        DATA_nf=[];
        DATA_f =[];
        Y_true =[''];
        return
    end 
    K = 0;
    TS = [];
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
    disp(['iSeg=' num2str(options.iSeg) ' iShank=' num2str(options.iShank) ...
          ' npar=' num2str(npar) ' num_features=' num2str(num_features)]);
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
    TS    = ftData.featStruct.seg(sg).shank(sh).ts;     
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

function [LL, LLf, Y, Q, TSM, PI, MU, SIGMA, UC, status] = ...
          compute_likelihoods(PI, MU, SIGMA, TS, DATA_nf, DATA_f, ...
                              Y_true, Alphabet, my_options)
    K = size(MU, 1);
    Y = [];
    Q = [];
    TSM = [];

    % Variables: MU, SIGMA^2, X, Y, Z for each cell 

    % Upper and lower bounds
    [lb, ub] = get_bound_constraints(my_options.M, K);

    % Equality constraints
    [Aeq, beq] = get_equality_constraints(my_options.M, K);
 
    % Inequality constraints
    if my_options.M > 1
        [Aineq, bineq] = get_inequality_constraints(my_options.M, K);
    else
        Aineq = [];
        bineq = [];
    end

    if strcmp(my_options.method, 'GEM')
        % Only bound constraints are allowed
        disp('Applying Expectation-Maximization algorithm ...');
        disp(' ');
        my_options.TolFun = 0.1;

        % Train the model/estimate parameters PI, MU, and SIGMA
        % on the non-fold data
        [LL, Y, Q, TSM, UC, PI, MU, SIGMA, status] = ...
            ExpectationMaximization(PI, MU, SIGMA, TS, DATA_nf, ...
            Y_true, Alphabet, my_options);     
          
        disp(['LL=' num2str(LL)]);
        if status < 0 || isinf(-LL)
            LLf = -Inf;
        elseif my_options.iFold == 0 || my_options.numFolds == 1
            LLf = LL;
        else
            % Apply the training results to the fold data:
            % LLf = LL{DATA_f | DATA \ Data_nf}
            LLf = TNC_MM_GaussMixLogLikelihood(PI, MU, SIGMA, DATA_f);
        end
        disp(['my_options.iFold=' num2str(my_options.iFold) ' size(DATA_nf)=' num2str(size(DATA_nf,1)) ...
              ' size(DATA_f)=' num2str(size(DATA_f,1)) ' LL=' num2str(LL) ' LLf=' num2str(LLf)]);
        disp(' ');
    elseif strcmp(my_options.method, 'KK') % Simulated Annealing (only bound constraints!)
        % Only bound constraints are allowed
        disp('Applying KlustaKwik algorithm ...');
        disp(['size(DATA_nf)=' num2str(size(DATA_nf))]);
        PI = [];
        MU = [];
        SIGMA = [];
        GAMMA = [];
        LLf = 0;
        if size(DATA_nf, 2) == 4
            [LL, M, Y] = KlustaKwik(DATA_nf, my_options);
            disp(' ');
            disp(['num_clusters=' num2str(M)]);
        elseif size(DATA_nf, 2) == 8
            DATA1 = DATA_nf(:,[1 3 5 7]);
            [LL, M1, Y] = KlustaKwik(DATA1, my_options);
            DATA2 = DATA_nf(:,[2 4 6 8]);
            [LL, M2, Y] = KlustaKwik(DATA2, my_options);
            disp(' ');
            disp(['num_clusters1=' num2str(M1) ' num_clusters2=' num2str(M2)]);
        end
    else
        disp([ 'Incorrectly specified name of global opt method: ' method ]);
        return;
    end

%-------------------------------------------------------------------------------

function  [lb, ub] = get_bound_constraints(M, K)

    lb = ones(1, 6*M);
    ub = ones(1, 6*M);
    
    lb_one = [0 -2.0  1.e-6  -50 -50   0]; % PI, MU, SIGMA, X, Y, Z
    ub_one = [1  0   10       50 250  50]; % PI, MU, SIGMA, X, Y, Z
    for m=1:M
        lb((1+(m-1)*6):(m*6)) = lb_one;
        ub((1+(m-1)*6):(m*6)) = ub_one;
    end                              

%-------------------------------------------------------------------------------

function [Aeq, beq] = get_equality_constraints(M, K)
    Aeq = zeros(1,M*6);
    Aeq(1,1+(0:(M-1))*6) = 1; % sum of PIs = 1, or A(1,1)*X(1)+A(1,7)*x(7)+...= 1
    beq = 1;

%-------------------------------------------------------------------------------

function [Aineq, bineq] = get_inequality_constraints(M, K)
    Aineq = zeros(M-1,M*6);
    for m=1:(M-1)
        Aineq(m,1+(m-1)*6) =  1;  
        Aineq(m,1+ m   *6) = -1; % MU_m <= MU_m+1 
    end
    bineq = -zeros(1, M-1);

%-------------------------------------------------------------------------------

function output_mm_data(LL,LLf,K,Y_true,Y,Q,TSM,classError,PI,MU,SIGMA, ...
                        UC,targetName,featStruct,options)
    M = length(PI);
    iSeg                    = options.iSeg;
    iShank                  = options.iShank;
    N                       = size(featStruct.seg(iSeg).shank(iShank).params,1);
    mmStruct.program        = 'TNC_MM_SortSpikes.m';
    mmStruct.paramNames     = {'PI', 'MU' 'SIGMA' 'XC' 'YC' 'ZC'};
    mmStruct.num_components = M;
    try
        mmStruct.thr        = featStruct.thr;
    catch
        mmStruct.thr        = 0;
    end
    mmStruct.method         = options.method;
    mmStruct.iFold          = options.iFold;
    mmStruct.numFolds       = options.numFolds;
    mmStruct.K              = K;
    mmStruct.shank(iShank).seg(iSeg).N     = N;                                                   
    mmStruct.shank(iShank).seg(iSeg).Y_true= Y_true;
    mmStruct.shank(iShank).seg(iSeg).Y     = Y;
    mmStruct.shank(iShank).seg(iSeg).Q     = Q;
    mmStruct.shank(iShank).seg(iSeg).TSM   = TSM; % matrix of time stamps, where rows = clusters
    mmStruct.shank(iShank).seg(iSeg).classError = classError;
    mmStruct.shank(iShank).seg(iSeg).PI    = PI;
    mmStruct.shank(iShank).seg(iSeg).MU    = MU;
    mmStruct.shank(iShank).seg(iSeg).SIGMA = SIGMA;
    mmStruct.shank(iShank).seg(iSeg).UC    = UC
    mmStruct.shank(iShank).seg(iSeg).LL    = LL;
    mmStruct.shank(iShank).seg(iSeg).LLf   = LLf;
    if options.fakeData
        mmStruct.coord.x_neurons = featStruct.x_neurons;
        mmStruct.coord.y_neurons = featStruct.y_neurons;
        mmStruct.coord.z_neurons = featStruct.z_neurons;
    end
    disp(['targetName=' targetName]);
    if options.iFold >= 0 && options.numFolds > 1
        output_name = strcat(targetName,'_shank',num2str(options.iShank),...
                             '_seg',num2str(options.iSeg), '_mc',num2str(M),...
                             '_mm_', num2str(options.iFold));
    else
        output_name = strcat(targetName, '_mm');
    end
    disp(['save ' output_name ' mmStruct']);
    eval(['save ' output_name ' mmStruct']);

%-------------------------------------------------------------------------------

function [Y, Q, TSM] = evaluate_parameters(PI, MU, SIGMA, DATA, GAMMA, TS, ...
                            Y_true, Alphabet, options)
   % Initialization
   [N, K] = size(DATA);
    M      = (length(PI));
    Y     = char('-'*ones(1, N));
    Q     = 100*ones(1, N);
    TSM = NaN * ones(M, 1);

    % Loop through all data
    for n=1:N
        Xn = DATA(n,:);
        Num = zeros(1,M);
        for m=1:M
            MU_m    = MU(:,m)';
            SIGMA_m = SIGMA(:,(K*(m-1)+1):(K*m));
            try
                Num(m) = PI(m)*mvnpdf(Xn, MU_m, SIGMA_m);
            catch
                Num(m) = 0;
            end
        end
        Denom = sum(Num);
        if ~isnan(Denom)
            for m=1:M
                if ~isnan(Num(m)/Denom)
                    [val ind] = max(Num);
                    GAMMA(n,m) = Num(m)/Denom;
                end
            end
        end
        desc_GAMMA = sort(GAMMA(n,:), 'descend');
        sum_GAMMA = sum(GAMMA(n,:));
        best_m    = min(find(GAMMA(n,:) == desc_GAMMA(1)));
        try
            Y(n)      = Alphabet(best_m);
        except
            disp(['best_m=' num2str(best_m) ' GAMMA(n,:)=' num2str(GAMMA(n,:))]);
        end

        if sum(~isnan(TSM(best_m,:))) < size(TSM,2)
            ind = min(find(isnan(TSM(best_m,:))));
            TSM(best_m,ind) = TS(n);
        else
            TSM = [ TSM NaN * ones(M, 1) ];
            ind = size(TSM, 2);
            TSM(best_m,ind) = TS(n);
        end
        if M == 1
            Q(n)  = 100;
        else
            Q(n)  = round(100*((desc_GAMMA(1)-desc_GAMMA(2))/sum_GAMMA));
        end
    end


%-------------------------------------------------------------------------------

function [LL, Y, Q, TSM, UC, PI, MU, SIGMA, status] = ...
    ExpectationMaximization(PI, MU, SIGMA, TS, DATA, Y_true, Alphabet, options)      

    K = int32(size(MU,1));           
    M = int32(options.M);            
    N = size(DATA, 1); 
    UC= [];   % unit coordinates
    Y = zeros(1, N);
    Q = zeros(1, N);
    TSM = [];
    LL      = -Inf;
    prev_LL = -Inf;
    status = 0;    % 0 = keep iterations; -1 = failure, LL = -Inf; +1 = success

    % Deterministic annealing parameters
    beta = 0.01;
    eps  = 1.e-2; % in Takekawa et al. (2010) it is = 1e-6
    t = 0;

    % Perform iterations
    disp('Performing iterations...');
    while beta < 1 || ~status
        % Expectation
        t = t + 1;
        if beta < 1
            beta = beta * 1.01;  % in Takekawa et al. (2010) the multiplier is 1.01
            if beta > 1
                beta = 1;
            end 
        end

        [PI, MU, SIGMA] = sort_components(PI, MU, SIGMA);

        if beta < 1
            disp(['t=' num2str(t) ' beta=' num2str(beta) ' PI=' num2str(PI)]);
        end
        [GAMMA, status] = Expectation_step(PI, MU, SIGMA, DATA, beta);

        if status < 0
            LL  = -Inf;
            return;
        end

        % Maximization
        [PI, MU, SIGMA] = Maximization_step(GAMMA, DATA);

        % Evaluation 
        if beta >= 1
            prev_LL = LL;
            LL      = TNC_MM_GaussMixLogLikelihood(PI, MU, SIGMA, DATA);
            if isinf(-LL) || abs(LL - prev_LL)/N < eps
                status = 1;
                break;
            end
            disp(['t=' num2str(t) ' beta=' num2str(beta) ' PI=' num2str(PI) ...
                  ' LL=' num2str(LL) ' incr=' num2str(abs(LL - prev_LL)/N) ...
                  ' eps=' num2str(eps)]);
        end
    end

    % Evaluate parameters
    if ~isinf(-LL)
        [Y, Q, TSM] = evaluate_parameters(PI, MU, SIGMA, DATA, GAMMA, TS, ...
                                          Y_true, Alphabet, options);
        UC     = get_unit_coordinates(PI, MU, SIGMA, options);
    end

    disp(['Final LL=' num2str(LL) ' (M=' num2str(M) ', node=' num2str(options.iRepl) ')']);
    if options.verbose
        disp(['      PI  =' num2str(PI)]);
        disp( '      MU  =' );
        MU
        disp( '      SIGMA=');
        SIGMA
        disp( '      UC=');
        UC
    end
    disp(' ');

%-------------------------------------------------------------------------------

% Compute "responsibilities" of different Gaussian components
%                            for each data point
% GAMMA is N x M matrix of "responsibilities"
%                       of different components for a given data point
% NOTE: use formula (9.13) from:
%       Christopher M. Bishop, "Pattern Recognition and Machine Learning"
%       Springer, 2006
%

function [GAMMA, status] = Expectation_step(PI, MU, SIGMA, DATA, beta)
%                    N x M                      1xM KxM KxK*M  NxK
    [N, K] = size(DATA);
    M      = (length(PI));
    status = 0;

    GAMMA  = zeros(N, M);
    Gauss  = zeros(N, M);
    MSIGMA = zeros(K,K,N);                            % K x K x N
    for m=1:M
        SIGMA_m = SIGMA(:, ((m-1)*K+1):(m*K));
%       disp(['   m=' num2str(m) ' det(SIGMA_m)=' num2str(det(SIGMA_m))]);
%       SIGMA_m
    end
    for m=1:M
        MMU     = MU(:,m) * ones(1,N);                % K x N
        %          1 x N       K x 1
        SIGMA_m = SIGMA(:,(K*(m-1)+1):(K*m));         % K x K

        for k=1:K
            MSIGMA(:,k,:) = SIGMA_m(:,k) * ones(1,N); % K x K x N
        end
%       disp(['   m=' num2str(m) ' det(MSIGMA(:,:,1)=' num2str(det(MSIGMA(:,:,1)))]);
%       MSIGMA(:,:,1)
        try
            Gauss(:, m) = mvnpdf(DATA, MMU', MSIGMA)';    % N x 1
            %                    NxK   NxK   KxKxN
        except
            status = -1;
            return;
        end
    end
    WeightedSumGauss = Gauss.^beta * PI'.^beta;   % column of size N
    for m=1:M
        GAMMA(:,m) = (Gauss(:, m) * PI(m)).^beta ./ WeightedSumGauss;
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
% NOTE: use, but in a vectorized form,
%       formulas (9.17), (9.18), (9.19) and (9.22) from:
%       Christopher M. Bishop, "Pattern Recognition and Machine Learning"
%       Springer, 2006

function [PI, MU, SIGMA]  = Maximization_step(GAMMA, DATA)
%                                             N x M  N x K
    N = (size(DATA,  1));
    M = (size(GAMMA, 2));
    K = (size(DATA,  2));

    Nm = sum(GAMMA);     % 1 x M

    % Compute PI
    PI = Nm / sum(Nm);   % 1 x M

    % Compute MU
    MU = DATA' * GAMMA;  % K x M 
    for m=1:M
        MU(:,m) = MU(:,m)/Nm(m);
    end

    % Compute covariance matrix
    SIGMA     = zeros(K, K*M);
%   disp(['size(DATA)=' num2str(size(DATA)) ' size(GAMMA)=' num2str(size(GAMMA)) ' size(MU)=' num2str(size(MU))]);
    for m=1:M
        MMU = (ones(N,1) * MU(:,m)');                    % N x K
        %          N x 1               1 x K
        SIGMA_m = (DATA' - MMU') * diag(GAMMA(:,m)') * (DATA - MMU) / Nm(m);
        % K x K        K x N            N x N              N x K       1 x 1
        %                           \___________  _______________/
        %                                       \/
        %                                     1 x K
        SIGMA(:,(K*(m-1)+1):(K*m)) = SIGMA_m;
    end

%-------------------------------------------------------------------------------

% Sort Gaussian components in the reverse order of PI

function[new_PI,new_MU,new_SIGMA] = sort_components(PI,    MU,    SIGMA)
    % Sort PI and MU in reverse order of PI
    K = size(MU, 1);
    M = length(PI);
    [new_PI, I]   = sort(PI, 'descend');
    new_MU        = zeros(K, M);
    new_SIGMA     = zeros(K, K*M);
    for m=1:length(I)
        new_MU(:,m) = MU(:, I(m));
        new_SIGMA(    :,(K*(m-1)+1):(K*m)) = SIGMA(    :,(K*(I(m)-1)+1):(K*I(m)));
    end

%-------------------------------------------------------------------------------

function UC = get_unit_coordinates(PI, MU, SIGMA, options)

    M  = size(MU, 2);
    K  = size(MU, 1);
    UC = zeros(M, 3);
    
    XC = zeros(1, M);
    YC = zeros(1, M);
    ZC = zeros(1, M);

    % Compute electrode coords 
    if K == 8
        [XE, YE, ZE] = TNC_MM_ElectrodeCoordinates('Buzsaki64', [1:K]);
    else
        [XE, YE, ZE] = TNC_MM_ElectrodeCoordinates('tetrode', [1:K]);
    end
    for m=1:M
        if K == 4
            tetrode = reshape([XE, YE, ZE], K, 3);
            [XC(m), YC(m), ZC(m)] = TNC_MM_Triangulation(MU(:,m)', tetrode', options);
            if options.verbose
                disp(['m=' num2str(m) ' MU=' num2str(MU(:,m)') ' XC(m), YC(m), ZC(m)=']);
                [XC(m), YC(m), ZC(m)] 
            end
        elseif K == 8
            tetrode1 = reshape([XE([1 3 5 7]), YE([1 3 5 7]), ZE([1 3 5 7])], 4, 3);
            tetrode2 = reshape([XE([2 4 6 8]), YE([2 4 6 8]), ZE([2 4 6 8])], 4, 3);
            [XC1, YC1, ZC1] = TNC_MM_Triangulation(MU([1 3 5 7],m)', tetrode1', options);
            [XC2, YC2, ZC2] = TNC_MM_Triangulation(MU([1 3 5 7],m)', tetrode1', options);
            XC(m) = (XC1 + XC2)/2;
            YC(m) = (YC1 + YC2)/2;
            ZC(m) = (ZC1 + ZC2)/2;
        else
            disp('TNC_MM_SortSpikes: K must be a multiple of 4');
            return
        end
        if options.verbose
            disp(['m=' num2str(m) ' MU=' num2str(MU(:,m)') ...
                  ' power=' num2str(options.power) ...
                  ' X,Y,Z=' num2str([XC(m), YC(m), ZC(m)])]);
        end
    end
    UC(:, 1) = XC';
    UC(:, 2) = YC';
    UC(:, 3) = ZC';

%-------------------------------------------------------------------------------

function estimate_unit_coords_accuracy(K, M, UC, options)
    xc_true = options.coord.x_neurons;
    yc_true = options.coord.y_neurons;
    zc_true = options.coord.z_neurons;
    for m=1:M
        xc = UC(m, 1);  
        yc = UC(m, 2);   
        zc = UC(m, 3);
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

function [LL, M, Y] = KlustaKwik(DATA, options)
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
    output_kk_mm_data(M, options);
    disp(' ');
    LL   = 0;

%-------------------------------------------------------------------------------

function output_kk_mm_data(M, options)
    mmStruct.program = 'TNC_MM_SortSpikes.m';
    mmStruct.best_mc = M;
    mmStruct.method  = 'KK';
    output_name = strcat(options.prefixName,'_kk_mm.', num2str(options.iRepl),'.mat');
    disp(['outputName=' output_name]);
    disp(['save ' output_name ' mmStruct']);
    eval(['save ' output_name ' mmStruct']);

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
            Error = TNC_MM_AlignClassLabels(temp_Y_comp, temp_Y_true, 0);
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
    classError = TNC_MM_AlignClassLabels(Y_best, Y_true, 1);

