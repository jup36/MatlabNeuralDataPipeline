% This function compute initial guess based on one of 2 options:
% 1) km: using k-means algorithm
%    - divide dataset into K subsets using K-means algorithm
%    - for each subset, determine the mean and covariance matrix
%    - use the means to determine the locations of point sources
% 2) randomized:
%    - randomly select the location of point sources
%    - for each point source, determine  vector ALPHA
%    - for each datapoint, determine vector AKPHA1
%    - classify datapoints by choosing the sourec with closest ALPHA

function [ X0, MU0, SIGMA0, DET0 ] = TNC_GM_InitialGuess(DATA, iter, options)
    %
    disp(' ');
    
    method  = options.method;
    iSeg    = options.iSeg;
    iShank  = options.iShank;
    verbose = options.verbose;
    iRepl   = options.iRepl;

    if strcmp(options.guess_method, 'km')     
        disp('Using K-means algorithm to compute initial guess...');
        [ X0, MU0, SIGMA0, DET0 ] = get_kmeans_initial_guess(DATA, iter, options);
    else         
        disp('Computing randomized initial guess...'); 
        [ X0, MU0, SIGMA0, DET0 ] = get_random_cell_initial_guess(DATA, iter, options);
    end

    K = size(DATA, 2);
    if options.verbose
        disp('Initial guess=');
        for i=1:((length(X0)-(options.model-1)*K)/6)
            disp(['           ' num2str(X0(int32(1:6)+(i-1)*6))]);
        end
        disp(' ')
        MU0
        DET0
        if options.model == 2
            disp(['           ' num2str(X0((length(X0)-K+1):length(X0)))]);
        end
    end

% ------------------------------------------------------------------------

function [ X0, MU, SIGMA, SIGMA_DET] = get_kmeans_initial_guess(DATA, iter, options)

    M = options.M;
    K = size(DATA, 2);
    MaxIter = 500;

    km_options = statset('Display', 'off', 'MaxIter', 400, 'UseParallel', false, 'UseSubstreams', false);

    % 'Normalize' data for more accurate estimation of Y
    N            = size(DATA, 1);
    DATA0        = zeros(N, K);
    for n=1:N
        DATA0(n, :) = DATA(n, :)/min(DATA(n, :));
    end

    disp(['M=' num2str(M) ' size(DATA0)=' num2str(size(DATA0))]);
    for i=1:(iter*int16(options.iRepl))
        r = rand(1,M);  % M-vector of random numbers between 0 and 1
    end
    N
    cols = round(N*r);
    seeds = DATA0(cols, :);
    cols
    seeds
    Y = kmeans(DATA0, M, 'start', seeds);
%   disp(['Y=' num2str(Y')]); 
    if options.verbose
        Y_unique = unique(Y');
        counts = zeros(1, length(Y_unique));
        for i=1:length(Y_unique)
            counts(i) = length(find(Y == Y_unique(i)));
        end
        disp(['Y_unique=', num2str(Y_unique) ' counts=' num2str(counts)]);
    end

    MU        = zeros(K, M);
    SIGMA     = zeros(K, K*M);
    SIGMA_DET = zeros(1, M);
    for m= 1:M
        DATA1 = DATA(find(Y == m),:);
        MU(:, m) = mean(DATA1)';
        SIGMA(:, ((m-1)*K+1):(m*K)) = cov(DATA1);
        SIGMA_DET(m) = det(SIGMA(:, ((m-1)*K+1):(m*K)));
    end

    XC = zeros(1, M);
    YC = zeros(1, M);
    ZC = zeros(1, M);
    if K == 8
        [XE, YE, ZE] = TNC_GM_ElectrodeCoordinates('Buzsaki64', [1:K]);
    else
        [XE, YE, ZE] = TNC_GM_ElectrodeCoordinates('tetrode', [1:K]);
    end
    for m=1:M
        if K == 4
            tetrode = reshape([XE, YE, ZE], K, 3);
            [XC(m), YC(m), ZC(m)] = TNC_GM_Triangulation(MU(:,m)', tetrode', options);
            if options.verbose >= 2
                disp(['m=' num2str(m) ' initial XC, YC, ZC=' num2str([XC(m), YC(m), ZC(m)])]);
            end
        elseif K == 8
            tetrode1 = reshape([XE([1 3 5 7]), YE([1 3 5 7]), ZE([1 3 5 7])], 4, 3);
            tetrode2 = reshape([XE([2 4 6 8]), YE([2 4 6 8]), ZE([2 4 6 8])], 4, 3);
            [XC1, YC1, ZC1] = TNC_GM_Triangulation(MU([1 3 5 7],m)', tetrode1', options);
            [XC2, YC2, ZC2] = TNC_GM_Triangulation(MU([2 4 6 8],m)', tetrode1', options);
            coord = ([XC1, YC1, ZC1] + [XC2, YC2, ZC2])/2; 
            XC(m) = coord(1);
            YC(m) = coord(2);
            ZC(m) = coord(3);
        else
            if options.verbose >= 2
                disp('TNC_GM_InitialGuess: K must be a multiple of 4');
            end
            return
        end
    end
    if options.model ~= 2
        X0 = ones(1, M*6);
    else
        X0 = ones(1, M*6 + K);
    end
    X0(1+(0:(M-1))*6) = 1/double(M);
    X0(2+(0:(M-1))*6) = mean(MU, 1); % vector of means for each row
    for m=1:M
        X0(3+(0:(M-1))*6) = mean(mean(SIGMA(:,((m-1)*K+1):(m*K))));
    end
    X0(4+(0:(M-1))*6) = XC; % cell X-coordinate
    X0(5+(0:(M-1))*6) = YC; % cell Y-coordinate
    X0(6+(0:(M-1))*6) = ZC; % cell Z-coordinate
   
%-------------------------------------------------------------------------------

function [flag] = standalone_detect()
    try
        dummy=which(mfilename); % this does not work at runtime in standalone mode
        flag=0; % i.e. running under matlab
    catch
        flag=1; % i.e. running in standalone mode
    end
 
% ------------------------------------------------------------------------

function [ X0, MU, SIGMA, SIGMA_DET ] = get_random_cell_initial_guess(DATA, iter, options)
    M = options.M;
    N = size(DATA, 1);
    K = size(DATA, 2);
    MU = zeros(K, M);
    SIGMA = zeros(K, K*M);
    SIGMA_DET = zeros(1, M);
    Xc = zeros(1, M);
    Yc = zeros(1, M);
    Zc = zeros(1, M);

    iRepl = options.iRepl;
    PI = 1/M;
    ALPHA = zeros(K, M);

    [dimX, dimY, dimZ, Xp, Yp, Zp] = get_mesh(K);
    num_knots = dimX*dimY*dimZ;

    if K == 8
        [XE, YE, ZE]  = TNC_GM_ElectrodeCoordinates('Buzsaki64', (1:8));
    else
        [XE, YE, ZE]  = TNC_GM_ElectrodeCoordinates('tetrode', (1:4));
    end

    % Determine coordinates and ALPHAs for M random points
    for m=1:M
        if standalone_detect()
            reset(RandStream.getGlobalStream,sum(100*clock));
        end
        knot = mod(round(rand(1)*(num_knots*iRepl-1)) + 1, num_knots);
        [iX, iY, iZ] = get_mesh_indices(knot, dimX, dimY, dimZ);
        if options.verbose
            disp(['knot=' num2str(knot) ' num_knots=' num2str(num_knots) ...
                  ' dimX,Y,Z=' num2str([dimX, dimY, dimZ]) ]);
            disp(['iX, iY, iZ=' num2str([iX, iY, iZ]) ]);
        end
        Xc(m) = Xp(iX);
        Yc(m) = Yp(iY);
        Zc(m) = Zp(iZ);
        ALPHA(:,m) = get_alpha(XE, YE, ZE, Xp(iX), Yp(iY), Zp(iZ), options.power);
        coord = [Xc(m), Yc(m), Zc(m)];
        coord
    end
    ALPHA

    % Classify datapoints using ALPHAs
    ALPHA0 = zeros(K, M); 
    for m=1:M
        ALPHA0(:,m) = sort(ALPHA(:,m)/min(ALPHA(:,m)));
    end
    DATA0  = zeros(N, K);
    for n=1:N
        DATA0(n, :) = DATA(n, :)/min(DATA(n, :));
    end
    Y = zeros(1, N);
    for n=1:N
        Xn = DATA0(n, :);
        ALPHAn = sort(Xn/min(Xn));
        distn = get_alpha_distances(ALPHAn, ALPHA);
%       disp(['Xn=' num2str(Xn) ' distn=' num2str(distn)]); 
        Y(n) = find(distn == min(distn)); 
    end
    Y

    % Setermine MU and SIGMA 
    for m=1:M
        IND = find(Y == m);
        if length(IND) > 0
            DATA1 = DATA(IND, :);
            MU(:,m) = mean(DATA1);
            SIGMA(:,(K*(m-1)+1):(K*m)) = cov(DATA1);
        else
            SIGMA(:,(K*(m-1)+1):(K*m)) = diag(1.e-8*ones(1,K));
        end
        SIGMA_DET(m) = det(SIGMA(:,(K*(m-1)+1):(K*m)));
    end
    MU
    SIGMA

    % Populate vector X0 
    if options.model ~= 2
        X0 = ones(1, M*6);
    else
        X0 = ones(1, M*6 + K);
    end
    X0(1+(0:(M-1))*6) = 1/M;
    X0(2+(0:(M-1))*6) = mean(MU, 1); % vector of means for each row
    for m=1:M
        X0(3+(0:(M-1))*6) = mean(mean(SIGMA(:,((m-1)*K+1):(m*K))));
    end
    X0(4+(0:(M-1))*6) = Xc; % cell X-coordinate
    X0(5+(0:(M-1))*6) = Yc; % cell Y-coordinate
    X0(6+(0:(M-1))*6) = Zc; % cell Z-coordinate

%-------------------------------------------------------------------------------

function dist = get_alpha_distances(An, ALPHA)
    M = size(ALPHA, 2);
    K = size(ALPHA, 1);
    An_sorted = sort(An);
    dist = zeros(1, M);
 
    for m=1:M
        ALPHA_sorted = sort(ALPHA(:,m));
        for k=1:K
            dist(m) = dist(m) + (ALPHA_sorted(k)-An_sorted(k))^2;
        end
        dist(m) = sqrt(dist(m)); 
    end    

%-------------------------------------------------------------------------------

function ALPHA = get_alpha(XE, YE, ZE, Xc, Yc, Zc, beta) 
    K = size(XE, 2);
    ALPHA = zeros(1,K);
   
    rev_dist_beta = zeros(1, K);
    for k=1:K
        rev_dist_beta(k) = 1/(sqrt((Xc - XE(k))^2 + (Yc - YE(k))^2 + (Zc - ZE(k))^2))^beta;
    end
    mean_rev_dist_beta = mean(rev_dist_beta);
    for k=1:K
        ALPHA(k) = rev_dist_beta(k)/mean_rev_dist_beta;
    end 

%-------------------------------------------------------------------------------

 % knot must satisfy:
 %
 % ind_knot = iX + (iY - 1)*dimX + (iZ - 1)*dimX*dimY

function[iX, iY, iZ] = get_mesh_indices(knot, dimX, dimY, dimZ)
    knot_XY = mod(knot, dimX*dimY);
    if knot_XY == 0
        iX = dimX;
        iY = dimY;
        iZ = knot/dimX/dimY;
    else
        iX    = mod(knot_XY, dimX);
        if iX == 0
            iX = dimX;
            iY = knot_XY/dimX;
            iZ = 1 + (knot    - iX - (iY-1)*dimX)/dimX/dimY;
        else
            iY = 1 + (knot_XY - iX              )/dimX;
            iZ = 1 + (knot    - iX - (iY-1)*dimX)/dimX/dimY;
        end
    end

%-------------------------------------------------------------------------------

function[iE1, iE2, iE3] = get_closest_electrodes(XE, YE, ZE, X, Y, Z)
    iE1 = 0;
    iE2 = 0;
    iE3 = 0;
    min_dist  = inf;
    min_dist2 = inf;
    min_dist3 = inf;
    for i=1:length(XE)
        dist = sqrt((XE(i)-X)^2 + (YE(i)-Y)^2 + (ZE(i)-Z)^2);
        if  min_dist  >= dist
            min_dist3 = min_dist2;
            min_dist2 = min_dist;
            min_dist  = dist;
            iE3       = iE2;
            iE2       = iE1;
            iE1       = i;
        elseif min_dist < dist && min_dist2 >= dist
            min_dist3 = min_dist2;
            min_dist2 = dist;
            iE3       = iE2;
            iE2       = i;
        elseif min_dist < dist && min_dist2 <  dist && min_dist3 >= dist
            min_dist3 = dist;
            iE3       = i;
        end
    end

%-------------------------------------------------------------------------------

function [dimX, dimY, dimZ, Xp, Yp, Zp] = get_mesh(K);
    if K == 8
        dimX = 6;
        dimY = 30;
        dimZ = 5;
        [XE, YE, ZE]  = TNC_GM_ElectrodeCoordinates('Buzsaki64', (1:8));

        Xp   = min(XE) - 20 + (max(XE)-min(XE) + 40)/(dimX-1)*(0:(dimX-1));
        Yp   = min(YE) - 20 + (max(YE)-min(YE) + 40)/(dimY-1)*(0:(dimY-1));
        Zp   = 5*ones(1,dimZ);

    else % K == 4
        dimX = 10;
        dimY = 10;
        dimZ = 5;
        [XE, YE, ZE]  = TNC_GM_ElectrodeCoordinates('tetrode', 1:4);     

        Xp   = min(XE) - 20 + (max(XE)-min(XE) + 40)/(dimX-1)*(0:(dimX-1));
        Yp   = min(YE) - 20 + (max(YE)-min(YE) + 40)/(dimY-1)*(0:(dimY-1));
        Zp   = 5*ones(1,dimZ);
    end

