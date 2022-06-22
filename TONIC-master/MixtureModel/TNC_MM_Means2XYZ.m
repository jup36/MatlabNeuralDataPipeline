function [ X, Y, Z ] = TNC_MM_Means2XYZ(MU1)
    % Given a set of mean values for clusters detected on each electrode, 
    % compute the most likely posiition X, Y, Z of the source
    global MU;
    global start2;

    MU = MU1;
    start2 = 1;
%   disp(['size(DATA)=' num2str(N)]);
%   options = optimset('LargeScale','off');
    options = optimoptions(@fminunc, 'TolFun', 0.00001, 'TolX', 0.00001, 'Display', 'off');
    [xmin, ~] = fminunc(@objective_fun1, [0 0 1], options);

    X = xmin(1);
    Y = xmin(2);
    Z = xmin(3);

%-------------------------------------------------------------------------------
function[fmin] = objective_fun1(X)
    % Compute the Gussian mixture likelihood function
    % Data is supposed to be stored in inputFile
    % X is a vector of parameters of the model:
    % X = [ MU_VECT, SIGMA_VECT, XC_VECT, YC_VECT, ZC_VECT ]
    persistent num_iterations;
    persistent previous_fmin;
    persistent best_fmin;
    persistent best_xmin;
    global MU;
    global start2;

    debug = 0;

    if start2 == 1
        num_iterations = 0;
        previous_fmin = 0.;
        best_fmin = inf;
        start2 = 0;
    end

    K = length(MU);
    if debug == 0
        [XE, YE, ZE]  = TNC_MM_ElectrodeCoordinates('Buzsaki64', (1:K));
    else
        [XE, YE, ZE]  = get_test_electrode_coordinates();
    end

    dist = [];
    for i=1:length(XE)
        dist1 = sqrt((X(1)-XE(i))^2 + (X(2)-YE(i))^2 + (X(3)-ZE(i))^2);
        dist = [ dist  dist1 ];
    end
    fmin = 0;
    for i=2:K                
       fmin = fmin + (dist(i)*MU(i) - dist(i-1)*MU(i-1))^2;    
    end
%   disp(['MU=' num2str(MU)]);
    if 0
        disp(['X=' num2str(X) ' dist=' num2str(dist) ' fmin1=' num2str(fmin)]);
    end
    if  best_fmin > fmin
        best_fmin = fmin;
        best_xmin = X;
    end
    if 0
        if fmin ~= previous_fmin
            num_iterations = num_iterations + 1;
            previous_fmin = fmin;
            disp(' ');
            disp(['Iteration# ' num2str(num_iterations) ]);
            disp(['Fmin=' num2str(fmin) ' best Fmin=' num2str(best_fmin)]);
            disp(['     X= ' num2str(X)]);
            disp(['best X= ' num2str(best_xmin)]);
            disp(' ');
        end
    end

%-------------------------------------------------------------------------------

function [XE, YE, ZE] = get_test_electrode_coordinates()
    iE = 1:3;
    X  = [0  5 10];
    Y  = [0  10 20];
    Z  = [0  20 -30];
    XE = X(iE);
    YE = Y(iE);
    ZE = Z(iE);

