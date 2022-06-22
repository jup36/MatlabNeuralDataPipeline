function [XC, YC, ZC] = TNC_GM_Triangulation(MU, tetrode, options)
    % Given a set of mean values for clusters detected on each electrode, 
    % compute the most likely posiition X, Y, Z of the source
    % MU is 4-vector of ithe valuesaof a feature
    % tetrode is 3 x 4 matrix of 3D coordinates of electrode tips
    % options.power   is the power (usually = 1) such that components of MU
    %      are propotional to the pw   powwer of theareciprocal distance
    %      from (point) source to the corresponding electrode
 
    K = size(MU, 2);
  
    [MU_sorted, IDX] = sort(abs(MU), 'descend');
    tetrode_sorted = tetrode(:, IDX);
 
    % Place 1st electrode to the origin of a coordinate frame
    tetrode_shifted = tetrode_sorted - [tetrode_sorted(:, 1) tetrode_sorted(:, 1) tetrode_sorted(:, 1) tetrode_sorted(:, 1)];

    % Rotate around z axis to make y-coordinate of the 2nd electrode zero
    ro_shifted = sqrt(tetrode_shifted(1, 2)^2 + tetrode_shifted(2, 2)^2);
    cos_theta = tetrode_shifted(1, 2)/ro_shifted;
    sin_theta = tetrode_shifted(2, 2)/ro_shifted;
    if options.debug
        disp(['cos_theta=' num2str(cos_theta) ' sin_theta=' num2str(sin_theta)]);
    end
    Rotation_xy_matrix = [cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];
    tetrode_xy_rotated = Rotation_xy_matrix * tetrode_shifted;

    % Rotate around y axis to make both y- and z-coordinate of the 2nd electrode zero
    ro_xy_rotated = sqrt(tetrode_xy_rotated(1, 2)^2 + tetrode_xy_rotated(3, 2)^2);
    cos_alpha = tetrode_xy_rotated(1, 2)/ro_xy_rotated;
    sin_alpha = tetrode_xy_rotated(3, 2)/ro_xy_rotated;
    if options.debug
        disp(['cos_alpha=' num2str(cos_alpha) ' sin_alpha=' num2str(sin_alpha)]);
    end
    Rotation_xz_matrix = [cos_alpha 0 sin_alpha; -sin_alpha 0 cos_alpha; 0 1 0];
    tetrode_xz_rotated = Rotation_xz_matrix * tetrode_xy_rotated;

    % Rotate around x axis to make z-coordinate of the 3rd electrode zero
    ro_xz_rotated = sqrt(tetrode_xy_rotated(2, 3)^2 + tetrode_xy_rotated(3, 3)^2);
    cos_psi   = tetrode_xz_rotated(2, 3)/ro_xz_rotated;
    sin_psi   = tetrode_xz_rotated(3, 3)/ro_xz_rotated;
    if options.debug
        disp(['cos_psi  =' num2str(cos_psi  ) ' sin_psi  =' num2str(sin_psi  )]);
    end
    Rotation_yz_matrix = [1 0 0; 0 cos_psi    sin_psi  ; 0 -sin_psi   cos_psi  ];
    tetrode_yz_rotated = Rotation_yz_matrix * tetrode_xz_rotated;

    ratios = zeros(1, K-1);
    for k=1:length(ratios)
        ratios(k) = MU_sorted(k)^(-options.power)/ MU_sorted(k+1)^(-options.power); % ratios of distance_1 : distance_1+1
    end
    if options.debug
        disp(['MU=' num2str(MU) ' MU_sorted=' num2str(MU_sorted) ' ratios=' num2str(ratios)]);
    end
    eps_step = 1.e-6;
    eps_abs  = 1.e-6;
    N = 1000;
    a = tetrode_yz_rotated(1, 2);
    coord_3 = tetrode_yz_rotated(:, 3)';
    coord_4 = tetrode_yz_rotated(:, 4)';
    if options.debug
        tetrode_yz_rotated
    end
    alpha_min = -1;
    alpha_max =  1;
    if ratios(1) > 1
        alpha_min = -sqrt(1-1/ratios(1));
        alpha_max =  sqrt(1-1/ratios(1));
    end
    [alpha, f3_min]  = bisection(@(x) f3(x, a, ratios, coord_3, coord_4, N, eps_step, eps_abs, options), ...
                                 alpha_min, alpha_max, N, eps_step, eps_abs);
    if options.debug
        disp(['In Main a=' num2str(a) ' ratios=' num2str(ratios) ...
              ' coord_3=' num2str(coord_3) ' coord_4=' num2str(coord_4) ...
              ' alpha=' num2str(alpha)]);
    end
    [a1, h, b, c]  = f1(alpha, a, ratios(1), options);
    [beta, f4_min] = bisection(@(x) f4(x, alpha, a, ratios, coord_3, coord_4, options), -1, 1, N, eps_step, eps_abs);
    if options.debug
        f3val = f3(alpha, a, ratios, coord_3, coord_4, N, eps_step, eps_abs, options);
        f4val = f4(beta, alpha, a, ratios, coord_3, coord_4, options);
        disp(['In Main alpha=' num2str(alpha) ' beta=' num2str(beta) ' a1=' num2str(a1) ' h=' num2str(h)]);
        disp(['Final f3=' num2str(f3val) ' f4=' num2str(f4val)]);
    end
    xc = a1;
    yc = h * beta;
    zc = h * sqrt(1-beta^2);
    if options.verbose 
        get_computed_ratios(xc, yc, zc, tetrode_yz_rotated, ratios)
    end
    Cell_coord_yz_rotated = [xc yc zc]';
    Cell_coord_xz_rotated = Rotation_yz_matrix' * Cell_coord_yz_rotated;
    Cell_coord_xy_rotated = Rotation_xz_matrix' * Cell_coord_xz_rotated;
    Cell_coord_shifted    = Rotation_xy_matrix' * Cell_coord_xy_rotated;
    Cell_coord            = Cell_coord_shifted  + tetrode_sorted(:, 1);
    XC = Cell_coord(1);
    YC = Cell_coord(2);
    ZC = Cell_coord(3);
   
%-------------------------------------------------------------------------------

function get_computed_ratios(xc, yc, zc, tetrode, ratios)
    dist = zeros(1,4);
    for i=1:4
        dist(i) = sqrt((tetrode(1,i)-xc)^2+(tetrode(2,i)-yc)^2+(tetrode(3,i)-zc)^2);
    end
    comp_ratios = zeros(1,3);
    for j=1:3
        comp_ratios(j)=dist(j)/dist(j+1);
    end
    disp(' ');
    disp(['Computed ratios=' num2str(comp_ratios)]); 
    disp(['True     ratios=' num2str(     ratios)]);

%-------------------------------------------------------------------------------

function [a1, h, b, c] = f1(alpha, a, r1, options)
    % Compute parameters of a triangle with sides a, b, c
    % Two vertices of the triangle are (0,0) and (a, 0)
    % coordinates of the 3rd vertice are (a1, h), where a1 is signed and h1 > 0
    % alpha = cos of the angle adjacent to vertice (0, 0)
    % r1 = b/c
    root_expression = 1 - r1^2*(1-alpha^2);
%   disp(['alpha=' num2str(alpha) ' r1=' num2str(r1) ' root_expression=' num2str(root_expression)]);
    if r1 < 1
        c     = a*(-r1*alpha + sqrt(root_expression))/(1-r1^2);
    elseif r1 == 1
        c = a/(2*alpha);
    else % r1 > 1 => alpha must satisfy:  -sqrt(1-1/r1) <= alpha <= sqrt(1-1/r1)
        if alpha^2 >= 1-1/r1
            c = a*(-r1*alpha - sqrt(root_expression))/(1-r1^2);
        else
            if options.verbose == 2
                disp('alpha must satisfy: alpha^2 >= 1 - 1/r1');
            end
            c = a*(-r1*alpha)/(1-r1^2);
        end
    end
    b  = r1 * c;
    h  = b  * sqrt(1 - alpha^2);
    a1 = b  * alpha;

%-------------------------------------------------------------------------------

function f = f3(x, a, ratios, coord_3, coord_4, N, eps_step, eps_abs, options)
    alpha = x;
    if options.debug
        disp(['    In f3: x=' num2str(x)]);
    end
    [a1, h, b, c]  = f1(alpha, a, ratios(1), options);
    [beta, f4_min] = bisection(@(x) f4(x, alpha, a, ratios, coord_3, coord_4, options),-1, 1, N, eps_step, eps_abs);
    xc = a1;
    yc = h * beta;
    zc = h * sqrt(1-beta^2);
    dist3 = sqrt((xc - coord_3(1))^2 + (yc - coord_3(2))^2 + zc^2);
    f = dist3 - b/(ratios(1)*ratios(2)); 
    if options.debug
        disp(['    In f3: a1=' num2str(a1) ' a=' num2str(a) ' alpha=' num2str(alpha) ...
              ' b=' num2str(b) ' c=' num2str(c) ' beta=' num2str(beta) ...
              ' ratios=' num2str(ratios) ' f=' num2str(f)]);
    end

%-------------------------------------------------------------------------------

function f = f4(x, alpha, a, ratios, coord_3, coord_4, options)
    beta = x;    
    [a1, h, b, c]  = f1(alpha, a, ratios(1), options);
    xc = a1;
    yc = h * beta;
    zc = h * sqrt(1-beta^2);
    dist4 = sqrt((xc - coord_4(1))^2 + (yc - coord_4(2))^2 + zc^2);
    f = dist4 - b/(ratios(1)*ratios(2)*ratios(3));
    if options.debug
        disp(['        In f4: a1=' num2str(a1) ' a=' num2str(a) ' alpha=' num2str(alpha) ...
              ' b=' num2str(b) ' c=' num2str(c) ' beta=' num2str(beta) ...
              ' ratios=' num2str(ratios) ' zc=' num2str(zc) ' alpha=' num2str(alpha) ' f= ' num2str(f)]);
    end

%-------------------------------------------------------------------------------

% f is a function handle
% Example: f = @(x) x^2 - 0.5;

function [x, w] = bisection(f, A, B, N, eps_step, eps_abs)
    u = f(A);
    v = f(B);
    e = B-A;
    if (u > 0 && v > 0) || (u < 0 && v < 0)
        x = A;
        w = u;   
        return;
    end;
    for k = 1:N
        e = e/2;
        C = A + e;
        w = f(C);
        if (abs(e) < eps_step || abs(w) < eps_abs)
            x = A;
            return;
        end
        if (w < 0 && u > 0) || (w > 0 && u < 0)
            B = C;
            v = w;
        else
            A = C;
            u = w;
        end
    end
    disp(['Iterations did not converge in ' num2str(N) ' steps']);
