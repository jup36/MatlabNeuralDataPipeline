function utemp = dca_one(Xtemp, Xijtemp, R_combined, u_0, column_indices, p)
    % performs distance covariance analysis for one dataset and one given re-centered distance matrix
    %       uses projected gradient descent
    % 
    %   X: (N x T), data in which we want to find the N x 1 dca dimension, where
    %       N is the number of variables, and T is the number of samples
    %   Xij: Xij = X_i - X_j (num_neurons x num_samples^2) for each dataset    
    %   R_combined: (T x T), combined re-centered distance matrix of the other sets of variables
    %   u_0: (N x 1), initial guess for the dca dimension
    %   p: (1 x 1), inputParser object which contains user constraints, such as number of iterations
    %   column_indices: think of it as a mapping between each element of
    %                   a pairwise matrix (e.g. 1000-by-1000) onto a long array containing the flattened version of the pairwise matrix (e.g. 1-by-1000^2) 
    %
    %  returns:
    %   u: (N x 1), the dimension of greatest distance covariance between D_X and R_combined
    
    % 12/6: replace Xtemp back to X, replace utemp back to u, replace Xijtemp back to Xij
    
    %%% PRE-PROCESSING
        N = size(Xtemp,1);  % number of neurons
        T = size(Xtemp,2);  % number of timepoints
        
        if (sum(var(Xtemp')) < 1e-10) % X has little variability left
            u = randn(N,1);
            u = u / norm(u);
            return;
        end

        utemp = u_0;    % set u (the dimension of greatest distance covariance between u'X and Y) to be initial guess
        
    %%% OPTIMIZATION  
        for istep = 1:p.Results.num_iters_per_dataset  % stop when num iters have been reached

            % COMPUTE GRAD DESCENT PARAMETERS
                D_uXij = get_D_uXij(utemp);  % get distance matrix of current u (u'X) - to be recentered in the next step (get_f)
                % dimension: sample x sample (distance matrix)
                
                f_val = get_f(D_uXij);       % compute current objective function value (Scalar, distance covariance between u'X and R_combined) and gradf for backtracking
                % dimension: scalar (objective function to be minimized)
                
                gradf = get_gradf(utemp, D_uXij);  % compute gradient of current solution (basically, re-centered and weighted Xij linearly combined with R_combined)
                % dimension: neuron x 1 (a vector whose dimension is the number of neurons(variables))
                
                t = 1;  % backtracking step size

            % BACKTRACKING LINE SEARCH
                % first check large intervals for t (so backtracking loop doesn't take forever)
                for candidate_power = 1:9
                    fprintf('.');
                    if (~backtrack_check(utemp, f_val, gradf, 10^-candidate_power))
                        break;
                    else
                        t = 10^-candidate_power;
                    end
                end

                % find more nuanced t
                while (backtrack_check(utemp, f_val, gradf, t) && t > 10^-9)
                    t = 0.7 * t;
                    fprintf('.')
                end 

            % PERFORM PROJECTED GRAD DESCENT
                u_unnorm = utemp - t * gradf; % gradient descent step

                norm_u = norm(u_unnorm); % project u_unnorm to the L2 unit ball
                if (norm_u > 1)
                    utemp = u_unnorm / norm_u;
                else
                    utemp = u_unnorm;   % allow solution to exist inside unit ball (for dca_one)
                end
        end

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NESTED DCA_ONE HELPER FUNCTIONS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function D = get_D_uXij(utemp) 
        % compute distance matrix of (X projected to u)
        
        D = squareform(pdist((utemp' * Xtemp)')); % D_uXij; a matrix whose dimension is sample x sample (e.g. 1000 x 1000)
    end

    function f = get_f(D_uXij)
        % compute objective function value (scalar)
        
        H = eye(T) - 1/T * ones(T);
        A = H * D_uXij * H;          % recenters distance matrix
        f = -R_combined(:)' * A(:);  % trying to minimize (objective function), so flip the sign!
    end

    function gradf = get_gradf(utemp, D_uXij)
        % computes the gradient for dca_one...note there are some tricky matrix operations in here!
        
        %%% weight Xij
            D_uXij(D_uXij == 0) = 1e-8; % add 1e-8 to avoid division by zero...other values
                                        % do not seem to matter
            % project Xij onto u
            XijT_u = Xijtemp' * utemp;  % a vector with dimension: length of Xij (sample^2) by 1 (e.g. 1000000x1)
                                        % XijT_u is closely related to D_uXij 
                                        
            % weight Xij by XijT_u ./ D_uXij(:)'
            Xij_weighted = bsxfun(@times, Xijtemp, (XijT_u ./ D_uXij(:))'); % just some sort of weighting, 'Xij_weighted' has the same dimension as 'Xij' (neuron x sample^2))
            % dimension: 
               % Xij:          neuron x (sample^2)
               % XijT_u:       (sample^2) x 1
               % D_uXij:       (sample^2) x 1
               % Xij_weighted: neuron x (sample^2)
            
        %%% subtract row, column, and matrix means (re-centering)
            Xij_row_means = blockproc(Xij_weighted, [N, T], @get_row_means); % blockproc(A,blockSize,fun) processes the image A by applying the function fun to each distinct block of A and concatenating the results into B
            Xij_col_means = Xij_row_means(:,column_indices); % to get Xij_col_means, just reorganize the Xij_row_means based on the column_indices
            Xij_matrix_mean = mean(Xij_weighted,2);

            Xij_weighted = bsxfun(@plus, Xij_weighted - Xij_row_means - Xij_col_means, Xij_matrix_mean); % re-centering

        %%% linearly combine with R_combined
            gradf = - Xij_weighted * R_combined(:);  % sign because we are minimizing negative dcov (a vector whose dimension is the number of variables)
            % a vector (dimension: neuron x 1) - still not sure (12/11) why this is the gradient for dca_one
            %                                  - gradf must be the derivative of the objective function f: D_uXij(normalized) * R_combined(:), thus
            %                                  specifically, Xij_weighted should be the derivative of D_uXij (a term in the objective function). 
           
                
    
    end

    function X_row = get_row_means(block_struct)
        % for blockproc, compute means along rows of distance matrix
        
        X_row = mean(block_struct.data,2);
        X_row = repmat(X_row, 1, size(block_struct.data,2));
    end

    function status = backtrack_check(utemp, f_next, gradf, t) 
        % check lecture 8 of ryan tibshirani opti class
        Gt = get_Gt(utemp, gradf, t); % t: backtracking stepsize, u: the dimension for each dataset that we are looking for maximizing the distance covariance 

        D_uXij_t = get_D_uXij(utemp - t * Gt);
        status = get_f(D_uXij_t) > f_next - t * gradf' * Gt + t/2 * Gt' * Gt;
    end

    function Gt = get_Gt(utemp, gradf, t)
        % vector used for backtracking check with projected gradient descent
        
        u_n = utemp - t * gradf; % gradient descent of u
        norm_u_n = norm(u_n);    % normalize u
        if (norm_u_n > 1)  % project to L2 unit ball
            u_norm = u_n / norm_u_n;
        else
            u_norm = u_n;
        end

        Gt = 1/t * (utemp - u_norm);
    end
end
