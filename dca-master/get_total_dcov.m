function total_dcov = get_total_dcov(R,D_given)
    % compute the total distance covariance across all datasets
    % R: (1 x num_datasets), re-centered matrices
    % D_given: (1 x num_given_datasets), combined re-centered matrix for given distance matrices

    R = [R D_given];
    
    Rtotal = 0;
    T = size(R{1},1);
    for iset = 1:length(R)
        for jset = (iset+1):length(R)
            Rtotal = Rtotal + sqrt(1/T^2 * R{iset}(:)' * R{jset}(:));  
        end
    end

    total_dcov = Rtotal / ((length(R)-1)*length(R)/2);
    
end

