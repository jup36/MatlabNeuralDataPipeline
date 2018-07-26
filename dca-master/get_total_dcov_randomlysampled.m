function total_dcov = get_total_dcov_randomlysampled(u, X, D_given, p)
 % computes dcov for a random subsample (for stochastic gradient descent)
 

    r = randperm(size(X{1},2));
    sample_indices = r(1:min(length(r), p.Results.num_samples_to_compute_stepwise_dcov));
    
    T = length(sample_indices);  % T = number of subsamples
    
    R = cell(1,length(X) + length(D_given));
    for iset = 1:length(X)
        R{iset} = get_recentered_matrix(u{iset}, X{iset}(:,sample_indices));
    end

    for iset = 1:length(D_given)
        R{iset + length(X)} = D_given{iset}(sample_indices, sample_indices);  % this is an approximation, we would really need to re-compute the re-centered distance matrix for each D_given
    end
    
    Rtotal = 0;
    for iset = 1:length(R)
        for jset = (iset+1):length(R)
            Rtotal = Rtotal + sqrt(1/T^2 * R{iset}(:)' * R{jset}(:));  
        end
    end

    total_dcov = Rtotal / ((length(R)-1)*length(R)/2);
end
