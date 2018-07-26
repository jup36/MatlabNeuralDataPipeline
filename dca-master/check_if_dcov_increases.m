function result = check_if_dcov_increases(p, total_dcov, total_dcov_old, itotal_dcov)
    % returns true if increase in dcov is greater than the percent threshold
    %   or if the number of iterations is less than iteration constraint
    % else returns false

    if (p.Results.num_stoch_batch_samples == 0) % full gradient descent
        percent_increase = abs(total_dcov - total_dcov_old)/abs(total_dcov_old);

        if (total_dcov - total_dcov_old < 0)  % if value goes down, stop
            result = false;
        elseif (percent_increase >= p.Results.percent_increase_criterion && ...
                    itotal_dcov <= p.Results.num_iters_foreach_dim)
            result = true;
        else
            result = false;
        end
    else  % stochastic gradient descent...just check number of iterations
        if (itotal_dcov <= p.Results.num_iters_foreach_dim)
            result = true;
        else
            result = false;
        end
    end
end