function [correct_idx2] = sort_trial_start_end1(idx1, idx2)
        % sort index 2 relative to index 1
        correct_idx2 = nan(1, length(idx1)); 
        for j = 1:length(idx1) 
            correct_idx = find(idx1(j)<idx2, 1, 'first'); 
            correct_idx2(1,j) = idx2(correct_idx); 
        end

end
