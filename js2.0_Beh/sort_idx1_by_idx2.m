function [correct_idx1] = sort_idx1_by_idx2(idx1, idx2)
        % sort index 1 relative to index 2
        correct_idx1 = nan(1, length(idx1)); 
        for j = 1:length(idx2) 
            correct_idx = find(idx1<idx2(j), 1, 'last'); 
            correct_idx1(1,j) = idx1(correct_idx); 
        end
end