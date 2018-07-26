function R_combined = get_recentered_combined(R, R_given)
    % compute the combined matrix of all re-centered distance matrices
    % returns a matrix, where each element is a pointwise-sum of all R
    % and R_given (remember, R_given is already re-centered)
    
    % initialize R_combined as a zero matrix
    if (~isempty(R))  
        R_combined = zeros(size(R{1}));  
    else
        R_combined = zeros(size(R_given));  % 
    end
    
    % iterate and add through R
    for iset = 1:length(R)
        R_combined = R_combined + R{iset}/length(R);
    end
    
    % incorporate R_given, if given
    if (~all(all(R_given==0))) % if R_given ~= 0, then D_given exists
        R_combined = (R_combined + R_given)/2;
    end

end
