function name = extractNameFromPath(path)
    % Use regular expression to match the part after the last slash
    % The pattern looks for any characters following the final slash in the path
    matches = regexp(path, '.*/(.*)$', 'tokens');
    
    % Check if there was a match
    if ~isempty(matches)
        % Extract the first (and only) match from the tokens
        name = matches{1}{1};
    else
        % If no match was found, return an empty string or handle as needed
        name = '';
    end
end
