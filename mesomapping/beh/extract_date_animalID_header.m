function header = extract_date_animalID_header(filepath)
    % Extract the portion of the file path that includes the preceding
    % string and the 6-digit date following the underscore.
    % Input:
    %   filepath - the input file path string
    % Output:
    %   header - the extracted string 'preceding_string_6digit_date'
    
    % Define the regex pattern: look for '\anyString_6digit' pattern
    pattern = '[\\]([a-zA-Z0-9]+_[0-9]{6})[\\]';
    
    % Use the regexp function to search and return the matched string
    match = regexp(filepath, pattern, 'tokens');
    
    % If a match is found, output the first match
    if ~isempty(match)
        header = match{1}{1};
    else
        header = ''; % Return an empty string if no match is found
    end
end