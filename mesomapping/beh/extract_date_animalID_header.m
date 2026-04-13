function header = extract_date_animalID_header(filepath)
    % Extracts 'animalID_6digitDate' or 'animalID_6digitDate-N'

    pattern = '([A-Za-z]+\d+_\d{6}(?:-\d+)?)'; %'(?:^|[\\/])([A-Za-z]+\d+_\d{6}(?:-\d+)?)(?=[\\/]|$)';

    match = regexp(filepath, pattern, 'tokens', 'once');

    if ~isempty(match)
        header = match{1};
    else
        header = '';
    end
end

% function header = extract_date_animalID_header(filepath)
%     % Extracts 'animalID_6digitDate' or 'animalID_6digitDate-N' from filepath
% 
%     % --- Pattern 1: standard form (no suffix)
%     pattern1 = '(?:^|[\\/])(\w+[0-9]+_[0-9]{6})(?=[\\/]|$)';
% 
%     % --- Pattern 2: with dash + single digit suffix
%     pattern2 = '(?:^|[\\/])(\w+[0-9]+_[0-9]{6}-\d)(?=[\\/]|$)';
% 
%     % Try first pattern
%     match = regexp(filepath, pattern1, 'tokens', 'once');
% 
%     if isempty(match)
%         % Try second pattern if first fails
%         match = regexp(filepath, pattern2, 'tokens', 'once');
%     end
% 
%     % Return result
%     if ~isempty(match)
%         header = match{1};
%     else
%         header = '';
%     end
% end

% function header = extract_date_animalID_header(filepath)
%     % Extract the portion of the file path that includes the preceding
%     % string and the 6-digit date following the underscore.
%     % Input:
%     %   filepath - the input file path string
%     % Output:
%     %   header - the extracted string 'preceding_string_6digit_date'
% 
%     % Define the regex pattern: look for either '\' or '/' followed by
%     % 'anyString_6digit', accounting for both Windows and Mac path separators.
%     pattern = '[\\/](\w+_[0-9]{6})[\\/]';
% 
%     % Use the regexp function to search and return the matched string
%     match = regexp(filepath, pattern, 'tokens');
% 
%     % If a match is found, output the first match
%     if ~isempty(match)
%         header = match{1}{1};
%     else
%         header = ''; % Return an empty string if no match is found
%     end
% end
% 
