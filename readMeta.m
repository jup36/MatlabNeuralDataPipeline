function meta = readMeta(binFile)
% Parse ini file into cell entries C{1}{i} = C{2}{i}
fid = fopen(replace(binFile, '.bin', '.meta'), 'r');
C = textscan(fid, '%[^=] = %[^\r\n]');
fclose(fid);

% New empty struct
meta = struct();

% Convert each cell entry into a struct entry
for i = 1:length(C{1})
    tag = C{1}{i};
    if tag(1) == '~'
        % remake tag excluding first character
        tag = sprintf('%s', tag(2:end));
        meta.(tag) = C{2}{i};
    else
        valueTemp = str2double(strsplit(C{2}{i}, ','));
        if isnan(valueTemp)
            meta.(tag) = C{2}{i};
        else
            meta.(tag) = valueTemp;
        end
    end
end