function option = readOption(binFile)
% Parse ini file into cell entries C{1}{i} = C{2}{i}
metaFile = replace(binFile, '.bin', '.meta');
fid = fopen(metaFile, 'r');
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
    end
    meta.(tag) = C{2}{i};
end
option = meta.imProbeOpt(1);


end