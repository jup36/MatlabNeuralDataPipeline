function stringout = readCsvAsString(fname)
%This is an option to read a CSV file as string with (multiple) headers as is.
stringout = string(readcell(fullfile(fname)));
end
