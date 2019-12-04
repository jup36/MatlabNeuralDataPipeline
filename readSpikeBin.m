function [data_array, meta] = readSpikeBin(bin_file, siteToRead)
if exist(bin_file, 'file') == 0; data_array = []; return; end


% BUFFER_SIZE = 32 * 1024^3; % 2^35
% if ispc
%     [~, sys] = memory;
%     current_memory = sys.PhysicalMemory.Available;
%     buffer_size = min([current_memory * 0.5, BUFFER_SIZE]);
% else
%     buffer_size = BUFFER_SIZE; 
% end

% load meta data
meta = readMeta(bin_file);

% file size check
data_info = dir(bin_file);
data_size = data_info.bytes;
assert(data_size == meta.fileSizeBytes, 'Corrected file');
assert(mod(data_size, 2 * meta.nSavedChans) == 0, 'Corrected file');

% set offset, read size
sample_number = data_size / (2 * meta.nSavedChans);
if strcmp(meta.typeThis, 'imec')
    %offset = meta.snsApLfSy(1);
    %n_sample = meta.snsApLfSy(3);
    data_type = '*uint16';
else
    %offset = meta.snsMnMaXaDw(1);
    %n_sample = meta.snsMnMaXaDw(2);
    data_type = '*int16';
end

sample_number_buffer = floor(buffer_size / meta.nSavedChans / 2);
n_read = ceil(sample_number / sample_number_buffer);
last_sample_number = mod(sample_number, sample_number_buffer);

% read file
fprintf('Loading %s ... ', bin_file);

data_array = zeros(1, sample_number, 'uint16');

tic;
fid = fopen(bin_file, 'rb');
fseek(fid, 0, 'bof');

for i = 0:meta.fileTimeSecs  
    
end




for i_read = 1:n_read
    fprintf('%d/%d ... ', i_read, n_read);
    if i_read < n_read
        read_n_sample = sample_number_buffer;
    else
        read_n_sample = last_sample_number;
    end
    
    data_array_temp = fread(fid, [meta.nSavedChans, read_n_sample], data_type);
    data_array(:, (1:read_n_sample) + (i_read - 1) * sample_number_buffer) = data_array_temp(siteToRead, :);
    clear data_array_temp
end
fprintf('\n');

data_array = data_array';

fclose(fid);
toc;