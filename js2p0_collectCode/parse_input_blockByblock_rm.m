function p = parse_input_blockByblock_rm(sbb0, file_path, vargs)
% sbb0: session-by-block behavioral data
default_avg_b = true;
default_avg_m = true;
default_b_pairs = {[1, 5], [2, 6], [3, 7], [4, 8]};

p = inputParser; % create parser object
addRequired(p, 'sbb0')
addRequired(p, 'file_path')
addParameter(p, 'avg_b', default_avg_b)
addParameter(p, 'avg_m', default_avg_m)
addParameter(p, 'b_pairs', default_b_pairs)

parse(p, sbb0, file_path, vargs{:})

end