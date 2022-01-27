    function p = parse_input_plot_handJsTrjXY(filePath, vargs)
        % parse input, and extract name-value pairs
        default_draw_blocks = []; % 
        default_draw_trials = []; %
        default_numb_to_draw = 10; 
        default_rewarded_only = true; 
        default_draw_only_to_pullstart = false;
        default_draw_trials_per_block = []; 
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addParameter(p,'draw_blocks', default_draw_blocks); 
        addParameter(p,'draw_trials', default_draw_trials); 
        addParameter(p,'numb_to_draw', default_numb_to_draw); 
        addParameter(p,'rewarded_only', default_rewarded_only); 
        addParameter(p,'draw_only_to_pullstart', default_draw_only_to_pullstart)
        addParameter(p,'draw_trials_per_block', default_draw_trials_per_block)
        
        parse(p, filePath, vargs{:})
    end