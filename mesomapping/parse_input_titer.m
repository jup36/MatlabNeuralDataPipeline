function p = parse_input_titer(varargin)
        % parse input, and extract name-value pairs
        default_aav_titer_ml = 1e+13;  % the original AAV titer (vg/ml)
        default_aav_aliquot_vol = 100; % the original AAV volume provided (ul)
        default_target_vg = 2e+11;  % the target vg in the injectate (genomes)
        default_inject_vol = 100; % injection volume (ul)
        % the volume (ul) lost in the dead space of the syringe hub and needle

        p = inputParser; % create parser object

        %addRequired(p,'filePath');
        addParameter(p, 'aav_titer_ml', default_aav_titer_ml);
        addParameter(p, 'aav_aliquot_vol', default_aav_aliquot_vol);
        addParameter(p, 'target_vg', default_target_vg);
        addParameter(p, 'inject_vol', default_inject_vol); 

        parse(p, varargin{:})
    end