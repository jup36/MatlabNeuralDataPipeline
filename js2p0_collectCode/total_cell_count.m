filePath = {'/Volumes/dudmanlab/junchol/js2p0/WR37_022119', ... % Dual recording without silencing, Trj checked (07/04/22) % GOOD TO BE USED
            '/Volumes/dudmanlab/junchol/js2p0/WR37_022219', ... % Dual recording without silencing, Trj checked            % BAD (Lots of Matrix singular or bad-scaled warnings)
            '/Volumes/dudmanlab/junchol/js2p0/WR37_022619', ... % Cg recording contra-Cg silencing, Trj checked            % GOOD TO BE USED
            '/Volumes/dudmanlab/junchol/js2p0/WR37_022719', ... % B Only contra-Cg delayed silencing, Trj checked          % BAD (Lots of Matrix singular or bad-scaled warnings)
            '/Volumes/dudmanlab/junchol/js2p0/WR38_052119', ... % Dual recording without silencing, Trj checked            % There was an issue with 'corr' after running iterations that needs to be revisited!
            '/Volumes/dudmanlab/junchol/js2p0/WR38_052219', ... % Dual recording without silencing, Trj checked            % GOOD
            '/Volumes/dudmanlab/junchol/js2p0/WR38_052319', ... % Cg recording contra-Cg silencing, Trj checked            % BAD (Crashed)
            '/Volumes/dudmanlab/junchol/js2p0/WR38_052419', ... % Corticostriatal recording M1 silencing, Trj checked       
            '/Volumes/dudmanlab/junchol/js2p0/WR39_091019', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR39_091119', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR39_100219', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/dudmanlab/junchol/js2p0/WR39_100319', ... % B Only contra-Cg delayed silencing, Trj checked      
            '/Volumes/dudmanlab/junchol/js2p0/WR40_081919', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/dudmanlab/junchol/js2p0/WR40_082019', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/dudmanlab/junchol/js2p0/WR44_031020'};    % Dual recording with contra Cg delayed silencing, Trj checked % GOOD

ctx_c = cell(length(filePath), 2); ctx_c(:,2) = deal({0});  
str_c = cell(length(filePath), 2); str_c(:,2) = deal({0});
cg_c = cell(length(filePath), 2); cg_c(:,2)  = deal({0});
        
for j = 1:length(filePath)
    cd(fullfile(filePath{j}, 'Matfiles'))
    wrI = strfind(filePath{j}, 'WR');
    
    m_name = filePath{j}(wrI:wrI+10);
    ctx_c{j, 1} = m_name;
    str_c{j, 1} = m_name;
    cg_c{j, 1} = m_name; 
    
    spkDir = dir(['binSpkCountSTRCTX', m_name, '.mat']);
    spkDir_cg = dir(['binSpkCountCg', m_name, '.mat']);
    
    if ~isempty(spkDir)
       load(fullfile(spkDir(1).folder,fullfile(spkDir(1).name)),'spkTimesCell'); 
       ctx_c{j, 2} = sum(~cell2mat(spkTimesCell(5, :)));
       str_c{j, 2} = sum(cell2mat(spkTimesCell(5, :)));
    end
    
    if ~isempty(spkDir_cg)
        load(fullfile(spkDir_cg(1).folder,fullfile(spkDir_cg(1).name)),'spkTimesCell'); 
        cg_c{j, 2} = size(spkTimesCell, 2); 
    end
    fprintf('processed file# %d\n', j)
end

ctx_sum = sum(cell2mat(ctx_c(:, 2)));
str_sum = sum(cell2mat(str_c(:, 2)));
cg_sum = sum(cell2mat(cg_c(:, 2)));

save(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData', 'js2p0_total_cell_count.mat'), 'ctx*', 'str*', 'cg*')

