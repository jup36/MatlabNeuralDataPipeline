
filePath = {'/Volumes/Extreme SSD/js2p0/WR37_022119', ... % Dual recording without silencing, Trj checked (07/04/22)
            '/Volumes/Extreme SSD/js2p0/WR37_022219', ... % Dual recording without silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR37_022619', ... % Cg recording contra-Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR37_022719', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR38_052119', ... % Dual recording without silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR38_052219', ... % Dual recording without silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR38_052319', ... % Cg recording contra-Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR38_052419', ... % Corticostriatal recording M1 silencing, Trj checked       
            '/Volumes/Extreme SSD/js2p0/WR39_091019', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR39_091119', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR39_100219', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR39_100319', ... % B Only contra-Cg delayed silencing, Trj checked      
            '/Volumes/Extreme SSD/js2p0/WR40_081919', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR40_082019', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR44_031020'};    % Dual recording with contra Cg delayed silencing, Trj checked

filePath_dPCA = '/Volumes/Extreme SSD/js2p0/dPCAForJunchol.mat'; 


%% load
dPCA = load(fullfile(filePath_dPCA)); 
dPCA = dPCA.("dPCAForJunchol"); 
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez.mat'), 'pull_force', 'rch_angle')
header_dPCA = fieldnames(dPCA); 
header_filePath = cell(length(filePath), 1); 
for f = 1:length(filePath)
    [~, header_filePath{f}] = fileparts(filePath{f}); 
end

for f = 1:length(header_dPCA)
   
    mId = cellfun(@(a) strcmpi(header_dPCA{f}, a), header_filePath); 

    %% dPCA separation riHi vs riLo
    riHi_riLo_M1{f, 1} = dPCA.(header_dPCA{f}).s_sep_M1(1); 
    riHi_riLo_M1{f, 2} = dPCA.(header_dPCA{f}).s_sep_M1(1)./sum(dPCA.(header_dPCA{f}).s_sep_M1)*100; 
    riHi_riLo_STR{f, 1} = dPCA.(header_dPCA{f}).s_sep_STR(1); 
    riHi_riLo_STR{f, 2} = dPCA.(header_dPCA{f}).s_sep_STR(1)./sum(dPCA.(header_dPCA{f}).s_sep_STR)*100; 
    
    %% dPCA separation riHi vs leHi
    riHi_leHi_M1{f, 1} = dPCA.(header_dPCA{f}).s_sep_M1(2); 
    riHi_leHi_M1{f, 2} = dPCA.(header_dPCA{f}).s_sep_M1(2)./sum(dPCA.(header_dPCA{f}).s_sep_M1)*100; 
    riHi_leHi_STR{f, 1} = dPCA.(header_dPCA{f}).s_sep_STR(2); 
    riHi_leHi_STR{f, 2} = dPCA.(header_dPCA{f}).s_sep_STR(2)./sum(dPCA.(header_dPCA{f}).s_sep_STR)*100; 

    %% dPCA separation riHi vs leLo
    riHi_leLo_M1{f, 1} = dPCA.(header_dPCA{f}).s_sep_M1(3); 
    riHi_leLo_M1{f, 2} = dPCA.(header_dPCA{f}).s_sep_M1(3)./sum(dPCA.(header_dPCA{f}).s_sep_M1)*100; 
    riHi_leLo_STR{f, 1} = dPCA.(header_dPCA{f}).s_sep_STR(3); 
    riHi_leLo_STR{f, 2} = dPCA.(header_dPCA{f}).s_sep_STR(3)./sum(dPCA.(header_dPCA{f}).s_sep_STR)*100; 
    
    %% dPCA separation riLo vs leHi
    riLo_leHi_M1{f, 1} = dPCA.(header_dPCA{f}).s_sep_M1(4); 
    riLo_leHi_M1{f, 2} = dPCA.(header_dPCA{f}).s_sep_M1(4)./sum(dPCA.(header_dPCA{f}).s_sep_M1)*100; 
    riLo_leHi_STR{f, 1} = dPCA.(header_dPCA{f}).s_sep_STR(4); 
    riLo_leHi_STR{f, 2} = dPCA.(header_dPCA{f}).s_sep_STR(4)./sum(dPCA.(header_dPCA{f}).s_sep_STR)*100; 

    %% dPCA separation riLo vs leLo
    riLo_leLo_M1{f, 1} = dPCA.(header_dPCA{f}).s_sep_M1(5); 
    riLo_leLo_M1{f, 2} = dPCA.(header_dPCA{f}).s_sep_M1(5)./sum(dPCA.(header_dPCA{f}).s_sep_M1)*100; 
    riLo_leLo_STR{f, 1} = dPCA.(header_dPCA{f}).s_sep_STR(5); 
    riLo_leLo_STR{f, 2} = dPCA.(header_dPCA{f}).s_sep_STR(5)./sum(dPCA.(header_dPCA{f}).s_sep_STR)*100; 

    %% dPCA separation leHi vs leLo
    leHi_leLo_M1{f, 1} = dPCA.(header_dPCA{f}).s_sep_M1(6); 
    leHi_leLo_M1{f, 2} = dPCA.(header_dPCA{f}).s_sep_M1(6)./sum(dPCA.(header_dPCA{f}).s_sep_M1)*100; 
    leHi_leLo_STR{f, 1} = dPCA.(header_dPCA{f}).s_sep_STR(6); 
    leHi_leLo_STR{f, 2} = dPCA.(header_dPCA{f}).s_sep_STR(6)./sum(dPCA.(header_dPCA{f}).s_sep_STR)*100; 
    
    %%
    % reach angle separation
    rch_angleC = rch_angle(mId, :);
    rch_angle_med = nanmedian(cell2mat(rch_angleC')); 
    rch_angle_medSubC = cellfun(@(a) a-rch_angle_med, rch_angleC, 'UniformOutput', false); 
    leLo_rchAngle{f} = cell2mat(rch_angle_medSubC(1, [1, 5])'); 
    leHi_rchAngle{f} = cell2mat(rch_angle_medSubC(1, [2, 6])'); 
    riLo_rchAngle{f} = cell2mat(rch_angle_medSubC(1, [3, 7])'); 
    riHi_rchAngle{f} = cell2mat(rch_angle_medSubC(1, [4, 8])'); 

    % pull force separation
    pull_forceC = pull_force(mId, :);
    leLo_pullForce{f} = cell2mat(pull_forceC(1, [1, 5])'); 
    leHi_pullForce{f} = cell2mat(pull_forceC(1, [2, 6])'); 
    riLo_pullForce{f} = cell2mat(pull_forceC(1, [3, 7])'); 
    riHi_pullForce{f} = cell2mat(pull_forceC(1, [4, 8])'); 

end

for f = 1:length(header_dPCA)  
    %% right high vs right low
    % right high vs right low pull Force
    riHi_riLo_dPCA_beh{f, 1} = riHi_riLo_M1{f, 1}; 
    riHi_riLo_dPCA_beh{f, 2} = riHi_riLo_M1{f, 2}; 
    riHi_riLo_dPCA_beh{f, 3} = riHi_riLo_STR{f, 1}; 
    riHi_riLo_dPCA_beh{f, 4} = riHi_riLo_STR{f, 2}; 
    riHi_riLo_dPCA_beh{f, 5} = abs(nanmean(riHi_rchAngle{f})-nanmean(riLo_rchAngle{f})); 
    riHi_riLo_dPCA_beh{f, 6} = abs(nanmean(riHi_pullForce{f})-nanmean(riLo_pullForce{f})); 

    %% right high vs left high
    % right high vs left high reach Angle
    riHi_leHi_dPCA_beh{f, 1} = riHi_leHi_M1{f, 1}; 
    riHi_leHi_dPCA_beh{f, 2} = riHi_leHi_M1{f, 2}; 
    riHi_leHi_dPCA_beh{f, 3} = riHi_leHi_STR{f, 1}; 
    riHi_leHi_dPCA_beh{f, 4} = riHi_leHi_STR{f, 2}; 
    riHi_leHi_dPCA_beh{f, 5} = abs(nanmean(riHi_rchAngle{f})-nanmean(leHi_rchAngle{f})); 
    riHi_leHi_dPCA_beh{f, 6} = abs(nanmean(riHi_pullForce{f})-nanmean(leHi_pullForce{f})); 

    %% right high vs left low
    riHi_leLo_dPCA_beh{f, 1} = riHi_leLo_M1{f, 1}; 
    riHi_leLo_dPCA_beh{f, 2} = riHi_leLo_M1{f, 2}; 
    riHi_leLo_dPCA_beh{f, 3} = riHi_leLo_STR{f, 1}; 
    riHi_leLo_dPCA_beh{f, 4} = riHi_leLo_STR{f, 2}; 
    riHi_leLo_dPCA_beh{f, 5} = abs(nanmean(riHi_rchAngle{f})-nanmean(leLo_rchAngle{f})); 
    riHi_leLo_dPCA_beh{f, 6} = abs(nanmean(riHi_pullForce{f})-nanmean(leLo_pullForce{f})); 

    %% right low vs left high
    riLo_leHi_dPCA_beh{f, 1} = riLo_leHi_M1{f, 1}; 
    riLo_leHi_dPCA_beh{f, 2} = riLo_leHi_M1{f, 2}; 
    riLo_leHi_dPCA_beh{f, 3} = riLo_leHi_STR{f, 1}; 
    riLo_leHi_dPCA_beh{f, 4} = riLo_leHi_STR{f, 2}; 
    riLo_leHi_dPCA_beh{f, 5} = abs(nanmean(riLo_rchAngle{f})-nanmean(leHi_rchAngle{f})); 
    riLo_leHi_dPCA_beh{f, 6} = abs(nanmean(riLo_pullForce{f})-nanmean(leHi_pullForce{f})); 

    %% right low vs left low
    riLo_leLo_dPCA_beh{f, 1} = riLo_leLo_M1{f, 1}; 
    riLo_leLo_dPCA_beh{f, 2} = riLo_leLo_M1{f, 2}; 
    riLo_leLo_dPCA_beh{f, 3} = riLo_leLo_STR{f, 1}; 
    riLo_leLo_dPCA_beh{f, 4} = riLo_leLo_STR{f, 2}; 
    riLo_leLo_dPCA_beh{f, 5} = abs(nanmean(riLo_rchAngle{f})-nanmean(leLo_rchAngle{f})); 
    riLo_leLo_dPCA_beh{f, 6} = abs(nanmean(riLo_pullForce{f})-nanmean(leLo_pullForce{f})); 

    %% left high vs left low
    leHi_leLo_dPCA_beh{f, 1} = leHi_leLo_M1{f, 1}; 
    leHi_leLo_dPCA_beh{f, 2} = leHi_leLo_M1{f, 2}; 
    leHi_leLo_dPCA_beh{f, 3} = leHi_leLo_STR{f, 1}; 
    leHi_leLo_dPCA_beh{f, 4} = leHi_leLo_STR{f, 2}; 
    leHi_leLo_dPCA_beh{f, 5} = abs(nanmean(leHi_rchAngle{f})-nanmean(leLo_rchAngle{f})); 
    leHi_leLo_dPCA_beh{f, 6} = abs(nanmean(leHi_pullForce{f})-nanmean(leLo_pullForce{f})); 
end

%%
ri_le_dPCA_rchAng(:, 1) = cell2mat([riHi_leHi_dPCA_beh(:, 1); riHi_leLo_dPCA_beh(:, 1); riLo_leLo_dPCA_beh(:, 1); riLo_leHi_dPCA_beh(:, 1)]); 
ri_le_dPCA_rchAng(:, 2) = cell2mat([riHi_leHi_dPCA_beh(:, 5); riHi_leLo_dPCA_beh(:, 5); riLo_leLo_dPCA_beh(:, 5); riLo_leHi_dPCA_beh(:, 5)]); 

ri_le_forceMatched_dPCA_rchAng(:, 1) = cell2mat([riHi_leHi_dPCA_beh(:, 1); riLo_leLo_dPCA_beh(:, 1)]); 
ri_le_forceMatched_dPCA_rchAng(:, 2) = cell2mat([riHi_leHi_dPCA_beh(:, 5); riLo_leLo_dPCA_beh(:, 5)]); 

scatter(ri_le_dPCA_rchAng(:, 1), ri_le_dPCA_rchAng(:, 2)) 
[corr.r_rchAng, corr.p_rchAng] = corr(ri_le_dPCA_rchAng(:, 1), ri_le_dPCA_rchAng(:, 2)); 

scatter(ri_le_forceMatched_dPCA_rchAng(:, 1), ri_le_forceMatched_dPCA_rchAng(:, 2)) 
[corr.r_rchAng, corr.p_rchAng] = corrcoeff(ri_le_forceMatched_dPCA_rchAng(:, 1), ri_le_forceMatched_dPCA_rchAng(:, 2)); 

%%
hi_lo_dPCA_pullForce(:, 1) = cell2mat([riHi_leLo_dPCA_beh(:, 2); leHi_leLo_dPCA_beh(:, 2); riLo_leHi_dPCA_beh(:, 2); riLo_leLo_dPCA_beh(:, 2)]); 
hi_lo_dPCA_pullForce(:, 2) = cell2mat([riHi_leLo_dPCA_beh(:, 6); leHi_leLo_dPCA_beh(:, 6); riLo_leHi_dPCA_beh(:, 6); riLo_leLo_dPCA_beh(:, 6)]); 

scatter(hi_lo_dPCA_pullForce(:, 1), hi_lo_dPCA_pullForce(:, 2)) 
[corr.r_pullForce, corr.p_pullForce] = corrcoef(hi_lo_dPCA_pullForce(:, 1), hi_lo_dPCA_pullForce(:, 2)); 

hi_lo_rchAngMatched_dPCA_pullForce(:, 1) = cell2mat([leHi_leLo_dPCA_beh(:, 2); riHi_leLo_dPCA_beh(:, 2)]); 
hi_lo_rchAngMatched_dPCA_pullForce(:, 2) = cell2mat([leHi_leLo_dPCA_beh(:, 6); riHi_leLo_dPCA_beh(:, 6)]); 

scatter(hi_lo_rchAngMatched_dPCA_pullForce(:, 1), hi_lo_rchAngMatched_dPCA_pullForce(:, 2)) 
[corr.r_pullForce_rchAngleMatched, corr.p_pullForce_rchAngleMatched] = corrcoef(hi_lo_rchAngMatched_dPCA_pullForce(:, 1), hi_lo_rchAngMatched_dPCA_pullForce(:, 2)); 

%%



