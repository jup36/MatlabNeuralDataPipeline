filePath = {'/Volumes/Extreme SSD/js2p0/WR37_022619', ... % Cg recording contra-Cg silencing, Trj checked            % GOOD TO BE USED
            '/Volumes/Extreme SSD/js2p0/WR38_052319', ... % Cg recording contra-Cg silencing, Trj checked            % BAD (Crashed)
            '/Volumes/Extreme SSD/js2p0/WR38_052419', ... % Corticostriatal recording M1 silencing, Trj checked       
            '/Volumes/Extreme SSD/js2p0/WR39_100219', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/Extreme SSD/js2p0/WR40_081919', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/Extreme SSD/js2p0/WR40_082019', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/Extreme SSD/js2p0/WR44_031020'};    % Dual recording with contra Cg delayed silencing, Trj checked % GOOD

%% MOp WR38_052419 use unit 109 as an example silencing unit 
load('/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles/binSpkCountSTRCTXWR38_052419.mat', 'tagLaser')

WR38_052419_combined_unit1 = cellfun(@(a, b, c) sort([a, b, c]), tagLaser.SpkTimes{117, 1}, tagLaser.SpkTimes{118, 1}, tagLaser.SpkTimes{119, 1}, 'un', 0 ); 

% pre-tagging vs tagging
fig1 = spikeRasterGramm( [5e3 5e3], {'taggingTrials'}, [2e3 2e3], WR38_052419_combined_unit1);
pbaspect([1 1 1])
print(fig1, fullfile('/Volumes/Extreme SSD/js2p0/WR38_052419','Figure',strcat('tagLaser','WR38_052419_combined_unit1')),'-dpdf','-bestfit')

%% MOs
load('/Volumes/Extreme SSD/js2p0/WR37_022619/Matfiles/binSpkCountCgWR37_022619.mat', 'tagLaser')

%% MOp during Cg silencing
load('/Volumes/Extreme SSD/js2p0/WR39_100219/Matfiles/binSpkCountCTXWR39_100219.mat', 'tagLaser')