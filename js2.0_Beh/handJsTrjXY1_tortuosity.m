%% Helper functions
function handJsTrjXY1_tortuosity(filePath)

cd(filePath)
file = dir('**/*js2p0_tbytSpkHandJsTrjBin_*.mat');
load(fullfile(file.folder, file.name), 'ss', 'jkvt')

valTrs = cell2mat(cellfun(@(a) ~isempty(a), {ss(:).hTrjB}, 'un', 0));
rwdTrI = [jkvt(:).rewarded];
valRwdTrs = find(valTrs & rwdTrI); % valid rewarded trials

nBaseBins = 50;

for j = 1:length(valRwdTrs)
    hXY_to_p1 = ss(valRwdTrs(j)).hXY_to_p1;
    if ~isempty(hXY_to_p1)
        % get tortuosity = distance (path length) / displacement (endpoints)
        if size(hXY_to_p1,2)>nBaseBins
            reachTrj = hXY_to_p1(:,max(1,nBaseBins-3):end); % tightly (right before the reachStart) around the reach portion of the trajectory
        else
            reachTrj = hXY_to_p1(:,1:end);
        end
        displace = sqrt(sum((reachTrj(:,end)-reachTrj(:,1)).^2)); % displacement: equivalently, sum(sqrt(sum(diff(shortTrj).^2,2)))
        distance = sum(sqrt(sum(diff(reachTrj).^2,2))); % actual distance
        ss(valRwdTrs(j)).hXYtort = min(15,distance/displace); % tortuosity
    end
end

save(fullfile(file.folder, file.name), 'ss', '-append')

end