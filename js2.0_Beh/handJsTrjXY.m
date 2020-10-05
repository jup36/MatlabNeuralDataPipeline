function [hXYinit,hXYutP1,hXYutPe,hXYtort,forceTrjToPeak,jXY1,rwdTrI] = handJsTrjXY(handXyzC,handXY1med,jsXY1,forceTrjC,rwdTrI)

nBaseBins = 50; % 20msx50bins = 1000ms before the reachStart  
nBackFromMaxForce = 11; % go back e.g. 11 bins from the peak force point

% get valid hand and joystick trajectories and normalize them to the common reference point
handXY1med = -handXY1med;
valHandI = cell2mat(cellfun(@(a) ~isempty(a), handXyzC, 'un', 0)); 
valHandXYc = cellfun(@(a) -a(1:2,:), handXyzC(valHandI), 'un', 0); 
valHxyC = cellfun(@(a) a-repmat(handXY1med,[1,size(a,2)]), valHandXYc, 'un', 0); 
hXyC = cell(1,length(handXyzC)); 
hXyC(valHandI) = valHxyC;  

valCellI = cell2mat(cellfun(@(a) ~isempty(a), hXyC, 'un', 0)); 
localhXY1med = nanmedian(cell2mat(cellfun(@(a) a(:,1), hXyC(valCellI), 'un', 0)),2); 

jXY1 = -jsXY1-handXY1med; 

%distV = @(a,b) sqrt((a(:,1)-b(:,1)).^2+(a(:,2)-b(:,2)).^2); % calculate the distance between two equal-length vectors  

for j = 1:length(hXyC) % trials 
    %% hand trajectory 
    if ~isempty(hXyC{j})
        hXY = hXyC{j}; % hand trajectory (cut at the pull end)
        if rwdTrI(j) % success trial
            [~,p1h] = min(sum((hXY-repmat(jXY1,[1,size(hXY,2)])).^2),[],2);
            %hXYutP1{1,j} = hXY(:,1:max(nBaseBins,p1h)); % hand XY trajectory upto the pull start
            hXYutP1{1,j} = hXY(:,1:p1h); % hand XY trajectory upto the pull start
            hXYutPe{1,j} = hXY; % hand XY trajectory upto the pull end
            hXYinit{1,j} = nanmedian(hXY(:,1:max(nBaseBins-3,1)),2); % initial hand position 
            % get tortuosity = distance (path length) / displacement (endpoints)
            if size(hXYutP1{1,j},2)>nBaseBins
                %reachTrj = hXYutP1{1,j}(:,1:end); 
                reachTrj = hXYutP1{1,j}(:,max(1,nBaseBins-3):end); % tightly (right before the reachStart) around the reach portion of the trajectory           
            else
                reachTrj = hXYutP1{1,j}(:,1:end); 
            end
            displace = sqrt(sum((reachTrj(:,end)-reachTrj(:,1)).^2)); % displacement: equivalently, sum(sqrt(sum(diff(shortTrj).^2,2)))
            distance = sum(sqrt(sum(diff(reachTrj).^2,2))); % actual distance
            hXYtort{1,j} = min(15,distance/displace); % tortuosity
        else % failed trial
            hXYutP1{1,j} = hXY; 
            hXYutPe{1,j} = hXY; 
            hXYtort{1,j} = nan; 
            hXYinit{1,j} = nanmedian(hXY(:,1:max(nBaseBins-3,1)),2); % initial hand position 
        end   
    else % empty trial
        hXYutP1{1,j} = nan; 
        hXYutPe{1,j} = nan; 
        hXYtort{1,j} = nan; 
        hXYinit{1,j} = nan; % initial hand position 
    end
    %% Force trajectory aligned to the max force point
    if ~isempty(forceTrjC{j})
        forceTrj = forceTrjC{j}; % force trajectory (cut at the pull end)
        if rwdTrI(j) % success trial
           [~,minFI]=min(forceTrj); 
           forceTrjToPeak{1,j} = forceTrj(max(1,minFI-nBackFromMaxForce+1):minFI); 
        else % failed trial
            forceTrjToPeak{1,j} = nan;  
        end
    else % empty trial
        forceTrjToPeak{1,j} = nan; 
    end   
end
end