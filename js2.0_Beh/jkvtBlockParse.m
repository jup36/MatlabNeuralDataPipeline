% parse trial/block information
function jkvt = jkvtBlockParse(jkvt)
tqd = diff([jkvt(1).pull_torque, jkvt(:).pull_torque]'); % torque change
p1d = diff([jkvt(1).reachP1, jkvt(:).reachP1]'); % position 1 change

blNumb = 1;  % block number
blType = []; % block type
% detect and parse block shifts across trials
for t = 1:size(jkvt,2)
    % first trial
    if t==1
        jkvt(t).blNumber = blNumb;
        jkvt(t).blType = 'first';
        jkvt(t).blShiftLogic = true;
    end
    % second trial and onward 
    if t>=2
        if tqd(t)==0 && p1d(t)==0 % if there's no change
            jkvt(t).blNumber = blNumb; % remains the same
            jkvt(t).blShiftLogic = false; % not shifted
            jkvt(t).blType = []; % block shift type
        else % if there's any change, classify the block shift type
            blNumb = blNumb+1;
            jkvt(t).blNumber = blNumb; % remains the same
            jkvt(t).blShiftLogic = true; % shifted
            % parse torque shift
            if tqd(t)>0
                tqi = 'tqUp';
            elseif tqd(t)<0
                tqi = 'tqDown';
            elseif tqd(t)==0
                tqi = 'tqSame';
            end
            % parse position shift
            if p1d(t)>0
                p1i = 'rightward';
            elseif p1d(t)<0
                p1i = 'leftward';
            elseif p1d(t)==0
                if jkvt(t).reachP1==min([jkvt(:).reachP1])
                    p1i = 'sameleft';
                elseif jkvt(t).reachP1==max([jkvt(:).reachP1])
                    p1i = 'sameright';
                end
            end
            jkvt(t).blType = strcat(tqi,p1i); % block shift type
        end
    end
end
end
