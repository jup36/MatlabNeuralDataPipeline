function [ tqPosList, tqPosTypes ] = posTrqCombinations( jkvt )
%This function identifies the reach position and joystick torque combination for each trial
tqPosList = cell(size(jkvt,2),1);

posList = unique([jkvt(:).reachP1]);
tqList = unique([jkvt(:).pull_torque]);

for tt = 1:size(jkvt,2)
    tqPosList{tt} = sprintf('p%d_t%d',find(jkvt(tt).reachP1==posList),find(jkvt(tt).pull_torque==tqList));
end
tqPosTypes = unique(tqPosList);
end