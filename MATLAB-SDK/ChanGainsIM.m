% =========================================================
% Return gain arrays for imec channels.
%
% Index into these with original (acquired) channel IDs.
%
function [APgain,LFgain] = ChanGainsIM(meta)
    C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
        'EndOfLine', ')', 'HeaderLines', 1 );
    APgain = double(cell2mat(C(1)));
    LFgain = double(cell2mat(C(2)));
end % ChanGainsIM
