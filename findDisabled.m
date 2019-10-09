function [exChan] = findDisabled(meta)     
    % read in the shank map    
    C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
    enabled = double(cell2mat(C(4)));
    % There's an entry in the shank map for each saved channel.
    % Get the array of saved channels:
    chan = OriginalChans(meta);
    % Find out how many are non-SY chans
    [AP,~,~] = ChannelCountsIM(meta);
    exChan = [];
    for i = 1:AP
        if enabled(i) == 0
            exChan = [exChan, chan(i)];
        end
    end
end % findDisabled
