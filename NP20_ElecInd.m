function [elecInd, shankInd, bankMask, connected] = NP20_ElecInd(meta)     
    pType = str2num(meta.imDatPrb_type);
    if pType == 21
        % Single shank probe
        % imro table entries: (channel, bank, refType, electrode #)
        C = textscan(meta.imroTbl, '(%*s %d %*s %d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        elecInd = int32(cell2mat(C(2)));
        bankMask = int32(cell2mat(C(1)));
        shankInd = zeros(size(elecInd));
        connected = ones(size(elecInd));
        exChan = findDisabled(meta);        
        for i = 1:numel(exChan)        	
            connected(elecInd == exChan(i)) = 0;
        end
        
    else
        % 4 shank probe
        % imro table entries: (channel, shank, bank, refType, electrode #)
        C = textscan(meta.imroTbl, '(%d %d %d %*s %d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        chan = double(cell2mat(C(1)));
        elecInd = int32(cell2mat(C(4)));
        bankMask = int32(cell2mat(C(3)));
        shankInd = double(cell2mat(C(2)));
        connected = ones(size(chan));
        exChan = findDisabled(meta);
        %exChan = [127];
        for i = 1:numel(exChan)       	
            connected(chan == exChan(i)) = 0;
        end
    end
end % NP20_ElecInd
