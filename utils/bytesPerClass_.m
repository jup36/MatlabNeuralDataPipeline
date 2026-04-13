function b = bytesPerClass_(cls)
switch cls
    case {'double'}
        b = 8;
    case {'single'}
        b = 4;
    case {'logical','char','int8','uint8'}
        b = 1;
    case {'int16','uint16'}
        b = 2;
    case {'int32','uint32'}
        b = 4;
    case {'int64','uint64'}
        b = 8;
    otherwise
        b = 8; % conservative default
end
end