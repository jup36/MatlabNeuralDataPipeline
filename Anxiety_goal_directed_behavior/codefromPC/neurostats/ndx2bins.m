function b = ndx2bins(ndx,M)
% function b = ndx2bins(ndx,M)
%
% ndx is a cell array of indices into a vector of length M
%
% b(k,j) is the number of times index j occurs in ndx{k}
%
% indices outside of 1,...,M are ignored

binflag = true;

if iscell(ndx)

    T = numel(ndx);

    b = zeros(M,T);

    for t = 1:T

        ndxt = ndx{t}(:);
        ndxt = ndxt(ndxt >= 1 & ndxt <= M);

        numndxt = numel(ndxt);
        
        for k = 1:numndxt

            ndxtk = ndxt(k);
            
            b(ndxtk,t) = b(ndxtk,t) + 1;

        end

        % logical checking
        if binflag
            if numndxt < M/2
                for k = 1:numndxt
                    if b(ndxt(k),t) > 1
                        binflag = false;
                        break
                    end
                end
            else
                for k = 1:M
                    if b(k,t) > 1
                        binflag = false;
                        break
                    end
                end
            end
        end

    end

else

    ndx = ndx(ndx >= 1 & ndx <= M);
    
    n = numel(ndx);

    b = zeros(M,1);

    for k = 1:n

        ndxk = ndx(k);

        b(ndxk) = b(ndxk) + 1;

    end
    
    % logical checking
    for k = 1:M
        if n < M/2
            for k = 1:n
                if b(ndx(k)) > 1
                    binflag = false;
                    break
                end
            end
        else
            if b(k,t) > 1
                binflag = false;
                break
            end
        end
    end

end

% convert to logical when possible
if  binflag,  b = b == 1; end
