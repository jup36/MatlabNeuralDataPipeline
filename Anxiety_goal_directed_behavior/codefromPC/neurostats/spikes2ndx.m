function [ndx,M] = spikes2ndx(t,p1,p2,p3)
%function [ndx,M] = spikes2ndx(s,edges)
%function [ndx,M] = spikes2ndx(s,tmin,tstep,tmax)
%
% converts a cell array (or a single vector) of spikes times s
% into indices into the vector of time bins with edges specified by edges
% or by edges = tmin:tstep:tmax, if the second usage is used
%
% the second usage is faster and does not require explicit creation of the 
% edge vector, but may not give identical results to the
% first usage when spikes land exactly on an edge (rounding errors)
%
% the convention is that the kth bin is [edges(k),edges(k+1))
%
% spike times < edges(1) or >= edges(end) are removed
%
% M = numel(edges)-1;

if nargin == 4

    p21 = 1./p2;
    M = fix((p3-p1).*p21);
    
    if iscell(t)
    
        ndx = cell(size(t));
    
        for k = 1:numel(t)
        
            tk = t{k};
            tk = fix((tk-p1).*p21)+1;
            
            ndx{k} = tk(tk >= 1 & tk <= M);
        
        end
    
    else
    
        t = fix((t(:)-p1).*p21)+1;
        ndx = t(t >= 1 & t <= M);
    
    end
    
elseif nargin == 2
    
    p1 = p1(:);
    M = numel(p1)-1;
    
    if iscell(t)
        
        ndx = cell(size(t));
        
        for k = 1:numel(t)
            
            tk = t{k}(:);
            
            if ~isempty(tk)
                [tmp,tk] = histc(tk,p1);
            
                ndx{k} = tk(tk >= 1 & tk <= M);

            end
            
        end
        
    else
        
        if ~isempty(t)
            [tmp,t] = histc(t(:),p1);
        end
        
        ndx = t(t >= 0 & t <= M);
        
    end
    
end
