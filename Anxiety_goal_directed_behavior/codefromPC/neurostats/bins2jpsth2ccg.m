function [c,N] = bins2jpsth2ccg(x,y,str,w,cstr,fftflag,xpsth,ypsth,xpsthSD,ypsthSD)
% function [c,N] = bins2jpsth2ccg(x,y,str,w,cstr,fftflag,xpsth,ypsth,xpsthSD,ypsthSD)
%
% J = bins2jpsth(x,y,str,xpsth,ypsth,xpsthSD,ypsthSD);
% [c,N] = jpsth2ccg(J,w,cstr);
%
% the arguments following cstr are optional
%
% This can provide a shortcut for computing ccg's, especially when w is
% small, by avoiding explicit construction of the jpsth.
%
% The optional argument fftflag (default = false) can be set to true to
% compute ccg's using the fft.  Sometimes this can be faster, especially
% for many bins and few trials.  Changing w does not affect the speed of
% the fft algorithm.

%%%%%%%%%%%%%%%%%%%%%%
% jpsth preprocessing
%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(x);

if ~isequal(size(x),size(y)), error('x and y must be the same size'), end

if nargin < 3 || isempty(str)
    str = 'raw';
end

str = lower(str);

rawF = strcmp(str,'raw');
covF = strcmp(str,'cov');
covubF = strcmp(str,'covub');
corF = strcmp(str,'cor');
corubF = strcmp(str,'corub');
covNSF = strcmp(str,'covns');
covubNSF = strcmp(str,'covubns');
corNSF = strcmp(str,'corns');
corubNSF = strcmp(str,'corubns');

%%%%%%%%%%%%%%%%%%%%%%
% ccg preprocessing
%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4 || isempty(w)
    w = m;
end

if nargin < 5 || isempty(cstr)
    cstr = 'full';
end

cstr = lower(cstr);

fullF = strcmp(cstr,'full');
rowvalidF = strcmp(cstr,'rowvalid');
colvalidF = strcmp(cstr,'colvalid');

if ~(fullF || rowvalidF || colvalidF)
    error('unknown cstr option')
end

if w < 1 || w > m
    error('w must between 1 and M')
end

if ~fullF && 2*(w-1) >= m
    error('w is too large for the valid options')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ccg setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6 || isempty(fftflag)
    fftflag = false;
end
if numel(fftflag) ~= 1 || ~islogical(fftflag)
    error('fftflag must be a logical scalar')
end

d = (1-w:w-1).';

if fullF

    rmin = 1;
    rmax = m;
    cmin = 1;
    cmax = m;
    md = m-abs(d);

elseif rowvalidF

    rmin = w;
    rmax = m-w+1;
    cmin = 1;
    cmax = m;
    md = m-2*w+2;

    if fftflag 
        % if using fft, go ahead and remove excluded data
        x(1:rmin-1,:) = 0;
        x(rmax+1:m,:) = 0;
    end 
        
elseif colvalidF

    rmin = 1;
    rmax = m;
    cmin = w;
    cmax = m-w+1;
    md = m-2*w+2;

    if fftflag
        % if using fft, go ahead and remove excluded data
        y(1:cmin-1,:) = 0;
        y(cmax+1:m,:) = 0;
    end

else
    
    error('unrecognized cstr option')
    
end

if nargout > 1
    
    if numel(md) == 1
        N = repmat(md,size(d));
    else
        N = md;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw ccg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fftflag && rawF

    c = sum(ifft(fft([flipud(x) ; zeros(m-1,n)]).*fft([y ; zeros(m-1,n)])),2);
    c = c(d+m) ./ (n.*md);
    return
end

if (rawF || covF || covubF) && ~fftflag

    c = DiagSumX(x(:,1),y(:,1),d,rmin,rmax,cmin,cmax);

    for k = 2:n

        c = c + DiagSumX(x(:,k),y(:,k),d,rmin,rmax,cmin,cmax);

    end

    if rawF
        c = c ./ (n.*md);
        return
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 7 || isempty(xpsth)

    xpsth = bins2psth(x);

end

if nargin < 8 || isempty(ypsth)

    ypsth = bins2psth(y);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cov, covub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if covF

    if fftflag
        c = sum(ifft(fft([flipud(x-xpsth*ones(1,n)) ; zeros(m-1,n)]).*fft([y-ypsth*ones(1,n) ; zeros(m-1,n)])),2);
        c = c(d+m) ./ (n.*md);
    else
        c = (c./n-DiagSumX(xpsth,ypsth,d,rmin,rmax,cmin,cmax)) ./ md;
    end
    return
end

if covubF
    
    if fftflag
        c = sum(ifft(fft([flipud(x-xpsth*ones(1,n)) ; zeros(m-1,n)]).*fft([y-ypsth*ones(1,n) ; zeros(m-1,n)])),2);
        c = c(d+m) ./ ((n-1).*md);
    else
        c = (c./(n-1)-DiagSumX(xpsth,ypsth,d,rmin,rmax,cmin,cmax).*(n/(n-1))) ./ md;
    end
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psth stdard deviations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9 || isempty(xpsthSD)
    
    xpsthSD = [];
    
    if corF
        xpsthSD = bins2psthSD(x,'raw',xpsth);
    end
    
    if corubF
        xpsthSD = bins2psthSD(x,'ub',xpsth);
    end
        
end

if nargin < 10 || isempty(ypsthSD)

    ypsthSD = [];
    
    if corF
        ypsthSD = bins2psthSD(y,'raw',ypsth);
    end
    
    if corubF
        ypsthSD = bins2psthSD(y,'ub',ypsth);
    end

end

% deal with divide by zero
xpsthSD = xpsthSD + sqrt(eps(xpsthSD));
ypsthSD = ypsthSD + sqrt(eps(ypsthSD));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cov, covub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if corF || corubF

    if fftflag
        
        c = sum(ifft(fft([flipud((x-xpsth*ones(1,n))./(xpsthSD*ones(1,n))) ; zeros(m-1,n)]).*fft([(y-ypsth*ones(1,n))./(ypsthSD*ones(1,n)) ; zeros(m-1,n)])),2);

        if corF
            c = c(d+m) ./ (n.*md);
            return
        end

        if corubF
            c = c(d+m) ./ ((n-1).*md);
            return
        end
        
    end

    x = x./(xpsthSD*ones(1,n));
    y = y./(ypsthSD*ones(1,n));

    c = DiagSumX(x(:,1),y(:,1),d,rmin,rmax,cmin,cmax);

    for k = 2:n
        c = c + DiagSumX(x(:,k),y(:,k),d,rmin,rmax,cmin,cmax);
    end
    
    if corF
        c = (c./n-DiagSumX(xpsth./xpsthSD,ypsth./ypsthSD,d,rmin,rmax,cmin,cmax)) ./ md;
        return
    end
    
    if corubF
        c = (c./(n-1)-DiagSumX(xpsth./xpsthSD,ypsth./ypsthSD,d,rmin,rmax,cmin,cmax).*(n/(n-1))) ./ md;
        return
    end
    
end

error('unrecognized str option')
