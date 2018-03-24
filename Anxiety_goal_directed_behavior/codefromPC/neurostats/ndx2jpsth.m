function J = ndx2jpsth(ndx,Mx,ndy,My,str,xpsth,ypsth,xpsthSD,ypsthSD)
% function J = ndx2jpsth(ndx,Mx,ndy,My,str,xpsth,ypsth,xpsthSD,ypsthSD)
%
% J = bins2jpsth(ndx2bins(ndx,Mx),ndx2bins(ndy,My),str,...)

% cast to cell for compatibility with ndx2bins

if ~iscell(ndx), ndx = {ndx}; end
if ~iscell(ndy), ndy = {ndy}; end

% input processing

T = numel(ndx);

if numel(ndy) ~= T, error('ndx and ndy must have the same number of elements'), end

if nargin < 5 || isempty(str)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rawF
    
    J = zeros(Mx,My);
    d = 1./T;
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % psths
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6 || isempty(xpsth)

        xpsth = ndx2psth(ndx,Mx);

    end

    if nargin < 7 || isempty(ypsth)

        ypsth = ndx2psth(ndy,My);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cov init
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if covF
        
        d = 1./T;
        
        J = (-xpsth)*(ypsth.');
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % covub init
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if covubF
        
        d = 1./(T-1);
        
        J = (xpsth.*(-T./(T-1)))*(ypsth.');
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw, cov, covub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rawF || covF || covubF

    for t = 1:T

        xt = ndx{t}(:);
        xt = xt(xt >= 1 & xt <= Mx);
        nx = numel(xt);

        yt = ndy{t}(:);
        yt = yt(yt >= 1 & yt <= My);
        ny = numel(yt);

        for k = 1:ny

            ytk = yt(k);

            for j = 1:nx

                xtj = xt(j);

                J(xtj,ytk) = J(xtj,ytk) + d;
            end
        end
    end

    return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psth stdard deviations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6 || isempty(xpsthSD)
    
    xpsthSD = [];
    
    if corF
        xpsthSD = ndx2psthSD(ndx,Mx,'raw',xpsth);
    end
    
    if corubF
        xpsthSD = ndx2psthSD(ndx,Mx,'ub',xpsth);
    end
        
end

if nargin < 7 || isempty(ypsthSD)

    ypsthSD = [];
    
    if corF
        ypsthSD = ndx2psthSD(ndy,My,'raw',ypsth);
    end
    
    if corubF
        ypsthSD = ndx2psthSD(ndy,My,'ub',ypsth);
    end

end

% deal with divide by zero
xpsthSD = xpsthSD + sqrt(eps(xpsthSD));
ypsthSD = ypsthSD + sqrt(eps(ypsthSD));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cor init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if corF

    dx = 1./(T.*xpsthSD);
    dy = 1./ypsthSD;

    J = (-xpsth./xpsthSD)*((ypsth./ypsthSD).');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corub init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if corubF

    dx = 1./((T-1).*xpsthSD);
    dy = 1./ypsthSD;

    J = ((xpsth.*(-T./(T-1)))./xpsthSD)*((ypsth./ypsthSD).');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cor, corub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if corF || corubF

    for t = 1:T

        xt = ndx{t}(:);
        xt = xt(xt >= 1 & xt <= Mx);
        nx = numel(xt);

        yt = ndy{t}(:);
        yt = yt(yt >= 1 & yt <= My);
        ny = numel(yt);

        for k = 1:ny

            ytk = yt(k);
            dyk = dy(ytk);

            for j = 1:nx

                xtj = xt(j);
                d = dx(xtj).*dyk;

                J(xtj,ytk) = J(xtj,ytk) + d;
            end
        end
    end

    return

end