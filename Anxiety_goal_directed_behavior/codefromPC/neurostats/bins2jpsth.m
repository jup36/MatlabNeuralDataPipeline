function J = bins2jpsth(x,y,str,xpsth,ypsth,xpsthSD,ypsthSD)
% function J = bins2jpsth(x,y,str,xpsth,ypsth,xpsthSD,ypsthSD)
%
% computes the joint peri-stimulus time histogram of the binned spike
% trains x and y, where the kth trial is in the kth columns of x and y
%
% J(m,n) = average of the product of the firing rate of x at time m and y
% at time n
%
% str is a string that specifies what type of jpsth to compute
%
% the options are
% [] : uses default
% 'raw' : (default) the raw jpsth (this is the J decribed above)
% 'cov' : the shuffle corrected jpsth obtained by subtracting the
%         product of the psths; can use xpsth and ypsth
% 'covub' : T/(T-1) times the result of 'cov'; can use xpsth and ypsth
% 'cor': the result of 'cov' divided by the product of the psth standard
%        deviations; can use xpsth, ypsth, xpsthSD, ypsthSD
% 'corub': the result of 'covub' divided by the product of the psth standard
%        deviations (using the ub option); can use xpsth, ypsth, xpsthSD, ypsthSD
%        this should give the same result (up to rounding errors) as cor
%        (but cor and corub use different psthSD's)
% 'covNS' : amplitude non-stationary correction to 'cov'
% 'covubNS' : amplitude non-stationary correction to 'covub'
% 'corNS' : amplitude non-stationary correction to 'cor'
% 'corubNS' : amplitude non-stationary correction to 'covub'
%
% The arguments after str are optional.  The algorithm assumes that
%
% xpsth = bins2psth(x); ypsth = bins2psth(y);
% xpsthSD = bins2psthSD(x); ypsthSD = bins2psthSD(y); % cor option
% xpsthSD = bins2psthSD(x,'ub'); ypsthSD = bins2psthSD(y,'ub'); % corub option
%
% They psth's should be used for efficiency only and should not be used to 
% trick the algorithm into doing a different computation.  The reason is 
% that the algorithm chooses among different computational procedures, some
% of which make explicit assumptions about how the psth's were computed.  
% The psthSD's are always used in the same way and can be used for 
% efficiency or other purposes.

[nx,n] = size(x);
[ny,m] = size(y);

if n ~= m, error('x and y must have the same number of columns'), end

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

if covNSF || covubNSF || corNSF || corubNSF
    error('option unimplemented')
end
if ~(rawF || covF || covubF || corF || corubF || covNSF || covubNSF || corNSF || corubNSF)
    error('unrecognized option')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rawF
    J = XtimesYt(x,y);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4 || isempty(xpsth)

    xpsth = bins2psth(x);

end

if nargin < 5 || isempty(ypsth)

    ypsth = bins2psth(y);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psth stdard deviations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6 || isempty(xpsthSD)
    
    xpsthSD = [];
    
    if corF
        xpsthSD = bins2psthSD(x,'raw',xpsth);
    end
    
    if corubF
        xpsthSD = bins2psthSD(x,'ub',xpsth);
    end
        
end

if nargin < 7 || isempty(ypsthSD)

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
% cov, covub, cor, corub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if covF || covubF || corF || corubF

    % check to see which method is fastest (rough heuristics)
    % assume that logical will get evaluated by the fast XtimesYt algorithm

    if ((corF || corubF) && nx*n + ny*m < .25*(nx*ny)) || (~(islogical(x) && islogical(y)) && nx*n + ny*m < nx * ny)

        % subtract the means first
        if covF
            J = XtimesYt(x - xpsth*ones(1,n),y - ypsth*ones(1,n));
            return
        end
        if covubF
            J = XtimesYt(x - xpsth*ones(1,n),y - ypsth*ones(1,n),n-1);
            return
        end
        % subtract the means and divide by standard deviations first
        if corF
            J = XtimesYt((x - xpsth*ones(1,n))./(xpsthSD*ones(1,n)), ...
                (y - ypsth*ones(1,n))./(ypsthSD*ones(1,n))); 
            return
        end
        if corubF
            J = XtimesYt((x - xpsth*ones(1,n))./(xpsthSD*ones(1,n)), ...
                (y - ypsth*ones(1,n))./(ypsthSD*ones(1,n)),n-1); 
            return
        end
        
    end

    % subtract the means afterwards
    if covF
        J = XtimesYt(x,y) - xpsth*(ypsth.');
        return
    end
    if covubF
        J = XtimesYt(x,y,n-1) - (xpsth.*(n./(n-1)))*(ypsth.');
        return
    end

    % subtract the means and divide by standard deviations afterwards
    if corF
        J = (XtimesYt(x,y) - xpsth*(ypsth.'))./(xpsthSD*(ypsthSD.'));
        return
    end
    if corubF
        J = (XtimesYt(x,y,n-1) - (xpsth.*(n./(n-1)))*(ypsth.'))./(xpsthSD*(ypsthSD.'));
        return
    end

end

