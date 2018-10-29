function [jsState, jsStateDt] = readStepperEncoder(pinAstate, pinBstate)
%This helper function reads the quadrature encoded signals from the encoder
% and extracts the trajectory of the 1-d joystick movement.
% pinABdt: any changes (rise or fall) detected by pin A or pin B (1st row: time points of change, 2nd row: label for change; pin A rise/fall: 1/-1, pin B rise/fall: 2/-2)
% pinAstate: the entire binarized pin A state
% pinBstate: the entire binarized pin B state
% For more information about quadrature encoded signals, https://www.pjrc.com/teensy/td_libs_Encoder.html

pinArise = find(diff([pinAstate(1), pinAstate])==1);
pinArise(2,:) = 1;  % pin A rise
pinAfall = find(diff([pinAstate(1), pinAstate])==-1);
pinAfall(2,:) = -1; % pin A fall

pinBrise = find(diff([pinBstate(1), pinBstate])==1);
pinBrise(2,:) = 2;  % pin B rise
pinBfall = find(diff([pinBstate(1), pinBstate])==-1);
pinBfall(2,:) = -2; % pin B fall

pinABdt = sortrows([pinArise, pinAfall, pinBrise, pinBfall]',1)'; % 1st row: event time, 2nd row: event label 

jsState = zeros(1,length(pinAstate)); % jsState accumulated on the whole length of the sequence
jsStateDt = nan(1,length(pinABdt));   % jsState tracking changes only

for pts = 1:size(pinABdt,2)-1
    
    if pinABdt(2,pts)==1 % pin A rise
        if pinABdt(2,pts+1)==2      % pin B rise follows
            jsState(1,pinABdt(1,pts))=-1;   % pull -1
            if pts==length(pinABdt)-1
                jsState(1,pinABdt(1,pts+1))=-1; % pull -1
            end
        elseif pinABdt(2,pts+1)==-2 % pin B fall follows
            jsState(1,pinABdt(1,pts))=1;    % push +1
            if pts==length(pinABdt)-1
                jsState(1,pinABdt(1,pts+1))=1;  % push +1
            end
        elseif pinABdt(2,pts+1)==-1 % pin A fall back down follows
            if pinBstate(pinABdt(1,pts))==0 % this means pull then push
                jsState(1,pinABdt(1,pts))=-1;   % pull -1
                if pts==length(pinABdt)-1
                    jsState(1,pinABdt(1,pts+1))=1;  % push +1
                end
            elseif pinBstate(pinABdt(1,pts))==1 % this means push then pull
                jsState(1,pinABdt(1,pts))=1;    % push 1
                if pts==length(pinABdt)-1
                    jsState(1,pinABdt(1,pts+1))=1;  % pull -1
                end
            end
        end
    end
    
    if pinABdt(2,pts)==-1 % pin A fall
        if pinABdt(2,pts+1)==2      % pin B rise follows
            jsState(1,pinABdt(1,pts))=1;   % push +1
            if pts==length(pinABdt)-1
                jsState(1,pinABdt(1,pts+1))=1; % push +1
            end
        elseif pinABdt(2,pts+1)==-2 % pin B fall follows
            jsState(1,pinABdt(1,pts))=-1;    % pull -1
            if pts==length(pinABdt)-1
                jsState(1,pinABdt(1,pts+1))=-1;  % pull -1
            end
        elseif pinABdt(2,pts+1)==1 % pin A rise back up follows
            if pinBstate(pinABdt(1,pts))==0 % this means push then pull
                jsState(1,pinABdt(1,pts))=1;    % push 1
                if pts==length(pinABdt)-1
                    jsState(1,pinABdt(1,pts+1))=-1; % pull -1
                end
            elseif pinBstate(pinABdt(1,pts))==1 % this means pull then push
                jsState(1,pinABdt(1,pts))=-1;   % pull -1
                if pts==length(pinABdt)-1
                    jsState(1,pinABdt(1,pts+1))=1;  % push 1
                end
            end
        end
    end
    
    if pinABdt(2,pts)==2 % pin B rise
        if pinABdt(2,pts+1)==1      % pin A rise follows
            jsState(1,pinABdt(1,pts))=1;   % push +1
            if pts==length(pinABdt)-1
                jsState(1,pinABdt(1,pts+1))=1; % push +1
            end
        elseif pinABdt(2,pts+1)==-1 % pin A fall follows
            jsState(1,pinABdt(1,pts))=-1;    % pull -1
            if pts==length(pinABdt)-1
                jsState(1,pinABdt(1,pts+1))=-1;  % pull -1
            end
        elseif pinABdt(2,pts+1)==-2 % pin B fall back down follows
            if pinAstate(pinABdt(1,pts))==0 % this means push then pull
                jsState(1,pinABdt(1,pts))=1;    % push 1
                if pts==length(pinABdt)-1
                    jsState(1,pinABdt(1,pts+1))=-1; % pull -1
                end
            elseif pinAstate(pinABdt(1,pts))==1 % this means pull then push
                jsState(1,pinABdt(1,pts))=-1;   % pull -1
                if pts==length(pinABdt)-1
                    jsState(1,pinABdt(1,pts+1))=1;  % push 1
                end
            end
        end
    end
    
    if pinABdt(2,pts)==-2 % pin B fall
        if pinABdt(2,pts+1)==1      % pin A rise follows
            jsState(1,pinABdt(1,pts))=-1;   % push +1
            if pts==length(pinABdt)-1
                jsState(1,pinABdt(1,pts+1))=-1; % push +1
            end
        elseif pinABdt(2,pts+1)==-1 % pin A fall follows
            jsState(1,pinABdt(1,pts))=1;    % push +1
            if pts==length(pinABdt)-1
                jsState(1,pinABdt(1,pts+1))=1;  % push +1
            end
        elseif pinABdt(2,pts+1)==2  % pin B rise back up follows
            if pinAstate(pinABdt(1,pts))==0 % this means pull then push
                jsState(1,pinABdt(1,pts))=-1;  % pull -1
                if pts==length(pinABdt)-1
                    jsState(1,pinABdt(1,pts+1))=1; % push +1
                end
            elseif pinAstate(pinABdt(1,pts))==1 % this means push then pull
                jsState(1,pinABdt(1,pts))=1;    % push 1
                if pts==length(pinABdt)-1
                    jsState(1,pinABdt(1,pts+1))=-1; % pull -1
                end
            end
        end
    end
    
    if pts < length(pinABdt)-1
        jsStateDt(1,pts) = jsState(1,pinABdt(1,pts));
    elseif pts==length(pinABdt)-1
        jsStateDt(1,pts) = jsState(1,pinABdt(1,pts));
        jsStateDt(1,pts+1) = jsState(1,pinABdt(1,pts+1));
    end
end

end

