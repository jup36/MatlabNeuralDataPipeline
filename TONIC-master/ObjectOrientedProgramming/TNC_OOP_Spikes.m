%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

classdef TNC_OOP_Spikes
    properties 
        position;     
        value;
        channel;
        num_channels;
        neighbors;
    end

    methods
        function obj = TNC_OOP_Spikes(varargin)  % class constructor 
            if isempty(varargin) 
                return;
            end
            indCell  = varargin{1}; 
            valCell  = varargin{2};
            chanCell = varargin{3};
            nc       = varargin{4};
            N=length(indCell);
            for i=1:N
                obj(i).position = indCell{i};
                obj(i).value    = valCell{i};
                obj(i).channel  = chanCell{i};
                obj(i).num_channels = nc;
                obj(i).neighbors = zeros(1, nc);
            end
        end
        
        function [obj,idx] = sort(obj, varargin)
            %sort object array with respect to 'position' property

             [A,idx] = sort([obj.position],varargin{:}); 
             obj = obj(idx);
        end % function sort

        function [obj] = identify_neighbors(obj, dist)
            for n=1:numel(obj)
                i = n;
                while i >=1  && abs(obj(i).position - obj(n).position) < dist
                    ch = obj(i).channel;
                    if obj(n).neighbors(ch) == 0 
                        obj(n).neighbors(ch) = i;
                    else
                        i1 = obj(n).neighbors(ch);
                        if abs(obj(i1).position - obj(n).position) > abs(obj(i).position - obj(n).position)
                            obj(n).neighbors(ch) = i;
                        end
                    end       
                    i = i-1;   
                end  
                i = n;
                while i <= numel(obj) && abs(obj(i).position - obj(n).position) < dist
                    ch = obj(i).channel;
                    if obj(n).neighbors(ch) == 0
                        obj(n).neighbors(ch) = i;
                    else
                        i1 = obj(n).neighbors(ch);
                        if abs(obj(i1).position - obj(n).position) > abs(obj(i).position - obj(n).position)
                            obj(n).neighbors(ch) = i;
                        end
                    end
                    i = i+1;
                end
            end
        end % function identify_neighbors
        function [min_dist, max_dist] = compute_distances(obj)
            min_dist = -1*ones(1, numel(obj));
            max_dist = -1*ones(1, numel(obj));
            for n=2:(numel(min_dist)-1)
                min_dist(n) = min(abs(obj(n).position-obj(n-1).position),abs(obj(n).position-obj(n+1).position));
                max_dist(n) = max(abs(obj(n).position-obj(n-1).position),abs(obj(n).position-obj(n+1).position));
            end 
        end % compute_distance                 
        function len = length(obj)
            len = numel(obj);
        end
    end
end
