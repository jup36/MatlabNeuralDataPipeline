%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

classdef TNC_OOP_SpikeSortingInputParser < inputParser
    properties (Constant)
        defaultMethod      = 'em';
        validMethods       = {'em', 'kk'};
        defaultGuessMethod = 'km';
        validGuessMethods  = {'rc','km'};  % 'random cell', 'k-means'
        defaultParNames    = 'min';
    end
    properties
        defaultReplica     = 0;
        defaultVerbose     = 0;
        defaultDebug       = 0;
        defaultModel       = 1;
        defaultPower       = 1;
        defaultDisplay     = 0;
        defaultFakeData    = 0;
        defaultMaxIter     = 500;  
        defaultNumFolds    = 1;
        defaultFold        = 1;
    end

    methods
        function obj = TNC_OOP_SpikeSortingInputParser % class constructor 
            obj  = obj@inputParser;                    % call base-class constructor
        end
        function parse(obj,  ftFile, iSeg, iShank, M, varargin)
            compiled_code = 0;
            % If compiled code     
            if ~isnan(str2double(iSeg))
                compiled_code = 1;
                iSeg   = str2double(iSeg);
                iShank = str2double(iShank);
                M      = str2double(M);
            end

            % Required inputs:
            obj.addRequired( 'ftFile', @ischar);
            obj.addRequired( 'iSeg',   @isnumeric);
            obj.addRequired( 'iShank', @isnumeric);    
            obj.addRequired( 'M',      @isnumeric);

            % If varargin{1} is a string convertable to number (i.e. this is compiled code)
            if size(varargin,1) > 0 && ~isnan(str2double(varargin{1}))
                varargin{1} = int32(str2double(varargin{1}));
            end

            % Optional, but position specific inputs
            obj.addOptional( 'iRepl', obj.defaultReplica);

            % Optional parameters
            checkMethod      = @(x) any(validatestring(x,obj.validMethods));
            checkGuessMethod = @(x) any(validatestring(x,obj.validGuessMethods));

            if compiled_code
                obj.defaultReplica  = num2str(obj.defaultReplica);
                obj.defaultVerbose  = num2str(obj.defaultVerbose);
                obj.defaultDebug    = num2str(obj.defaultDebug);   
                obj.defaultModel    = num2str(obj.defaultModel);
                obj.defaultPower    = num2str(obj.defaultPower);
                obj.defaultDisplay  = num2str(obj.defaultDisplay);
                obj.defaultFakeData = num2str(obj.defaultFakeData);
                obj.defaultMaxIter  = num2str(obj.defaultMaxIter);
                obj.defaultNumFolds = num2str(obj.defaultNumFolds);
                obj.defaultFold     = num2str(obj.defaultFold);
            end

            obj.addParamValue('method',      obj.defaultMethod,     checkMethod);
            obj.addParamValue('guess_method',obj.defaultGuessMethod,checkGuessMethod);
            obj.addParamValue('model',       obj.defaultModel       );
            obj.addParamValue('verbose',     obj.defaultVerbose     );
            obj.addParamValue('debug',       obj.defaultDebug       );
            obj.addParamValue('power',       obj.defaultPower       );
            obj.addParamValue('dispOn',      obj.defaultDisplay     );
            obj.addParamValue('fakeData',    obj.defaultFakeData    );
            obj.addParamValue('maxIter',     obj.defaultMaxIter     );
            obj.addParamValue('numFolds',    obj.defaultNumFolds    );
            obj.addParamValue('iFold',       obj.defaultFold        );
            obj.addParamValue('parNames',    obj.defaultParNames    );

            obj.KeepUnmatched = true;

            parse@inputParser(obj, ftFile, iSeg, iShank, M, varargin{:});
        end
    end
end
