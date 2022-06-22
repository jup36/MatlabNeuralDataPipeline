classdef TNC_OOP_FeatureExtractionInputParser < inputParser
    properties
        defaultNode        = 0;
        defaultSeg         =-1;
        defaultShank       =-1;
        defaultVerbose     = 0;
        defaultDisplay     = 0;
        defaultFakeData    = 0;
        defaultRawData     = 1;
        defaultsnrThresh   =-1;
        defaultSubtrCoMean = 0;
    end

    methods
        function obj = TNC_OOP_FeatureExtractionInputParser % class constructor 
            obj  = obj@inputParser;     
        end
        function parse(obj, dataFileName, arrayType, numSegs, varargin)
            compiled_code = 0;
            if ~isnan(str2double(numSegs))
                numSegs = str2double(numSegs);
                compiled_code = 1;
            end

            % Required inputs
            obj.addRequired( 'dataFileName', @ischar);
            obj.addRequired( 'arrayType',    @ischar);    
            obj.addRequired( 'numSegs'     , @isnumeric);    

            % If varargin{1} is a string convertable to number (i.e. this is compiled code)
            if size(varargin,1) > 0 && ~isnan(str2double(varargin{1}))
                varargin{1} = int32(str2double(varargin{1}));
            end

            % Optional, but position specific inputs
            obj.addOptional( 'iNode',     obj.defaultNode);

            posdot = findstr(dataFileName, '.');
            defaultTarget = dataFileName(1:(posdot-1));

            % Optional parameters
            if compiled_code
                obj.defaultSeg         = num2str(obj.defaultSeg); 
                obj.defaultShank       = num2str(obj.defaultShank);
                obj.defaultVerbose     = num2str(obj.defaultVerbose);
                obj.defaultDisplay     = num2str(obj.defaultDisplay);
                obj.defaultFakeData    = num2str(obj.defaultFakeData);
                obj.defaultRawData     = num2str(obj.defaultRawData);
                obj.defaultsnrThresh   = num2str(obj.defaultsnrThresh);
                obj.defaultSubtrCoMean = num2str(obj.defaultSubtrCoMean);
            end

            % Optional parameters
            obj.addParamValue('targetName',     defaultTarget, @ischar);    
            obj.addParamValue('iSeg',       obj.defaultSeg);
            obj.addParamValue('iShank',     obj.defaultShank);
            obj.addParamValue('verbose',    obj.defaultVerbose);
            obj.addParamValue('dispOn',     obj.defaultDisplay);
            obj.addParamValue('fakeData',   obj.defaultFakeData);
            obj.addParamValue('rawData',    obj.defaultRawData);
            obj.addParamValue('snr',        obj.defaultsnrThresh);
            obj.addParamValue('scm',        obj.defaultSubtrCoMean);

            obj.KeepUnmatched = true;

            parse@inputParser(obj, dataFileName, arrayType, numSegs, varargin{:});
        end
    end
end
