function p = findArea(width,varargin)

   defaultUnits = 'inches';

   p = inputParser;

   addRequired(p,'width');
   addParameter(p,'units',defaultUnits);
   
   parse(p,width,varargin{:});
   
end