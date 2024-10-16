function ctpath(CT0,varargin)
%  CTPATH(CT,{CMODEL},{MARKERSIZE},{OPTIONS})
%  Plot the trajectory of a color table in a selected color space.
%  Create a datacursor to help inspection of the trajectory.
%
%  CT is a Mx3 color table of any standard image class
%    The values in CT must be scaled appropriately for the class of CT.
%  CMODEL specifies the color space to use when plotting (default 'rgb')
%    'rgb', 'lrgb', 'hsv', 'hsl', 'hsi', 
%    'ypbpr', 'luv', 'lab', 'srlab', 'oklab', 
%    'lchbr', 'lchuv', 'lchab', 'lchsr', 'lchok', 
%    'hsy', 'huslab', 'husluv', 'hsyp', 'huslpab', 'huslpuv'
%  MARKERSIZE optionally specifies the plot marker size (default 30)
%  OPTIONS includes the following keys:
%    'noline' omits the lines connecting the markers
%    'invert' inverts the marker colors for use on an inverted display.
%    '3D' plots the table in 3D (default)
%    '2D' plots the table in a simplified 2D fashion
%
%  See also: csview, ccmap, makect, rgbplot

bigtable = { % create as one big array for sake of compactness, readability, ease of maintenance
% csnames	plotaxlabel		xyz			plotrange
'rgb',		'R','G','B',	[1 2 3],	[0 1; 0 1; 0 1];
'lrgb',		'R','G','B',	[1 2 3],	[0 1; 0 1; 0 1];

'hsv',		'H','S','V',	[1 2 3],	[0 360; 0 1; 0 1];
'hsl',		'H','S','L',	[1 2 3],	[0 360; 0 1; 0 1];
'hsi',		'H','S','I',	[1 2 3],	[0 360; 0 1; 0 1];

'ypbpr',	'Pb','Pr','Y',	[2 3 1],	[-1 1; -1 1; 0 1];
'luv',		'U','V','L',	[2 3 1],	[-180 180; -180 180; 0 100];
'lab',		'A','B','L',	[2 3 1],	[-140 140; -140 140; 0 100];
'srlab',	'A','B','L',	[2 3 1],	[-110 110; -110 110; 0 100];
'oklab',	'A','B','L',	[2 3 1],	[-35 35; -35 35; 0 100];

'lchbr',	'L','C','H',	[1 2 3],	[0 1; 0 0.55; 0 360];
'lchuv',	'L','C','H',	[1 2 3],	[0 100; 0 180; 0 360];
'lchab',	'L','C','H',	[1 2 3],	[0 100; 0 135; 0 360];
'lchsr',	'L','C','H',	[1 2 3],	[0 100; 0 105; 0 360];
'lchok',	'L','C','H',	[1 2 3],	[0 100; 0 35; 0 360];

'hsy',		'H','S','Y',	[1 2 3],	[0 360; 0 1; 0 1];
'huslab',	'H','S','L',	[1 2 3],	[0 360; 0 100; 0 100];
'husluv',	'H','S','L',	[1 2 3],	[0 360; 0 100; 0 100];
'hsyp',		'H','S','Y',	[1 2 3],	[0 360; 0 1; 0 1];
'huslpab',	'H','S','L',	[1 2 3],	[0 360; 0 100; 0 100];
'huslpuv',	'H','S','L',	[1 2 3],	[0 360; 0 100; 0 100];
};


% defaults
spc = 'rgb';
invert = false;
markersize = 30;
noline = false;
is3D = true;

% process inputs
if nargin>1
	for k = 1:numel(varargin)
		thisarg = varargin{k};
		if ischar(thisarg)
			switch lower(thisarg)
				case '2d'
					is3D = false;
				case '3d'
					is3D = true;
				case 'invert'
					invert = true;
				case 'noline'
					noline = true;
				otherwise
					spc = thisarg;
			end
		elseif isnumeric(thisarg)
			markersize = thisarg(1);
		end
	end
end

% fetch space name, sanitize, validate
spc = lower(spc);
spc = spc(spc ~= ' ');
csnames = bigtable(:,1);
[validcsname,idx] = ismember(spc,csnames);
if ~validcsname
	error('CTPATH: unknown colorspace name %s \n',spc)
end
selectedcs = mod(idx-1,size(csnames,1))+1; % index of selected CS

% get associated params
plotaxlabels = bigtable(selectedcs,2:4);
xyz = bigtable{selectedcs,5};
plotrange = bigtable{selectedcs,6};

% convert CT and prepare
if ndims(CT0)>2 || size(CT0,2)~=3 %#ok<ISMAT>
	error('CTPATH: expected CT to be a Mx3 matrix')
end
CT0 = imcast(CT0,'double');
CT = ctflop(fromrgb(ctflop(CT0),spc));

if is3D
	CT = CT(:,xyz); % reorder axes
	if invert
		CT0 = 1-CT0;
	end

	% plot things
	hp = plot3(CT(:,1),CT(:,2),CT(:,3),':'); hold on; grid on
	set(hp,'color',[1 1 1]*0.5);
	if noline
		set(hp,'visible',false);
	end
	scatter3(CT(:,1),CT(:,2),CT(:,3),markersize,CT0,'filled')

	xlim(plotrange(1,:))
	ylim(plotrange(2,:))
	zlim(plotrange(3,:))

	xlabel(plotaxlabels{1})
	ylabel(plotaxlabels{2})
	zlabel(plotaxlabels{3})

	% i'm going to assume nobody else wants a dark bg on a noninverted display
	% if ~invert
	% 	set(gca,'color',[0 0 0],'gridcolor',[1 1 1]*0.85)
	% end

else
	% reorder labels, limits
	plotaxlabels(xyz) = plotaxlabels;
	plotrange(xyz,:) = plotrange;
	
	% curve colors
	if invert
		lc = {'c','m','y'};
	else
		lc = {'r','g','b'};
	end

	% plot first channel to set up main axes
	x = 1:size(CT,1);
	plot(x,CT(:,1),lc{1}); hold on;
	set(gca,'xgrid','on');
	set(gca,'ycolor',lc{1});
	ylim(plotrange(1,:))
	xlim([1 size(CT,1)])
	
	% plot other channels
	for kc = 2:3
		addaxis(x,CT(:,kc),plotrange(kc,:),lc{kc});
	end
	
	xlabel('CT index');
	hax = getaddaxisdata(gca,'axisdata');
	hax = [hax{1:2}];
	hax = [gca hax(1,:)];
	
	x0 = 0.1; % spacing between outermost rulers and figure edge
	xgap = 0.08; % spacing between rulers
	fontsize = 8; % base font size for all rulers/labels
	
	% get info about ruler locations
	nlr = nnz(strcmp({hax.YAxisLocation},'left'));
	nrr = numel(hax) - nlr;
	
	% calculate x-offsets for rulers, reorder to match axes order
	xposl = x0+(0:nlr-1)*xgap;
	xposr = 1-(x0+(nrr-1:-1:0)*xgap);
	xpos = zeros(1,numel(hax));
	xpos(1:2:nlr*2) = fliplr(xposl);
	xpos(2:2:nrr*2) = xposr;
	
	% width of primary axes
	wax = 1 - ((nlr-1)+(nrr-1))*xgap - 2*x0;
	
	% apply properties
	for k = 1:numel(hax)
		hax(k).Position(1) = xpos(k);
		hax(k).YLabel.String = plotaxlabels{k};
		hax(k).FontSize = fontsize;
	end
	hax(1).Position(3) = wax;
end

% create custom datacursor to view index
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@myupdatefcn);

% this is kind of flaky, but i don't care that much
function txt = myupdatefcn(~,event_obj)
	% get rid of old tips
	alltips = findall(gcf,'Type','hggroup');
	if numel(alltips)>2
		delete(alltips(3:end))
	end
	
	if is3D
		pos = get(event_obj,'Position');
		didx = get(event_obj, 'DataIndex');
		[~,xyzinv] = sort(xyz);
		txt = {sprintf('%s: %.4f',plotaxlabels{xyzinv(1)},pos(xyzinv(1))),...
			   sprintf('%s: %.4f',plotaxlabels{xyzinv(2)},pos(xyzinv(2))),...
			   sprintf('%s: %.4f',plotaxlabels{xyzinv(3)},pos(xyzinv(3))),...
			   ['Index: ',num2str(didx')]};
	else
		didx = get(event_obj, 'DataIndex');
		txt = {sprintf('%s: %.4f',plotaxlabels{1},CT(didx,1)),...
			   sprintf('%s: %.4f',plotaxlabels{2},CT(didx,2)),...
			   sprintf('%s: %.4f',plotaxlabels{3},CT(didx,3)),...
			   ['Index: ',num2str(didx')]};
	end
end

end % END MAIN SCOPE


function out = fromrgb(f,spc)
	switch spc
		case 'rgb'
			out = f;
		case 'lrgb'
			out = rgb2linear(f);
		case 'hsi'
			out = rgb2hsi(f);
		case 'hsl'
			out = rgb2hsl(f);
			H = out(:,:,1);
			H(isnan(H)) = 0;
			out(:,:,1) = H;
		case 'hsv'
			out = rgb2hsv(f);
			out(:,:,1) = out(:,:,1)*360;
		case {'hsy'}
			out = rgb2hsy(f);
		case {'hsyp'}
			out = rgb2hsy(f,'pastel');
		case {'huslab'}
			out = rgb2husl(f,'lab');
		case {'husluv'}
			out = rgb2husl(f,'luv');
		case {'huslpab'}
			out = rgb2husl(f,'labp');
		case {'huslpuv'}
			out = rgb2husl(f,'luvp');
		case 'lchbr'
			out = rgb2lch(f,'ypbpr');
		case 'lchab'
			out = rgb2lch(f,'lab');
		case 'lchuv'
			out = rgb2lch(f,'luv');
		case 'lchsr'
			out = rgb2lch(f,'srlab');
		case 'srlab'
			out = rgb2lch(f,'srlab');
			out = lch2lab(out);
		case 'lchok'
			out = rgb2lch(f,'oklab');
		case 'oklab'
			out = rgb2lch(f,'oklab');
			out = lch2lab(out);
		case 'lab'
			out = rgb2lch(f,'lab');
			out = lch2lab(out);
		case 'luv'
			out = rgb2lch(f,'luv');
			out = lch2lab(out);
		case {'ypbpr','ypp'}
			out = imappmat(f,gettfm(spc));
		otherwise
			% this shouldn't ever happen
			error('CTPATH: unknown colorspace %s',spc)
	end
end




