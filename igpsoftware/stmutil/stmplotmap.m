function h=stmplotmap(stmfile,varargin)
%stmplotmap   Plot map of the space time matrix points
%  STMPLOTMAP(STM) plot map of points in the space time matrix dataset
%  STM.
%
%  H=STMPLOTMAP(STM,'option',value,...) returns the plot handle H and
%  allows to the following option/value pairs
%
%    'background'     matfile or lat/lon with coastlines, lakes, borders 
%                     (default 'igp-background.mat') 
%    'defoborder'     matfile or lat/lon with deformation border (default none)
%    'roi'            matfile or lat/lon with roi (default none)
%
%  Example:
%
%      stmplotmap('../3_integrate/waddenzee_insarReducedWithDiag_LevDiag.mat');
%
%  See also STMDISP and STMPLOTPROJECTMAP.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:  12 April 2023 by Hans van der Marel
% Modified: 23 Oct 2023 by Hans van der Marel
%            - branched off from more complicated stmplotprojectmap and
%              stmplotvelmap

%% Check arguments and process options

if nargin < 1 
    error('stmplotprojectmap expects at least one argument')
end 

opt.background='igp-background.mat';   % matfile or lat/lon with map background (coastlines, lakes, borders) 
opt.defoborder='igp-unstablearea.mat'; % matfile or lat/lon with deformation border
opt.roi='';                            % matfile or lat/lon with region of interest
opt.datatips=true;                     % use extended datatips

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end


%% Load space time matrix

st=stmread(stmfile,'NODATA');

stmfile = st.datasetId;
pntCrd = st.pntCrd;
pntName = st.pntName;

pntMask=true(size(pntName));

%% Prepare map

% Check if datatips are supported

if opt.datatips && exist('dataTipTextRow','file') ~= 2
    fprintf('Warning: Datatips are not supported (by this Matlab version)')
    opt.datatips=false;
end

% Map labels and aspect ratio

if isfield(opt,'doplots_mapcrd') && opt.doplots_mapcrd
   xylabel = { 'East [km]', 'North [km]' };
   xyaspect = [1 1 1];
else
   xylabel = { 'Longitude [deg]', 'Latitude [deg]' };
   xyaspect = [ 1/cosd(mean(pntCrd(:,1))) 1 1];
end

% Get background, deformation border and roi

if isnumeric(opt.background)
    igpbackground=opt.background;
elseif exist(opt.background,'file')
    igpbackground=load(opt.background);
    if isstruct(igpbackground)
        igpbackground = [ igpbackground.Lat(:) igpbackground.Lon(:) ];
    end
else
    igpbackground=[];
end

if isnumeric(opt.defoborder)
    defoborder=opt.defoborder;
elseif exist(opt.defoborder,'file')
    defoborder=load(opt.defoborder);
    if isstruct(defoborder)
        defoborder = [ defoborder.Lat(:) defoborder.Lon(:) ];
    end
else
    defoborder=[];
end

if isnumeric(opt.roi)
    roi=opt.defoborder;
elseif exist(opt.roi,'file')
    roi=load(opt.defoborder);
else
    roi=[];
end

% Prepare plotting of names

showpointnames = numel(pntName) <= 50;

%% Plot map

h=figure;

if showpointnames
    ht=text(pntCrd(pntMask,2),pntCrd(pntMask,1),pntName(pntMask),'FontSize',7,'VerticalAlignment','top');
    hold on; 
end
hp=plot(pntCrd(pntMask,2),pntCrd(pntMask,1),'.','Markersize',6);
hold on;
xl=xlim();yl=ylim(); 
if opt.datatips
   datatiprows= [ ...
            dataTipTextRow('pntName',pntName(pntMask)), ...
            dataTipTextRow('Lon [deg]',pntCrd(pntMask,2)) , ...
            dataTipTextRow('Lat [deg]',pntCrd(pntMask,1)) ];
   hp.DataTipTemplate.DataTipRows = datatiprows;
end   
if ~isempty(igpbackground)
    plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
end
if ~isempty(defoborder)
    plot(defoborder(:,2),defoborder(:,1),':m')
end
if ~isempty(roi)
    plot(roi(:,2),roi(:,1),':m')
end
xlim(xl);ylim(yl);

daspect(xyaspect);
xlabel(xylabel{1});
ylabel(xylabel{2});

title(stmfile,'interpreter','none')


end
