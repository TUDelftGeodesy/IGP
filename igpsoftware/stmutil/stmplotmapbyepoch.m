function h=stmplotmapbyepoch(stmfile,varargin)
%stmplotmapbyepoch   Plot map of the space time matrix points by epoch.
%  STMPLOTMAPBYEPOCH(STM) plot map of points in the space time matrix dataset
%  STM, counting the number of epochs, and if feasible a subplot of the 
%  network for each campaign (epoch).
%
%  H=STMPLOTMAP(STM,'option',value,...) returns the plot handle H and
%  allows to the following option/value pairs
%
%    'background'     matfile or lat/lon with coastlines, lakes, borders 
%                     (default 'igp-background.mat') 
%    'defoborder'     matfile or lat/lon with deformation border (default none)
%    'roi'            matfile or lat/lon with roi (default none)
%    'saveplt'        Directory for pdf plots (default [] is none)
%    'maxepochs'      Maximum number of subplots with individual epochs (default is 30)
%
%  Example:
%
%      stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_ameland.mat')
%
%  See also STMDISP, STMPLOTMAP and STMPLOTPROJECTMAP.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023-2024

% Created:  12 April 2023 by Hans van der Marel
% Modified: 23 Oct 2023 by Hans van der Marel
%            - branched off from more complicated stmplotprojectmap and
%              stmplotvelmap
%           17 June 2024 by Hans van der Marel
%            - added saveplt option
%            - modified default plot name
%            - added option for maximum number of epochs
%            - support for cases with multiple observation types

%% Check arguments and process options

if nargin < 1 
    error('stmplotprojectmap expects at least one argument')
end 

opt.background='igp-background.mat';   % matfile or lat/lon with map background (coastlines, lakes, borders) 
opt.defoborder='igp-unstablearea.mat'; % matfile or lat/lon with deformation border
opt.roi='';                            % matfile or lat/lon with region of interest
opt.datatips=true;                     % use extended datatips
opt.saveplt=[];                        % directory for pdf plots (default [] is none)
opt.maxepochs=30;                      % Maximum number of subplots with individual epochs (default is 30)

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end


%% Load space time matrix

%st=stmread(stmfile,'NODATA');
st=stmread(stmfile,'NOSTOCHDATA');

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

if ~isempty(opt.saveplt)
   if ~exist(opt.saveplt,'dir')
       mkdir(opt.saveplt)
   end
   pdffile=fullfile(opt.saveplt,[stmfile '_network.pdf']);
   fprintf('The plots will be saved to %s\n',pdffile);
end

vis=~isnan(st.obsData);
if ndims(vis) == 3
   vis=any(vis,3);
end

pntCount=sum(vis,2);
epochMaskNonEmpty=~all(vis,1);
pntMaskAll=all(vis(:,epochMaskNonEmpty),2);
if any(pntMaskAll)
   fprintf('PointIds of points present in every campaign:\n')
   disp(pntName(pntMaskAll))
else
   pntMaskAll= pntCount == max(pntCount);
   fprintf('PointIds of points present in %d campaigns (out of %d):\n',max(pntCount),st.numEpochs);
   disp(pntName(pntMaskAll))
end

h=figure;

if showpointnames
    ht=text(pntCrd(pntMask,2),pntCrd(pntMask,1),pntName(pntMask),'FontSize',7,'VerticalAlignment','top');
    hold on; 
end
%hp=plot(pntCrd(pntMask,2),pntCrd(pntMask,1),'o','Markersize',6,'MarkerFaceColor','green');
hp=scatter(pntCrd(pntMask,2),pntCrd(pntMask,1),10+50*pntCount/st.numEpochs,pntCount,"filled");
hc=colorbar();
ylabel(hc,'#epochs');
colormap(turbo)
hold on;
plot(pntCrd(pntMaskAll,2),pntCrd(pntMaskAll,1),'s','Markersize',12);
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

title(stmfile ,'interpreter','none')

if ~isempty(opt.saveplt)
   exportgraphics(h, pdffile, 'Append', false);
end

%% Plot map for each epoch

if st.numEpochs <= opt.maxepochs

    [~,kidx]=sort(st.epochDyear);
    
    for kk=1:st.numEpochs
    
        k=kidx(kk);
        pntMask=~isnan(st.obsData(:,k));
    
        if ~any(pntMask), continue, end
    
        h=figure;
        
        if showpointnames
            ht=text(pntCrd(pntMask,2),pntCrd(pntMask,1),pntName(pntMask),'FontSize',7,'VerticalAlignment','top');
            hold on; 
        end
        hp=plot(pntCrd(pntMask,2),pntCrd(pntMask,1),'o','Markersize',6,'MarkerFaceColor','green');
        hold on;
        plot(pntCrd(pntMaskAll,2),pntCrd(pntMaskAll,1),'s','Markersize',12);
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
    
        if isfield(st.epochAttrib,'prjName')
           stmtitle=[stmfile ' ' st.epochAttrib.prjName{k} ' (' num2str(st.epochDyear(k)) ')' ];
        else
           stmtitle=[stmfile ' ' num2str(st.epochDyear(k)) ];
        end
        title(stmtitle ,'interpreter','none')
    
        if ~isempty(opt.saveplt)
           exportgraphics(h, pdffile, 'Append', true);
        end
    
    end

end

end
