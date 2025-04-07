function h=stmplotvel(stm,varargin)
%STMPLOTVEL   Estimate and plot velocities from space time matrix datasets
%  STMPLOTVEL(STM) estimate velocities the space time matrix STM using 
%  STMVELOCITY and plot. STM can be a single space time matrix, or else the 
%  name of a space time matrix dataset, or a character cell array with 
%  the names of multiple space time matrix datasets.
%
%  H=STMPLOTVEL(STM,'option',value,...) returns the plot handle(s) H 
%  and allows to the following options
%
%    'saveplt'            Directory for pdf plots (default [] is none)
%
%    'background'         Map background (Default 'igp-background.mat')
%    'unstablearea'       Polygon defining unstable area (Default 'igp-unstablearea.mat')
%
%    'refsystem'          Choice of reference system (default 'min-velocities'):
%                         - 'min-velocities'
%                              The mean velocity of the points outside 
%                              'unstablearea' is set to zero.
%                         - 'inherit'        
%                              Same as 'relative'=false 
%                         - 'minimum-norm'   
%                              Minimum norm solution (use only when offsets are close to zero)
%                         - 'min-refseries'  
%                              Minimize the reference series corrections.
%
%    'ROI'                Region of interest, as [latmin lonmin ; latmax lonmax] 
%                         bounding box, or lat/lon polygon (Default [], is all)
%    'ignoreStochModel'   If true, ignore stochastic model from stm, use 
%                         unit matrix instead (default false) 
%
%    'stmplotvelmap'      Cell array with options for stmplotvelmap
%
%  For each space time matrix dataset three plots are created:
%  - plot with the estimated velocity parameters and residuals
%  - plot of the covariance matrices
%  - map of the area with estimated velocities
%  
%  See also STMVELOCITY, STMPLOTVELPAR, STMPLOTVELCOV and STMPLOTVELMAP.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:   5 Oct 2023 by Hans van der Marel
% Modified: 23 Oct 2023 by Hans van der Marel
%            - refactored stmvelocity_test.m into stmplotvel.m function

% Check input arguments and process options

if nargin < 1
   error('This function expects at least one input argument.')
end

opt.saveplt=[];                          % directory for pdf plots (default [] is none)
opt.refsystem='min-velocities';          % refsystem option for stmvelocity
opt.ignoreStochModel=false;              % If true, ignore stochasticmodel from stm, use unit matrix instead 
opt.ROI=[];                              % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default none)
opt.background='igp-background.mat';     % map background
opt.unstablearea='igp-unstablearea.mat'; % polygon defining unstable area
opt.stmplotvelmap={};                    % cell array with options for stmplotvelmap

opt.append=false;                        % undocumented / hidden / use with care

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end

% Check what we have to do, if necessary, call recursively

if iscell(stm)
    % Cell array with space time matrix datasets, call this function recursively
    stmfiles=stm;
    h={};
    append=opt.append;
    for k=1:numel(stmfiles)
        stmfile=stmfiles{k};
        copt=opt2cell(opt);
        hs=stmplotvel(stmfile,copt{:});
        % Save plots to pdf (optional)
        if ~isempty(opt.saveplt)
           exportgraphics(hs(1), fullfile(opt.saveplt,'velpar.pdf'), 'ContentType', 'vector','Append',append);
           exportgraphics(hs(2), fullfile(opt.saveplt,'velcov.pdf'), 'ContentType', 'vector','Append',append);
           for i=3:numel(hs)
              exportgraphics(hs(i), fullfile(opt.saveplt,'velmap.pdf'), 'ContentType', 'vector','Append',append);
           end
           append=true;
        end
    end
    h=[h;hs];
    return
elseif ischar(stm) || isstring(stm)
   % Load space time matrix dataset
   stmfile=stm;
   stm = stmread(stmfile);
end

% Get polygon with unstable area

unstable_area = load(opt.unstablearea);
if isstruct(unstable_area)
   unstable_area = [ unstable_area.Lat(:) unstable_area.Lon(:) ];
end

% Estimate velocities

if strcmp(opt.refsystem,'min-velocities')
    % minimize velocities of selected points (pntRefMask)
    pntRefMask = ~inpolygon(stm.pntCrd(:,2),stm.pntCrd(:,1),unstable_area(:,2),unstable_area(:,1));
    fprintf('Number of stable points: %d\n',sum(pntRefMask))
    [vel,t0,offset,ref,omt,emat,dmat,qx,qe,qy]=stmvelocity(stm, ...
         'refsystem',pntRefMask,'ROI',opt.ROI,'ignoreStochModel',opt.ignoreStochModel);
else
    % default refsystem option (min-refseries) - stochmodel from stm (THIS IS THE DEFAULT)
    [vel,t0,offset,ref,omt,emat,dmat,qx,qe,qy]=stmvelocity(stm, ...
        'refsystem',opt.refsystem,'ROI',opt.ROI,'ignoreStochModel',opt.ignoreStochModel);
end

% Plot

h1=stmplotvelpar(stm,vel,t0,offset,ref,omt,emat,qx);
h2=stmplotvelcov(stm,vel,qx,qe,qy);
h3=stmplotvelmap(stm,vel,'background',opt.background,'defoborder',opt.unstablearea, ...
     opt.stmplotvelmap{:});

% Save plots to pdf (optional)

if ~isempty(opt.saveplt) 
   if ~exist(opt.saveplt,'dir')
       mkdir(opt.saveplt)
   end
   [~,stmbase]=fileparts(stmfile);
   pltfile=fullfile(opt.saveplt,[ stmbase '_velocity.pdf']);
   exportgraphics(h1, pltfile, 'ContentType', 'vector');
   exportgraphics(h2, pltfile, 'ContentType', 'vector','Append',true);
   for i=1:numel(h3)
      exportgraphics(h3(i), pltfile, 'ContentType', 'vector','Append',true);
   end 
end

% return plot handles

h=[h1 h2 h3];

end

function C=opt2cell(opt)
C=[fieldnames(opt).'; struct2cell(opt).'];
C=C(:).';
end

%%

function h=stmplotvelpar(st,vel,t0,offset,ref,omt,emat,qx)
%STMPLOTVELPAR Plot velocities, offsets and refseries for space time matrices.
%  H=STMPLOTVELPAR(ST,VEL,T0,OFFSET,REF,OMT,EMAT,QX) plots the velocities
%  VEL, offsets OFFSET and reference series REFSERIES estimated by STMVELOCITY
%  for the space time matrix ST. Other inpunts from STMVELOCITY are the 
%  reference epoch T0, overall model test values OMT, residual matrix EMAT
%  and covariance matrix QX of the estimated parameters.
% The function returns the plot handle H.
%
%  See also STMVELOCITY, STMPLOTVELOCITY, STMPLOTVELCOV and STMPLOTVELMAP.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:   2 Oct 2023 by Hans van der Marel
% Modified: 

stmfile=st.datasetId;

numObsTypes=size(vel,2);

% Get standard deviation values

[svel,soffset,sref]=qx2sigma(vel,offset,ref,qx);

% Check if extended datatips are supported

opt.datatips=true;
if opt.datatips && exist('dataTipTextRow','file') ~= 2
    fprintf('Warning: Datatips are not supported (by this Matlab version)')
    opt.datatips=false;
end

% Plot using tiled layout 

h=figure('units','normalized','outerposition',[0.05 0.05 .95 0.05+0.3.*numObsTypes]);

%tcl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
tcl = tiledlayout(numObsTypes,4,'TileSpacing','compact','Padding','compact');

for l=1:numObsTypes

   nexttile
   %plot(vel(:,l),'ob','MarkerFaceColor','b','MarkerSize',3); 
   h1=errorbar(vel(:,l),svel(:,l),'ob','MarkerFaceColor','b','MarkerSize',3); 
   ylabel(['Velocity ' st.obsTypes{l} ' [mm/y]'])
   text(0.95,0.95, ['omt=' num2str(omt(l),3) ],'Units','normalized','VerticalAlignment','top','HorizontalAlignment','right')
   if opt.datatips
        dtRows= [  dataTipTextRow('pntName',st.pntName) , ...
                   dataTipTextRow(['Vel ' st.obsTypes{l} ' [mm/y]' ],vel(:,l)) , ...
                   dataTipTextRow(['Std ' st.obsTypes{l} ' [mm/y]' ],svel(:,l)) , ...
                   dataTipTextRow('Lon [deg]',st.pntCrd(:,1)) , ...
                   dataTipTextRow('Lat [deg]',st.pntCrd(:,2)) ];
        h1.DataTipTemplate.DataTipRows = dtRows;
   end

   nexttile
   %plot(offset(:,l),'ob','MarkerFaceColor','b','MarkerSize',3)
   h2=errorbar(offset(:,l),soffset(:,l),'ob','MarkerFaceColor','b','MarkerSize',3);
   ylabel(['Offset ' st.obsTypes{l} ' [mm]'])
   if opt.datatips
        dtRows= [  dataTipTextRow('pntName',st.pntName) , ...
                   dataTipTextRow(['Offset ' st.obsTypes{l} ' [mm]' ],offset(:,l)) , ...
                   dataTipTextRow(['Std ' st.obsTypes{l} ' [mm]' ],soffset(:,l)) , ...
                   dataTipTextRow('Lon [deg]',st.pntCrd(:,1)) , ...
                   dataTipTextRow('Lat [deg]',st.pntCrd(:,2)) ];
        h2.DataTipTemplate.DataTipRows = dtRows;
   end

   nexttile
   %plot(st.epochDyear,ref(:,l),'ob','MarkerFaceColor','b','MarkerSize',3)
   h3=errorbar(st.epochDyear,ref(:,l),sref(:,l),'ob','MarkerFaceColor','b','MarkerSize',3);
   ylabel(['Refseries ' st.obsTypes{l} ' [mm]'])
   if opt.datatips
        dtRows= [  dataTipTextRow('dYear',st.epochDyear) , ...
                   dataTipTextRow(['Refseries ' st.obsTypes{l} ' [mm]' ],ref(:,l)) , ...
                   dataTipTextRow(['Std ' st.obsTypes{l} ' [mm]' ],sref(:,l)) ];
        h3.DataTipTemplate.DataTipRows = dtRows;
   end

   nexttile
   imagesc(emat(:,:,l),'AlphaData', 1-isnan(emat(:,:,l)))
   hc=colorbar();
   ylabel(hc,['Residual ' st.obsTypes{l} ' [mm]'])

end

title(tcl,[ stmfile ' (t0=' sprintf('%.2f',t0) ')'] ,'interpreter','none')

end

function [svel,soffset,sref]=qx2sigma(vel,offset,ref,qx)
%QX2SIGMA   Supporting function for STMPLOTVELPAR.

numObsTypes=size(vel,2);

svel=nan(size(vel));
soffset=nan(size(offset));
sref=nan(size(ref));

for l=1:numObsTypes
   pntMask=~isnan(vel(:,l));
   epochMask=~isnan(ref(:,l));
   np=sum(pntMask);
   ne=sum(epochMask);
   sigmas=sqrt(diag(qx{l}));
   svel(pntMask,l)=sigmas(1:np);
   soffset(pntMask,l)=sigmas(np+1:2*np);
   if length(sigmas) > 2*np
       sref(epochMask,l)=sigmas(2*np+1:end);
   else
       sref(epochMask,l)=zeros(ne,1);
   end

end

end


function h=stmplotvelcov(st,vel,qx,qe,qy)
%STMPLOTVELCOV  Plot covariance matrix of velocities and residuals.
%  H=STMPLOTVELCOV(ST,VEL,QX,QE,QY) plots the covariance matrix QX
%  for the estimated parameters (velocities, offsets, refseries) by 
%  STMVELOCITY(ST,...), the velocity part of QX, the covariance matrix of
%  the estimated residuals QE, and covariance matrix Qy of the observations.
%  The input VEL is needed to determine the size of the velocity part of
%  QC, but is otherwise not needed. The function returns the plot handle H.
%
%  See also STMVELOCITY, STMPLOTVELOCITY, STMPLOTVELPAR and STMPLOTVELMAP.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:   3 Oct 2023 by Hans van der Marel
% Modified: 26 Mar 2024 by Freek van Leijen
%           - inserted check on Matlab release since before 2022 the
%             clim function was called caxis.

numObsTypes=size(vel,2);
stmfile=st.datasetId;

h=figure('units','normalized','outerposition',[0.05 0.05 .95 0.05+0.3.*numObsTypes]);

%t=tiledlayout('flow');
t = tiledlayout(numObsTypes,5,'TileSpacing','compact','Padding','compact');

for i=1:size(qe,1)
   qxhat=qx{i};
   np=size(vel(~isnan(vel(:,i))),1);
   d=diag(qx{i});
   
   nexttile
   imagesc(qxhat)
   if max(d) > 2
    matlab_version = version('-release');
    if str2num(matlab_version(1:4))<2022
      cl=caxis;
      cl(2)=quantile(d,.95);
      caxis(cl);
    else
      cl=clim;
      cl(2)=quantile(d,.95);
      clim(cl);
    end
    colorbar
      title(['Qxhat (maxVar=' num2str(max(d)) ')'])
   else
      colorbar
      title('Qxhat')
   end

   nexttile
   imagesc(qxhat(1:np,1:np))
   colorbar
   title('Qxhat (velocities)')

   nexttile
   imagesc(qy{i})
   colorbar
   title('Qy')
   
   nexttile
   imagesc(qe{i})
   colorbar
   title('Qehat')
   
   nexttile
   sig=diag(qy{i});
   sig(sig <= 0)=0;
   sig=sqrt(sig);
   plot(sig,'b.')
   hold on
   sig=diag(qe{i});
   sig(sig <= 0)=0;
   sig=sqrt(sig);
   plot(sig,'r.')
   legend('\sigma_y','\sigma_e')
   title('sqrt(diag(Q))')
end
title(t, stmfile ,'interpreter','none')

end

function h=stmplotvelmap(st,vel,varargin)
%STMPLOTVELMAP   Map with estimated velocities for space time matrix.
%  STMPLOTVELMAP(ST,VEL) plots the velocities VEL, estimated by STMVELOCITY,
%  from the space time matrix ST. A separate plot is created for each
%  observation type, except for optional horizontal components which are
%  plotted as vectors. 
%
%  H=STMPLOTVELMAP(ST,VEL,'option',value,...) returns the plot handle(s) H 
%  and allows to the following options
%
%   'background'     matfile or lat/lon with coastlines, lakes, borders 
%                    (default 'igp-background.mat') 
%   'defoborder'     matfile or lat/lon with deformation border (default 
%                    'igp-unstablearea.mat')
%   'width'          width of plot in normalized coordinates (default 0.8),
%                    the height is computed respecting the aspect ratio
%   'zmax'           colorrange in up dimension (default 6 [mm/y])
%   'quiverStretch'  stretch factor for quiver part of plots (default 0.2)
%
%  See also STMVELOCITY, STMPLOTVELOCITY, STMPLOTVELPAR and STMPLOTVELCOV.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:   5 Oct 2023 by Hans van der Marel
% Modified: 

% Check input arguments and process options

if nargin < 2
   error('This function expects at least two input arguments.')
end

opt.background='igp-background.mat';   % matfile or lat/lon with map background (coastlines, lakes, borders) 
opt.defoborder='igp-unstablearea.mat'; % matfile or lat/lon with deformation border
opt.width=0.8;                         % width of plot in normalized coordinates
opt.zmax=6;                            % colorrange in z-dimension for colormaps
opt.quiverStretch=0.2;                 % stretch factor for quiver part of plots
opt.datatips=true;                     % use extended datatips

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end

% Check if extended datatips are supported

if opt.datatips && exist('dataTipTextRow','file') ~= 2
    fprintf('Warning: Datatips are not supported (by this Matlab version)')
    opt.datatips=false;
end

% Get background and deformation border

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

% get point mask (from velocities) and set plotsize

pntMask=any(~isnan(vel),2);

width=opt.width;
height=   ( max(st.pntCrd(pntMask,1)) - min(st.pntCrd(pntMask,1)) ) / ...
        ( ( max(st.pntCrd(pntMask,2)) - min(st.pntCrd(pntMask,2)) ) * cosd(mean(st.pntCrd(pntMask,1))) ) * width * 16/9;
if height > width
    tmp=height/width;
    height=height/tmp;
    width=width/tmp;
end

% Prepare observation types and plotting of names

obsTypes=st.obsTypes;
stmfile=st.datasetId;

isnorth=ismember(st.obsTypes,{'North'});
iseast=ismember(st.obsTypes,{'East'});

showpointnames = sum(any(~isnan(vel),2)) <= 50;

% Plot map

k=0;
for l=1:numel(obsTypes)

    if isnorth(l) || iseast(l), continue; end
    
    k=k+1;
    h(k)=figure('units','normalized','innerPosition',[0.1 0.1 width height]);

    hs=scatter(st.pntCrd(pntMask,2),st.pntCrd(pntMask,1),[],vel(pntMask,l),'filled');
    scolormap(opt.zmax);
    hc=colorbar;
    ylabel(hc,[ st.obsTypes{l} ' Velocity [mm/y]'])
    if showpointnames
        text(st.pntCrd(pntMask,2),st.pntCrd(pntMask,1),st.pntName(pntMask),'FontSize',7,'VerticalAlignment','top')
    end
    hold on; 
    xl=xlim();yl=ylim(); 
    if ~isempty(igpbackground)
        plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
    end
    if ~isempty(defoborder)
        plot(defoborder(:,2),defoborder(:,1),':m')
    end
    xlim(xl);ylim(yl);
    daspect([1/cosd(mean(yl)) 1 1])
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    if opt.datatips
        dtRows= [  dataTipTextRow('pntName',st.pntName(pntMask)) , ...
                   dataTipTextRow('Lon [deg]',st.pntCrd(pntMask,1)) , ...
                   dataTipTextRow('Lat [deg]',st.pntCrd(pntMask,2)) , ...
                   dataTipTextRow(['Vel ' st.obsTypes{l} ' [mm/y]' ],vel(pntMask,l)) ];
        if any(iseast)
            dtRows(end+1) = dataTipTextRow('Vel East [mm/y]',vel(pntMask,iseast));
        end
        if any(isnorth)
            dtRows(end+1) = dataTipTextRow('Vel North [mm/y]',vel(pntMask,isnorth));
        end
        hs.DataTipTemplate.DataTipRows = dtRows;
    end
    if any(iseast) && any(isnorth)
       hq1=quiver(st.pntCrd(pntMask,2),st.pntCrd(pntMask,1),vel(pntMask,iseast),vel(pntMask,isnorth)*cosd(mean(yl)),opt.quiverStretch);
       if opt.datatips
          hq1.DataTipTemplate.DataTipRows = dtRows;
       end
       hasonlyeast = pntMask & ~isnan(vel(:,iseast)) & isnan(vel(:,isnorth));
       if any(hasonlyeast) 
          hq2=quiver(st.pntCrd(hasonlyeast,2),st.pntCrd(hasonlyeast,1),vel(hasonlyeast,iseast),zeros(size(vel(hasonlyeast,isnorth)))*0,opt.quiverStretch);
          if opt.datatips
             dtRows= [  dataTipTextRow('pntName',st.pntName(hasonlyeast)) , ...
                        dataTipTextRow('Lon [deg]',st.pntCrd(hasonlyeast,1)) , ...
                        dataTipTextRow('Lat [deg]',st.pntCrd(hasonlyeast,2)) , ...
                        dataTipTextRow(['Vel ' st.obsTypes{l} ' [mm/y]' ],vel(hasonlyeast,l)) , ...
                        dataTipTextRow('Vel East [mm/y]',vel(hasonlyeast,iseast)) ];
             hq2.DataTipTemplate.DataTipRows = dtRows;
          end
       end
    end
    title([stmfile ' - ' st.obsTypes{l} ],'interpreter','none')

end

end

function scolormap(zmax,zint)
%SCOLORMAP  Divergent colormap for ground motion.
%   SCOLORMAP(ZMAX) sets colormap property to divergent RdYlBu colormap from
%   Colorbrewer and set color limit for the current axis to [-zmax zmax].
%
%   The divergent RdYlBu colormap is colorblind safe and suitable for
%   showing ground motion. 
%
%   The colormap is discrete. The number of classes N is computed from 
%   ZMAX such that (i) the number of classes is about 11 for each half of the
%   colormap and (ii) the class interval ZINT is integer or an integer fraction.  
%
%   SCOLORMAP(ZMAX,ZINT) overrides the default class interval ZINT.
%
%   See also colormap and clim (caxis in older matlab versions).
%
%   (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:  21 Oct 2023 by Hans van der Marel
% Modified: 26 Mar 2024 by Freek van Leijen
%           - inserted check on Matlab release since before 2022 the
%             clim function was called caxis.

NCLASS2=11;
if nargin < 2
    zint = 1/round(NCLASS2/zmax);
end
N=2*zmax/zint;

colormap(rdylbu(N))
matlab_version = version('-release');
if str2num(matlab_version(1:4))<2022
  caxis([-zmax zmax ])
else
  clim([-zmax zmax ])  
end

end

function map=rdylbu(N)
%RDYLBU  Divergent RdYlBu colormap from Colorbrewer.
%   MAP=RDYLBU(N) returns the three-column matrix MAP of RGB triplets, defining 
%   a Matlab colormap, with the divergent RdYlBu colormap from Colorbrewer.
%   N is the number of classes for the map (default is 11).
% 
%   The divergent RdYlBu colormap is colorblind safe.
%
%   (c) 2023, Delft Univerversity of Technology, Hans van der Marel

% Created:  21 Oct 2023 by Hans van der Marel
% Modified: 

% Define RdYlBu colormap ( https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=11 )

CB_div_RdYlBu_11 = [ ...
165   0  38  ; ...
215  48  39  ; ...
244 109  67  ; ...
253 174  97  ; ...
254 224 144  ; ...
255 255 191  ; ...
224 243 248  ; ...
171 217 233  ; ...
116 173 209  ; ...
 69 117 180  ; ...
 49  54 149  ] / 255;

num = size(CB_div_RdYlBu_11,1);

% Check input arguments

if nargin < 1
    N=num;
end

% Interpolate Colorbrewer map to desired number of classes

map = CB_div_RdYlBu_11;

if N ~= num
    ido = linspace(1,num,N);   
    mid = ceil(num/2);
    ida = 1:mid;
    idz = mid:11;
    map = [...
	    interp1(ida,map(ida,:),ido(ido<=mid),'pchip');...
	    interp1(idz,map(idz,:),ido(ido>mid),'pchip')];
end

% Reverse the map

map = map(end:-1:1,:);

end
