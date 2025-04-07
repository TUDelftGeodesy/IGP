function stmexport(inputFilename,options)
%stmexport   Export of Space-Time Matrix in desired format(s).
%  STMEXPORT(INPUTFILENAME,OPTIONS) Function to export an STM in
%  an output format. The STM is specified by the INPUTFILENAME.
%
%  STMEXPORT(INPUTFILENAME,'option',value,...) allows to specify options
%  for the export
%
%  'format'        array with output format(s) (can be multiple, but
%                  at least one). Options: 'csv' (default),'png',
%                  'geotiff','shp'.
%
%     Examples:
%       options.format = { 'csv' 'shp'};
%       options.format = { 'png' 'geotiff' };
%       options.format = { 'png' };
%       options.format = { 'geotiff' }:
%       options.format = { 'geotiff' 'csv' 'png' 'shp'};
%
%   'unit'         unit of output, 'm' or 'mm' (default)
%
%   'dpi'          string with resolution (in DPI) for png plots.
%                  Default is '300';
%   'cmap'         color map for plots. Default is 'defo'.
%   'visible'      visibility of plots, 'on' (default) or 'off'
%
%

%  (c) Freek van Leijen, Delft University of Technology, 2023.

% Created:   16 Nov 2023 by Freek van Leijen
% Modified: 
%         - 31 March 2024 by Freek van Leijen
%           Added runId input parameter for unique output directory.
%           Use of datasetId instead of techniqueId for output filename.
%
%  TODO: - finalize geotiff, shp output
%

% Check input arguments and process options

if nargin < 2
   error('This function expects at least two input arguments.')
end

progname = 'stmexport';

% Default options

opt.format = 'csv';                 % Output format specficication (can be multiple). Default is csv.
opt.dpi =  '300';                   % Plot output resolution (in DPI).
opt.cmap = 'defo';                  % Colormap.
opt.visible = 'on';                 % Visibility of plots.
opt.unit = 'mm';                    % Unit ([m] or [mm] (default)).
opt.background = '';                % Plot background

%% Catch errors and start timing
try

[~,inputFileroot]=fileparts(inputFilename);
runId = [ inputFileroot '_export_' datestr(now,30) ];
diary([ runId '.log' ])
  
fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Check the options and if necessary overwrite the default values
outputFilename = [];
flag = 'overwrite';
[inputFilename,outputFilename,opt]= ...
    stmcheckarguments({inputFilename},outputFilename,opt,options,flag);
inputFilename = inputFilename{1}; % Convert to character string

%% Read space time matrix dataset
stmIn = stmread(inputFilename);


fprintf('\n');
fprintf('Start export ...\n');

%% Create output
for v = 1:stmIn.techniqueAttrib.numLayers
  for w = 1:numel(opt.format)

    switch opt.format{w}
      case 'csv'
        exportcsv(stmIn,runId,opt,v);
      case 'png'
         exportpng(stmIn,runId,opt,v);
      case 'geotiff'
        error('Geotiff: not fully implemented yet.')
        exportgeotiff(stmIn,runId,opt,v);
        %exportgeotiff(stmIn,runId,opt,v,cmapdefo);
      case 'shp'
        error('Shp: not fully implemented yet.')
        exportshape(stmIn,runId,opt,v);
      otherwise
        error('You specified a unsupported output format.');
    end
  end
end

%% Finish the function

fprintf('\n');
fprintf('%s finished at %s  (elapsed time %.2f s)\n',progname,datestr(now),toc)
diary OFF

catch ME

   getReport(ME)
   fprintf('%s ABORTED with an error at %s (elapsed time %.2f s)\n',progname,datestr(now),toc)
   diary OFF
   rethrow(ME)
    
end

end


%% Sub-functions

function exportcsv(stm,runId,opt,layer)
%exportcsv   Export of Space-Time Matrix in csv format.
%  EXPORTCSV(STM,RUNID,OPT,LAYER) Function to export an STM in
%  csv format. The RUNID provides the ID for the output, used to
%  create an output directory, containing the original filename of
%  the STM and the timestamp of the export. The LAYER indicates
%  the index of the layer exported. OPT contains the options
%  as specified when calling the parent function stmexport. Here
%  used:
%
%  'unit'          unit of output, 'm' or 'mm' (default)
%
%  (c) Freek van Leijen, Delft University of Technology, 2023. 

% Created:  16 Nov 2023 by Freek van Leijen
% Modified: 
%         - 31 March 2024 by Freek van Leijen
%           Added runId input parameter for unique output directory.
%           Use of datasetId instead of techniqueId for output filename.
%
% TODO:
%         - Make the function compatible with all STM files.
%

%% Create output directory and file
csvDir = [runId '_csv'];
if ~exist(['./' csvDir],'dir')
  mkdir(csvDir);
end

csv_fid = fopen([csvDir '/' stm.datasetId '_' stm.obsTypes{layer} '_' ...
                   stm.techniqueAttrib.space{1} '_' stm.techniqueAttrib.crs '.csv'],'w');

%% Create header
if isempty(strmatch(upper(stm.techniqueAttrib.crs),{'WGS84','EPSG4326'}))
  attrib = { 'ID' 'Lon [deg]' 'Lat [deg]' 'X [m]' 'Y [m]'};
else
  attrib = { 'ID' 'Lon [deg]' 'Lat [deg]' };
end
numAttrib = length(attrib);

for v = 1:numel(stm.epochDyear)
  attrib{numAttrib+v} = [num2str(stm.epochDyear(v)) ' [' opt.unit ']'];
end
numAttrib = length(attrib);

for v = 1:numel(stm.epochDyear)
  attrib{numAttrib+v} = ['std ' num2str(stm.epochDyear(v)) ' [' opt.unit ']'];
end
numAttrib = length(attrib);


for v = 1:numAttrib
  if v==numAttrib
    fprintf(csv_fid,'%s\n',attrib{v});
  else
    fprintf(csv_fid,'%s,',attrib{v});
  end
end

%% Write data to file
for v = 1:stm.numPoints
  if any(stm.obsData(v,:,layer)) % skip points without results
  
    if isempty(strmatch(upper(stm.techniqueAttrib.crs),{'WGS84','EPSG4326'}))
      if strcmp(opt.unit,'mm')
        fprintf(csv_fid,['%s,%12.8f,%12.8f,%12.2f,%12.2f,' repmat('%8.2f,',1,2*stm.numEpochs-1) '%8.2f\n'],...
                stm.pntName{v},stm.pntCrd(v,[2 1]),stm.pntAttrib.pntCrdOut(v,1:2), ...
                stm.obsData(v,:,layer),sqrt(stm.stochData(v,:,layer)));
      elseif strcmp(opt.unit,'m')
        fprintf(csv_fid,['%s,%12.8f,%12.8f,%12.2f,%12.2f,' repmat('%8.5f,',1,2*stm.numEpochs-1) '%8.5f\n'],...
                stm.pntName{v},stm.pntCrd(v,[2 1]),stm.pntAttrib.pntCrdOut(v,1:2), ...
                1e-3*stm.obsData(v,:,layer),1e-3*sqrt(stm.stochData(v,:,layer)));    
      else
        error('The unit you requested is not supported.');
      end
    else  
      if strcmp(opt.unit,'mm')
        fprintf(csv_fid,['%s,%12.8f,%12.8f,' repmat('%8.2f,',1,2*stm.numEpochs-1) '%8.2f\n'],...
                stm.pntName{v},stm.pntCrd(v,[2 1]),stm.obsData(v,:,layer),sqrt(stm.stochData(v,:,layer)));
      elseif strcmp(opt.unit,'m')
        fprintf(csv_fid,['%s,%12.8f,%12.8f,' repmat('%8.5f,',1,2*stm.numEpochs-1) '%8.5f\n'],...
                stm.pntName{v},stm.pntCrd(v,[2 1]),1e-3*stm.obsData(v,:,layer),1e-3*sqrt(stm.stochData(v,:,layer)));
      else
        error('The unit you requested is not supported.');
      end
    end
  end
end
%% Close file
fclose(csv_fid);

end


function exportshape(stm,runId,opt,layer)
% EXPORTSHAPE Export of Space-Time Matrix in shape format.
%  EXPORTSHAPE(STM,RUNID,OPT,LAYER) Function to export an STM in
%  shape format. The RUNID provides the ID for the output, used to
%  create an output directory, containing the original filename of
%  the STM and the timestamp of the export. The LAYER indicates
%  the index of the layer exported. OPT contains the options
%  as specified when calling the parent function stmexport. Here
%  used:
%
%  'unit'          unit of output, 'm' or 'mm' (default)
%
%  NOTE: this function is not fully implemented yet!
%
%  Alternative export functions:
%    - EXPORTCSV       .csv file
%    - EXPORTGEOTIFF   geotiff file
%    - EXPORTPNG       .png file
%
%  (c) Freek van Leijen, Delft University of Technology, 2023. 

% Created:  16 Nov 2023 by Freek van Leijen
% Modified: 
%         - 31 March 2024 by Freek van Leijen
%           Added runId input parameter for unique output directory.
%           Use of datasetId instead of techniqueId for output filename.
%
% TODO:
%         - Finish the script.
%

%% Create output directory
shapeDir = [runId '_shp'];
if ~exist(['./' shapeDir],'dir')
  mkdir(shapeDir);
end

shapeFile = [shapeDir '/' stm.datasetId '_' stm.par{layer} '_' ...
                stm.space{1} '_' stm.ref{1} '.shp'];
prjFile = [shapeDir '/' stm.datasetId '_' stm.par{layer} '_' ...
                   stm.space{1} '_' stm.ref{1} '.prj']

%% Get coordinate system
switch stm.ref{1}
  case 'RD'
    prj = 'PROJCS["RD_New",GEOGCS["GCS_Amersfoort",DATUM["D_Amersfoort",SPHEROID["Bessel_1841",6377397.155,299.1528128]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Double_Stereographic"],PARAMETER["False_Easting",155000.0],PARAMETER["False_Northing",463000.0],PARAMETER["Central_Meridian",5.38763888888889],PARAMETER["Scale_Factor",0.9999079],PARAMETER["Latitude_Of_Origin",52.15616055555555],UNIT["Meter",1.0]]';
  case 'WGS84'
    prj = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]';
  otherwise
    error('The specified reference/projection is not supported yet.');
end

%% Create header
attrib = { 'Geometry' 'Lat' 'Lon' 'ID' };
numAttrib = length(attrib);

for v = 1:numel(stm.epochDyear)
  attrib{numAttrib+v} = ['y' num2str(stm.epochDyear(v))];
end
attrib = matlab.lang.makeValidName(attrib);
numAttrib = length(attrib);

%% Collect data
data = zeros( stm.numPoints, numAttrib );
if strcmp(opt.unit,'mm')
  data(:,[2:3 5:end]) = [stm.pntCrd(:,1) stm.pntCrd(:,2) stm.obsData(:,:,layer)];
elseif strcmp(opt.unit,'m')
  data(:,[2:3 5:end]) = [stm.pntCrd(:,1) stm.pntCrd(:,2) 1e-3*stm.obsData(:,:,layer)];
else
  error('The unit you requested is not supported.');
end
data = num2cell(data);
data(:,1) = cellstr(repmat('Point',stm.numPoints,1));
data(:,4) = stm.pntName;

%% Write shape file
S = cell2struct( data, attrib, 2 );
shapewrite(S,shapeFile);
prj_fid = fopen(prjFile,'w');
fprintf(prj_fid,'%s',prj);
fclose(prj_fid);

end


function exportgeotiff(stm,runId,opt,layer,cmap)
% EXPORTGEOTIFF Export of Space-Time Matrix in geotiff format.
%  EXPORTGEOTIFF(STM,RUNID,OPT,LAYER) Function to export an STM in
%  geotiff format. The RUNID provides the ID for the output, used to
%  create an output directory, containing the original filename of
%  the STM and the timestamp of the export. The LAYER indicates
%  the index of the layer exported. OPT contains the options
%  as specified when calling the parent function stmexport. Here
%  used:
%
%  'unit'          unit of output, 'm' or 'mm' (default)
%
%  NOTE: this function is not fully implemented yet!
%
%  Alternative export functions:
%    - EXPORTCSV       .csv file
%    - EXPORTSHAPE     shape file
%    - EXPORTPNG       .png file
%
%  (c) Freek van Leijen, Delft University of Technology, 2023. 

% Created:  16 Nov 2023 by Freek van Leijen
% Modified: 
%         - 31 March 2024 by Freek van Leijen
%           Added runId input parameter for unique output directory.
%           Use of datasetId instead of techniqueId for output filename.
%
% TODO:
%         - Finish the script.
%

% Create output directory and file
geotiffDir = [runId '_geotiff'];
if ~exist(['./' geotiffDir],'dir')
  mkdir(geotiffDir);
end

outputFilename = [geotiffDir '/' stm.datasetId '_' stm.par{layer} '_' ...
                    stm.space{1} '_' stm.ref{1} '.tif']

%% Collect data
if strcmp(opt.unit,'mm')
  output = stm.obsData(:,:,layer);
elseif strcmp(opt.unit,'m')
  output = 1e-3*stm.obsData(:,:,layer);
else
  error('The unit you requested is not supported.');
end

% Create geotiff
[Nx,Ny,Nz] = size(stm);

R = maprasterref();
R.XLimWorld = [xmin xmax];
R.YLimWorld = [ymin ymax];
R.RasterSize = [Ny Nx];
R.ColumnsStartFrom = 'south';
R.RowsStartFrom = 'west';
%R.RasterInterpretation = 'cells';
%R.RasterInterpretation = 'postings';
%R.DeltaX = dx_out;
%R.DeltaY = dx_out;
%R.RasterWidthInWorld = R.XLimWorld(2)-R.XLimWorld(1);
%R.RasterHeightInWorld = R.YLimWorld(2)-R.YLimWorld(1);
%R.XLimIntrinsic = [0.5000 Nx+0.5];
%R.YLimIntrinsic = [0.5000 Ny+0.5];
%R.TransformationType = 'rectilinear';
%R.CoordinateSystemType = 'planar';

%% Write geotiff
if ~isempty(cmap)
  geotiffwrite(outputFilename, output, cmap, R, 'CoordRefSysCode', coordRefSysCode);
else
  geotiffwrite(outputFilename, output, R, 'CoordRefSysCode', coordRefSysCode);
end

%%Problem: nodata value
%%[status, result] = system(sprintf('gdal_translate -a_nodata %g -of GTiff %s %s', nodata, tmpfilename, filename));

end


function exportpng(stm,runId,opt,layer)
% EXPORTPNG Export of dataset in .png image format.
%  EXPORTPNG(STM,RUNID,OPT,LAYER) exports an STM file in .png image format,
%  for a certain LAYER in the dataset. The RUNID is used to create
%  a unique output directory based on a datestring. The output image can 
%  have the form of points, a grid, or contours (TO BE IMPLEMENTED).
%
%  NOTE: currently, this function can only be called for Space-Time matrices
%  created by the STMPREDICT function. The objective is to enable 
%  calling this function for any Space Time matrix dataset.
%
%  Alternative export functions:
%    - EXPORTCSV       .csv file
%    - EXPORTSHAPE     shape file
%    - EXPORTGEOTIFF   geotiff file
%
%  (c) Freek van Leijen, Delft University of Technology, 2023. 

% Created:  16 Nov 2023 by Freek van Leijen
% Modified: 
%         - 31 March 2024 by Freek van Leijen
%           Added runId input parameter for unique output directory.
%           Use of datasetId instead of techniqueId for output filename.
%         - 14 November 2024 by Freek van Leijen
%           Added background, and some minor bug fixes.
%
% TODO:
%         - Make the function compatible with all STM files.
%         - Implement contour plot option.
%
 
pngDir = [runId '_png'];
if ~exist(['./' pngDir],'dir')
  mkdir(pngDir);
end

% Get background and deformation border

if isnumeric(opt.background)
    igpbackground = opt.background;
elseif exist(opt.background,'file')
    igpbackground = load(opt.background);
    if isstruct(igpbackground)
        if isfield(igpbackground,'Lat')
          igpbackground = [ igpbackground.Lat(:) igpbackground.Lon(:) ];
        elseif isfield(igpbackground,'X');
          igpbackground = [ igpbackground.X(:) igpbackground.Y(:) ];
        else
          error('The opt.background does not have the correct format.');
        end
    end
else
    igpbackground=[];
end

if ~isfield(stm,'temporalRef')
  epochRef = stm.epochDyear(1);
else
  epochRef = stm.temporalRef;
end
    
if strmatch(opt.cmap,'defo')
  load('cmap_defo');
  cmap = cmap_defo;
else
  cmap = parula;
end

cmap2 = flipud(hot);

clab = 'mm';
clab2 = 'mm^2';

markersize = 10;
fontsize = 14;

cmin = nanmin(nanmin(stm.obsData(:,:,layer)));
cmax = nanmax(nanmax(stm.obsData(:,:,layer)));
cmin2 = nanmin(nanmin(stm.stochData(:,:,layer)));
cmax2 = nanmax(nanmax(stm.stochData(:,:,layer)));

switch stm.techniqueAttrib.space{1}

  % Points
  case {'original','list'}

   for v = 1:stm.numEpochs

     % Plot displacement
     figure;hold on;
     if strcmp(opt.visible,'off')
       set(gcf,'visible','off');
     end
     set(gcf,'PaperPositionMode','auto','InvertHardcopy','off');
     set(gcf,'color','w');
     if strmatch('deg',stm.techniqueAttrib.space{3})
       scatter(stm.pntAttrib.pntCrdOut(:,2),stm.pntAttrib.pntCrdOut(:,1), ...
         markersize,stm.obsData(:,v,layer),'filled');caxis([cmin cmax]);colormap(cmap);k=colorbar;
       set(get(k,'title'),'string',clab,'fontweight','bold','fontsize',fontsize-2);
       xl=xlim();yl=ylim(); 
       if ~isempty(igpbackground)
         plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
       end
       xlim(xl);ylim(yl);
       daspect([1/cosd(mean(yl)) 1 1])
       xlabel(['Lon [deg] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
       ylabel(['Lat [deg] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
     elseif strmatch('m',stm.techniqueAttrib.space{3})
       scatter(0.001*stm.pntAttrib.pntCrdOut(:,1),0.001*stm.pntAttrib.pntCrdOut(:,2), ...
         markersize,stm.obsData(:,v,layer),'filled');caxis([cmin cmax]);colormap(cmap);k=colorbar;
       set(get(k,'title'),'string',clab,'fontweight','bold','fontsize',fontsize-2);
       xl=xlim();yl=ylim(); 
       if ~isempty(igpbackground)
         plot(0.001*igpbackground(:,1),0.001*igpbackground(:,2),'Color',[.5 .5 .5])
       end
       xlim(xl);ylim(yl);
       axis equal
       xlabel(['X [km] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
       ylabel(['Y [km] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
     end
     title([stm.obsTypes{layer} ' ' num2str(epochRef,'%6.2f') ' - ' ...
            num2str(stm.epochDyear(v),'%6.2f')])

     print([pngDir '/' stm.datasetId '_' stm.obsTypes{layer} '_' ...
       num2str(epochRef,'%6.2f') ' - ' num2str(stm.epochDyear(v),'%6.2f') '.png'], '-dpng');

     % Plot variance
     figure;hold on;
     if strcmp(opt.visible,'off')
       set(gcf,'visible','off');
     end
     set(gcf,'PaperPositionMode','auto','InvertHardcopy','off');
     set(gcf,'color','w');
     if strmatch('deg',stm.techniqueAttrib.space{3})
       scatter(stm.pntAttrib.pntCrdOut(:,2),stm.pntAttrib.pntCrdOut(:,1), ...
         markersize,stm.stochData(:,v,layer),'filled');caxis([cmin2 cmax2]);colormap(cmap2);k=colorbar;
       set(get(k,'title'),'string',clab2,'fontweight','bold','fontsize',fontsize-2);
       xl=xlim();yl=ylim(); 
       if ~isempty(igpbackground)
         plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
       end
       xlim(xl);ylim(yl);
       daspect([1/cosd(mean(yl)) 1 1])
       xlabel(['Lon [deg] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
       ylabel(['Lat [deg] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
     elseif strmatch('m',stm.techniqueAttrib.space{3})
       scatter(0.001*stm.pntAttrib.pntCrdOut(:,1),0.001*stm.pntAttrib.pntCrdOut(:,2), ...
         markersize,stm.stochData(:,v,layer),'filled');caxis([cmin2 cmax2]);colormap(cmap2);k=colorbar;
       set(get(k,'title'),'string',clab2,'fontweight','bold','fontsize',fontsize-2);
       xl=xlim();yl=ylim(); 
       if ~isempty(igpbackground)
         plot(0.001*igpbackground(:,1),0.001*igpbackground(:,2),'Color',[.5 .5 .5])
       end
       xlim(xl);ylim(yl);
       axis equal
       xlabel(['X [km] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
       ylabel(['Y [km] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
     end
     title([stm.obsTypes{layer} ' error variance ' num2str(epochRef,'%6.2f') ' - ' ...
            num2str(stm.epochDyear(v),'%6.2f')])

     print([pngDir '/' stm.datasetId '_' stm.obsTypes{layer} '_error_variance_' ...
       num2str(epochRef,'%6.2f') ' - ' num2str(stm.epochDyear(v),'%6.2f') '.png'], '-dpng');

   end

  % Grid
  case {'grid'}

   for v = 1:stm.numEpochs

     % Plot displacement
     figure;hold on;
     if strcmp(opt.visible,'off')
       set(gcf,'visible','off');
     end
     set(gcf,'PaperPositionMode','auto','InvertHardcopy','off');
     set(gcf,'color','w');
     if strmatch('deg',stm.techniqueAttrib.space{3})
       datalayer = reshape(stm.obsData(:,v,layer),stm.techniqueAttrib.grid.numY, ...
                           stm.techniqueAttrib.grid.numX);
       h = imagesc(stm.techniqueAttrib.grid.X, stm.techniqueAttrib.grid.Y, ...
                 datalayer);caxis([cmin cmax]);colormap(cmap);k=colorbar;
       set(h,'alphadata',~isnan(datalayer));
       set(get(k,'title'),'string',clab,'fontweight','bold','fontsize',fontsize-2);
       xl=xlim();yl=ylim(); 
       if ~isempty(igpbackground)
         plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
       end
       xlim(xl);ylim(yl);
       daspect([1/cosd(mean(yl)) 1 1])
       xlabel(['Lon [deg] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
       ylabel(['Lat [deg] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
     elseif strmatch('m',stm.techniqueAttrib.space{3})
       datalayer = reshape(stm.obsData(:,v,layer),stm.techniqueAttrib.grid.numY, ...
                           stm.techniqueAttrib.grid.numX);
       h = imagesc(0.001*stm.techniqueAttrib.grid.X, 0.001*stm.techniqueAttrib.grid.Y, ...
                 datalayer);caxis([cmin cmax]);colormap(cmap);k=colorbar;
       set(h,'alphadata',~isnan(datalayer));
       set(get(k,'title'),'string',clab,'fontweight','bold','fontsize',fontsize-2);
       xl=xlim();yl=ylim(); 
       if ~isempty(igpbackground)
         plot(0.001*igpbackground(:,1),0.001*igpbackground(:,2),'Color',[.5 .5 .5])
       end
       xlim(xl);ylim(yl);
       axis equal
       xlabel(['X [km] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
       ylabel(['Y [km] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
     end
     title([stm.obsTypes{layer} ' ' num2str(epochRef,'%6.2f') ' - ' ...
            num2str(stm.epochDyear(v),'%6.2f')])

     print([pngDir '/' stm.datasetId '_' stm.obsTypes{layer} '_' ...
       num2str(epochRef,'%6.2f') ' - ' num2str(stm.epochDyear(v),'%6.2f') '.png'], '-dpng');

     % Plot variance
     figure;hold on;
     if strcmp(opt.visible,'off')
       set(gcf,'visible','off');
     end
     set(gcf,'PaperPositionMode','auto','InvertHardcopy','off');
     set(gcf,'color','w');
     if strmatch('deg',stm.techniqueAttrib.space{3})
       datalayer = reshape(stm.stochData(:,v,layer),stm.techniqueAttrib.grid.numY, ...
                           stm.techniqueAttrib.grid.numX);
       h = imagesc(stm.techniqueAttrib.grid.X, stm.techniqueAttrib.grid.Y, ...
                 datalayer);caxis([cmin2 cmax2]);colormap(cmap2);k=colorbar;
       set(h,'alphadata',~isnan(datalayer));
       set(get(k,'title'),'string',clab2,'fontweight','bold','fontsize',fontsize-2);
       xl=xlim();yl=ylim(); 
       if ~isempty(igpbackground)
         plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
       end
       xlim(xl);ylim(yl);
       daspect([1/cosd(mean(yl)) 1 1])
       xlabel(['Lon [deg] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
       ylabel(['Lat [deg] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
     elseif strmatch('m',stm.techniqueAttrib.space{3})
       datalayer = reshape(stm.stochData(:,v,layer),stm.techniqueAttrib.grid.numY, ...
                           stm.techniqueAttrib.grid.numX);
       h = imagesc(0.001*stm.techniqueAttrib.grid.X, 0.001*stm.techniqueAttrib.grid.Y, ...
                 datalayer);caxis([cmin2 cmax2]);colormap(cmap2);k=colorbar;
       set(h,'alphadata',~isnan(datalayer));
       set(get(k,'title'),'string',clab2,'fontweight','bold','fontsize',fontsize-2);
       xl=xlim();yl=ylim(); 
       if ~isempty(igpbackground)
         plot(0.001*igpbackground(:,1),0.001*igpbackground(:,2),'Color',[.5 .5 .5])
       end
       xlim(xl);ylim(yl);
       axis equal
       xlabel(['X [km] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
       ylabel(['Y [km] (' stm.techniqueAttrib.crs ')'],'fontsize',fontsize);
     end
     title([stm.obsTypes{layer} ' error variance ' num2str(epochRef,'%6.2f') ' - ' ...
            num2str(stm.epochDyear(v),'%6.2f')])

     print([pngDir '/' stm.datasetId '_' stm.obsTypes{layer} '_error_variance_' ...
       num2str(epochRef,'%6.2f') ' - ' num2str(stm.epochDyear(v),'%6.2f') '.png'], '-dpng');

   end

  % Contour
  case {'contour'}
   %Not implemented yet.
   %contour(x,y,vx);
end

end

