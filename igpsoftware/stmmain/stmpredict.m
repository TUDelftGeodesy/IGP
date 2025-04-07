function stmpredict(inputFilename,outputFilename,options,flag)
%stmpredict   Prediction based on Space Time Matrix.
%   STMPREDICT(INPUTFILENAME,OUTPUTFILENAME,OPTIONS) predicts
%   observations based on the input Space Time Matrix, given by INPUTFILENAME. 
%   The results are writen to the output Space Time Matrix OUTPUTFILENAME.
%   OPTIONS is a cell array or structure with the processing options, if 
%   empty, or when fields are missing, default values are used.
%
%   The output is generated for a certain Region of Interest (ROI) and 
%   Period of Interest (POI) (full spatio-temporal coverage in case the 
%   input OPTIONS are empty/missing). Also, a desired Coordinate Reference
%   System (CRS) for the output product can be indicated (default WGS84).
%
%   Examples:
%      stmpredict('integratedstm.mat','predictedstm.mat',options)      
%
%   Examples options:
%
%   % Set coordinate transformation toolbox, RDNAPTRANS or PROJ4
%   options.crstoolbox = 'RDNAPTRANS';
%   options.crstoolbox = 'rdnaptrans';
%   options.crstoolbox = 'proj4';
%   options.crstoolbox = 'PROJ4';
%
%   % Set Region of interest, as [lonmin latmin ; lonmax latmax]
%   % or [xmin ymin; xmax ymax] bounding box, or 
%   % polygon (.kml/.shp) (Default: none, resulting in full 
%   % extent of input stm).
%   options.ROI = 'area.shp'; % .kml, .shp
%   options.ROI = 'pointlist.txt'; % .txt, .csv, .xlsx
%   options.ROI = [52.3434 6.34343; 52.2534 6.34342];
%   options.ROI = [230000 570000; 250000 590000];
%   options.ROI = 'groningen_test_roi.kml';
%   options.ROI = [];
%
%   % Set Coordinate Reference System (CRS) of ROI.
%   options.crsROI = { 'UTM' };
%   options.crsROI = { 'WGS84' };
%   options.crsROI = { 'RD' };
%   options.crsROI = { 'EPSG28992' };
%   options.crsROI = [];
%
%   % Set Period of Interest (in dates or decimal years)
%   % periods {{start,stop},{start,stop}, ...}, epoch list (cellarray,.txt, .csv, .xlsx)
%   options.POI = { {'1984-01-01','2021-01-01'} }; % Date format should be 'YYYY-MM-DD' (for now)
%   options.POI = { {'1984-01-01','1990-01-01'} {'2000-01-01','2021-01-01'} };
%   options.POI = [ 1984.00 2021.00 ];
%   options.POI = [ 1984.00 1990.00; 2000.00 2021.00 ];
%   options.POI = [];
%
%   % Set Spatial Reference
%   options.spatialRef = [];
%   options.spatialRef = {'unstablearea.coo' 'EPSG28992' 'reverse'}; %
%   options.spatialRef = {'filename.kml' 'wgs84' 'normal'};
%
%   % Set Temporal Reference
%   options.temporalRef = [];
%   options.temporalRef = 1984;
%
%   % Set Coordinate Reference System (CRS) of the output
%   options.crsOutput = 'UTM';
%   options.crsOutput = 'WGS84';
%   options.crsOutput = 'RD';
%   options.crsOutput = 'EPSG28992';
%
%   % Set spatial output, the orignal points (default), a grid, or a list of coordinates
%   options.spaceOutput = { 'original' };
%   options.spaceOutput = { 'grid' 500 'm' };
%   options.spaceOutput = { 'grid' 0.001 'deg' };
%   options.spaceOutput = { 'list' 'pointlist.txt'};
%
%   % Set temporal output, the original epochs (default), a sampling, or a list of epochs
%   options.timeOutput = { 'original' };
%   options.timeOutput = { 'sampling' 1 'year' }; % e.g., 0.2 year
%   options.timeOutput = { 'list' 'epochlist.txt' };
%
%   % Set use of stochastic model for velocity estimation, true or false
%   options.ignoreStochModel = true;
%
%   % Set reference system (see stmvelocity.m)
%   options.refsystem = = 'min-velocities' (default here), using options.spatialRef as reference
%
%   % Set covariance function of signal ('estimated' or one model per par)
%   % Use double cells, can be multiple models
%   % e.g., options.signalModel = { {'[modelstring]' '[modelstring]'} {'[modelstring]'} };
%   options.signalModel = { 'estimated' }; %NOT IMPLEMENTED YET
%   options.signalModel = {{'tudpred1(s20=NaN,s2t=3,s2s=3,Rt=1,Rs=5)'} ...
%                          {'tudpred1(s20=NaN,s2t=2,s2s=2,Rt=1,Rs=3)'}};
%   options.signalModel = {{'tudpred1(s20=NaN,s2t=4,s2s=3,Rt=1,Rs=4)'} ...
%                          {'tudpred1(s20=NaN,s2t=2,s2s=2,Rt=1,Rs=3)'}};
%   options.signalModel = {{'sure(s2z=0.8,L=4,q=0.85,tref=2003)'} ...
%                          {'sure(s2z=0.8,L=4,q=0.85,tref=2003)'} ...
%                          {'sure(s2z=0.8,L=4,q=0.85,tref=2003)'} ...
%                         };
%   options.signalModel = {{'sure(s2z=0.8,L=4,q=0.85,tlag=0)'} ... % tlag=1 means tref becomes 1 year earlier, hence 1984 -> 1983
%                          {'sure(s2z=0.8,L=4,q=0.85,tlag=0)'} ...
%                          {'sure(s2z=0.8,L=4,q=0.85,tlag=0)'} ...
%                         };
%   options.signalModel = {{'sure(s2z=0.8,L=4.0,q=0.85,tref=2010)'} ...
%                          {'houtenbos(s2t=0.16,pt=1.24)'}}; 
%
%   % Set convergence threshold for iterative prediction (sum of
%   % absolute values) [mm]
%   options.predConvThres = 10000;
%
%   % Set maximum number of iterations for prediction
%   options.predMaxIter = 10;
%
%   % Set interpolation method (for trend interpolation), 'linear' for triangulation-based 
%   % linear interpolation or 'idw' for inverse distance weighting.
%   options.interpMethod = 'linear';
%   options.interpMethod = 'idw';
%
%   % Set power for inverse distance weighting
%   options.idwPower = 2;
%   options.idwPower = 4;
%
%   % Set background for plotting (default [] is none)
%   options.background = 'igp-background.mat';
%
%   % Set plot mode, 'last' for last epoch only (default) or 'timeseries'
%   options.plotmode = 'last';
%   options.plotmode = 'timeseries';
%
%   Directory for pdf plots (default [] is none), automatically generated
%   using 'generated'
%   options.saveplt = [];
%   options.saveplt = 'generated';
%

%
%  (c) Freek van Leijen, Delft University of Technology, 2020.

% Created:   19 October 2020 by Freek van Leijen
% Modified:  02 January 2024 by Freek van Leijen
%            - restructured in separate remove, compute, restore steps
%            - documented
%            21 March 2024 by Freek van Leijen
%            - removed the use of the project STM
%            15 October 2024 by Freek van Leijen
%            - added clustering of points
%            30 October 2024 by Freek van Leijen
%            - added a check on the spatial overlap between the input
%              and output point locations.
%            - converted all crs descriptions to EPSG codes for consistency.
%

%% Check the input arguments and options

% For testing purposes we don't call this as function, but hardcode the inputs.

if nargin < 3
    error('This function expects at least three input arguments.')
end

progname = 'stmpredict';

% Default options

opt.crstoolbox = 'proj4';          % Coordinate transformation method
opt.ROI = [];                      % Region of interest, as [lonmin latmin ; lonmax latmax]
                                   % or [xmin ymin; xmax ymax] bounding box, or 
                                   % polygon (.kml/.shp) (Default: none, resulting in full 
                                   % extent of input stm).
opt.crsROI = 'WGS84';              % Coordinate Reference System (CRS) of ROI.
opt.POI = [];                      % Period of interest [ dYearStart dYearEnd ] (Default none).
opt.spatialRef = [];               % Spatial reference.
opt.temporalRef = [];              % Temporal reference.

opt.crsOutput = 'WGS84';           % Coordinate Reference System (CRS) of output product.
opt.spaceOutput = 'original';      % Spatial output specification. Default is original points.
opt.timeOutput = 'original';       % Temporal output specification. Default is original epochs.

opt.ignoreStochModel = false;      % If true, ignore stochastic model from stm, use 
                                   % unit matrix instead (default false).
opt.refsystem = 'min-velocities';  % Choice of reference system (default 'min-velocities').

opt.signalModel = [];              % Covariance model of signal.

opt.clusteringMethod = 'distance'; % clustering method, {'benchmark','distance'}
opt.spatialTolerance = 1;          % [m], maximum distance between benchmarks
opt.doplots = 0;                   % make plots

opt.predConvThres = 10000;         % Convergence threshold for iterative prediction (sum of
                                   % absolute values) [mm].
opt.predMaxIter = 10;              % Maximum number of iterations for iterative prediction.

opt.interpMethod = 'linear';       % Interpolation method, 'linear' for triangulation-based 
                                   % linear interpolation or 'idw' for inverse distance weighting.
opt.idwPower = 2;                  % Power for inverse distance weighting (if used).


opt.background = [];               % Map background (default [] is none).
opt.plotmode = 'last';             % Plot mode, 'last' for last epoch only (default) or 'timeseries'.
opt.saveplt = [];                  % Directory for pdf plots (default [] is none).


%% Catch errors and start timing
%try

[~,outputFileroot]=fileparts(outputFilename);
runId = [ outputFileroot '_' datestr(now,30) ];
diary([ runId '.log' ])
  
fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Check the options and if necessary overwrite the default values

[inputFilename,outputFilename,opt]= ...
    stmcheckarguments({inputFilename},outputFilename,opt,options,flag);
if isempty(outputFilename)
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    return;
end
inputFilename = inputFilename{1}; % Convert to character string

if strcmp(opt.saveplt,'generated')
  opt.saveplt = [runId '_figures'];
end

% Convert crs to EPSG codes for consistency
if strmatch(upper(opt.crsOutput),{'RD','RDNAP'})
  opt.crsOutput = 'EPSG28992';
elseif strmatch(upper(opt.crsOutput),{'WGS84','ETRS','ITRS'})
  opt.crsOutput = 'EPSG4326';
end

if strmatch(upper(opt.crsROI),{'RD','RDNAP'})
  opt.crsROI = 'EPSG28992';
elseif strmatch(upper(opt.crsOutput),{'WGS84','ETRS','ITRS'})
  opt.crsROI = 'EPSG4326';
end

%% Read space time matrix dataset

stmIn = stmread(inputFilename);
projectId = stmIn.datasetAttrib.projectId;

if ~strcmp('estimate',opt.signalModel) & numel(stmIn.obsTypes)~=numel(opt.signalModel)
  error('The number of observation types and covariance models does not match.');
end

% Get pntNeu and pntCrdRef
if isfield(stmIn.pntAttrib,'pntCrdRef')
  % Get from project file
  plh0 = stmIn.techniqueAttrib.pntCrdRef;
  if ~isfield(stmIn.pntAttrib,'pntNeu') %TODO check whether this is correct
    stmIn.pntAttrib.pntNeu = stmIn.pntAttrib.pntNeu;
  end
else
  % Convert latitude/longitude into local topocentric coordinates, and
  % keep plh0 for later reuse for consistency
  [stmIn.pntAttrib.pntNeu,plh0] = plh2neusp(stmIn.pntCrd);     % deg/deg/m -> m/m/m
  stmIn.pntAttrib.pntNeu(:,1:2) = stmIn.pntAttrib.pntNeu(:,1:2)./1000;   % m/m/m -> km/km/m
end

%% Create output stm
stmOut = stm(projectId,'predicted','displ');
stmOut.techniqueAttrib.crs = opt.crsOutput;
stmOut.techniqueAttrib.signalModel = opt.signalModel;

%% Set spatial output

% TODO: separate ROI and sampling (as for time), so also ROI for original points
% 1 Get ROI
% 2 Transform to desired output format
% 3 Create grid
% 4 Transform to lat,lon
% 5 Transform to NEU (for prediction)

stmOut.techniqueAttrib.space = opt.spaceOutput;
switch opt.spaceOutput{1}

  case 'original'
    stmOut.numPoints = stmIn.numPoints;
    stmOut.pntName = stmIn.pntName;
    stmOut.pntCrd = stmIn.pntCrd;
    stmOut.pntAttrib.pntIdx = 1:stmOut.numPoints; % Update after ROI

    % Create output coordinates
    if strmatch(upper(stmOut.techniqueAttrib.crs),'EPSG4326')
      stmOut.pntAttrib.mapCrd = stmOut.pntCrd;
    else
      if isfield(stmIn,'mapCrd') 
        stmOut.pntAttrib.mapCrd = stmIn.mapCrd; % Use the original coordinates
      else
        stmOut.pntAttrib.mapCrd = crstrans(stmOut.pntCrd,'EPSG4326',stmOut.techniqueAttrib.crs,opt.crstoolbox);
      end
    end

  case {'grid','contour'}
    if isempty(opt.ROI)
      % Always in lat, lon BBox = [latmin,lonmin;latmax,lonmax]
      ROI.BoundingBox = [min(stmIn.pntCrd(:,1)) min(stmIn.pntCrd(:,2)); ... 
                         max(stmIn.pntCrd(:,1)) max(stmIn.pntCrd(:,2))];
      opt.crsROI = 'EPSG4326'; % Overwrite potential other input
    elseif ischar(opt.ROI)
      ROI = roi2poly(opt.ROI,'geostruct'); % get polygon from file
    elseif isnumeric(opt.ROI)
      if (size(opt.ROI,1)==2) & (size(opt.ROI,2)==2)
        ROI.BoundingBox = opt.ROI;
      else
        error('The ROI bounding box you specified has a wrong dimension, should be 2x2.');
      end
    else
      error('You specified an unsupported ROI.');
    end

    % Transform ROI in desired output CRS
    if isempty(strmatch(upper(opt.crsROI),upper(stmOut.techniqueAttrib.crs),'exact')) % Different CRS, transform polygon/Bounding Box

      if isfield(ROI,'Lon') % Transform polygon, output always as X,Y (even for lon, lat)

        nanIdx = ~isnan(ROI.Lon); % Take care of NaNs, e.g. for multi-polygons
        stmOut.techniqueAttrib.ROI.X = NaN(size(ROI.Lon));
        stmOut.techniqueAttrib.ROI.Y = NaN(size(ROI.Lon));

        ROImap = crstrans([ROI.Lat(nanIdx) ROI.Lon(nanIdx)], ...
                          opt.crsROI,stmOut.techniqueAttrib.crs,opt.crstoolbox);
        [stmOut.techniqueAttrib.ROI.X(nanIdx),stmOut.techniqueAttrib.ROI.Y(nanIdx)] = ROImap(:,1:2);

        % Get associated Bounding Box
        stmOut.techniqueAttrib.ROI.BoundingBox = [min(stmOut.techniqueAttrib.ROI.X(nanIdx)) ...
                                                min(stmOut.techniqueAttrib.ROI.Y(nanIdx)); ...
                                                max(stmOut.techniqueAttrib.ROI.X(nanIdx)) ...
                                                max(stmOut.techniqueAttrib.ROI.Y(nanIdx));];

      elseif isfield(ROI,'X') % Transform polygon, output always as X,Y (even for lon, lat)

        nanIdx = ~isnan(ROI.X); % Take care of NaNs, e.g. for multi-polygons
        stmOut.techniqueAttrib.ROI.X = NaN(size(ROI.X));
        stmOut.techniqueAttrib.ROI.Y = NaN(size(ROI.X));

        ROImap = crstrans([ROI.X(nanIdx) ROI.Y(nanIdx)], ...
                          opt.crsROI,stmOut.techniqueAttrib.crs,opt.crstoolbox);
        [stmOut.techniqueAttrib.ROI.X(nanIdx),stmOut.techniqueAttrib.ROI.Y(nanIdx)] = ROImap(:,1:2);
        

        % Get associated Bounding Box
        stmOut.techniqueAttrib.ROI.BoundingBox = [min(stmOut.techniqueAttrib.ROI.X(nanIdx)) ...
                                                min(stmOut.techniqueAttrib.ROI.Y(nanIdx)); ...
                                                max(stmOut.techniqueAttrib.ROI.X(nanIdx)) ...
                                                max(stmOut.techniqueAttrib.ROI.Y(nanIdx));];

      else % Bounding Box only

        ROImap = crstrans([ROI.BoundingBox(:,1),ROI.BoundingBox(:,2)],opt.crsROI, ...
                          stmOut.techniqueAttrib.crs,opt.crstoolbox);
        stmOut.techniqueAttrib.ROI.BoundingBox = ROImap(:,1:2);

      end

    else % Copy the roi
      stmOut.techniqueAttrib.ROI = ROI;
    end


    % Get the Bounding Box consistent with the grid size (may be slightly larger)

    dx = stmOut.techniqueAttrib.space{2};
    stmOut.techniqueAttrib.ROI.BoundingBoxNew = ...
         [floor(stmOut.techniqueAttrib.ROI.BoundingBox(1,1)/dx)*dx ...
           ceil(stmOut.techniqueAttrib.ROI.BoundingBox(1,2)/dx)*dx; ...
          floor(stmOut.techniqueAttrib.ROI.BoundingBox(2,1)/dx)*dx ...
           ceil(stmOut.techniqueAttrib.ROI.BoundingBox(2,2)/dx)*dx];

    % Create grid

    if strmatch(upper(stmOut.techniqueAttrib.crs),'EPSG4326')
      lonRange = stmOut.techniqueAttrib.ROI.BoundingBoxNew(1,1):dx: ...
                 stmOut.techniqueAttrib.ROI.BoundingBoxNew(2,1);
      latRange = stmOut.techniqueAttrib.ROI.BoundingBoxNew(1,2):dx: ...
                 stmOut.techniqueAttrib.ROI.BoundingBoxNew(2,2);
      stmOut.techniqueAttrib.grid.X = lonRange;
      stmOut.techniqueAttrib.grid.Y = latRange;
      stmOut.techniqueAttrib.grid.numX = numel(lonRange);
      stmOut.techniqueAttrib.grid.numY = numel(latRange);

      [latGrid,lonGrid] = meshgrid(latRange,lonRange);
      lat = latGrid(:);
      lon = lonGrid(:);

      stmOut.pntAttrib.pntCrdOut = [lat lon zeros(size(lat))];
      stmOut.pntCrd = stmOut.pntAttrib.pntCrdOut;

    else % Can also be lon,lat in other CRS, but stil transformation needed
      xRange = stmOut.techniqueAttrib.ROI.BoundingBoxNew(1,1):dx: ...
               stmOut.techniqueAttrib.ROI.BoundingBoxNew(2,1);
      yRange = stmOut.techniqueAttrib.ROI.BoundingBoxNew(1,2):dx: ...
               stmOut.techniqueAttrib.ROI.BoundingBoxNew(2,2);
      stmOut.techniqueAttrib.grid.X = xRange;
      stmOut.techniqueAttrib.grid.Y = yRange;
      stmOut.techniqueAttrib.grid.numX = numel(xRange);
      stmOut.techniqueAttrib.grid.numY = numel(yRange);

      [xGrid,yGrid] = meshgrid(xRange,yRange);
      stmOut.pntAttrib.pntCrdOut = [xGrid(:) yGrid(:) zeros(size(xGrid(:)))];

      % Transform to WGS84
      stmOut.pntCrd = crstrans([xGrid(:),yGrid(:)],stmOut.techniqueAttrib.crs, ...
                         'EPSG4326',opt.crstoolbox);

    end

  case 'list'
    listCrd = load(opt.spaceOutput{2});
    stmOut.pntAttrib.pntCrdOut = [listCrd(:,1) listCrd(:,2) zeros(size(listCrd(:,1)))];

    if isempty(strmatch(upper(stmOut.techniqueAttrib.crs),'EPSG4326'))
      % Transform to WGS84
      stmOut.pntCrd = crstrans([stmOut.pntAttrib.pntCrdOut(:,1),stmOut.pntAttrib.pntCrdOut(:,2)], ...
                               stmOut.techniqueAttrib.crs,'EPSG4326',opt.crstoolbox);
    else
      stmOut.pntCrd = stmOut.pntAttrib.pntCrdOut;
    end

  otherwise
    error('You indicated an unsupported spatial output format.');

end

stmOut.numPoints = size(stmOut.pntCrd,1);


% Convert latitude/longitude into local topocentric coordinates, using
% plh0 for consistency
stmOut.pntAttrib.pntNeu = plh2neusp(stmOut.pntCrd,plh0);     % deg/deg/m -> m/m/m
stmOut.pntAttrib.pntNeu(:,1:2) = stmOut.pntAttrib.pntNeu(:,1:2)./1000;   % m/m/m -> km/km/m

% Check spatial overlap between input and output points
chull = convhull(stmIn.pntAttrib.pntNeu(:,2),stmIn.pntAttrib.pntNeu(:,1));
inPoly = inpolygon(stmOut.pntAttrib.pntNeu(:,2),stmOut.pntAttrib.pntNeu(:,1), ...
                     stmIn.pntAttrib.pntNeu(chull,2),stmIn.pntAttrib.pntNeu(chull,1));

if isempty(find(inPoly))
  error('All output points lie outside the convex hull of the input points. No proper prediction possible.');
end


%% Set temporal output

stmOut.techniqueAttrib.time = opt.timeOutput;
if isfield(opt,'POI')
  if iscell(opt.POI)
    numPeriods = numel(opt.POI); % number of periods
    stmOut.techniqueAttrib.POI = NaN(numPeriods,2);
    for v = 1:numPeriods
      stmOut.techniqueAttrib.POI(v,1) = date2dyear(opt.POI{v}{1});
      stmOut.techniqueAttrib.POI(v,2) = date2dyear(opt.POI{v}{2});
    end
  elseif isnumeric(opt.POI) %already in decimal years
    if isempty(opt.POI)
      stmOut.techniqueAttrib.POI = [min(stmIn.epochDyear) max(stmIn.epochDyear)];
    else
      stmOut.techniqueAttrib.POI = opt.POI;
    end
  else
    error('Something is wrong with the specified POI.');
  end
else
  stmOut.techniqueAttrib.POI = [min(stmIn.epochDyear) max(stmIn.epochDyear)];
end

%% Get output epochs
stmOut.epochDyear = [];
switch opt.timeOutput{1}
  case 'original' % original epochs
    stmOut.epochAttrib.epochIdx = [];
    for v = 1:size(stmOut.techniqueAttrib.POI,1)
      epochIdx = find(stmIn.epochDyear>=stmOut.techniqueAttrib.POI(v,1) & ...
                      stmIn.epochDyear<=stmOut.techniqueAttrib.POI(v,2));
      stmOut.epochAttrib.epochIdx = [stmOut.epochAttrib.epochIdx epochIdx];
      stmOut.epochDyear = [stmOut.epochDyear stmIn.epochDyear(epochIdx)];
    end

  case 'sampling' % sampling within observed period
    % First, determine possible epoch samples within dataset
    epochDyearSamples = [ceil(min(stmIn.epochDyear)/opt.timeOutput{2})*opt.timeOutput{2} : ...
             opt.timeOutput{2}:floor(max(stmIn.epochDyear)/opt.timeOutput{2})*opt.timeOutput{2}];
    % Second, select epochs within POI
    for v = 1:size(stmOut.techniqueAttrib.POI,1)
      epochIdx = find(epochDyearSamples>=stmOut.techniqueAttrib.POI(v,1) & ...
                      epochDyearSamples<=stmOut.techniqueAttrib.POI(v,2));
      stmOut.epochDyear = [stmOut.epochDyear epochDyearSamples(epochIdx)];
    end

  case 'list' % epochs from list, still apply POI (for easy subsets based on same list)
    epochDyearSamples = load(opt.timeOutput{2});
    for v = 1:size(stmOut.techniqueAttrib.POI,1)
      epochIdx = find(epochDyearSamples>=stmOut.techniqueAttrib.POI(v,1) & ...
                      epochDyearSamples<=stmOut.techniqueAttrib.POI(v,2));
      stmOut.epochDyear = [stmOut.epochDyear epochDyearSamples(epochIdx)];
    end
    stmOut.epochDyear = sort(stmOut.epochDyear);
end

stmOut.numEpochs = numel(stmOut.epochDyear);

%% Set output layers
stmOut.obsTypes = stmIn.obsTypes;
stmOut.techniqueAttrib.numLayers = numel(stmOut.obsTypes);
stmOut.techniqueAttrib.obsIdx = [1:stmOut.techniqueAttrib.numLayers];

%% Set spatial reference
if isempty(opt.spatialRef{1})
  stmOut.spatialRef = [];
elseif ischar(opt.spatialRef{1})
  stmOut.spatialRef.crs = opt.spatialRef{2};
  stmOut.spatialRef.direction = opt.spatialRef{3};
  stmOut.spatialRef.poly = roi2poly(opt.spatialRef{1},'geostruct'); % get polygon from file

  if strmatch(upper(stmOut.spatialRef.crs),'EPSG4326')
    stmOut.spatialRef.polyWGS84 = stmOut.spatialRef.poly;
  else
    % Transform to WGS84
    if isfield(stmOut.spatialRef.poly,'Lon')
      for v = 1:numel(stmOut.spatialRef.poly)
        temp = crstrans([stmOut.spatialRef.poly(v).Lat' stmOut.spatialRef.poly(v).Lon'],stmOut.spatialRef.crs, ...
                          'EPSG4326',opt.crstoolbox);
        stmOut.spatialRef.polyWGS84(v).Geometry = stmOut.spatialRef.poly(v).Geometry;
        stmOut.spatialRef.polyWGS84(v).Lat = temp(:,1)';
        stmOut.spatialRef.polyWGS84(v).Lon = temp(:,2)';
        stmOut.spatialRef.polyWGS84(v).BoundingBox = [min(temp(:,1)) min(temp(:,2));max(temp(:,1)) max(temp(:,2))];
      end
    elseif isfield(stmOut.spatialRef.poly,'X')
      for v = 1:numel(stmOut.spatialRef.poly)
        temp = crstrans([stmOut.spatialRef.poly(v).X' stmOut.spatialRef.poly(v).Y'],stmOut.spatialRef.crs, ...
                          'EPSG4326',opt.crstoolbox);
        stmOut.spatialRef.polyWGS84(v).Geometry = stmOut.spatialRef.poly(v).Geometry;
        stmOut.spatialRef.polyWGS84(v).Lat = temp(:,1)';
        stmOut.spatialRef.polyWGS84(v).Lon = temp(:,2)';      
        stmOut.spatialRef.polyWGS84(v).BoundingBox = [min(temp(:,1)) min(temp(:,2));max(temp(:,1)) max(temp(:,2))];
      end
    else
      error('The polygon of the spatial reference has a wrong format.');
    end

  end
  
  % Convert latitude/longitude into local topocentric coordinates, using
  % plh0 for consistency
  
  stmOut.spatialRef.polyNeuVector.X = [];
  stmOut.spatialRef.polyNeuVector.Y = [];  
  stmOut.spatialRef.polyWGS84Vector.Lat = [];
  stmOut.spatialRef.polyWGS84Vector.Lon = [];  
  for v = 1:numel(stmOut.spatialRef.polyWGS84)
    temp = plh2neusp([stmOut.spatialRef.polyWGS84(v).Lat' ...
                      stmOut.spatialRef.polyWGS84(v).Lon' ...
                      zeros(size(stmOut.spatialRef.polyWGS84(v).Lon'))],plh0); % deg/deg/m -> m/m/m
    temp(:,1:2) = temp(:,1:2)./1000; % m/m/m -> km/km/m
    stmOut.spatialRef.polyNeu(v).Geometry = stmOut.spatialRef.poly(v).Geometry;
    stmOut.spatialRef.polyNeu(v).X = temp(:,2)';
    stmOut.spatialRef.polyNeu(v).Y = temp(:,1)';
    stmOut.spatialRef.polyNeu(v).BoundingBox = [min(temp(:,2)) min(temp(:,1));max(temp(:,2)) max(temp(:,1))];
    stmOut.spatialRef.polyNeuVector.X = [stmOut.spatialRef.polyNeuVector.X [stmOut.spatialRef.polyNeu(v).X NaN]];
    stmOut.spatialRef.polyNeuVector.Y = [stmOut.spatialRef.polyNeuVector.Y [stmOut.spatialRef.polyNeu(v).Y NaN]];

    stmOut.spatialRef.polyWGS84Vector.Lat = [stmOut.spatialRef.polyWGS84Vector.Lat [stmOut.spatialRef.polyWGS84(v).Lat NaN]];
    stmOut.spatialRef.polyWGS84Vector.Lon = [stmOut.spatialRef.polyWGS84Vector.Lon [stmOut.spatialRef.polyWGS84(v).Lon NaN]];
  end
else
  error('The specification of the spatial reference is wrong.');
end

if ~isempty(stmOut.spatialRef)

  % Select points in spatial reference area

  inPoly = inpolygon(stmIn.pntAttrib.pntNeu(:,2),stmIn.pntAttrib.pntNeu(:,1), ...
                     stmOut.spatialRef.polyNeuVector.X,stmOut.spatialRef.polyNeuVector.Y);
                       
  if strmatch(stmOut.spatialRef.direction,'normal')
    pntRefMask = inPoly;
  elseif strmatch(stmOut.spatialRef.direction,'reverse')
    pntRefMask = ~inPoly;
  else
    error('You indicated an invalid direction (normal/reverse) of the spatial reference.');
  end
else
  pntRefMask = [];
  %% Alternative: 1:numPoints. Thereby, the mean will be set as reference
end
fprintf('Number of stable points: %d\n',sum(pntRefMask));


%% Set temporal reference
if isempty(opt.temporalRef)
  % Set temporal reference equal to start POI (could be first epoch in STM)
  stmOut.temporalRef = stmOut.techniqueAttrib.POI(1);
elseif isnumeric(opt.temporalRef)
  stmOut.temporalRef = opt.temporalRef;
else
  error('The specification of the temporal reference is wrong.');
end
  

%% Remove-compute-restore

if strmatch(lower(opt.spaceOutput{1}),'original') & strmatch(lower(opt.timeOutput{1}),'original')
  % No action needed, possible only subset
  stmOut.obsData = stmIn.obsData(stmOut.pntAttrib.pntIdx,stmOut.epochAttrib.epochIdx,:);
  stmOut.stochData = []; % To be updated
else

  %% Select unique points in input STM
  [uniqueMask,commonPoints] = stmclustering(stmIn,'clusteringMethod',opt.clusteringMethod, ...
      'spatialTolerance',opt.spatialTolerance,'doplots',opt.doplots);
  
  %% Remove
  if strcmp(opt.refsystem,'min-velocities')
    % minimize velocities of selected points (pntRefMask)
    [vel,t0,offset,ref,omt,emat,dmat,qx,qe,qy,Atrend] = stmvelocity(stmIn,'relative',true, ...
       'ignoreStochModel',opt.ignoreStochModel,'refsystem',pntRefMask,'t0',stmOut.temporalRef,...
       'pntMask',uniqueMask);
  else
    % default refsystem option (min-refseries) - stochmodel from stm (THIS IS THE DEFAULT)
    [vel,t0,offset,ref,omt,emat,dmat,qx,qe,qy,Atrend] = stmvelocity(stmIn,'relative',true, ...
       'ignoreStochModel',opt.ignoreStochModel,'refsystem',opt.refsystem,'t0',stmOut.temporalRef,...
       'pntMask',uniqueMask);
  end
  
  % Note: the final output of stmvelocity is based on the total 'pntMask', created within the function.
  % This is the combination of the 'pntMask' given as input to the function, and masks created based
  % on additional criteria, such as ROI and the minimum number of displacement epochs ('minObsDisplPoint').
  % Hence, for the stmprediction and stmrestore step, we use emat and vel, respectively, to get the 
  % final mask.
  
  %% Prediction
  [predSignal,predSignalErrorVar,QzyQyyinvAFull] = stmprediction(stmIn,stmOut,emat,qy,Atrend, ...
         'convThres',opt.predConvThres,'maxIter',opt.predMaxIter);

  %% Restore
  [predTrend,predTrendErrorVar,predTotal,predTotalErrorVar,vel_int] = ...
                    stmrestore(stmIn,stmOut,vel,qx,predSignal,predSignalErrorVar,QzyQyyinvAFull,...
                    'interpMethod',opt.interpMethod,'idwPower',opt.idwPower,'t0',stmOut.temporalRef);

  %% Plot
  fprintf('\n');
  fprintf('Start plots ...\n');

  h1 = stmplotvelpar(stmIn,vel,t0,offset,ref,omt,emat,qx);
  h2 = stmplotvelcov(stmIn,vel,qx,qe,qy);
  h3 = stmplotvelmap(stmIn,vel,'background',opt.background,'defoborder',...
         [stmOut.spatialRef.polyWGS84Vector.Lat' stmOut.spatialRef.polyWGS84Vector.Lon']);
         
  h4 = stmplotpred(stmIn,stmOut,emat,predSignal,predSignalErrorVar,'plotmode',opt.plotmode, ...
         'background',opt.background,'defoborder',...
         [stmOut.spatialRef.polyWGS84Vector.Lat' stmOut.spatialRef.polyWGS84Vector.Lon']);

  h5 = stmplotrestore(stmIn,stmOut,vel,vel_int,predTrend,predSignal,predTotal, ...
         'plotmode',opt.plotmode, 'background',opt.background,'defoborder',...
         [stmOut.spatialRef.polyWGS84Vector.Lat' stmOut.spatialRef.polyWGS84Vector.Lon']);

  h6 = stmplotrestorevar(stmIn,stmOut,vel,vel_int,predTrendErrorVar,...
         predSignalErrorVar,predTotalErrorVar,'plotmode',opt.plotmode, ...
         'background',opt.background,'defoborder',...
         [stmOut.spatialRef.polyWGS84Vector.Lat' stmOut.spatialRef.polyWGS84Vector.Lon']);

  % Save plots to pdf (optional)
  if ~isempty(opt.saveplt) 
     if ~exist(opt.saveplt,'dir')
         mkdir(opt.saveplt)
     end
     [~,stmbase]=fileparts(outputFilename);
     pltfile=fullfile(opt.saveplt,[ stmbase '.pdf']);
     exportgraphics(h1, pltfile, 'ContentType', 'vector');
     exportgraphics(h2, pltfile, 'ContentType', 'vector','Append',true);
     for v = 1:numel(h3)
        exportgraphics(h3(v), pltfile, 'ContentType', 'vector','Append',true);
     end 
     for v = 1:numel(h4)
        exportgraphics(h4(v), pltfile, 'ContentType', 'vector','Append',true);
     end 
     for v = 1:numel(h5)
        exportgraphics(h5(v), pltfile, 'ContentType', 'vector','Append',true);
     end 
     for v = 1:numel(h6)
        exportgraphics(h6(v), pltfile, 'ContentType', 'vector','Append',true);
     end 
  
  end

  % return plot handles
  h=[h1 h2 h3 h4 h5 h6];
                              
end

%% Fill STM

% Point name
stmOut.pntName=cellstr(olcencode(stmOut.pntCrd(:,1),stmOut.pntCrd(:,2),3));

% Observation data
stmOut.obsData = predTotal;

stmOut.sensitivityMatrix = NaN(stmOut.numPoints,3,3);
stmOut.sensitivityMatrix(:,:,1) = [ ones(stmOut.numPoints,1) zeros(stmOut.numPoints,1) zeros(stmOut.numPoints,1) ];
stmOut.sensitivityMatrix(:,:,2) = [ zeros(stmOut.numPoints,1) ones(stmOut.numPoints,1) zeros(stmOut.numPoints,1) ];
stmOut.sensitivityMatrix(:,:,3) = [ zeros(stmOut.numPoints,1) zeros(stmOut.numPoints,1) ones(stmOut.numPoints,1) ];     

% Stochastic model
% Only the variances of predicted values are stored, in the same format as obsData, hence
% different compared to other stochData

stmOut.stochModel{1} = {'Variances'};
stmOut.stochData = predTotalErrorVar;

% Auxiliary data (predicted signal and trend separately, including variances)
stmOut.auxData = cat(3,predSignal,predTrend,predSignalErrorVar,predTrendErrorVar);
stmOut.auxTypes = {'Signal-North' 'Signal-East' 'Signal-Up' ...
                   'Trend-North' 'Trend-East' 'Trend-Up' ...
                   'SignalVar-North' 'SignalVar-East' 'SignalVar-Up' ...
                   'TrendVar-North' 'TrendVar-East' 'TrendVar-Up'};

% Dataset history
inputDatasets = stmOut.inputDatasets;
inputDatasets(1).datasetId = stmIn.datasetId;
inputDatasets(1).techniqueId = stmIn.techniqueId;
inputDatasets(1).datasetAttrib = stmIn.datasetAttrib;
inputDatasets(1).numPoints = stmIn.numPoints;
inputDatasets(1).numEpochs = stmIn.numEpochs;
inputDatasets(1).inputDatasets = stmIn.inputDatasets;
stmOut.inputDatasets = inputDatasets;


%% Write STM
stmwrite(stmOut,outputFilename);


%% Finish the function

fprintf('\n');
fprintf('%s finished at %s  (elapsed time %.2f s)\n',progname,datestr(now),toc)
diary OFF

%catch ME
%
%   getReport(ME)
%   fprintf('%s ABORTED with an error at %s (elapsed time %.2f s)\n',progname,datestr(now),toc)
%   diary OFF
%   rethrow(ME)
%    
%end

end


%% Sub-functions

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
   plot(sqrt(diag(qy{i})),'b.')
   hold on
   plot(sqrt(diag(qe{i})),'r.')
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


function h=stmplotpred(stmIn,stmOut,emat,predSignal,predSignalErrorVar,varargin);
%STMPLOTPRED   Map with signal predictions in space time matrix.
%  STMPLOTPRED(STMIN,STMOUT,EMAT,PREDSIGNAL,PREDSIGNALERRORVAR) plots 
%  1) the original residuals wrt the trend EMAT at the observation points and 
%  epochs (first column in plot), 2) the predicted signal PREDSIGNAL at the 
%  evaluation points and epochs (second column in plot), and 3) the associated
%  prediction error variance PREDSIGNALERRORVAR (third column in plot). A 
%  separate series of plots is created for each observation type (rows). 
%
%  H=STMPLOTPRED(STMIN,STMOUT,EMAT,PREDSIGNAL,PREDSIGNALERRORVAR,'option',value,...)
%  returns the plot handle(s) H and allows to the following options
%
%   'background'     matfile or lat/lon with coastlines, lakes, borders 
%                    (default 'igp-background.mat') 
%   'defoborder'     matfile or lat/lon with deformation border (default 
%                    'igp-unstablearea.mat')
%   'width'          width of plot in normalized coordinates (default 0.8),
%                    the height is computed respecting the aspect ratio
%   'plotmode'       the plotmode: 'last' (default) for only the last predicted
%                    epoch or 'timeseries' for a the full timeseries
%
%  See also STMPLOTRESTORE, STMPLOTRESTOREVAR.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:   12 Dec 2023 by Freek van Leijen
% Modified: 

% Check input arguments and process options

if nargin < 5
  error('This function expects at least five input arguments.')
end

opt.background='igp-background.mat';   % matfile or lat/lon with map background (coastlines, lakes, borders) 
opt.defoborder='igp-unstablearea.mat'; % matfile or lat/lon with deformation border
opt.width=0.8;                         % width of plot in normalized coordinates
opt.plotmode = 'last';                 % 'last' (default) or 'timeseries'

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
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

pntMask=logical(sum(any(~isnan(predSignal),2),3));

width=opt.width;
height=   ( max(stmOut.pntCrd(pntMask,1)) - min(stmOut.pntCrd(pntMask,1)) ) / ...
        ( ( max(stmOut.pntCrd(pntMask,2)) - min(stmOut.pntCrd(pntMask,2)) ) * ...
        cosd(mean(stmOut.pntCrd(pntMask,1))) ) * width * 16/9;
if height > width
    tmp=height/width;
    height=height/tmp;
    width=width/tmp;
end

% Prepare observation types and plotting of names
isnorth=ismember(stmOut.obsTypes,{'North'});
iseast=ismember(stmOut.obsTypes,{'East'});

showpointnames = sum(any(~isnan(emat),2)) <= 50;

numObsTypes=numel(stmOut.obsTypes);

% Set plotmode
switch opt.plotmode
  case 'last'
    plotEpochs = numel(stmOut.epochDyear);
  case 'timeseries'
    plotEpochs = 1:numel(stmOut.epochDyear);
  otherwise
    error('You specified an unsupported plot mode.')
end

k = 0;    
for v = plotEpochs

    k = k+1;
    h(k) = figure('units','normalized','outerposition',[0.05 0.05 .95 0.05+0.3.*numObsTypes]);
   
    %tcl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
    tcl = tiledlayout(numObsTypes,3,'TileSpacing','compact','Padding','compact');

    [~,closestEpoch] = min(abs(stmIn.epochDyear-stmOut.epochDyear(v)));
    
    for l=1:numObsTypes

        %plot last epochs of stmIn and stmOut
        nexttile;
        scatter(stmIn.pntCrd(:,2),stmIn.pntCrd(:,1),30,emat(:,closestEpoch,l),'filled');
        caxis([min(min(emat(:,:,l))) max(max(emat(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmIn.obsTypes{l} ' Residuals closest observation epoch [mm]']);
        if showpointnames
            text(stmIn.pntCrd(:,2),stmIn.pntCrd(:,1),stmIn.pntName,'FontSize',7,'VerticalAlignment','top');
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
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');        
        
        nexttile;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),[],predSignal(:,v,l),'filled');
        caxis([min(min(predSignal(:,:,l))) max(max(predSignal(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmOut.obsTypes{l} ' Predicted signal [mm]']);
        hold on; 
        xl=xlim();yl=ylim(); 
        if ~isempty(igpbackground)
          plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
        end
        if ~isempty(defoborder)
          plot(defoborder(:,2),defoborder(:,1),':m')
        end
        xlim(xl);ylim(yl);
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');

        nexttile;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),[],predSignalErrorVar(:,v,l),'filled');
        caxis([min(min(predSignalErrorVar(:,:,l))) max(max(predSignalErrorVar(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmOut.obsTypes{l} ' Predicted signal variance [mm^2]']);
        hold on; 
        xl=xlim();yl=ylim(); 
        if ~isempty(igpbackground)
          plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
        end
        if ~isempty(defoborder)
          plot(defoborder(:,2),defoborder(:,1),':m')
        end
        xlim(xl);ylim(yl);
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');

    end
        
end

end

function h=stmplotrestore(stmIn,stmOut,vel,vel_int,predTrend,predSignal,predTotal,varargin);
%STMPLOTRESTORE   Map with trend, signal and total predictions in space time matrix.
%  STMPLOTRESTORE(STMIN,STMOUT,VEL,VEL_INT,PREDTREND,PREDSIGNAL,PREDTOTAL) plots 
%  1) the velocities (trend) VEL at the observation points and the velocities 
%  VEL_INT at the prediction points (first column in plot), 2) the predicted trend
%  PREDTREND at the evaluation points and epochs (second column in plot), the 
%  predicted signal PREDSIGNAL at the evaluation points and epochs (third column
%  in plot)and 4) the total predictions PREDTOTAL (forth column in plot). A 
%  separate series of plots is created for each observation type (rows). 
%
%  H=STMPLOTPRED(STMIN,STMOUT,VEL,VEL_INT,PREDTREND,PREDSIGNAL,PREDTOTAL,'option',value,...)
%  returns the plot handle(s) H and allows to the following options
%
%   'background'     matfile or lat/lon with coastlines, lakes, borders 
%                    (default 'igp-background.mat') 
%   'defoborder'     matfile or lat/lon with deformation border (default 
%                    'igp-unstablearea.mat')
%   'width'          width of plot in normalized coordinates (default 0.8),
%                    the height is computed respecting the aspect ratio
%   'plotmode'       the plotmode: 'last' (default) for only the last predicted
%                    epoch or 'timeseries' for a the full timeseries
%
%  See also STMPLOTPRED, STMPLOTRESTOREVAR.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:   12 Dec 2023 by Freek van Leijen
% Modified: 

% Check input arguments and process options

if nargin < 7
  error('This function expects at least seven input arguments.')
end

opt.background='igp-background.mat';   % matfile or lat/lon with map background (coastlines, lakes, borders) 
opt.defoborder='igp-unstablearea.mat'; % matfile or lat/lon with deformation border
opt.width=0.8;                         % width of plot in normalized coordinates
opt.plotmode = 'last';                 % 'last' (default) or 'timeseries'

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
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

pntMask=logical(sum(any(~isnan(predSignal),2),3));

width=opt.width;
height=   ( max(stmOut.pntCrd(pntMask,1)) - min(stmOut.pntCrd(pntMask,1)) ) / ...
        ( ( max(stmOut.pntCrd(pntMask,2)) - min(stmOut.pntCrd(pntMask,2)) ) * ...
        cosd(mean(stmOut.pntCrd(pntMask,1))) ) * width * 16/9;
if height > width
    tmp=height/width;
    height=height/tmp;
    width=width/tmp;
end

% Prepare observation types and plotting of names
isnorth=ismember(stmOut.obsTypes,{'North'});
iseast=ismember(stmOut.obsTypes,{'East'});

showpointnames = sum(any(~isnan(vel),2)) <= 50;

numObsTypes=numel(stmOut.obsTypes);

% Set plotmode
switch opt.plotmode
  case 'last'
    plotEpochs = numel(stmOut.epochDyear);
  case 'timeseries'
    plotEpochs = 1:numel(stmOut.epochDyear);
  otherwise
    error('You specified an unsupported plot mode.')
end

k = 0;    
for v = plotEpochs

    k = k+1;
    h(k) = figure('units','normalized','outerposition',[0.05 0.05 .95 0.05+0.3.*numObsTypes]);
   
    %tcl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
    tcl = tiledlayout(numObsTypes,4,'TileSpacing','compact','Padding','compact');

    for l=1:numObsTypes

        nexttile;
        scatter(stmIn.pntCrd(:,2),stmIn.pntCrd(:,1),30,vel(:,l),'filled');
        hold on;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),10,vel_int(:,l),'filled');
        hc=colorbar;
        ylabel(hc,[ stmIn.obsTypes{l} ' Velocity [mm/y]']);
        if showpointnames
            text(stmIn.pntCrd(:,2),stmIn.pntCrd(:,1),stmIn.pntName,'FontSize',7,'VerticalAlignment','top');
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
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');        
        
        nexttile;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),[],predTrend(:,v,l),'filled');
        caxis([min(min(predTrend(:,:,l))) max(max(predTrend(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmOut.obsTypes{l} ' Predicted trend [mm]']);
        hold on; 
        xl=xlim();yl=ylim(); 
        if ~isempty(igpbackground)
          plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
        end
        if ~isempty(defoborder)
          plot(defoborder(:,2),defoborder(:,1),':m')
        end
        xlim(xl);ylim(yl);
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');

        nexttile;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),[],predSignal(:,v,l),'filled');
        caxis([min(min(predSignal(:,:,l))) max(max(predSignal(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmOut.obsTypes{l} ' Predicted signal [mm]']);
        hold on; 
        xl=xlim();yl=ylim(); 
        if ~isempty(igpbackground)
          plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
        end
        if ~isempty(defoborder)
          plot(defoborder(:,2),defoborder(:,1),':m')
        end
        xlim(xl);ylim(yl);
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');

        nexttile;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),[],predTotal(:,v,l),'filled');
        caxis([min(min(predTotal(:,:,l))) max(max(predTotal(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmOut.obsTypes{l} ' Predicted total [mm]']);
        hold on; 
        xl=xlim();yl=ylim(); 
        if ~isempty(igpbackground)
          plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
        end
        if ~isempty(defoborder)
          plot(defoborder(:,2),defoborder(:,1),':m')
        end
        xlim(xl);ylim(yl);
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');

    end
        
end

end

function h=stmplotrestorevar(stmIn,stmOut,vel,vel_int,predTrendErrorVar,...
predSignalErrorVar,predTotalErrorVar,varargin);
%STMPLOTRESTOREVAR   Map with trend, signal and total prediction error variances
%  in space time matrix.
%  STMPLOTRESTOREVAR(STMIN,STMOUT,VEL,VEL_INT,PREDTRENDERRORVAR,PREDSIGNALERRORVAR,
%  PREDTOTALERRORVAR) plots 1) the velocities (trend) VEL at the observation points 
%  and the velocities VEL_INT at the prediction points (first column in plot), 2) 
%  the predicted trend error variances PREDTRENDERRORVAR at the evaluation points 
%  and epochs (second column in plot), the predicted signal error variances 
%  PREDSIGNALERRORVAR at the evaluation points and epochs (third column in plot) 
%  and 4) the total prediction error variances PREDTOTALERRORVAR (forth column in
%  plot). A separate series of plots is created for each observation type (rows). 
%
%  H=STMPLOTPRED(STMIN,STMOUT,VEL,VEL_INT,PREDTRENDERRORVAR,PREDSIGNALERRORVAR, ...
%  PREDTOTALERRORVAR,'option',value,...) returns the plot handle(s) H and allows 
%  to the following options
%
%   'background'     matfile or lat/lon with coastlines, lakes, borders 
%                    (default 'igp-background.mat') 
%   'defoborder'     matfile or lat/lon with deformation border (default 
%                    'igp-unstablearea.mat')
%   'width'          width of plot in normalized coordinates (default 0.8),
%                    the height is computed respecting the aspect ratio
%   'plotmode'       the plotmode: 'last' (default) for only the last predicted
%                    epoch or 'timeseries' for a the full timeseries
%
%  See also STMPLOTPRED, STMPLOTRESTORE.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:   12 Dec 2023 by Freek van Leijen
% Modified: 

% Check input arguments and process options

if nargin < 7
  error('This function expects at least seven input arguments.')
end

opt.background='igp-background.mat';   % matfile or lat/lon with map background (coastlines, lakes, borders) 
opt.defoborder='igp-unstablearea.mat'; % matfile or lat/lon with deformation border
opt.width=0.8;                         % width of plot in normalized coordinates
opt.plotmode = 'last';                 % 'last' (default) or 'timeseries'

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
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

pntMask=logical(sum(any(~isnan(predSignalErrorVar),2),3));

width=opt.width;
height=   ( max(stmOut.pntCrd(pntMask,1)) - min(stmOut.pntCrd(pntMask,1)) ) / ...
        ( ( max(stmOut.pntCrd(pntMask,2)) - min(stmOut.pntCrd(pntMask,2)) ) * ...
        cosd(mean(stmOut.pntCrd(pntMask,1))) ) * width * 16/9;
if height > width
    tmp=height/width;
    height=height/tmp;
    width=width/tmp;
end

% Prepare observation types and plotting of names
isnorth=ismember(stmOut.obsTypes,{'North'});
iseast=ismember(stmOut.obsTypes,{'East'});

showpointnames = sum(any(~isnan(vel),2)) <= 50;

numObsTypes=numel(stmOut.obsTypes);

% Set plotmode
switch opt.plotmode
  case 'last'
    plotEpochs = numel(stmOut.epochDyear);
  case 'timeseries'
    plotEpochs = 1:numel(stmOut.epochDyear);
  otherwise
    error('You specified an unsupported plot mode.')
end

k = 0;    
for v = plotEpochs

    k = k+1;
    h(k) = figure('units','normalized','outerposition',[0.05 0.05 .95 0.05+0.3.*numObsTypes]);
   
    %tcl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
    tcl = tiledlayout(numObsTypes,4,'TileSpacing','compact','Padding','compact');

    for l=1:numObsTypes

        nexttile;
        scatter(stmIn.pntCrd(:,2),stmIn.pntCrd(:,1),30,vel(:,l),'filled');
        hold on;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),10,vel_int(:,l),'filled');
        hc=colorbar;
        ylabel(hc,[ stmIn.obsTypes{l} ' Velocity [mm/y]']);
        if showpointnames
            text(stmIn.pntCrd(:,2),stmIn.pntCrd(:,1),stmIn.pntName,'FontSize',7,'VerticalAlignment','top');
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
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');        
        
        nexttile;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),[],predTrendErrorVar(:,v,l),'filled');
        caxis([min(min(predTrendErrorVar(:,:,l))) max(max(predTrendErrorVar(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmOut.obsTypes{l} ' Predicted trend error variance [mm^2]']);
        hold on; 
        xl=xlim();yl=ylim(); 
        if ~isempty(igpbackground)
          plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
        end
        if ~isempty(defoborder)
          plot(defoborder(:,2),defoborder(:,1),':m')
        end
        xlim(xl);ylim(yl);
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');

        nexttile;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),[],predSignalErrorVar(:,v,l),'filled');
        caxis([min(min(predSignalErrorVar(:,:,l))) max(max(predSignalErrorVar(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmOut.obsTypes{l} ' Predicted signal error variance [mm^2]']);
        hold on; 
        xl=xlim();yl=ylim(); 
        if ~isempty(igpbackground)
          plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
        end
        if ~isempty(defoborder)
          plot(defoborder(:,2),defoborder(:,1),':m')
        end
        xlim(xl);ylim(yl);
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');

        nexttile;
        scatter(stmOut.pntCrd(:,2),stmOut.pntCrd(:,1),[],predTotalErrorVar(:,v,l),'filled');
        caxis([min(min(predTotalErrorVar(:,:,l))) max(max(predTotalErrorVar(:,:,l)))]);
        hc=colorbar;
        ylabel(hc,[ stmOut.obsTypes{l} ' Predicted total error variance [mm^2]']);
        hold on; 
        xl=xlim();yl=ylim(); 
        if ~isempty(igpbackground)
          plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
        end
        if ~isempty(defoborder)
          plot(defoborder(:,2),defoborder(:,1),':m')
        end
        xlim(xl);ylim(yl);
        daspect([1/cosd(mean(yl)) 1 1]);
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');

        title([stmOut.obsTypes{l} ' - ' num2str(stmOut.temporalRef,'%.2f') '-' ...
                   num2str(stmOut.epochDyear(v),'%.2f')],'interpreter','none');

    end
        
end

end

