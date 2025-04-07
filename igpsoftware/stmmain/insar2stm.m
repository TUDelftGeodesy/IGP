function varargout=insar2stm(inputfilenames,outputfilename,varargin)
%insar2stm   Import InSAR data files into space time matrix dataset
%   INSAR2STM(INPUTFILENAME,OUTPUTFILENAME,OPTIONS) imports InSAR data file
%   given by INPUTFILENAME into a single Space Time Matrix dataset with 
%   the name OUTPUTFILENAME. INPUTFILENAME must be the name of an file 
%   containing the input filename, a character string with wildcards, 
%   or a cell array with the input filenames. OUTPUTFILENAME is a character 
%   string with the name of the output Space Time Matrix dataset. OPTIONS 
%   is a cell array or structure with the processing options, if empty, 
%   or when fields are missing, default values are used.
%
%   INSAR2STM(...,UPDATEFLAG) affects the processing in case OUTPUTFILENAME
%   exists. The following values of UPDATEFLAG are accepted
%
%     'create'     abort processing if the output dataset exists
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
%
%   STAT=GNSS2STM(...) returns a status code STAT. A status code of 0
%   indicates success.
%
%   OPTIONS is a cell array or structure with the processing options, 
%   if empty, or when fields are missing, default values are used. Valid 
%   OPTIONS are
%
%      verbose=0                Verbosity level, higher is more output, 0 is almost nothing
%
%      inputDir=''              Directory with the inputfiles
%
%      terrainFilename=''       TerrainFilename (Default none)
%
%      inputFormat='Shell'      Input format (Default 'Shell' format [1])
%      headingAngle=<float>     Heading angle in [degrees] (MANDATORY for Shell)
%      system=<string>          System identifier (MANDATORY for Shell)
%      systemMode=<string>      System mode ('sm') (MANDATORY for Shell)
%      mode=''                  Ascending (Asc) or Decending (Dsc) mode
%      polarization=''          Polarization
%
%      stochModelParameterFile  Name of the file with stochastic model parameters
%                               (Default 'insarStochasticModelParameters.txt')
%
%      highLowThres=3           Threshold for highLow classification [m] (Default 3 m)
%      classificationOrder={'highLow'} Classification order {'highLow','deep','total'} (Default {'highLow'})
%
%      ROI=[]                   Region of interest, as [latmin lonmin ; latmax lonmax] bounding 
%                               box, or lat/lon polygon, or kml/shape file(Default all) [2]
%      POI=[-Inf +Inf]          Period of interest [ dYearStart dYearEnd ] or 
%                               [ dYearStart dYearEnd ; dYearStart dYearEnd ;...  ]  (Default all) [2]
%
%      globalAttrib=[]          Struct with global attributes for output dataset (Empty by default) 
%
%   [1] InSAR Deformation Monitoring, Shell Geomatics Technical
%   Specification. Document number: ON-TS-105. Revision: 06, 27 Sep 2017.
%   [2] ROI and POI are NOT implemented at the moment...
%
%   Some of the options are mandatory for the Shell format as their values
%   cannot be deduced from the file itself. 
%
%   Within the function also the stochastic model is inserted in the Space
%   Time Matrix structure. This model is based on the library of stochastic
%   models specified in the file: insarStochasticModelParameters.txt . This
%   file should exist in the PATH, and the stochastic model for the various
%   InSAR datasets used should have been specified in this file.

% 
%   Examples:
%      insar2stm('tsx/*.csv','tsx.mat',options)      
%      insar2stm('sentinel1/*.csv','sentinel1.mat',options,'update')      
%
%   See also STM, STMCHECKARGUMENTS, STMREAD and STMDIFF.
%
%  (c) Freek van Leijen, Delft University of Technology, 2020, 2021. 

% Created:  18 August 2020 by Freek van Leijen
% Modified: 24 August 2020 by Hans van der Marel
%              - New call syntax
%              - Adapted to the new datastructure format and STM functions
%              - Removed loop processing
%           25 Sep 2020 by Freek van Leijen
%              - pntCrd in double precison
%              - implementation stochastic model
%           09 Oct 2020 by Freek van Leijen
%              - made attribute names compatible with stm.m
%           21 Nov 2020 by Freek van Leijen
%              - inserted log
%           27 Dec 2021 by Hans van der Marel
%              - use the same syntax as other modules
%              - dataset history and global attributes
%              - set highLow and deformationRegime attributes
%              - possibility to force deformationRegime to [deep|total]
%              - save RD coordinates as mapCrd and set CRS
%              - various new options (e.g. for stochasticModelParams)
%           06 Jan 2024 by Freek van Leijen
%              - added comments on the classificationOrder options

% TO DO:
%              - implement ROI and POI (use getPointMask and getEpochMask)

%% Check the input arguments and options

if nargin < 2
   error('This function expects at least two input arguments.')
end

progname='insar2stm';

% Default options

opt.verbose=0;                                  % Default verbosity level, higher is more output, 0 is almost nothing

opt.inputDir='';                                % Directory with the inputfiles
opt.inputFormat='Shell';                        % Input format

opt.terrainFilename='';                         % TerrainFilename (Default none)
opt.headingAngle='';                            % Heading angle in [degrees] (MANDATORY for Shell)
opt.system='';                                  % System identifier (MANDATORY for Shell)
opt.systemMode='';                              % System mode ('sm') (MANDATORY for Shell)
opt.mode='';                                    % Ascending (Asc) or Decending (Dsc) mode
opt.polarization='';                            % Polarization

opt.stochModelParameterFile='insarStochasticModelParameters.txt';  % Name of the file with stochastic model parameters

opt.highLowThres=3;                             % High-low threshold [m] for high-low classification (Default 3 m)
opt.classificationOrder={'highLow'};            % Classification order {'highLow','deep','total'} (Default {'highLow'})

%opt.ROI=[];                                     % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon, kml/shape file (Default none)
%opt.POI=[-Inf +Inf];                            % Period of interest [ dYearStart dYearEnd ] or [ dYearStart dYearEnd ; dYearStart dYearEnd ; ... ] (Default none)

opt.globalAttrib=struct([]);                    % Struct, if empty (default) globalAttributes in stm is empty

opt.projectId='';                               % Project Id (default directory of the outputfile, or empty if no directory is specified)

% Duplicate output to file, and catch errors, and start timing

try

[~,outputfileroot]=fileparts(outputfilename);
diary([ outputfileroot '_' datestr(now,30) '.log' ])

fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Check the options and if necessary overwrite the default values

[inputfilenames,outputfilename,opt]= ...
    stmcheckarguments(inputfilenames,outputfilename,opt,varargin{:});
if isempty(outputfilename)
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    return;
end

% Check input file format and mandatory parameters 

switch lower(opt.inputFormat)
    case 'shell'
        % Check
        if isempty(opt.headingAngle); error('Heading angle is not specified properly'); end
        if isempty(opt.system); error('System is not specified properly'); end
        if isempty(opt.systemMode); error('System mode is not specified properly'); end
        if isempty(opt.mode)
            if mod(opt.headingAngle-90,360) > 180
               opt.mode='Asc';
            else
               opt.mode='Dsc';
            end
        end
        if numel(inputfilenames) > 1
            error('Shell format supports only single input file')
        end
        inputFilename=char(fullfile(opt.inputDir,inputfilenames{1}));
        if ~exist(inputFilename,'file')
            fprintf('Shell input file %s not found\n',inputFilename)
            fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
            return;
        end
    otherwise
        fprintf('Unsupported input file format %s\n',opt.inputFormat)
        fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
        return;
end

% Check system and systemMode options

switch opt.system
    case {'RadarSAT-2','radarsat-2','Rsat-2','rsat-2','rsat2'}
        opt.system = 'RadarSAT-2';
    case {'TerraSAR-X','terrasar-x','TSX','tsx'}
        opt.system = 'TerraSAR';
    case {'Sentinel-1','sentinel-1','S1','s1'}
        opt.system = 'Sentinel-1';
    otherwise
        error('You specified an unsupported system.');
end

switch opt.systemMode
    case {'Standard','standard','S3','s3'}
        opt.systemMode = 'Standard';
    case {'XF','xf'}
        opt.systemMode = 'XF';
    case {'Stripmap','stripmap','SM','sm'}
        opt.systemMode = 'Stripmap';
    case {'IW','iw'}
        opt.systemMode = 'IW';
    otherwise
        error('You specified an unsupported system.');
end

% Get InSAR stochastic model parameters

if ~exist(opt.stochModelParameterFile,'file')
    error('File with stochastic model parameters does not exist.');
end    

fid = fopen(opt.stochModelParameterFile,'r');
modelParam = textscan(fid,['%s%s%s%f32%f32%f32%f32%f32'], ...
     'Delimiter',',','Headerlines',1);
fclose(fid);

idx1 = strmatch(opt.system,modelParam{1});
idx2 = strmatch(opt.systemMode,modelParam{2}(idx1));
idx = idx1(idx2);
if isempty(idx)
    error('The stochastic model for this InSAR dataset is not specified.');
else
    stochModelParam.model = modelParam{3}{idx};
    stochModelParam.s20 = modelParam{4}(idx);
    stochModelParam.s2t = modelParam{5}(idx);
    stochModelParam.s2s = modelParam{6}(idx);
    stochModelParam.Rt = modelParam{7}(idx);
    stochModelParam.Rs = modelParam{8}(idx);
end

%% Build the output space time matrix structure

[projectId,datasetId,~]=fileparts(outputfilename);
if numel(opt.projectId) > 1
   projectId=opt.projectId;
end
techniqueId='insar';

dataset = stm(projectId,datasetId,techniqueId,'STRUCT',outputfilename);

% Technique attributes

techniqueAttrib = [];
techniqueAttrib.system = opt.system;
techniqueAttrib.systemMode = opt.systemMode;
techniqueAttrib.mode = opt.mode;
techniqueAttrib.polarization = opt.polarization;
techniqueAttrib.mapCrs='RD';
dataset.techniqueAttrib = techniqueAttrib;

% Dataset attributes

datasetAttrib=dataset.datasetAttrib;
datasetAttrib.softwareOptions=opt;
dataset.datasetAttrib=datasetAttrib;

dataset.parTypes = {'North' 'East' 'Up'};
dataset.obsTypes = {'los'};
dataset.auxTypes = [];


%% Read InSAR deformation header

disp(['Reading ' inputFilename ' ...']);

fid1 = fopen(inputFilename,'r');
header = textscan(fid1,'%s',1);
header = strsplit(char(header{1}),',')';
fclose(fid1);

dates = char(header(20:end)); % Fixed Shell format
dates = cellstr([dates(:,5:8) dates(:,3:4) dates(:,1:2)]);
numEpochs = size(dates,1);
epochDyear = date2dyear(dates,'yyyymmdd');
epochId = cellstr(num2str(epochDyear,'%6.2f'));

dataset.numEpochs = numEpochs;
dataset.epochDyear = epochDyear';

epochAttrib = [];
epochAttrib.epochId = epochId';
dataset.epochAttrib = epochAttrib;

%% Read InSAR deformation data

fid1 = fopen(inputFilename,'r');
data = textscan(fid1,['%s' repmat('%f32',1,4) '%f64%f64' repmat('%f32',1,5) ...
     '%f64' repmat('%f32',1,6) repmat('%*f32',1,numEpochs)], ...
     'Delimiter',',','Headerlines',1,'CollectOutput',1);
fclose(fid1);

numPoints = size(data{1},1);

dataset.numPoints = numPoints;
dataset.pntName = data{1}(:,1); % Unique identifier for each measurement point.

hEll = nap2etrs(data{3}(:,1),data{3}(:,2),data{5}(:,1)); % NAP or EGM86 to Ellipsoid, to be improved!
dataset.pntCrd = [data{3}(:,1) data{3}(:,2) hEll]; % Lat [Deg], Lon [Deg], h [m]
clear hEll

% Store pntAttrib first in temporary array, before moving it to dataset (objects only support one level of indexing)
pntAttrib = [];
pntAttrib.pntClass = data{2}(:,1); % REF(1)=reference point/point inside reference area; 
                                   % PS(2)=independent permanent scatterer;
                                   % SL(3)=Sidelobe; 
                                   % DS(4)=distributed scatterer; 
                                   % CR(5)=corner reflector; 
                                   % (0)=other.
pntAttrib.area = data{2}(:,2); % [m2] effective area covered by the scatterer (for DS, 0 for PS and CR).
pntAttrib.az = data{2}(:,3); % azimuth (along-track) sub-pixel coordinate [pixels] in the original SLC resolution stack geometry 
pntAttrib.rg = data{2}(:,4); % slant-range (across-track) sub-pixel coordinate [pixels] in the original SLC resolution stack geometry
pntAttrib.x = data{4}(:,1); % Easting [m]
pntAttrib.y = data{4}(:,2); % Northing [m]
pntAttrib.mapCrd = [ data{4}(:,1) data{4}(:,2)];  % New style map coordinates (can replace x and y)
pntAttrib.incAngle = data{4}(:,3); % incidence angle [degrees] using the ellipsoidal model of the Earth
pntAttrib.azAngle = repmat(single(opt.headingAngle-90),numPoints,1); % azimuth angle [degrees] using the ellipsoidal model of the Earth
pntAttrib.coh = data{4}(:,4); % ensemble phase coherence
pntAttrib.hDem = data{4}(:,5); % DEM height [m] (if topographic phase contribution is subtracted in PSI processing)
pntAttrib.sigH = data{6}(:,1); % precision absolute height (height in pntCRD) [m]
pntAttrib.lin = data{6}(:,2); % average deformation rate in line-of-sight (full time series) [mm/y]
pntAttrib.sigLin = data{6}(:,3); % precision average deformation rate in line-of-sight (full time series) [mm/y]
pntAttrib.linYear = data{6}(:,4); % average deformation rate in line-of-sight (last year) [mm/y]
pntAttrib.sigLinYear = data{6}(:,5); % precision average deformation rate in line-of-sight (last year) [mm/y]
pntAttrib.sigDefo = data{6}(:,6); % precision deformation along line-of-sight (standard deviation single deformation estimates assuming the applied deformation model) [mm]

clear data

%% InSAR terrain file, high-low classification and deformation regime

% initialize highLow Classification and deformation regime

pntAttrib.highLow = zeros(numPoints,1,'int16'); %0=unclassified, 1=high, 2=low
pntAttrib.defoRegime = zeros(numPoints,1,'int16'); %0=unclassified, 1=deep, 2=total

% Read and save InSAR Terrain .csv file

if ~isempty(char(opt.terrainFilename))
   fid2 = fopen(fullfile(opt.inputDir,char(opt.terrainFilename)),'r');
   data2 = textscan(fid2,'%s%f32','Delimiter',',','Headerlines',1, ...
        'CollectOutput',1);
   fclose(fid2);
   
   pntName2 = data2{1}(:,1);
   hTerrain = data2{2}(:,1);
   hTerrain(hTerrain==-327.7) = NaN; % Fixed Shell format
   clear data2

   % Get index (to account for missing values)
   [~,idx,~] = intersect(dataset.pntName,pntName2,'stable');

   % Get AHN terrain height wrt ellipsoid [m]
   hTerrainEll = nap2etrs(dataset.pntCrd(idx,1),dataset.pntCrd(idx,2),hTerrain);

   pntAttrib.hTerrain = NaN(dataset.numPoints,1);
   pntAttrib.hTerrain(idx) = hTerrainEll;
   
   % highLow classification
   
   dHeight = dataset.pntCrd(:,3) - pntAttrib.hTerrain; % heights wrt ellipsoid

   pntAttrib.highLow(dHeight>=opt.highLowThres) = 1;
   pntAttrib.highLow(dHeight<opt.highLowThres) = 2;
        
   clear hTerrain hTerrainEll pntName2 idx dHeight
end

% A-priori deformation regime %0=unclassified, 1=deep, 2=total

for k = numel(opt.classificationOrder):-1:1 % start with last option, overwrite when better method available
    switch opt.classificationOrder{k}
        case 'highLow'
            pntAttrib.defoRegime(pntAttrib.highLow==1) = 1; % High points -> sensitive to deep deformation only
            pntAttrib.defoRegime(pntAttrib.highLow==2) = 2; % Low points -> sentitive to total deformation
        case 'deep'
            pntAttrib.defoRegime = ones(numPoints,1); % If all points are sensitive to deep deformation, e.g., IGRS
        case 'total'
            pntAttrib.defoRegime = repmat(2,numPoints,1); % If all points are sensitive to total deformation, e.g., shallow founded corner reflector or transponder
        case {'defoModels','histogram','auxiliary'}           
            warning('You specified an advanced classification method which still needs to be implemented.');
        otherwise
            error('You specified an unsupported classification method.');
    end
end

% Save point attributes to dataset

dataset.pntAttrib = pntAttrib;


%% Calculate sensitivity

dataset.sensitivityMatrix = [-sind(pntAttrib.incAngle)*cosd(opt.headingAngle-270) ... % North
                           -sind(pntAttrib.incAngle)*sind(opt.headingAngle-270) ... % East
                            cosd(pntAttrib.incAngle)]; % Up

%% Read InSAR deformation time series

fid3 = fopen(inputFilename,'r');
dataset.obsData = cell2mat(textscan(fid3,['%*s' repmat('%*f32',1,18) repmat('%f32',1,numEpochs)], ...
                                       'Delimiter',',','Headerlines',1,'CollectOutput',0));
fclose(fid3);


%% Set stochastic model
dataset.stochModel = {[stochModelParam.model ...
                    '(s20=' num2str(stochModelParam.s20) ...
                    ',s2t=' num2str(stochModelParam.s2t) ...
                    ',s2s=' num2str(stochModelParam.s2s) ...
                    ',Rt=' num2str(stochModelParam.Rt) ...
                    ',Rs=' num2str(stochModelParam.Rs) ')']};
dataset.stochData = [];


%% Input dataset history and global attributes

inputDatasets(1).datasetId=datasetId;
inputDatasets(1).techniqueId=techniqueId;
datasetAttrib=[];
datasetAttrib.softwareName='';
datasetAttrib.softwareOptions=[];
datasetAttrib.fileFormat=opt.inputFormat;
datasetAttrib.fileName=inputFilename;
datasetAttrib.projectId='';
inputDatasets(1).datasetAttrib=datasetAttrib;
inputDatasets(1).numPoints=numPoints;
inputDatasets(1).numEpochs=numEpochs;
inputDatasets(1).inputDatasets=[];
dataset.inputDatasets=inputDatasets;

dataset.globalAttrib=opt.globalAttrib;
 

%% Write dataset file

stmwrite(dataset);
clear dataset;

% Finish the function

fprintf('%s finished at %s  (elapsed time %.2f s)\n',progname,datestr(now),toc)
diary OFF

catch ME

   getReport(ME)
   fprintf('%s ABORTED with an error at %s (elapsed time %.2f s)\n',progname,datestr(now),toc)
   diary OFF
   rethrow(ME)
    
end

end

