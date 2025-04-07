function varargout=cupido2stm(inputFilename,outputFilename,varargin)
%CUPIDO2STM    CUPiDO v1 to Space Time Matrix format conversion.
%   CUPIDO2STM(INPUTFILENAME,OUTPUTFILENAME,OPTIONS) converges a CUPiDO v1
%   dataset into a Space Time Matrix. The content can be based on any
%   kind of measurement technique, e.g., levelling, GNSS, ... 
%   INPUTFILENAME is a CUPiDO NetCDF (.nc) file. OUTPUTFILENAME is the
%   STM output filename.  OPTIONS is a cell array or structure with the 
%   processing options, if empty, or when fields are missing, default 
%   values are used.
%
%   CUPIDO2STM(...,UPDATEFLAG) specifies conditional processing in case the
%   output file OUTPUTFILENAME exists. Possible values for UPDATEFLAG are
%
%     'create'     abort processing if the output dataset exists (the default)
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
%
%   SUCCESSFLAG=CUPIDO2STM(...) returns a successflag when the function is
%   run successfully.
% 
%   OPTIONS is a cell array or structure with the processing options, 
%   if empty, or when fields are missing, default values are used. Valid 
%   OPTIONS are
%
%      verbose=0             Verbosity level, higher is more output, 0 is almost nothing
%      inputDir=''           Directory with the inputfile
%      heightOnly=false      Only import the height component to the STM, ignore everything else
%      splitNetwork='warn'   Split network {'warn','delete','keeplargest','split'}; warning in case 
%                            there are multiple networks in an epoch, delete the epoch, keep only 
%                            the largest network, split in separate epochs/campaigns.
%      minObs=2              Minimum number of observations per point
%      ROI=[]                Region of interest, as [latmin lonmin ; latmax lonmax] bounding 
%                            box, or lat/lon polygon, or kml/shape file (Default all)
%      POI=[-Inf +Inf]       Period of interest [ dYearStart dYearEnd ] or 
%                            [ dYearStart dYearEnd ; dYearStart dYearEnd ;...  ]  (Default all)
%      opt.includePrjName={} Cell array with patterns of project names to include (if empty, all 
%                            are excluded), this rule is applied after POI selection.
%      opt.excludePrjName={} Cell array with patterns of project names to exclude (if empty, none 
%                            are excluded), this rule is applied after POI and include rule selection
%      globalAttrib=[]       Struct, if empty (default), take only globalAttributes from Cupido NetCDF
%      projectId=''          Project Id (default directory of the outputfile, or empty if no 
%                            directory is specified)
%
%   [Note: splitNetwork options {'keeplargest','split'} are NOT implemented at the moment...]
%
%   Examples (standard form):
%
%      cupido2stm(infile,outfile,options) 
%      cupido2stm(infile,outfile,options,'update')
%      successflag=cupido2stm(...)
%
%  (c) Freek van Leijen, Hans van der Marel, Delft University of Technology, 2020. 

% Created:  12 Aug 2020 by Freek van Leijen
% Modified: 25 Sep 2020 by Freek van Leijen
%              - pntCrd in double precison
%              - implementation stochastic model
%           21 Nov 2020 by Freek van Leijen
%              - inserted log
%           20 Oct 2021 by Hans van der Marel
%              - include cupido_read_netcdf.m as local function
%              - 'GPS and 'GNSS' both accepted as techniques (req. Hermann)
%              - added options to input arguments for consistency with
%                other functions in stmmain (make it futere proof)
%              - check input/output parsing and final cleanup
%              - actually produce the advertised successflag
%              - added global attributes and file history to output stm
%              - include coordinates in mapprojection (pntCrdProj / projCrs)
%           24 Jan 2022 by Freek van Leijen
%              - fixed bug in conversion covariance matrix 
%           11 Aug 2022 by Hans van der Marel
%              - fixed boolean test for techniques
%           22 Aug 2022 by Hans van der Marel
%              - major update to handle Cupido datasets with contains more
%                than only the height component (ndim gives the dimension) 
%              - added option to supress anything other than heigth  
%              - implemented check whether there are separate networks at
%                single epoch, implement three actions: a) delete the epoch
%                b) keep only largest subnetwork, c) split in separate 
%                epochs/campaigns. This is determined by opt.splitNetwork 
%                = {'warn','delete','keeplargest','split'}. Only action
%                a) is implemented at the moment.
%              - remove synthetic benchmarks from output space time matrix 
%              - enabled selection of POI and ROI
%              - enabled verbose option
%           06 Oct 2022 by Hans van der Marel
%              - added options to include and exclude project names
%              - added output of epoch table
%           03 Nov 2022 by Hans van der Marel
%              - correct units for stochastic model (cupodo m -> stm mm)
%           16 Feb 2023 by Hans van der Marel
%              - fixed bug in cupido_remove_points
%              - remove points which have no observations left (e.g. after
%                selecting epochs), implemented in cupido_remove_points_noobs
%           14 Mar 2023 by Hans van der Marel
%              - added option to include selected integer sdobsflag values
%            6 Oct 2023 by Hans van der Marel
%              - fix scaling bug in output space time matrix
%           11 Apr 2024 by Hans van der Marel
%              - add from points to space time matrix as zero (one per campaign)
%              - added from points to epochAttrib in stm
%              - remove epochs with zero observations
%              - update counts in case splitNetwork='delete'
%              - reindex observation table in case of synthetic Benchmarks
%            6 June 2024 by Hans van der Marel 
%              - remove points with less than opt.minObs(=2) observations

%% Check the input arguments and options

if nargin < 3
    error('This function expects at least two input arguments.')
end

progname='cupido2stm';

% Default options

opt=[];
opt.verbose=0;                   % Default verbosity level, higher is more output, 0 is almost nothing
opt.inputDir='';                 % Directory with the inputfiles
opt.ROI=[];                      % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon, kml/shape file (Default none)
opt.POI=[-Inf +Inf];             % Period of interest [ dYearStart dYearEnd ] or [ dYearStart dYearEnd ; dYearStart dYearEnd ; ... ] (Default none)
opt.includePrjName={};           % Patterns of project names to include (if empty, all are excluded), this rule is applied after POI selection
opt.excludePrjName={};           % Patterns of project names to exclude (if empty, noe are excluded), this rule is applied after POI and include rule selection
opt.heightOnly=false;            % If true, then only import the height component to the STM, ignore everything else
opt.splitNetwork='warn';         % Split network {'warn','delete','keeplargest','split'}; warning in case there are multiple networks in an epoch,
                                 % delete the epoch, keep only the largest network, split in separate epochs/campaigns.
opt.sdObsFlag=[];                % select observations with these integer sdobsflag values, if empty, all are taken
opt.minObs=2;                    % Minimum number of observations per point
opt.projectId='';                % Project Id (default directory of the outputfile, or empty if no directory is specified)
opt.globalAttrib=struct([]);     % Struct, if empty (default), take only globalAttributes from Cupido NetCDF

% Duplicate output to file, and catch errors, and start timing

try

[~,outputFileRoot]=fileparts(outputFilename);
diary([ outputFileRoot '_' datestr(now,30) '.log' ])

fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Check the options and if necessary overwrite the default values

[inputFilename,outputFilename,opt]= ...
    stmcheckarguments(inputFilename,outputFilename,opt,varargin{:});
if isempty(outputFilename)
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    return;
end

%% Read Cupido NetCDF file

[pntname,pntcrd,pntclass,prjname,prjepoch,prjclass, ...
   obstable,sdobs,sdcov,sdobsflag,sensitivity,epoch, ...
   finfo] = cupido_read_netcdf(fullfile(opt.inputDir,char(inputFilename)));
sdcov=single(sdcov); % remove once cupido_read_netcdf is updated
%maybe give pntcrd as double, so that conversion later on is not needed

%% Conversion preliminaries 

[projectId,datasetId,~]=fileparts(outputFilename);
if numel(opt.projectId) > 1
   projectId=opt.projectId;
end

techIdx=strmatch('technique',{finfo.Attributes.Name});
if strmatch('Levelling',finfo.Attributes(techIdx).Value)
  techniqueId='lev';
elseif strmatch('GNSS',finfo.Attributes(techIdx).Value) 
  techniqueId='gnss';
elseif strmatch('GPS',finfo.Attributes(techIdx).Value)
  techniqueId='gnss';
else
  error('The technique in the input file is not supported.');
end

numObs = size(obstable,1);
numPoints = size(pntname,1);
numEpochs = size(prjname,1);

% Create epochs in dyears    
prjepochStr = datestr(prjepoch,'yyyymmdd');
epochDyear = date2dyear(prjepochStr,'yyyymmdd');

% Create lat,lon coordinates
% (for now ETRS89, change?)
mapCrs='RD';
mapCrd= double(pntcrd);
pntCrd = rdnap2etrs([mapCrd zeros(numPoints,1)],'PLH'); %double for double precision output
pntCrd(:,1:2) = pntCrd(:,1:2)*180/pi;


% Select period of interest (POI)

fprintf('Number of epochs in input space time dataset %d\n',numEpochs);
if isempty(opt.POI) 
   epochMask=true(1,numEpochs);
else
   epochMask=getepochmask(epochDyear,opt.POI);
   fprintf('Number of epochs inside POI %d (out of %d)\n',sum(epochMask),numEpochs);
end

% Update epochMask based on project names to include and exclude

if ~isempty(opt.includePrjName)
   tf=contains(cellstr(prjname),opt.includePrjName);
   fprintf('Option includePrjMask caused %d epochs (out of %d) to be excluded\n',sum(~tf & epochMask),sum(epochMask));
   epochMask = epochMask & tf; 
end
if ~isempty(opt.excludePrjName)
   tf=~contains(cellstr(prjname),opt.excludePrjName);
   fprintf('Option excludePrjMask caused %d epochs (out of %d) to be excluded\n',sum(~tf & epochMask),sum(epochMask));
   epochMask = epochMask & tf; 
end

% remove epochs

fprintf('Number of selected epochs %d (out of %d)\n',sum(epochMask),numEpochs);
if any(~epochMask)
   numEpochsRemove=sum(~epochMask);
   cupido_remove_epochs;
   fprintf("Removed %d epochs out of %d from observation table, %d observations out of %d are kept.\n",numEpochsRemove,numEpochs,size(obstable,1),numObs);
   numEpochs = size(prjname,1);  
   numObs = size(obstable,1);  
end

% Select regions of interest  (ROI) 
% - must be done after adjusting POI
% - there is a catch: from points must be inside the ROI, if not, observations will be removed
% - synthetic benchmarks must be dealt with separately (later)

fprintf('Number of points in input space time dataset %d\n',numPoints);
if isempty(opt.ROI)
   pntMask=true(numPoints,1);
else
   % There is a potential pitfall here, if a from points in a campaign is outside the mask, the whole campagin is removed
   pntMask=getpntmask(pntCrd,opt.ROI);
   fprintf('Number of points inside ROI %d (out of %d)\n',sum(pntMask),numPoints);
   if any(~pntMask)
       numPointsRemove=sum(~pntMask);
       cupido_remove_points;
       fprintf("Removed %d points out of %d from observation table, %d observations out of %d are kept.\n",numPointsRemove,numPoints,size(obstable,1),numObs);
       numPoints = size(pntname,1);  
       numObs = size(obstable,1);
   end
end

% remove observations based on sdobsflag

if ~isempty(opt.sdObsFlag)
   fprintf('Keep observations with sdobsflag value(s): '); for k=1:numel(opt.sdObsFlag), fprintf('%d ',opt.sdObsFlag(k)); end; fprintf('\n');
   itmp=unique(sdobsflag);
   fprintf('sdobsflag  numobs\n')
   for k=1:numel(itmp)
       ktmp=sum(ismember(sdobsflag,itmp(k)));
       if ismember(itmp(k),opt.sdObsFlag), ctmp='keep';, else ctmp='remove';, end;
       fprintf('   %2d    %6d     %s\n',itmp(k),ktmp,ctmp)
   end
   removeObs = ~ismember(sdobsflag, opt.sdObsFlag);
   cupido_remove_obs;
   fprintf('Number of removed observations: %d\n',numObs-size(obstable,1));
   numObs = size(obstable,1);
end

% check that points have observations

pntHasObs=false(numPoints,1);
pntHasObs(unique(obstable(:,1)))=true;
pntHasObs(unique(obstable(:,2)))=true;
if any(~pntHasObs)
   numObs = size(obstable,1);
   cupido_remove_points_noobs;
   fprintf("Removed %d points out of %d because they have no observations left, rehashed observation table (check %d =? %d)\n",sum(~pntHasObs),numPoints,size(obstable,1),numObs);
   numPoints = size(pntname,1);  
   numObs = size(obstable,1);
end


%% Conversion of observation data

% Determine/check dimension of the dataset

totalSensitivity = sum(sensitivity);
switch techniqueId
    case 'lev'
      if any(totalSensitivity(1:2) > 0)
         error("Levelling Cupido dataset has North and/or East components, this cannot be...")
      end
      obsTypes = {'Up'};
      ndim = 1;
    case 'gnss'
      if any(totalSensitivity(1:2) > 0)
          fprintf("Cupido GNSS dataset has North and/or East components.\n")
          obsTypes = {'North' 'East' 'Up'};
          ndim = 3;
      else
          fprintf("Cupido GNSS dataset has only Up components.\n")
          obsTypes = {'Up'};
          ndim = 1;
      end
    otherwise
        % e.g. InSAR, dim=1, but no one's in the sensitivity matrix  
        obsTypes = {'los'};
        ndim = 1;
end

% Optionally keep only height components

if ndim > 1 && opt.heightOnly
    obsTypes = {'Up'};
    ndim = 1;
    removeObs = any( sensitivity(:,1:2) > 0, 2 );  
    obstable(removeObs,:)=[];
    sdobs(removeObs)=[];
    sdobsflag(removeObs)=[];
    sensitivity(removeObs,:)=[];
    epoch(removeObs)=[];
    sdcov(removeObs,:)=[]; 
    sdcov(:,removeObs)=[]; 
    fprintf("Observation dimension reduced from 3 to 1; out of %d observations only %d up observations are kept.\n",numObs,size(obstable,1));
    numObs = size(obstable,1);  
end

% handle multi-dimensional cases (not yet implemented)

if ndim > 1
    error("Support for multi-dimensional Cupido datasets is not yet implemented, use opt.heightOnly=true to continue.")
end

% Space-Time matrix format is less flexible than Cupido format, first check if 
% conversion is possible (see also opt.splitNetwork)
% - only one from point can be handled every epoch
% - duplicates within the same epoch are not used
% We check this by counting the from and to observations

fromSpaceTimeCount = zeros(numPoints,numEpochs);
toSpaceTimeCount = zeros(numPoints,numEpochs);
for w = 1:numObs
  fromSpaceTimeCount(obstable(w,1),obstable(w,3)) = fromSpaceTimeCount(obstable(w,1),obstable(w,3)) + 1;
  toSpaceTimeCount(obstable(w,2),obstable(w,3)) = toSpaceTimeCount(obstable(w,2),obstable(w,3)) + 1;
end
fromCountByEpoch = sum(fromSpaceTimeCount > 0, 1);
if any(fromCountByEpoch < 1)                          % Added 11 April 2024
   % remove epochs with zero observations
   epochMask(fromCountByEpoch < 1) = false ;
   cupido_remove_epochs;
   fprintf("Removed %d epochs without observations from observation table, out of %d observations %d are kept.\n",sum(fromCountByEpoch < 1),numObs,size(obstable,1));
   numEpochs = size(prjname,1);  
   numObs = size(obstable,1);  
   % Adjust toSpaceTimeCount and fromSpaceTimeCount
   toSpaceTimeCount(:,fromCountByEpoch < 1) = [];
   fromSpaceTimeCount(:,fromCountByEpoch < 1) = [];
   fromCountByEpoch(fromCountByEpoch < 1) = [];
end
if any(fromCountByEpoch ~= 1)
    if opt.verbose > 1
       fprintf('fromCountByEpoch');
       disp(fromCountByEpoch)
    end
    for k=1:numEpochs
        if fromCountByEpoch(k) > 1
            fprintf('Epoch: %d (%s), numFromPoints: %d\n', k, prjepochStr(k,:), fromCountByEpoch(k))
            for l=1:numPoints
               if fromSpaceTimeCount(l,k) > 0
                    fprintf('    From point: %d (%s), numFromObs: %d\n', l, pntname(l,:), fromSpaceTimeCount(l,k))
               end
            end
        end
    end
    switch lower(opt.splitNetwork) 
        case 'warn'
            error("Option splitNetwork=warn -> Cupido dataset cannot be translated into a space-time matrix; there is more than one from point per epoch.")
        case 'delete'
            % remove epochs
            epochMask(fromCountByEpoch > 1) = false ;
            cupido_remove_epochs;
            fprintf("Option splitNetwork=delete -> removed %d epochs from observation table, out of %d observations %d are kept.\n",sum(fromCountByEpoch > 1),numObs,size(obstable,1));
            numEpochs = size(prjname,1);  
            numObs = size(obstable,1);  
            % Adjust toSpaceTimeCount and fromSpaceTimeCount
            toSpaceTimeCount(:,fromCountByEpoch > 1) = [];
            fromSpaceTimeCount(:,fromCountByEpoch > 1) = [];
        case 'keeplargest'
            error("Unimplemented option splitNetwork=keeplargest.")
        case 'split'
            error("Unimplemented option splitNetwork=split.")
        otherwise
            error("Invalid option splitNetwork, must be [warn|keeplargest|split].")
    end
end
if any(toSpaceTimeCount(:) > 1)
    fprintf('Maximum toSpaceTimeCount %d\n\n',max(toSpaceTimeCount(:)))
    if opt.verbose > 1
       fprintf('toSpaceTimeCount');
       disp(toSpaceTimeCount)
    end
    toCountByEpoch = sum(toSpaceTimeCount > 1);
    for k=1:numEpochs
        if toCountByEpoch(k) > 0 
            fprintf('Epoch: %d (%s), numToPoints exceeded: %d\n', k, prjepochStr(k,:), toCountByEpoch(k))
            for l=1:numPoints
               if toSpaceTimeCount(l,k) > 1
                    fprintf('    To point: %d (%s), numToObs: %d\n', l, pntname(l,:), toSpaceTimeCount(l,k))
               end
            end
        end
    end
    error("Cupido dataset cannot be translated into a space-time matrix; there are duplicate observations in an epoch.")
end

% Remove points with less than opt.minObs(=2) observations - Added 6 June 2024 by HM

fromSpaceTimeCount = zeros(numPoints,numEpochs);
toSpaceTimeCount = zeros(numPoints,numEpochs);
for w = 1:numObs
  fromSpaceTimeCount(obstable(w,1),obstable(w,3)) = fromSpaceTimeCount(obstable(w,1),obstable(w,3)) + 1;
  toSpaceTimeCount(obstable(w,2),obstable(w,3)) = toSpaceTimeCount(obstable(w,2),obstable(w,3)) + 1;
end
% Do the test, but we must leave from always in (even if used in only one campaign)
pntMask=sum( fromSpaceTimeCount+toSpaceTimeCount > 0 ,2) >= opt.minObs | sum(fromSpaceTimeCount > 0, 2) > 0;
fprintf('Number of points with less than %d observations: %d (out of %d)\n',opt.minObs,sum(~pntMask),numPoints);
if any(~pntMask)
     numPointsRemove=sum(~pntMask);
     toSpaceTimeCount(~pntMask,:) = [];
     fromSpaceTimeCount(~pntMask,:) = [];
     cupido_remove_points;
     fprintf("Removed %d points out of %d from observation table, %d observations out of %d are kept.\n",numPointsRemove,numPoints,size(obstable,1),numObs);
     numPoints = size(pntname,1);  
     numObs = size(obstable,1);
end

% Remove points without coordinates (e.g. synthetic benchmarks in case of GNSS)

synBench =  isnan(pntCrd(:,1)) | isnan(pntCrd(:,2)) | strcmpi('SYN_BM',pntname);
if any(synBench)
   % remove points
   pntname(synBench,:) = [];
   mapCrd(synBench,:) = [];
   pntCrd(synBench,:) = [];
   pntclass(synBench,:) = [];
   pntMask(synBench) = []; 
   % update arrays with counts
   toSpaceTimeCount(synBench,:) = [];
   fromSpaceTimeCount(synBench,:) = [];
   % reindex observation table (both from and to points) - New 11 April 2024
   pntidx = nan(size(synBench));
   pntidx(~synBench) = 1:sum(~synBench);
   obstable(:,1) = pntidx(obstable(:,1));
   obstable(:,2) = pntidx(obstable(:,2));
   fprintf("Removed %d synthetic benchmark from point list, out of %d points %d are kept.\n",sum(synBench),numPoints,size(pntname,1));
   numPoints = size(pntname,1);  
   numObs = size(obstable,1);
end

% Create Space-Time matrix
spaceTimeMatrix = NaN(numPoints,numEpochs,'single');
for w = 1:numObs
  %spaceTimeMatrix(obstable(w,1),obstable(w,3)) = 0; 
  spaceTimeMatrix(obstable(w,2),obstable(w,3)) = sdobs(w) * 1e3;
end

% Add from points (one per campaign - with zero observation) to space-time matrix - New 11 April 2024
idxFromPoints = nan(1,numEpochs);
fromPoints=cell(1,numEpochs);
for k = 1:numEpochs
  if ~isempty(find(fromSpaceTimeCount(:,k) > 0))
     idxFromPoints(k) = find(fromSpaceTimeCount(:,k) > 0);
     if isnan(spaceTimeMatrix(idxFromPoints(k),k)) 
        spaceTimeMatrix(idxFromPoints(k),k) = 0; 
     else
        error('Something is wrong with the fromPoints, this should never happen')
     end
     fromPoints{k}=strtrim(pntname(idxFromPoints(k),:));
  end
end

% Create Stochastic Model
stochData=zeros(numPoints,numPoints,numEpochs,'single');
for k=1:numEpochs
   idx=find(obstable(:,3)==k);
   stochData(obstable(idx,2),obstable(idx,2),k)=sdcov(idx,idx) * 1e6;
end

% Create Sensitivity matrix (shape and size is different from Cupido)
switch techniqueId
  case 'lev'
    sensitivityMatrix = [zeros(numPoints,2,'single') ones(numPoints,1,'single')];
  case 'gnss'
    if ndim == 1
       sensitivityMatrix = [zeros(numPoints,2,'single') ones(numPoints,1,'single')];
    else
       sensitivityMatrix = zeros(numPoints,3,ndim,'single');
       sensitivityMatrix(:,1,1) = ones(numPoints,1,'single');
       sensitivityMatrix(:,2,2) = ones(numPoints,1,'single');
       sensitivityMatrix(:,3,3) = ones(numPoints,1,'single');
    end
end 

% Create flag matrix (if needed)
auxTypes={};
flagMatrix=[];
if ~isempty(sdobsflag) || sum(sdobsflag)==0
  flagMatrix=NaN(numPoints,numEpochs,2,'single');
  for w = 1:numObs
    flagMatrix(obstable(w,2),obstable(w,3),1) = sdobsflag(w);
    flagMatrix(obstable(w,2),obstable(w,3),2) = epoch(w);
  end
  auxTypes={'obsflag' 'epoch'};
end

%% Print epoch table

printEpochTable(cellstr(prjname),epochDyear,spaceTimeMatrix)

%% Create stm dataset structure 

stout = stm(projectId,datasetId,techniqueId);

% Dataset attributes
datasetAttrib=stout.datasetAttrib;
datasetAttrib.softwareOptions=opt;
stout.datasetAttrib=datasetAttrib;

% Technique attributes
techniqueAttrib=stout.techniqueAttrib;
techniqueAttrib.mapCrs=mapCrs;             % Name of coordinate reference system for mapCrd
techniqueAttrib.mode='Campaign';
stout.techniqueAttrib=techniqueAttrib;

% Epochs
stout.numEpochs = numEpochs;
stout.epochDyear = epochDyear';

% Store pntAttrib first in temporary array, before moving it to dataset (objects only support one level of indexing)
epochAttrib=[];
epochAttrib.prjName = cellstr(prjname)';
epochAttrib.prjClass = cellstr(prjclass)';
epochAttrib.fromPoints = fromPoints;            % New 11 April 2024
stout.epochAttrib = epochAttrib;

% Points
stout.numPoints = numPoints;
stout.pntName = cellstr(pntname);
stout.pntCrd = pntCrd;

% Store pntAttrib first in temporary array, before moving it to dataset (objects only support one level of indexing)
pntAttrib=[];
pntAttrib.pntClass = cellstr(pntclass);
pntAttrib.mapCrd = mapCrd;
stout.pntAttrib=pntAttrib;

% Types
stout.parTypes = {'North' 'East' 'Up'};
stout.obsTypes = obsTypes;
stout.auxTypes = auxTypes;

% Matrices
stout.sensitivityMatrix = sensitivityMatrix;
stout.obsData = spaceTimeMatrix;
stout.auxData = flagMatrix;

% Stochastic model
stout.stochModel={['covmatrix(format=blkdiag,nd=' num2str(numEpochs) ')']};
stout.stochData=stochData;

% Input dataset history
inputDatasets(1).datasetId=datasetId;
inputDatasets(1).techniqueId=techniqueId;
datasetAttrib=[];
datasetAttrib.softwareName='';
datasetAttrib.softwareOptions=[];
datasetAttrib.fileFormat='Cupido NetCDF';
datasetAttrib.fileName=inputFilename;
datasetAttrib.projectId='';
inputDatasets(1).datasetAttrib=datasetAttrib;
inputDatasets(1).numPoints=numPoints;
inputDatasets(1).numEpochs=numEpochs;
inputDatasets(1).inputDatasets=[];
stout.inputDatasets=inputDatasets;

% Global attributes
globalAttrib=[];
for k=1:numel(finfo.Attributes)
  globalAttrib.(finfo.Attributes(k).Name)=finfo.Attributes(k).Value;
end
allfields=fieldnames(struct(opt.globalAttrib));
for k=1:numel(allfields)
  globalAttrib.(allfields{k})=opt.globalAttrib.(allfields{k});
end
stout.globalAttrib=globalAttrib;

%% Write the stm outputfile

stmwrite(stout,outputFilename);

% Finish the function

if nargout > 0
   varargout{1}=0;
end

fprintf('%s finished at %s  (elapsed time %.2f s)\n',progname,datestr(now),toc)
diary OFF

catch ME

   getReport(ME)
   fprintf('%s ABORTED with an error at %s (elapsed time %.2f s)\n',progname,datestr(now),toc)
   diary OFF
   rethrow(ME)
    
end

%% Nested functions (share variables with calling function; only possible when used within a function)

    function cupido_remove_epochs       
        if any(~epochMask)
            % remove epochs
            prjname = prjname(epochMask,:);
            prjclass = prjclass(epochMask,:);
            prjepochStr = prjepochStr(epochMask,:);
            epochDyear = epochDyear(epochMask);
            % remove observations belonging to deleted epochs
            removeObs = ismember( obstable(:,3) , find(~epochMask));  
            obstable(removeObs,:) = [];
            sdobs(removeObs) = [];
            sdobsflag(removeObs) = [];
            sensitivity(removeObs,:) = [];
            epoch(removeObs) = [];
            sdcov(removeObs,:) = []; 
            sdcov(:,removeObs) = []; 
            % reindex observation table
            epochidx = nan(size(epochMask));
            epochidx(epochMask) = 1:sum(epochMask);
            obstable(:,3) = epochidx(obstable(:,3));
            % adjust epochMask
            epochMask = epochMask(epochMask); 
        end
    end

    function cupido_remove_points
        if any(~pntMask)
            % remove points
            pntname = pntname(pntMask,:);
            mapCrd = mapCrd(pntMask,:);
            pntCrd = pntCrd(pntMask,:);
            pntclass = pntclass(pntMask,:);
            % remove observations belonging to deleted from_ and to_points 
            removeObs = ismember( obstable(:,1) , find(~pntMask)) | ismember( obstable(:,2) , find(~pntMask));  
            obstable(removeObs,:) = [];
            sdobs(removeObs) = [];
            sdobsflag(removeObs) = [];
            sensitivity(removeObs,:) = [];
            epoch(removeObs) = [];
            sdcov(removeObs,:) = []; 
            sdcov(:,removeObs) = []; 
            % reindex observation table (both from and to points)
            pntidx = nan(size(pntMask));
            pntidx(pntMask) = 1:sum(pntMask);
            obstable(:,1) = pntidx(obstable(:,1));
            obstable(:,2) = pntidx(obstable(:,2));
            % adjust pntMask
            pntMask = pntMask(pntMask); 
        end
    end

    function cupido_remove_obs
        if any(removeObs)
            % remove observations belonging to deleted from_ and to_points 
            obstable(removeObs,:) = [];
            sdobs(removeObs) = [];
            sdobsflag(removeObs) = [];
            sensitivity(removeObs,:) = [];
            epoch(removeObs) = [];
            sdcov(removeObs,:) = []; 
            sdcov(:,removeObs) = []; 
        end
    end

    function cupido_remove_points_noobs
        if any(~pntHasObs)
            % remove points
            pntname = pntname(pntHasObs,:);
            mapCrd = mapCrd(pntHasObs,:);
            pntCrd = pntCrd(pntHasObs,:);
            pntclass = pntclass(pntHasObs,:);
            pntMask = pntMask(pntHasObs); 
            % reindex observation table (both from and to points)
            pntidx = nan(size(pntHasObs));
            pntidx(pntHasObs) = 1:sum(pntHasObs);
            obstable(:,1) = pntidx(obstable(:,1));
            obstable(:,2) = pntidx(obstable(:,2));
        end
    end

end

%% Local functions (have their own workspace)

function varargout = cupido_read_netcdf(netcdf_file)
%CUPIDO_READ_NETCDF  Read CUPiDO Netcdf file 
%   CUPIDO_READ_NETCDF(NETCDF_FILE) reads the CUPiDO Netcdf file NETCDF_FILE
%   and prints a summary of the content.
%
%   [pntname,pntcrd,pntclass, prjname,prjepoch,prjclass,obstable,sdobs,...
%   sdcov,sdobsflag,sensitivity,epoch, finfo] = cupido_read_netcdf(netcdf_file)
%   reads the CUPiDO Netcdf file NETCDF_FILE and outputs the content
%   as variables.
%
%   Example:
%
%      cupido_read_netcdf('cupido_gps.nc');
%
%      [pntname,pntcrd,pntclass, prjname,prjepoch,prjclass, ...
%              obstable,sdobs,sdcov,sdobsflag,sensitivity,epoch, ...
%              finfo] = cupido_read_netcdf('cupido_gps.nc');  
%
%   See also cupido_write_netcdf and cupido_merge_netcdf.
%
%   (c) Hans van der Marel, Delft University of Technology, 2016. 

% Created:  12 Apr 2016 by Hans van der Marel
% Modified: 12 Apr 2016 by Hans van der Marel
%              - Initial version
%           24 Aug 2016 by Hans van der Marel
%              - New structure
%           14 Sep 2016 by Hans van der Marel
%              - Additional plots 
%           10 Oct 2016 by Hans van der Marel
%              - converted script to function

%% Define netcdf file name 

if nargin ~=1, error('This function requires one input argument.');, end


%% Display NetCDF file schema in command window

% ncdisp(netcdf_file)

%% Get information about NetCDF file into structure finfo and print attributes

fprintf('Netcdf file: %s\n\n',netcdf_file)

finfo=ncinfo(netcdf_file);

fprintf('Global Attributes:\n\n')
for k=1:numel(finfo.Attributes)
  fprintf('%s: %s\n',finfo.Attributes(k).Name,finfo.Attributes(k).Value)
end
fprintf('\n\n')

%% Read netcdf file


% Read point data

pntname=ncread(netcdf_file,'station_name');
pntcrd(:,1)=ncread(netcdf_file,'x');
pntcrd(:,2)=ncread(netcdf_file,'y');
pntclass=ncread(netcdf_file,'station_class');

% Project data

prjname=ncread(netcdf_file,'project_name');
prjepoch=ncread(netcdf_file,'project_epoch')+datenum('1-Jan-1970 00:00:00');
prjclass=ncread(netcdf_file,'project_class');

% Observations

%  sdobstable     observation table with index to from_point, to_point and project 
%  epoch          array with epoch (Matlab date number)
%  sdobs          array with the observed height difference [m]
%  sdcov          covariance matrix [m]
%  sdobsflag      integer observation flag (default 0)
%  sensitivity    sensitivity matrix [0-1]
 
stationFromIndex=ncread(netcdf_file,'stationFromIndex');
stationToIndex=ncread(netcdf_file,'stationToIndex');
projectIndex=ncread(netcdf_file,'projectIndex');

sdobs=ncread(netcdf_file,'sdObs');
sdcov=ncread(netcdf_file,'sdCov');

epoch=ncread(netcdf_file,'epoch')+datenum('1-Jan-1970 00:00:00');
sdobsflag=ncread(netcdf_file,'sdObsFlag');
sensitivity=ncread(netcdf_file,'sensitivity');

sdobstable=[stationFromIndex stationToIndex projectIndex ];

%% 

if nargout >=1 
  varargout={  pntname,pntcrd,pntclass, ...
               prjname,prjepoch,prjclass, ...
               sdobstable,sdobs,sdcov,sdobsflag,sensitivity,epoch,finfo }; 
  return;
end

%% Print point data

numpnt=size(pntname,1);

fprintf('\nBenchmarks (%d points):\n',numpnt)

fprintf('\nPNTNAME             X_RD         Y_RD   CLASS\n\n')
for k=1:numpnt
   fprintf('%-10s  %12.3f %12.3f   %s\n',pntname(k,:),pntcrd(k,:),pntclass(k,:))
end
fprintf('\n');

%% Print project data

numprj=size(prjname,1);

fprintf('\nProjects (%d projects):\n',numprj)

fprintf('\nPRJNAME      MEAN_EPOCH   CLASS\n\n')
dfmt='yyyy-mm-dd';
for k=1:numprj
   fprintf('%-10s   %s   %s\n',prjname(k,:), ...
       datestr(prjepoch(k),dfmt),prjclass(k,:));
end
fprintf('\n');

%% Print observation data

numobs=size(sdobstable,1);

fprintf('\nObservations (%d observations):\n',numobs)

fprintf('\nFROM       TO         PROJECT      OBS [m] STDEV [mm]  FLAG   SENSITIVITY   EPOCH\n\n')
dfmt='yyyy-mm-dd HH:MM';
for k=1:numobs
   fprintf('%-10s %-10s %-10s%10.3f %10.3f  %4d   %3.1f %3.1f %3.1f   %s\n', ...
       pntname(sdobstable(k,1),:),pntname(sdobstable(k,2),:),prjname(sdobstable(k,3),:), ...
       sdobs(k),sqrt(sdcov(k,k))*1000,sdobsflag(k),sensitivity(k,:),datestr(epoch(k),dfmt));
end
fprintf('\n');


%% Done

end

%% Function to print epoch table

function printEpochTable(prjName,epochDyear,obsData,epochMask)

if nargin < 4
    epochMask=true(size(epochDyear));
end

npoints = sum(~isnan(obsData(:,epochMask,end)));

[t,idx]=sort(epochDyear(epochMask));
tmpname=prjName(epochMask);
tmpname=tmpname(idx);
npoints=npoints(idx);

fprintf('\nprjName               dYear     yyyy-mmm-dd  dDays  #Pnts\n')
fprintf('--------------------  --------  -----------  -----  -----\n')
last=dyear2date(t(1));
for k=1:numel(tmpname)
    current=dyear2date(t(k)); 
    fprintf('%-20s  %8.3f  %11s  %5d %6d\n' ,tmpname{k},t(k),datestr(current,'yyyy-mmm-dd'),current-last,npoints(k))
    last=current;
end
fprintf('\n\n')

end
