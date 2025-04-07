function varargout=stmintegrate(inputfilenames,outputfilename,options,flag)
%stmintegrate   Integrated processing of geodetic displacement observations.
%   STMINTEGRATE(INPUTFILENAMES,OUTPUTFILENAME,OPTIONS) integrates multiple 
%   input Space Time Matrix datasets, given by INPUTFILENAMES, into a single 
%   integrated output Space Time Matrix dataset, with the name OUTPUTFILENAME,
%   using least-squares adjustment and testing on common points between
%   the input datasets. INPUTFILENAMES must be the name of an file containing 
%   the input filenames, or a cell array with the input filenames. 
%   OUTPUTFILENAME is a character strin with the name of the output 
%   Space Time Matrix dataset, with the displacement time series for the
%   combined points, the transformation parameters for the input datasets,
%   and the result of the statistical testing. OPTIONS is a cell array or
%   structure with the processing options, if empty, or when fields are missing, 
%   default values are used.
%
%   STMINTEGRATE(...,UPDATEFLAG) affects the processing in case OUTPUTFILENAME
%   exists. The following values of UPDATEFLAG are accepted
%
%     'create'     abort processing if the output dataset exists (default)
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
%
%   STAT=STMINTEGRATE(...) returns a status code STAT. A status code of 0
%   indicates success.
% 
%   OPTIONS is a cell array or structure with the processing options, 
%   if empty, or when fields are missing, default values are used. Valid 
%   OPTIONS are in the following tables
%
%   Main (with default)         Description of main input options
%   --------------------------  -------------------------------------------------------
%   selectTechnique={'lev' 'gnss' 'insar'}   Cell array with techniques to include (default is all supported techniques)
%   selectDataset= ''           Cell array with datasetIds to select, you can give a pattern to select, default is all
%
%   ROI=[]                      Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default all)
%   POI=[-Inf +Inf]             Period of interest [ dYearStart dYearEnd ] (Default all)
%
%   excludePointIds=[]          Cell array with pointIds to be excluded from the integration (Default none)
%   excludeEpochIds=[]          Cell array with epochIds to be excluded from the integration (Default none)
%
%   aprioriDispl=''             Name of dataset with a-priori displacement (Default none)
%   filterLos=false             If true, skip displacements which have only ascending, or only decending, InSAR observations (los qualifier)
%   filterNear=false            If true, skip displacements which have only InSAR observations (near qualifier)
%   filterPrioriEast=false      If true, skip displacements which depend on a-priori data for the East component
%   filterPrioriNorth=false     If true, skip displacements which depend on a-priori data for the North component
%
%   overrideStochModel=struct() Structure with new stochastic models (opt), see updStochModel for help
%   addStochModel=struct()      Structure with stochastic models to add to existing (opt), see updStochModel for help
%
%   Output (with default)       Description of output options
%   --------------------------  -------------------------------------------------------
%   verbose=0                   Verbosity level, higher is more output, 0 is almost nothing
%   doplots=0                   Plot level, higher is more plots, 0 is no plots
%
%   saveResiduals=true          If true, save residuals and transformation parameters (default) to a file OUTPUTFILENAME_res
%   saveIntermediate=false      If true, save intermediate least-squares results (optional) to a file OUTPUTFILENAME_lsq
%
%   covmatrix='byobstype'       Ouput stochastic model  {'chol', 'full', 'byobstype'}
%
%   Testing (with default)      Description of options for statistical testing
%   --------------------------  -------------------------------------------------------
%   alpha0=0.001                Probability of false alarm
%
%   critWtest=nan               Criterion for printing test statistics, if nan, then
%   critDisplOmt=nan            the criterion is computed from the probability of false alarm
%   critPointOmt=nan            alpha and the redundancy 
%   critEpochOmt=nan       
%
%   maxWtest=50                 Maximum of test statistics to print
%   maxDisplOmt=20
%   maxPointOmt=10
%
%   Internals (with default)    Description of internal options
%   --------------------------  -------------------------------------------------------
%   minEpochs=2                 Minumum number of epochs in combination, if below, stop with error (Default is 2)
%   minPoints=2                 Minimum number of points in combination, if below, stop with error (Default is 2)
%
%   minDatasetsEpoch=2          Minumum number of datasets per epoch, if below, epoch is removed (Default is 2)
%   minDatasetsPoint=2          Minimum number of datasets per point, if below, point is removed (Default is 2)
%   forceDatasetPoints=''       DatasetIds for which points should be included, even if not redundant (Default is '')
%
%   minDisplEpoch=2             Minimum number of displacements per epoch, if below, epoch is removed (Default is 2)
%   minDisplPoint=2             Minimum number of displacements per point, if below, point is removed (Default is 2)
%   minDispl=4                  Minimum number of displacements in total for a component, if below, component is removed (Default is minDisplEpoch*minDisplPoint) 
%
%   minObsDispl=1               Minimum number of observations in an input space time matrix row, if below, row is ignored
%
%   onlyRedundantEpochs=true    Keep only epochs with redundant data, ie. with connecting points between datasets
%   onlyRedundantPoints=true    Keep only points with redundant data  
%
%   pntCrdType='deg/m'          Coordinate type in <stm>.pntCrd {'deg/m', 'km/m'}, default 'deg/m' (Latitude, Longitude, Height)
%   tolEpoch=.3                 Tolerance for identical epochs in case epochs haven't been harmonized
%   tolPointDist=.001           Tolerance for distance between identical points (in km) in case points haven't been harmonized
%   ignoreProjectCheck=false    If true, project data is ignored in the check for harmonized points and epochs (default is false)
%
%   cbaseOffset='firstEpoch'    Computing base for ref epoch { 'firstEpoch' , 'lastEpoch', 'firstCommonEpoch' , 'lastCommonEpoch' , <decimalYear> }
%   cbaseTransform='firstPoint' Computing base for ref time series { 'firstPoint' , 'lastPoint', 'firstCommonPoint' , 'lastCommonPoint' , <pointId> }
%                               - In case of '[first|last][Epoch|Point]' always the first or last epoch (point) is chosen with all components observed. 
%                                 This is not necessarily the same point for each epoch, or vice versa.
%                               - In case of ''[first|last]Common[Epoch|Point]' the first or last epoch (point) is chosen with all 
%                                 components and points (epochs) observed. When this is not possible, the epoch (point) is selected with
%                                 most observed points (epochs).
%                               When a specific epoch (point) is selected it must be an epoch with all components available and
%                               all points (epochs) observed for that epoch (point).
%
%   sreg=1e-3                   Regularization parameter for output cov-matrices
%
%   Examples:
%      stmintegrate('simtrue1_filenames.txt',simtrue1_combined.mat')      
%      stmintegrate('simtrue1_filenames.txt',simtrue1_combined.mat',options,'update')      
%      stmintegrate('simtrue1_filenames.txt',simtrue1_combined.mat')      
%
%   See also STM, STMCHECKARGUMENTS, STMREAD and STMDIFF.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020-2023.

% Created:   8 July 2020 by Hans van der Marel
% Modified: 20 July 2020 by Hans van der Marel
%            - split development version in ispEstDefo2.m and ispLsqDefo2.m
%            - ispEstDefo2.m contains the symbolic analysis of the observation structure
%            - ispLsqDefo2.m contains only bare minumum for Lsq solution
%           21 August 2020 by Hans van der Marel
%            - ispLsqDefo2.m renamed to ispLsqDefo2b.m
%            - changed to new data structure format and functions from stmutil
%            - compatible with ispSimDefo2b.m
%            - merged code for solvability analysis, call to ispSolvability,
%              works now for all techniques
%           20 October 2020 by Hans van der Marel
%            - first alpha release of the integration module (ready for testing on real data)
%            - complete overhaul of the software since the 21 August version:
%              Meer uitgebreid testen van alle opties en alternatieven met “underobserved” 
%              displacement parameters (nearUp/nearEast/losAsc/locDsc). De verschillende methodes zijn:
%                1. met a-priori deformatie: deze wordt ingelezen als stm en vervolgens geïnterpoleerd naar de punten 
%                   en tijdstippen in de gemeenschappelijke dataset (gereed, muv van temporele interpolatie). 
%                   Ook weer twee opties
%                   a. als keiharde constraint (momenteel geïmplementeerd)
%                   b. als pseudo waarneming met zijn eigen standaard afwijking (nog niet geïmplementeerd, maar eenvoudig te doen)
%                2. zonder a-priori deformatie: de parameters krijgen een vlag (ook in bovenstaande optie), echter als 
%                   verschillen soorten parameters voor een en hetzelfde punt (bv GPS start later) voorkomt moeten 
%                   verschillende transformatie parameters geschat worden.
%            - many improvements to the output
%            - tested on simulated data
%           28 October 2020 by Hans van der Marel
%            - second alpha release (following real data testing)
%            - selection of ROI and POI implemented
%            - selection of observations and parameters based on options
%            2 November 2020 by Hans van der Marel
%            - check on minimum number of points and epochs in datasets after filtering
%            - problem with non-positive definite cov matrix in tudinsar4rd
%            - implemented overwrite functions for stochastic model
%           13 February 2021 by Hans van der Marel
%            - alpha release of statistical testing
%            - output of intermediate least squares results to mat file
%           02 March 2021 by Hans van der Marel
%            - saved transformation and offset parameters to residual file
%            - wrap up of statistical testing
%            - first beta release
%           19 October 2021 by Hans van der Marel
%            - remove points and epochs when there are no connections,
%              depending on input option (default is yes)
%            - remove points from dataset when all observations are Nan's in obsData
%            - halt the processing when not enough epochs/points left for one dataset
%            - reduced some long points list for lower verbosity levels
%           02 June 2022 by Hans van der Marel
%            - newly defined pntName (not anymore identical to pntIds)
%            - use pntName in print outputs 
%            - added pntName and pntCrd to residual and intermediate files
%           26 October 2022 by Hans van der Marel
%            - fixed bug in call to updStochModel (only affected print output)
%            - remove duplicates (maintaining the order) in pntName
%            - increased length of pntName in print statements (28 -> 38)
%            - compute and print condition numbers of Qyy, U22 and U11
%            - option to plot diagonal of Choleski factors (opt.plotCond)
%           10 November 2022 by Hans van der Marel
%            - New option opt.minObsDispl and related code 
%           28 Sep 2023 by Hans van der Marel
%            - fixed bug epochAttrib assignment for stochModel comp
%           14 Oct 2023 by Hans van der Marel
%            - New option opt.forceDatasetPoints giving DatasetIds for which 
%              points should be included, even if not redundant
%            - Modified computation of pntMaskDs and epochMaskDs, updating
%              pntMaskDs whenever epochMaskDs has been changed TODO: this
%              task is currently done for each dataset, but should be done
%              for each layer as it involves actual data
%            - Fixed bug in application of Ry factor (needed transpose) and
%              updated comment lines associated with this change
%            - Fixed epochAttrib.numDisplPerEpoch in output stm (forgot transpose) 
%           29 Oct 2023 by Hans van der Marel
%            - Updates to code changes of 14 Oct 2023
%           30 Jan 2024 by Hans van der Marel
%            - Output of covariance matrix instead of cholesky (chol2stmcov)
%            - Removed stochIndex from stm ouput (now part of stochData)
%            - New options associated with the change
%            7 Mar 2024 by Hans van der Marel
%            - Minor changes (debugging output under opt.verbose)
%
%           xx xxxxx xxxx
%            - optimize processing in case of block diagonal covariance matrix?
%            - extend options for computing base and connection to stable benchmarks?
%            - rigorous computation of statistical test statistics (using Qe instead
%              of approximated redundancy numbers)?


%% Check the input arguments and options

% For testing purposes we don't call this as function, but hardcode the inputs.
% When called as script
% - comment out the function statetement (line 1)
% - comment out the next three lines
% - comment out the end statement at the end of the file, above local functions

if nargin < 3
     error('This function expects at least two input arguments.')
end

progname='stmintegrate';

% Default options 

opt.verbose=0;                                  % Default verbosity level, higher is more output, 0 is almost nothing
opt.doplots=0;                                  % Default plot level, higher is more plots, 0 is no plots
opt.plotCond=false;                             % Plot condition numbers (diagonal of Choleski factor)

opt.outputTarget='create';                      % Default mode for output STM {'create','overwrite','update'}
opt.saveResiduals=true;                         % If true, save residuals and transformation parameters (default) to a file OUTPUTFILENAME_res
opt.saveIntermediate=false;                     % If true, save intermediate least-squares results (optional) to a file OUTPUTFILENAME_lsq

opt.selectTechnique= { 'lev' 'gnss' 'insar' };  % Default is to include all supported techniques
opt.selectDataset= '';                          % Default is to not select on datasetIds, but you can give a pattern to select datasets

opt.ROI=[];                                     % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default none)
opt.POI=[-Inf +Inf];                            % Period of interest [ dYearStart dYearEnd ] (Default none)

opt.excludePointIds={};                         % Cell array with pointIds to be excluded from the integration
opt.excludeEpochIds={};                         % Cell array with epochIds to be excluded from the integration

opt.aprioriDispl='';                            % Name of dataset with a-priori deformation
opt.filterLos=false;                            % If true, skip displacements which have only ascending, or only decending, InSAR observations (los qualifier)
opt.filterNear=false;                           % If true, skip displacements which have only InSAR observations (near qualifier)
opt.filterPrioriEast=false;                     % If true, skip displacements which depend on a-priori data for the East component
opt.filterPrioriNorth=false;                    % If true, skip displacements which depend on a-priori data for the North component

opt.overrideStochModel=struct();                % Structure with new stochastic models (opt), see updStochModel for help
opt.addStochModel=struct();                     % Structure with stochastic models to add to existing (opt), see updStochModel for help

opt.minEpochs=2;                                % Minumum number of epochs in combination, if below, stop with error
opt.minPoints=2;                                % Minimum number of points in combination, if below, stop with error

opt.minDatasetsEpoch=2;                         % Minumum number of datasets per epoch, if below, epoch is removed
opt.minDatasetsPoint=2;                         % Minimum number of datasets per point, if below, point is removed
opt.forceDatasetPoints='';                      % DatasetIds for which points should be included, even if not redundant

opt.minDisplEpoch=2;                            % Minimum number of displacements per epoch, if below, epoch is removed
opt.minDisplPoint=2;                            % Minimum number of displacements per point, if below, point is removed
%opt.minDispl=4;                                % Minimum number of displacements in total for a component, if below, component is removed 
opt.minDispl=opt.minDisplEpoch*opt.minDisplPoint; 

opt.minObsDispl=1;                              % Minimum number of observed displacements, if below, input space time matrix row is ignored 

opt.onlyRedundantEpochs=true;                   % Keep only epochs with redundant data, ie. with connecting points between datasets
opt.onlyRedundantPoints=true;                   % Keep only points with redundant data  

opt.pntCrdType='deg/m';                         % Coordinate type in <stm>.pntCrd {'deg/m', 'km/m'}, default 'deg/m' (Latitude, Longitude, Height)
opt.tolEpoch=.3;                                % Tolerance for identical epochs in case epochs haven't been harmonized
opt.tolPointDist=.001;                          % Tolerance for distance between identical points (in km) in case points haven't been harmonized
opt.ignoreProjectCheck=false;                   % If true, project data is ignored in the check for harmonized points and epochs (default is false)

opt.cbaseOffset='firstEpoch';                   % Computing base for ref epoch { 'firstEpoch' , 'lastEpoch', 'firstCommonEpoch' , 'lastCommonEpoch' , <decimalYear> }  1)
opt.cbaseTransform='firstPoint';                % Computing base for ref time series { 'firstPoint' , 'lastPoint', 'firstCommonPoint' , 'lastCommonPoint' , <pointId> }

opt.alpha0=0.001;                               % probability of false alarm

opt.critWtest=nan;                              % criterion for printing test statistics,
opt.critDisplOmt=nan;                           % if nan, then the criterion is computed
opt.critPointOmt=nan;                           % from the probability of false alarm
opt.critEpochOmt=nan;                           % opt.alpha and the redundancy

opt.maxWtest=50;                                % maximum of test statistics to print
opt.maxDisplOmt=20;
opt.maxPointOmt=10;

opt.covmatrix='byobstype';                      % output stoch model {'chol', 'full', 'byobstype'}
opt.sreg=1e-3;                                  % regularization parameter for cov-matrices

% Notes:
% 
% 1) In case of '[first|last][Epoch|Point]' always the first or last epoch (point) is chosen with all components observed. 
%    This is not necessarily the same point for each epoch, or vice versa.
%    In case of ''[first|last]Common[Epoch|Point]' the first or last epoch (point) is chosen with all 
%    components and points (epochs) observed. When this is not possible, the epoch (point) is selected with
%    most observed points (epochs).
%    When a specific epoch (point) is selected it must be an epoch with all components available and
%    all points (epochs) observed for that epoch (point).

% Duplicate output to file, and catch errors, and start timing

try

[~,outputfileroot]=fileparts(outputfilename);
diary([ outputfileroot '_' datestr(now,30) '.log' ])

fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Check the options and if necessary overwrite the default values

%[inputfilenames,outputfilename,opt]= ...
%    stmcheckarguments(inputfilenames,outputfilename,opt,varargin);
 [inputfilenames,outputfilename,opt]= ...
    stmcheckarguments(inputfilenames,outputfilename,opt,options,flag);
if isempty(outputfilename)
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    return;
end

if ~strcmpi(opt.covmatrix,{'chol';'byobstype';'full'})
    error(['Illegal value for opt.covmatrix ' opt.covmatrix ])
end

%% Load space time matrix datasets 

% Load space time datasets 

datasets=[];
for k=1:numel(inputfilenames)
   datasets{k}=stmread(inputfilenames{k});
end

% Select the techniques to process

techniqueIds=cellfun(@(x) x.techniqueId,datasets,'UniformOutput',false);
datasetIds=cellfun(@(x) x.datasetId,datasets,'UniformOutput',false);
selected=ismember(techniqueIds,opt.selectTechnique) & contains(datasetIds,opt.selectDataset);

datasets(~selected)=[];

numDatasets=numel(datasets);

%% Print overview of datasets, with datasetId, TechniqueId and projectInfo
    
fprintf('\nDatasetId                              TechId createdBy       creationDate    softwareName    projectId       projectFile (projectFileDate)\n')
fprintf(  '-------------------------------------  ------ --------------- --------------- --------------- --------------- -------------------------------\n')
for k=1:numDatasets
   datasetId=datasets{k}.datasetId;
   if length(datasetId) > 38
      datasetId=[datasetId(1:35) '...'];
   end
   datasetAttrib=datasets{k}.datasetAttrib;
   fprintf('%-38s %-5s  %-15s %-15s %-15s %-15s %s (%s)\n',datasetId,datasets{k}.techniqueId, ...
          datasetAttrib.createdBy,datasetAttrib.creationDate,datasetAttrib.softwareName,datasetAttrib.projectId,datasetAttrib.projectFile,datasetAttrib.projectFileDate);
end
fprintf('\n')

%%  Make arrays pntIds and epochIds with the output point and epoch identifiers
%
% Each input dataset may contains its own point coordinates and epoch
% times, which we need to link. The first step in this process is to make
% array with unique point and epoch Ids, point coordinates and epoch
% times:
%
%   pntIds       - cell column array with the point id's
%   pntCrd       - numeric matrix with the geographic point coordinates (deg/m)
%   pntNeu       - numeric matrix with the topocentric coordinates (km/m)
%   epochIds     - cell row array with the epoch id's 
%   epochDyears  - row vector with the decimal years 
%
% The pntIds and/or epochIds are build from the coordinates and/or epoch
% times, unless the input datasets contain attributes pntId and/or
% epochId, which alreay contain the harmonized point and epoch Ids that
% must agree (be the same) between the individual input datasets. This
% condition must be checked, and the result is given by two logical called
% haveHarmonizedPoints and haveHarmonizedEpochs.

% Loop over the datasets and collect information for building the point and epoch Id arrays 

pntIds=[];
pntCrd=[];
epochIds=[];
epochDyears=[];
inputDatasetAttrib=[];

haveHarmonizedPoints=true;
haveHarmonizedEpochs=true;

for k=1:numDatasets  
    % build concatenated list of points
    pntAttrib=datasets{k}.pntAttrib;
    if isfield(pntAttrib,'pntId')
       pntIds=[ pntIds ; pntAttrib.pntId ]; 
    else
       pntIds=[ pntIds ; datasets{k}.pntName ]; 
       haveHarmonizedPoints=false;
    end
    pntCrd=[ pntCrd ; datasets{k}.pntCrd ];
    % build concatenated list of epochs
    epochAttrib=datasets{k}.epochAttrib;
    if isfield(epochAttrib,'epochId')
       epochIds=[ epochIds ; epochAttrib.epochId(:) ]; 
    else
       haveHarmonizedEpochs=false;
    end
    epochDyears=[epochDyears ; datasets{k}.epochDyear(:)]; 
    % collect input dataset attributes (needed for harmonization integrity)
    inputDatasetAttrib{k}=datasets{k}.datasetAttrib;
end

% check if input pntId and pntEpoch attributes come from the same source

projectId=unique(cellfun(@(x) x.projectId,inputDatasetAttrib,'UniformOutput',false));
projectFile=unique(cellfun(@(x) x.projectFile,inputDatasetAttrib,'UniformOutput',false));
projectFileDate=unique(cellfun(@(x) x.projectFileDate,inputDatasetAttrib,'UniformOutput',false));

if numel(projectId) ~= 1 || numel(projectFile) ~= 1 || numel(projectFileDate) ~=1 
   warning('Input datasets come from different projects, proceed with caution...')
   if haveHarmonizedPoints && ~opt.ignoreProjectCheck
       warning('PntId attributes come from different projects, will be recomputed, proceed with caution...')
       haveHarmonizedPoints=false;
   end
   if haveHarmonizedEpochs && ~opt.ignoreProjectCheck
       warning('EpochId attributes come from different projects, will be recomputed, proceed with caution...')
       haveHarmonizedEpochs=false;
   end
end
projectId=char(projectId{1});
projectFile=char(projectFile{1});
projectFileDate=char(projectFileDate{1});

% Convert geographic coordinates in pntCrd into topocentric coordinates in pntNeu 

switch lower(opt.pntCrdType)
    case 'deg/m'
       % convert latitude/longitude into local topocentric coordinates
       [pntNeuDs,plh0] = plh2neusp(pntCrd);     % deg/deg/m -> m/m/m
       pntNeuDs(:,1:2)=pntNeuDs(:,1:2)./1000;   % m/m/m -> km/km/m
    case 'km/m'
       % coordinates are already in the right units, this is exceptional,
       % and only happens for simulations, just copy
       pntNeuDs=pntCrd;
    otherwise
       error('unknown pntCrdType option')        
end

% Assemble pntId, pntCrd pntCrdTopo arrays for the integrated dataset

if haveHarmonizedPoints 
    %pntIds=unique(pntIds); 
    [pntCrd]=grpstats(pntCrd,pntIds,{'mean'});
    [pntNeu,pntIds]=grpstats(pntNeuDs,pntIds,{'mean' 'gname'});
else
    % find unique points between datasets, but the algorithm is simple and will fail if points in a single dataset are with opt.tolPointDist
    warning('Point in input datasets have not been harmonized, use coordinates to find matching points, but you are warned...')
    clusterIdx=clusterdata(pntNeuDs(:,1:2),'cutoff',opt.tolPointDist,'Criterion','distance');
    [pntCrd]=grpstats(pntCrd,clusterIdx,{'mean'});
    [pntNeu,pntRange,PntCount]=grpstats(pntNeuDs,clusterIdx,{'mean' 'range' 'numel'});
    % assign pntId from one of the datasets, should be replaced by Plus Location identifier
    pntIds2=cell(size(pntNeu,1),1);
    pntIds2(clusterIdx)=pntIds;      
    pntIds=pntIds2;
end
% sort the pntIds
[pntIds,idx]=sort(pntIds);
pntCrd=pntCrd(idx,:);
pntNeu=pntNeu(idx,:);

% Assemble epochId and epochDyears arrays for the integrated dataset

if haveHarmonizedEpochs
    %epochIds=unique(epochIds);
    [epochDyears,epochIds]=grpstats(epochDyears,epochIds,{'mean' 'gname'});
    epochDyears=epochDyears';
    epochIds=epochIds';
else
    % epochId is missing, use epochDyear to find out common epochs
    warning('Epochs in input datasets have not been harmonized, we will give it a try now ...')
    dyears=sort(epochDyears);
    idx=[ find(diff(dyears) > opt.tolEpoch); numel(dyears) ];
    epochDyears=nan(1,numel(idx));
    epochIds=cell(1,numel(idx));
    l2=0;
    for l=1:numel(idx)
       l1=l2+1;
       l2=idx(l);
       epochDyears(l)=mean(dyears(l1:l2));
       epochIds{l}=sprintf('%8.3f',epochDyears(l));
    end
end
% sort the epochs
[epochDyears,idx]=sort(epochDyears);
epochIds=epochIds(idx);

% Count the number of points and epochs

numPoints=numel(pntIds);
numEpochs=numel(epochIds);

% % Print overview of datasets, number of points and epochs
% 
% fprintf('\nDatasetId                              TechId  Pnts Epochs  Dim\n')
% fprintf(  '-------------------------------------  ------ ----- ------ ----\n')
% for k=1:numDatasets
%    fprintf('%-38s %-5s %6d %6d  %3d\n',datasets{k}.datasetId,datasets{k}.techniqueId, ...
%      datasets{k}.numPoints,datasets{k}.numEpochs,numel(datasets{k}.obsTypes));
% end
% 
% fprintf('                                             ------  -----\n')
% fprintf('Distinct                                     %6d %6d\n\n',numPoints,numEpochs)

%% Create intermediate data holding, add parameter index arrays and count observations
%
% In this step the points and epochs in the input datasets are linked to
% the arrays pntId and epochId created in the previous section. This is
% done by creating two index arrays and adding them to the datasets.
%
% We don't store intermediate product in datasets{k}, as these can be objects, 
% and we don't want to update those. Instead all intermediate results for the 
% integration, including the index arrays, are stored in dslocal. With
%
%    dslocal{k}.datasetId      datasetId
%    dslocal{k}.techniqueId    techniqueId
%    dslocal{k}.numPoints      number of points in dataset
%    dslocal{k}.numEpochs      number of epochs in dataset
%    dslocal{k}.numObsTypes    number of observation types in dataset
%    dslocal{k}.numObs         number of observations in dataset (not including missing observations)
%    dslocal{k}.pntNeu         point coordinates in local topocentric system (km/m)
%    dslocal{k}.idxPnt         index array to pntId and pntCrd
%    dslocal{k}.idxEpoch       index array to epochId and epochDyear 
%    dslocal{k}.obsTypes       observation types 
%    dslocal{k}.stochModel     stochastic model type 
%
% dslocal{k}.idxPnt and dslocal{k}.idxEpoch have the dimension of the observation 
% dataset, and do not contain any zeros, and provide pointers to pntIds  and
% epochDyears . The observation types in dslocal{k} are slightly different
% from the name in the input datasets, instead of 'los', we distinquish
% between ascending and descending tracks ('asclos' and 'dsclos'). This is
% something which we also could move into the input datasets (discuss with freek)
%
% dslocal does not contain the pntAttrib, epochAttrib, obsData, stochData
% and sensitivityMatrix. These have to be retrieved from datasets.
%
% dslocal with also contain a reference to intermediate results from the 
% least squares processing, and residuals, which will be added later.
%
% At the same time we will build a matrix SO that hold information on the
% type of observations for each point and epoch. SO is a matrix with dimension 
% [numPoints,numEpochs,5] and contains counts of the observations in East, 
% North, Up, Asc and Dsc directions, with names set in expObsTypes.
% The matrix SO with observation counts is used later to classify points
% and displacement parameters.
%
% We will also create obsTypes cell array with actually used observation types.

dslocal=cell(numDatasets, 1);

expObsTypes={'North' 'East' 'Up' 'asclos' 'dsclos' };
SO=zeros(numPoints,numEpochs,numel(expObsTypes));
obsTypes=[];

k2=0;
for k=1:numDatasets
    k1=k2+1;
    k2=k2+datasets{k}.numPoints;
    techniqueId=datasets{k}.techniqueId;             % techniqueId
    datasetId=datasets{k}.datasetId;                 % datasetId
    if length(datasetId) > 38
        datasetId=[datasetId(1:35) '...'];
    end
    dslocal{k}.datasetId=datasetId;
    dslocal{k}.techniqueId=datasets{k}.techniqueId;  % techniqueId
    dslocal{k}.numPoints=datasets{k}.numPoints;      % Number of points in dataset
    dslocal{k}.numEpochs=datasets{k}.numEpochs;      % Number of points in dataset
    dslocal{k}.epochDyear=datasets{k}.epochDyear;    % Epoch time in decimal years
    dslocal{k}.pntNeu=pntNeuDs(k1:k2,1:2);           % Point coordinates in local topocentric system
    pntAttrib=datasets{k}.pntAttrib;
    epochAttrib=datasets{k}.epochAttrib;
    % If not harmonized, build pntId for the dataset, and add to pntAttrib
    if ~haveHarmonizedPoints
       [~,idx]=pdist2(pntNeu(:,1:2),pntNeuDs(k1:k2,1:2),'euclidean','smallest',1);
       pntAttrib.pntId=pntIds(idx);
    end
    % Index points
    [~,ia,ib]=intersect(pntAttrib.pntId,pntIds);
    ic=ia;ic(ia)=ib;
    dslocal{k}.idxPnt=ic;
    % If not harmonized, build epochId for the dataset, and add to epochAttrib
    if ~haveHarmonizedEpochs
       [~,idx]=pdist2(epochDyears(:),datasets{k}.epochDyear(:),'euclidean','smallest',1);
       epochAttrib.epochId=epochIds(idx);
    end
    % Index epochs
    [~,ja,jb]=intersect(epochAttrib.epochId,epochIds);
    jc=ja;jc(ja)=jb;
    dslocal{k}.idxEpoch=jc;
    % Observation types in dataset 
    if ( strcmpi(datasets{k}.techniqueId,'insar') && all(mod(pntAttrib.azAngle,360) > 180) )
       dslocal{k}.obsTypes={'asclos'}; 
       obsTypes=[ obsTypes {'asclos'} ]; 
    elseif ( strcmpi(datasets{k}.techniqueId,'insar') && all(mod(pntAttrib.azAngle,360) < 180) )
       dslocal{k}.obsTypes={'dsclos'}; 
       obsTypes=[ obsTypes {'dsclos'} ]; 
    else
       dslocal{k}.obsTypes=datasets{k}.obsTypes; 
       obsTypes=[ obsTypes datasets{k}.obsTypes ]; 
    end
    dslocal{k}.numObsTypes=numel(dslocal{k}.obsTypes);    
    % Optionally modify (overwrite/add) the stochastic models
    stochModel=datasets{k}.stochModel;
    [stochModel, modified] = updStochModel(stochModel,techniqueId,datasetId,dslocal{k}.obsTypes,opt);
    dslocal{k}.stochModel=stochModel;
    % Observation type count matrix (handles missing observations)
    numObsDs=zeros(1,numel(expObsTypes));
    for ll=1:numel(dslocal{k}.obsTypes)
       l=find(ismember(expObsTypes,dslocal{k}.obsTypes{ll}));
       numObsDs(l)=numel(find(~isnan(datasets{k}.obsData(:,:,ll))));
       SO(ic,jc,l)=SO(ic,jc,l)+ ~isnan(datasets{k}.obsData(:,:,ll));
    end    
    dslocal{k}.numObs=numObsDs;
end
obsTypes=unique(obsTypes);

%% Build point names (pntName)
%
% The point name consists of the original marker names and pntId joined
% into a single string seperated by a forward slash

pntName=pntIds;
for k=1:numDatasets
   ic=dslocal{k}.idxPnt;
   tf = strcmp(datasets{k}.pntName,pntIds(ic));
   if any(~tf)
      pntName(ic(~tf))=join([datasets{k}.pntName(~tf) pntName(ic(~tf))],'/');
   end
end

% Remove duplicates (maintaining the order)
for k=1:numel(pntName)
   pntName{k}=strjoin(unique(strsplit(pntName{k},'/'),'stable'),'/');
end

if opt.verbose > 1
   fprintf('Point names:\n\n')
   pntName
end

%% Print overview of datasets, number of points and epochs

numObsTot=zeros(1,numel(expObsTypes));

fprintf('\n                                                                       Observation count per Type\n\n')
fprintf(  'DatasetId                              TechId   Pnts Epochs  Dim    North   East     Up ascLos dscLos    total\n')
fprintf(  '-------------------------------------  ------ ------ ------ ----   ------ ------ ------ ------ ------   ------\n')
for k=1:numDatasets
   fprintf('%-38s %-5s  %6d %6d  %3d   %6d %6d %6d %6d %6d   %6d\n',dslocal{k}.datasetId,dslocal{k}.techniqueId, ...
     dslocal{k}.numPoints,dslocal{k}.numEpochs,numel(dslocal{k}.obsTypes),dslocal{k}.numObs,sum(dslocal{k}.numObs));
   numObsTot=numObsTot+dslocal{k}.numObs;
end

fprintf(  '                                              ------ ------        ------ ------ ------ ------ ------   ------\n')
fprintf('Distinct / total                              %6d %6d        %6d %6d %6d %6d %6d   %6d\n\n',numPoints,numEpochs,numObsTot,sum(numObsTot))


%% Print stochastic models

if opt.verbose > 0
   fprintf('Stochastic models:\n\n')
   for k=1:numel(dslocal)
      printStochModel(dslocal{k})
   end
end

%% Ascii graphic / plot of the time line 

% Ascii graphic of time line

fprintf('\n')                                                     
ctmp=char(epochIds)';
tmpline=blanks(numEpochs*2);
for l=1:size(ctmp,1)-1
   tmpline(1:2:numEpochs*2)=ctmp(l,:);
   fprintf(  '                                       %s\n',tmpline);
end
tmpline(1:2:numEpochs*2)=ctmp(end,:);
tmpcount=zeros(numEpochs,1);

fprintf(  'DatasetId                              %s\n',tmpline);
fprintf(  '-------------------------------------  %s\n',repmat('- ',[1,numEpochs]))
for k=1:numDatasets
   tmpline=blanks(numEpochs*2);
   tmpline(dslocal{k}.idxEpoch*2-1)='*';
   tmpcount(dslocal{k}.idxEpoch)=tmpcount(dslocal{k}.idxEpoch)+1;
   fprintf('%-38s %s\n',dslocal{k}.datasetId,tmpline)
end
fprintf(  '                                       %s\n',repmat('- ',[1,numEpochs]))
ctmp=num2str(tmpcount,'%-d')';
for l=1:size(ctmp,1)
   tmpline(1:2:numEpochs*2)=ctmp(l,:);
   fprintf(  '                                       %s\n',tmpline);
end
fprintf('\n')

% Ascii Histogram of number of common epochs/points in datasets 

tmpcounte=zeros(numDatasets,numEpochs);
tmpcountp=zeros(numDatasets,numPoints);
for k=1:numDatasets
   tmpcounte(k,dslocal{k}.idxEpoch)=tmpcounte(k,dslocal{k}.idxEpoch)+1;
   tmpcountp(k,dslocal{k}.idxPnt)=tmpcountp(k,dslocal{k}.idxPnt)+1;
end

hce=histcounts(sum(tmpcounte,1),0.5:numDatasets+0.5);
hcp=histcounts(sum(tmpcountp,1),0.5:numDatasets+0.5);

fprintf('Number of occurences of epochs/points in datasets:\n\n')
fprintf('##  Epochs Points   (## is number of datasets)\n')
fprintf('--  ------ ------\n')
for k=1:numDatasets
   if ( hce(k) > 0 || hcp(k) > 0 )
      fprintf('%2d  %6d %6d\n',k,hce(k),hcp(k)) 
   end
end
fprintf('\n')

% Plot the time line and point distribution (optional)

if opt.doplots > 0

    datasetIds=cellfun(@(x) x.datasetId,dslocal,'UniformOutput',false);

    figure('name','Timeline');
    spy(tmpcounte)
    h=gca;
    h.TickLabelInterpreter='none';
    h.YTick=1:numDatasets;
    h.YTickLabel=datasetIds;
    h.XTick=1:numEpochs;
    h.XTickLabel=epochIds;
    h.XTickLabelRotation=90;
    hold on
    stairs(0.5:numEpochs+1,numDatasets-[sum(tmpcounte) sum(tmpcounte(:,end)) ] +1 );
    title('Time line')

    figure('name','Point Distribution');
    spy(tmpcountp')
    h=gca;
    h.TickLabelInterpreter='none';
    h.DataAspectRatio=[ 1  10^floor(log10(numPoints/5))*2 1];
    h.XTick=1:numDatasets;
    h.XTickLabel=datasetIds;
    h.XTickLabelRotation=45;
    h.YTick=1:10^floor(log10(numPoints/5)):numPoints;
    h.YTickLabel=pntIds(h.YTick);
    hold on
    stairs([sum(tmpcountp) sum(tmpcountp(:,end)) ] ,0.5:numPoints+1);
    title('Point Distribution')

end

%% Analyze the observation structure and determine the parameter types
%
% SO is an integer space-time matrix with the observation count for the 
% expected observation types (expObsTypes) {'North' 'East' 'Up' 'asclos'
% 'dsclos' }, that was created in the previous section. We use it here
% to classify the displacement parameters and points.
%
% SP.<partype> is derived from SO, but converted to extended parameters types 
% {'North' 'East' 'Up' 'nearEast', 'nearUp', 'losAsc' 'losDsc'}. Each
% SP.<partype> is a two dimensional space time matrix. It contains
% the number of observations on each extended parameter type. These counts
% are not necessarily integer, as some observations have to be "divided"
% over multiple extended parameters.
%
% SPmat has the same info as SP.<partype>, but it is a 3-dimensional space
% time matrix with in the 3rd dimension the <partype> in the order given by
% extParTypes ({'North' 'East' 'Up' 'nearEast', 'nearUp', 'losAsc'
% 'losDsc'}). SPmat is a regular 3 dimensional Matlab matrix.
%
% The matrix SC is a numPoints x numEpochs x numNomParTypes space time matrix 
% with a displacement parameter code 
%
%   0 - Not observed
%   1 - Not observed, but necessary for the other components (must have a-priori value)
%   2 - losAsc parameter, depends on a-priori values for two other components
%   3 - losDsc parameter, depends on a-priori values for two other components
%   4 - nearUp or nearEast parameter, depends on a-priori values for one other component
%   7 - Fully estimable parameter
%
% for the nominal parameter types {'North' 'East' 'Up'}.

% Initialize the SP.<partype> counts with the observation counts, they will
% be updated subsequently

SP.North=SO(:,:,1);        
SP.East=SO(:,:,2);
SP.Up=SO(:,:,3);
SP.nearEast=zeros(numPoints,numEpochs);
SP.nearUp=zeros(numPoints,numEpochs);
SP.losAsc=SO(:,:,4);
SP.losDsc=SO(:,:,5);

% Find out how many nearEast and nearUp parameters we might have from
% the InSar observations. We need at least to have ascending AND descending
% observations, if so, divide the InSAR observations over nearEast and NearUp.

bSP=SP.losAsc > 0 & SP.losDsc > 0;                    % true if both ascending and descending observations
SP.nearEast(bSP)=min(SP.losAsc(bSP),SP.losDsc(bSP));  % divide over nearEast and nearUp  
SP.nearUp(bSP)=SP.losAsc(bSP)+SP.losDsc(bSP)-SP.nearEast(bSP);
SP.losAsc(bSP)=0;                                     % remove them from the los counts
SP.losDsc(bSP)=0;         

% When the East and/or Up component has been observed the (remaining) line of 
% sight observations can be also be counted as East, nearEast, Up and/or nearUp.

bSP=SP.East > 0 & SP.Up > 0;           
SP.nearEast(bSP)=SP.nearEast(bSP)+(SP.losAsc(bSP)+SP.losDsc(bSP))*1/3; % Divide observations
SP.nearUp(bSP)=SP.nearUp(bSP)+(SP.losAsc(bSP)+SP.losDsc(bSP))*2/3;
SP.losAsc(bSP)=0;                                     % remove them from the los counts
SP.losDsc(bSP)=0;         

bSP=SP.East > 0 & SP.Up == 0;           
SP.nearUp(bSP)=SP.nearUp(bSP)+(SP.losAsc(bSP)+SP.losDsc(bSP));% Count as up observations 
SP.losAsc(bSP)=0;                                     % remove them from the los counts
SP.losDsc(bSP)=0;         

bSP= SP.East == 0 & SP.Up > 0;           
SP.nearEast(bSP)=SP.nearEast(bSP)+(SP.losAsc(bSP)+SP.losDsc(bSP));% Count as East observations 
SP.losAsc(bSP)=0;                                     % remove them from the los counts
SP.losDsc(bSP)=0;         

% When (also) the North component has been observed, then nearEast and nearUp can be
% counted as true East and Up components.

bSP=SP.North > 0;           
SP.East(bSP)=SP.East(bSP)+SP.nearEast(bSP);
SP.Up(bSP)=SP.Up(bSP)+SP.nearUp(bSP);
SP.nearEast(bSP)=0;
SP.nearUp(bSP)=0;

% When both the North and East component has been observed, then the line of sight
% observations can be counted partly as true East and partly as Up. The
% division is based on 30 deg incidence angle.

bSP=SP.North > 0 & SP.East > 0;      
SP.East(bSP)=SP.East(bSP)+(SP.losAsc(bSP)+SP.losDsc(bSP))*1/3; % Divide 
SP.Up(bSP)=SP.Up(bSP)+(SP.losAsc(bSP)+SP.losDsc(bSP))*2/3;
SP.losAsc(bSP)=0;                                     % remove them from the los counts
SP.losDsc(bSP)=0;         

% When the East and nearEast, or Up and nearUp, have been observed simultaneously,
% then nearEast and nearUp can be counted also as East and/or Up, and North
% can be set as observed (indirectly)

bSP= SP.Up > 0;
SP.Up(bSP)=SP.Up(bSP)+SP.nearUp(bSP);
SP.nearUp(bSP)=0;    

bSP=SP.East > 0;
SP.East(bSP)=SP.East(bSP)+SP.nearEast(bSP);
SP.nearEast(bSP)=0;    

% Collect the information of the SP structure into a numPoints x numEpochs x numExtParTypes matrix SPmat

extParTypes={'North' 'East' 'Up' 'nearEast', 'nearUp', 'losAsc' 'losDsc'};

SPmat=zeros(numPoints,numEpochs,numel(extParTypes));
SPmat(:,:,1)=SP.North;
SPmat(:,:,2)=SP.East;
SPmat(:,:,3)=SP.Up;
SPmat(:,:,4)=SP.nearEast;
SPmat(:,:,5)=SP.nearUp;
SPmat(:,:,6)=SP.losAsc;
SPmat(:,:,7)=SP.losDsc;

% Make numPoints x numEpochs x numNomParTypes matrix SC with a displacement parameter code
%
%   0 - Not observed
%   1 - Not observed, but necessary for the other components (must have a-priori value)
%   2 - losAsc parameter, depends on a-priori values for two other components
%   3 - losDsc parameter, depends on a-priori values for two other components
%   4 - nearUp or nearEast parameter, depends on a-priori values for one other component
%   7 - Fully estimable parameter

nomParTypes={'North' 'East' 'Up'};
numNomParTypes=numel(nomParTypes);

SC=zeros(numPoints,numEpochs,numel(nomParTypes));
SC(:,:,3)=2*(SP.losAsc > 0) + 3*(SP.losDsc > 0) + 4*(SP.nearUp > 0) + 7*(SP.Up > 0);
SC(:,:,2)=1*( SP.losAsc > 0 | SP.losDsc > 0) + 4*(SP.nearEast > 0) + 7*(SP.East > 0);
SC(:,:,1)=1*( SP.losAsc > 0 | SP.losDsc > 0 | SP.nearEast > 0 | SP.nearUp > 0) + 7*(SP.North > 0);

%% Analyze the point capability

% Aggregate the number of observations and parameters over time (one row per point)
%
%   pntObs...  are numPoints x numObsTypes matrices
%   pntPar...  are numPoints x numExtParTypes matrices

pntObsMax=squeeze(max(SO,[],2));       % Maximum number of observations of this type in a single epoch              
pntObsCount=squeeze(sum(SO>0,2));      % Number of epochs with observations of this type
pntObsFlag=pntObsCount > 0;            % True if there are observations of this type      

pntParMaxSP=squeeze(max(SPmat,[],2));  % Maximum number of extended parameter type contributions in a single epoch
pntParCount=squeeze(sum(SPmat>0,2));   % Number of epochs with extended parameter contributions of this type
pntParFlag=pntParCount > 0;            % True if there are extended parameter contributions of this type

% Aggregate the number of observations and parameters over points (one value per epoch)
%
%   epochObs...  are numPoints x numObsTypes matrices
%   epochPar...  are numPoints x numExtParTypes matrices

epochObsMax=squeeze(max(SO,[],1));     % Maximum number of observations of this type for a single point        
epochObsCount=squeeze(sum(SO>0,1));    % Number of epochs with observations of type               
epochObsFlag=epochObsCount > 0;        % True if there are observations of this type           

epochParMaxSP=squeeze(max(SPmat,[],1));% Maximum number of extended parameter type contributions for a single point
epochParCount=squeeze(sum(SPmat>0,1)); % Number of points with extended parameter contributions of this type
epochParFlag=epochParCount > 0;        % True if there are extended parameter contributions of this type

% Determine point and epoch capabilities (binary coded)
%
%   pntSC      is a numPoints x numNomParType matrix derived from SC
%   epochSC    is a numEpochs x numNomParType matrix derived from SC
%
%     bit     set when at least one epoch is       
%      1      not observed
%      2      not observed, but necessary for the other components (must have a-priori value)
%      3      losAsc parameter, depends on a-priori values for two other components
%      4      losDsc parameter, depends on a-priori values for two other components
%      5      nearUp or nearEast parameter, depends on a-priori values for one other component
%      6      fully estimable parameter
%
% To evaluate the bit values use bitget(pntSC,bit) or dec2bin(pntSC)

pntSC =    1*(squeeze(any(SC == 1,2))) + 2*(squeeze(any(SC == 1,2))) + ...
           3*(squeeze(any(SC == 2,2))) + 8*(squeeze(any(SC == 3,2))) + ...
          16*(squeeze(any(SC == 4,2)))+ 32*(squeeze(any(SC == 7,2)));
epochSC =  1*(squeeze(any(SC == 1,1))) + 2*(squeeze(any(SC == 1,1))) + ...
           3*(squeeze(any(SC == 2,1))) + 8*(squeeze(any(SC == 3,1))) + ...
          16*(squeeze(any(SC == 4,1)))+ 32*(squeeze(any(SC == 7,1)));

% Make cell string arrays with point capabilities      
      
defCapability = {'x' 'no' 'under' 'asc' 'dsc' 'near' 'full'};
      
pntNorthCapability=join(defCapability([ 1+bitget(pntSC(:,1),1) 1+2*bitget(pntSC(:,1),2) 1+3*bitget(pntSC(:,1),3) 1+4*bitget(pntSC(:,1),4) 1+5*bitget(pntSC(:,1),5) 1+6*bitget(pntSC(:,1),6)]),'/');
pntNorthCapability=strrep(pntNorthCapability,'x/','');
pntNorthCapability=join([ pntNorthCapability, repmat({'North'},[numPoints,1])],'');

pntEastCapability=join(defCapability([ 1+bitget(pntSC(:,2),1) 1+2*bitget(pntSC(:,2),2) 1+3*bitget(pntSC(:,3),2) 1+4*bitget(pntSC(:,2),4) 1+5*bitget(pntSC(:,2),5) 1+6*bitget(pntSC(:,2),6)]),'/');
pntEastCapability=strrep(pntEastCapability,'x/','');
pntEastCapability=join([ pntEastCapability, repmat({'East'},[numPoints,1])],'');

pntUpCapability=join(defCapability([ 1+bitget(pntSC(:,3),1) 1+2*bitget(pntSC(:,3),2) 1+3*bitget(pntSC(:,3),3) 1+4*bitget(pntSC(:,3),4) 1+5*bitget(pntSC(:,3),5) 1+6*bitget(pntSC(:,3),6)]),'/');
pntUpCapability=strrep(pntUpCapability,'x/','');
pntUpCapability=join([ pntUpCapability, repmat({'Up'},[numPoints,1])],'');
      
pntCapability = join([ pntNorthCapability pntEastCapability pntUpCapability ],','); 
      
% Compute the point dimension (counts if bits 3, 4 or higher are set)
      
pntDim=sum(pntSC > 2^2,2);

% Print the results

if opt.verbose > 2 || ( opt.verbose > 1 && numPoints <= 200 ) || ( opt.verbose > 0 && numPoints <= 50 )

   fprintf('Maximum number of observations in a single epoch (for points):\n\n')
   fprintf('pntName                                  N    E    U  asc  dsc  ->     N     E     U   nearE nearU    asc  dsc   dim   Point capability\n')
   for i=1:numPoints
      fprintf('%-38s %3d  %3d  %3d  %3d  %3d     %5.2g %5.2g %5.2g   %5.2g %5.2g    %3d  %3d     %1d   %s\n', ...
          pntName{i}, ...
          pntObsMax(i,:), ...
          pntParMaxSP(i,:), ...
          pntDim(i), ...
          pntCapability{i})
   end
   fprintf('                                       ---  ---  ---  ---  ---       ---   ---   ---     ---   ---    ---  ---\n') 
   fprintf('Average                              %5.1f%5.1f%5.1f%5.1f%5.1f    %6.1f%6.1f%6.1f  %6.1f%6.1f  %5.1f%5.1f\n',sum(pntObsMax)./numPoints,sum(pntParMaxSP)./numPoints)    
   fprintf('Point count                          %5d%5d%5d%5d%5d     %5d %5d %5d   %5d %5d  %5d%5d\n\n', sum(pntObsFlag),sum(pntParFlag))

   fprintf('Number of epochs with observations (for points):\n\n')
   fprintf('pntName                                  N    E    U  asc  dsc  ->     N     E     U   nearE nearU    asc  dsc   dim   Point capability\n')
   for i=1:numPoints
      fprintf('%-38s %3d  %3d  %3d  %3d  %3d     %5.2g %5.2g %5.2g   %5.2g %5.2g    %3d  %3d     %1d   %s\n', ...
          pntName{i}, ...
          pntObsCount(i,:), ...
          pntParCount(i,:), ...
          pntDim(i), ...
          pntCapability{i})
   end
   fprintf('                                       ---  ---  ---  ---  ---       ---   ---   ---     ---   ---    ---  ---\n') 
   fprintf('Average                              %5.1f%5.1f%5.1f%5.1f%5.1f    %6.1f%6.1f%6.1f  %6.1f%6.1f  %5.1f%5.1f\n',sum(pntObsCount)./numPoints,sum(pntParCount)./numPoints)    
   fprintf('Point count                          %5d%5d%5d%5d%5d     %5d %5d %5d   %5d %5d  %5d%5d\n\n', sum(pntObsFlag),sum(pntParFlag))

end

%%
fprintf('Number of points with:\n')
fprintf('- North capability    %6d     fullNorth=%d (%.1f%%)\n',...
    sum(pntParFlag(:,1)),...
    sum(pntParFlag(:,1)),100*sum(pntParCount(:,1))./(numPoints*numEpochs));
fprintf('- East capability     %6d     fullEast=%d (%.1f%%)\n', ...
    sum(pntParFlag(:,2)|pntParFlag(:,4)),...
    sum(pntParFlag(:,2)&~pntParFlag(:,4)),100*sum(pntParCount(pntParFlag(:,2)&~pntParFlag(:,4),2))./(numPoints*numEpochs) )
fprintf('                                 nearEast=%d (%.1f%%), mixed full/nearEast=%d (%.1f%%/%.1f%%)\n',...
    sum(~pntParFlag(:,2)&pntParFlag(:,4)),100*sum(pntParCount(~pntParFlag(:,2)&pntParFlag(:,4),4))./(numPoints*numEpochs), ...
    sum(pntParFlag(:,2)&pntParFlag(:,4)),100*sum(pntParCount(pntParFlag(:,2)&pntParFlag(:,4),2))./(numPoints*numEpochs),100*sum(pntParCount(pntParFlag(:,2)&pntParFlag(:,4),4))./(numPoints*numEpochs))
fprintf('- Vertical capability %6d     fullUp=%d (%.1f%%)\n', ...
    sum(pntParFlag(:,3)|pntParFlag(:,5)|pntParFlag(:,6)|pntParFlag(:,7)),...
    sum(pntParFlag(:,3)&~pntParFlag(:,5)),100*sum(pntParCount(pntParFlag(:,3)&~pntParFlag(:,5),3))./(numPoints*numEpochs) )
fprintf('                                 nearUp=%d (%.1f%%), mixed full/nearUp=%d (%.1f%%/%.1f%%),\n',...
    sum(~pntParFlag(:,3)&pntParFlag(:,5)),100*sum(pntParCount(~pntParFlag(:,3)&pntParFlag(:,5),5))./(numPoints*numEpochs), ...
    sum(pntParFlag(:,3)&pntParFlag(:,5)),100*sum(pntParCount(pntParFlag(:,3)&pntParFlag(:,5),3))./(numPoints*numEpochs),100*sum(pntParCount(pntParFlag(:,3)&pntParFlag(:,5),5))./(numPoints*numEpochs))
fprintf('                                 ascLos=%d (%.1f%%), mixed full/ascLos=%d (%.1f%%/%.1f%%), mixed nearUp/ascLos=%d (%.1f%%/%.1f%%),\n',...
    sum(~pntParFlag(:,3)&~pntParFlag(:,5)& pntParFlag(:,6)&~pntParFlag(:,7)),100*sum(pntParCount(~pntParFlag(:,3)&~pntParFlag(:,5)& pntParFlag(:,6)&~pntParFlag(:,7),6))./(numPoints*numEpochs), ...
    sum( pntParFlag(:,3)&~pntParFlag(:,5)& pntParFlag(:,6)&~pntParFlag(:,7)),100*sum(pntParCount( pntParFlag(:,3)&~pntParFlag(:,5)& pntParFlag(:,6)&~pntParFlag(:,7),3))./(numPoints*numEpochs),100*sum(pntParCount( pntParFlag(:,3)&~pntParFlag(:,5)& pntParFlag(:,6)&~pntParFlag(:,7),6))./(numPoints*numEpochs), ...
    sum(~pntParFlag(:,3)& pntParFlag(:,5)& pntParFlag(:,6)&~pntParFlag(:,7)),100*sum(pntParCount(~pntParFlag(:,3)& pntParFlag(:,5)& pntParFlag(:,6)&~pntParFlag(:,7),5))./(numPoints*numEpochs),100*sum(pntParCount(~pntParFlag(:,3)& pntParFlag(:,5)& pntParFlag(:,6)&~pntParFlag(:,7),6))./(numPoints*numEpochs) ) 
fprintf('                                 dscLos=%d (%.1f%%), mixed full/dscLos=%d (%.1f%%/%.1f%%), mixed nearUp/dscLos=%d (%.1f%%/%.1f%%), mixed ascLos/dscLos=%d (%.1f%%/%.1f%%)\n',...
    sum(~pntParFlag(:,3)&~pntParFlag(:,5)&~pntParFlag(:,6)& pntParFlag(:,7)),100*sum(pntParCount(~pntParFlag(:,3)&~pntParFlag(:,5)&~pntParFlag(:,6)& pntParFlag(:,7),7))./(numPoints*numEpochs), ...
    sum( pntParFlag(:,3)&~pntParFlag(:,5)&~pntParFlag(:,6)& pntParFlag(:,7)),100*sum(pntParCount( pntParFlag(:,3)&~pntParFlag(:,5)&~pntParFlag(:,6)& pntParFlag(:,7),3))./(numPoints*numEpochs),100*sum(pntParCount( pntParFlag(:,3)&~pntParFlag(:,5)&~pntParFlag(:,6)& pntParFlag(:,7),7))./(numPoints*numEpochs), ...
    sum(~pntParFlag(:,3)& pntParFlag(:,5)&~pntParFlag(:,6)& pntParFlag(:,7)),100*sum(pntParCount(~pntParFlag(:,3)& pntParFlag(:,5)&~pntParFlag(:,6)& pntParFlag(:,7),5))./(numPoints*numEpochs),100*sum(pntParCount(~pntParFlag(:,3)& pntParFlag(:,5)&~pntParFlag(:,6)& pntParFlag(:,7),7))./(numPoints*numEpochs), ... 
    sum(~pntParFlag(:,3)&~pntParFlag(:,5)& pntParFlag(:,6)& pntParFlag(:,7)),100*sum(pntParCount(~pntParFlag(:,3)&~pntParFlag(:,5)& pntParFlag(:,6)& pntParFlag(:,7),6))./(numPoints*numEpochs),100*sum(pntParCount(~pntParFlag(:,3)&~pntParFlag(:,5)& pntParFlag(:,6)& pntParFlag(:,7),7))./(numPoints*numEpochs)); 
fprintf('\n\n');

   
% %pntParPerc=100*sum(pntParCount)./sum(pntParFlag)./numEpochs;
% pntParPerc=100*sum(pntParCount)./(numPoints*numEpochs);
%
% fprintf('Number of points with:\n')
% fprintf('- North capability    %6d     fullNorth=%d (%.1f%%)\n',...
%     sum(pntParFlag(:,1)),...
%     sum(pntParFlag(:,1)),pntParPerc(:,1) );
% fprintf('- East capability     %6d     fullEast=%d (%.1f%%), nearEast=%d (%.1f%%)\n',...
%     sum(pntParFlag(:,2)|pntParFlag(:,4)),...
%     sum(pntParFlag(:,2)),pntParPerc(:,2), ...
%     sum(pntParFlag(:,4)),pntParPerc(:,4) );
% fprintf('- Vertical capability %6d     fullUp=%d (%.1f%%), nearUp=%d (%.1f%%), ascLos=%d (%.1f%%), dscLos=%d (%.1f%%)\n',...
%     sum(pntParFlag(:,3)|pntParFlag(:,5)|pntParFlag(:,6)|pntParFlag(:,7)),...
%     sum(pntParFlag(:,3)),pntParPerc(:,3), ...
%     sum(pntParFlag(:,5)),pntParPerc(:,5),...
%     sum(pntParFlag(:,6)),pntParPerc(:,6), ...
%     sum(pntParFlag(:,7)),pntParPerc(:,7) );
% fprintf('\n\n');
%
% f=zeros(6,6,3);
% for k=1:6
%    for l=1:k
%        for i=1:3
%            f(k,l,i)=sum(bitget(pntSC(:,i),k) & bitget(pntSC(:,i),l));
%        end
%    end
% end
% 
% fprintf('Number of points with:\n')
% fprintf('- North capability    %6d\n',sum(pntParFlag(:,1)))
% fprintf('- East capability     %6d     (East=%d, nearEast=%d)\n',sum(pntParFlag(:,2))+sum(pntParFlag(:,4)),sum(pntParFlag(:,2)),sum(pntParFlag(:,4)))
% fprintf('- Vertical capability %6d     (Up=%d, nearUp=%d, ascLos=%d, dscLos=%d)\n',sum(pntParFlag(:,3))+sum(pntParFlag(:,5))+sum(pntParFlag(:,6))+sum(pntParFlag(:,7)),...
%     sum(pntParFlag(:,3)),sum(pntParFlag(:,5)),sum(pntParFlag(:,6)),sum(pntParFlag(:,7)))
% fprintf('\n\n');
% 
% fprintf('Number of points with:\n')
% fprintf('- Full 3D capability                %6d\n',sum( pntParFlag(:,1) & pntParFlag(:,2) & pntParFlag(:,3) ))
% fprintf('- 2D horizontal only capability     %6d\n',sum( pntParFlag(:,1) & pntParFlag(:,2) & ~( pntParFlag(:,3) | pntParFlag(:,5) | pntParFlag(:,6) | pntParFlag(:,7) ) ))
% fprintf('- 2D nearEast/Up only capability    %6d\n',sum( pntParFlag(:,4) & pntParFlag(:,3) ));
% fprintf('- 2D nearEast/nearUp only capability%6d\n',sum( pntParFlag(:,4) & pntParFlag(:,5) ));
% fprintf('- 1D Up only capability             %6d\n',sum( pntParFlag(:,3) & ~( pntParFlag(:,1) | pntParFlag(:,2) | pntParFlag(:,4) | pntParFlag(:,5) ) ));
% fprintf('- 1D Los only capability            %6d\n',sum( pntParFlag(:,6) | pntParFlag(:,7) ));
% fprintf('\n\n');

%% Determine the output dimension and parameters types
%
% There is still something to think over here, what do we do if we have
% for instance both 'Up' and 'nearUp' types for the same point, do
% we introduce extra parameter types?
%
% The present structure is quite complicated. Maybe better to use fixed
% positions ( North [East|nearEast] [Up|nearUp|losUp] ), and use a flag
% for each point to determine the subtype (true[|near[|los]] ).
%
% Or use six standard positions ( North East Up nearEast nearUp losUp ), 
% and reduce later if some components are unused. This solution would
% allow East|nearEast and Up|nearUp|losUp for the same point (e.g. when
% some observations are missing at certains epochs (e.g. GPS starts later).

%ndim=max(pntDim);

pnt3D = pntParFlag(:,1) & pntParFlag(:,2) & pntParFlag(:,3);
pnt2D = pntParFlag(:,1) & pntParFlag(:,2) & ~( pntParFlag(:,3) | pntParFlag(:,5) | pntParFlag(:,6) | pntParFlag(:,7) );
pnt2DnearE = pntParFlag(:,4) & pntParFlag(:,3);
pnt2DnearEU = pntParFlag(:,4) & pntParFlag(:,5);
pnt1D = pntParFlag(:,3) & ~( pntParFlag(:,1) | pntParFlag(:,2) | pntParFlag(:,4) | pntParFlag(:,5) );
pnt1Dlos = pntParFlag(:,6) | pntParFlag(:,7);

if any(pnt3D)
   ndim=3;
   parTypes={'North','East','Up'};
elseif any(pnt2D) && ~(any(pnt2DnearE) || any(pnt2DnearEU) )
   ndim=2;
   parTypes={'North','East'};
elseif any(pnt2DnearE) && ~any(pnt2DnearEU)
   ndim=2;
   parTypes={'nearEast','Up'};
elseif any(pnt2DnearEU)
   ndim=2;
   parTypes={'nearEast','nearUp'};
elseif any(pnt1D)
   ndim=1;
   parTypes={'Up'};
elseif any(pnt1Dlos)
   ndim=1;
   parTypes={'losUp'};
else
   error('Hey, there must be something wrong with determining the parameter types.')    
end

fprintf('Output spacetime matrix dimension   [ %d x %d x %d ]\n\n',numPoints,numEpochs,ndim)

fprintf('Parameter types: %s\n\n',strjoin(parTypes));

%% Count number of observations, parameters and rank-defects
%
% Another structure that is created is dsstat, this contains counts of the
% observations in array format:
%
%    dsstat.numPoints    numDatasets x 1 matrix
%    dsstat.numEpochs    numDatasets x 1 matrix
%    dsstat.numObsTypes  numDatasets x 1 matrix
%    dsstat.numObs       numDatasets x 1 matrix

dsstat.numPoints=cellfun(@(x) x.numPoints,dslocal,'UniformOutput',true);
dsstat.numEpochs=cellfun(@(x) x.numEpochs,dslocal,'UniformOutput',true);
dsstat.numObsTypes=cellfun(@(x) x.numObsTypes,dslocal,'UniformOutput',true);
dsstat.numObs=cellfun(@(x) sum(x.numObs),dslocal,'UniformOutput',true);

% Number of observations

numObs=sum(dsstat.numObs);
maxObs=sum(dsstat.numPoints.*dsstat.numEpochs.*dsstat.numObsTypes);

% Number of displacement, transformation and offset parameters

hasdispl= ( SC > 1 );
numDisplPerEpoch=squeeze(sum(hasdispl,1));
numDisplPerPoint=squeeze(sum(hasdispl,2));

numParDispl=sum(sum(pntParCount));
numParTransf=sum(dsstat.numEpochs.*dsstat.numObsTypes);
numParOffset=sum(dsstat.numPoints.*dsstat.numObsTypes);
numRankDefect1=sum(dsstat.numObsTypes);

maxParDispl=numPoints*numEpochs*ndim;

numPar=numParDispl+numParTransf+numParOffset;

% Number of rank defect in displacements (computing base) and redundancy
%numRankDefect2=numEpochs*ndim+sum(pntDim)-ndim;
numoffset=sum(numDisplPerPoint(:)>0);
numtransform=sum(numDisplPerEpoch(:)>0);
numRankDefect2=numtransform+numoffset-ndim;
numRedundancy=numObs-numPar+numRankDefect1+numRankDefect2;

% Print

fprintf('\nObservation and parameter count\n')
fprintf('                                               actual     max\n')
fprintf('Observations                                  %7d %7d\n',numObs,maxObs)
fprintf('Parameters                 actual     max     %7d\n',numPar)
fprintf('  Displacements           %7d %7d\n',numParDispl,maxParDispl)
fprintf('  Transformations         %7d\n',numParTransf)
fprintf('  Offsets                 %7d\n',numParOffset)
fprintf('Rank-defect Transf&Offsets                    %7d\n',numRankDefect1)
fprintf('Rank-defect Displacements                     %7d\n',numRankDefect2)
fprintf('                                              -------\n')
fprintf('Reduncancy                                    %7d\n\n',numRedundancy)


%% Analyze the redundancy (1)

RP=SPmat-1;    
RP(RP<0)=0;

% Aggregate the number of observations and parameters over time (one value per point)

pntRedMax=max(RP,[],2);                % Maximum redunancy for a extended parameter type in a single epoch
pntRedCount=squeeze(sum(RP>0,2));      % Number of epochs with redundancy for an extended parameter type
pntRedFlag=any(RP,2);                  % True if there is redundancy of this type 

% Print the redundancy

if opt.verbose > 2 || ( opt.verbose > 1 && numPoints <= 200 ) || ( opt.verbose > 0 && numPoints <= 50 )

   fprintf('Maximum possible redundancy in points (in a single epoch):\n\n')
   fprintf('pntName                                    N     E     U  nearE  nearU    asc  dsc\n')
   for i=1:numPoints
      if ( sum(pntRedMax(i,:)) == 0 ), continue; end
      fprintf('%-38s %5.2g %5.2g %5.2g   %5.2g %5.2g    %3d  %3d\n', ...
          pntName{i}, pntRedMax(i,:))
   end
   fprintf('                                          ---   ---   ---     ---   ---    ---  ---\n') 
   fprintf('Average                                %6.1f%6.1f%6.1f  %6.1f%6.1f  %5d%5d\n',sum(pntRedMax)./numPoints)    
   fprintf('Point count                            %5d %5d %5d   %5d %5d  %5d%5d\n\n', sum(pntRedFlag))
   fprintf('\n\n')
   
   fprintf('Number of epochs with redundancy in points:\n\n')
   fprintf('pntName                                     N     E     U  nearE  nearU    asc  dsc\n')
   for i=1:numPoints
      if ( sum(pntRedCount(i,:)) == 0 ), continue; end
      fprintf('%-38s  %5.2g %5.2g %5.2g   %5.2g %5.2g    %3d  %3d\n', ...
          pntName{i}, pntRedCount(i,:))
   end
   fprintf('                                          ---   ---   ---     ---   ---    ---  ---\n') 
   fprintf('Average                                %6.1f%6.1f%6.1f  %6.1f%6.1f  %5d%5d\n',sum(pntRedCount)./numPoints)    
   fprintf('Point count                            %5d %5d %5d   %5d %5d  %5d%5d\n\n', sum(pntRedFlag))
   fprintf('\n\n')

end

%% Analyze the redundancy (2)

redPnt=squeeze(sum(RP,2));           % redundancy per point
redEpo=squeeze(sum(RP,1))';          % redundacy per epoch
idxRedPnt=find(any(redPnt,2));
idxRedEpo=find(any(redEpo,1));

redTrans=squeeze(max(RP,[],1))';     % max redundancy per epoch
redOffset=squeeze(max(RP,[],2));     % max redundancy per point

% Print

if opt.verbose > 2 || ( opt.verbose > 1 && numPoints <= 200 ) || ( opt.verbose > 0 && numPoints <= 50 )

    fprintf('Number of conditions per point:\n\n')
    fprintf('pntNum pntName                                 ___________ Number of conditions _________  _____ Max conditions per point/epoch _____\n')
    fprintf('                                                    N     E     U nearE nearU   asc   dsc       N     E     U nearE nearU   asc   dsc\n')
    for k=1:numel(idxRedPnt)
      fprintf('%6d %-38s %6.6g%6.6g%6.6g%6.6g%6.6g%6.6g%6.6g  %6.6g%6.6g%6.6g%6.6g%6.6g%6.6g%6.6g\n', ...
          idxRedPnt(k),pntName{idxRedPnt(k)},round(redPnt(idxRedPnt(k),:),1),round(redOffset(idxRedPnt(k),:),1))
    end
    fprintf('\n\n')

end

if opt.verbose > 0

    fprintf('Number of conditions per epoch:\n\n')
    fprintf('epoNum   dyear  ___________ Number of conditions _________  _____ Max conditions per point/epoch _____\n')
    fprintf('                     N     E     U nearE nearU   asc   dsc       N     E     U nearE nearU   asc   dsc\n')
    for k=1:numel(idxRedEpo)
      fprintf('%6d  %6.1f  %6.6g%6.6g%6.6g%6.6g%6.6g%6.6g%6.6g  %6.6g%6.6g%6.6g%6.6g%6.6g%6.6g%6.6g\n', ...
          idxRedEpo(k),epochDyears(idxRedEpo(k)),round(redEpo(:,idxRedEpo(k)),1),round(redTrans(:,idxRedEpo(k)),1))
    end
    fprintf('\n\n')

end

fprintf('Number of ...\n\n')
fprintf('Points with redundant observations %7d    (%d)\n',size(idxRedPnt,1),numPoints)
fprintf('Epochs with redundant observations %7d    (%d)\n\n',size(idxRedEpo,2),numEpochs)

fprintf('                                 N     E     U nearE nearU   asc   dsc\n')
fprintf('                            ------ ----- ----- ----- ----- ----- -----\n')
fprintf('Number of conditions        %6.6g%6.6g%6.6g%6.6g%6.6g%6.6g%6.6g\n',round(squeeze(sum(sum(RP))),1))
fprintf('Average per epoch           %6.6g%6.6g%6.6g%6.6g%6.6g%6.6g%6.6g\n',round(sum(redTrans,2)./numEpochs,1))
fprintf('Average per point           %6.6g%6.6g%6.6g%6.6g%6.6g%6.6g%6.6g\n',round(sum(redOffset)./numPoints,1))


%% Normal matrix analysis of redundancy


%% Create space time matrix with a-priori displacements

aprioriDispl=zeros(numPoints,numEpochs,3);
if ~isempty(opt.aprioriDispl) 
   if exist(opt.aprioriDispl,'file')

      % read space time matrix dataset with a-priori displacements
      
      apriori=stmread(opt.aprioriDispl);
      switch lower(opt.pntCrdType)
         case 'deg/m'
            % convert latitude/longitude into local topocentric coordinates
            tmpNeu = plh2neusp(apriori.pntCrd,plh0);  % deg/deg/m -> m/m/m
            tmpNeu(:,1:2)=tmpNeu(:,1:2)./1000;        % m/m/m -> km/km/m
         case 'km/m'
            % coordinates are already in the right units, this is exceptional,
            % and only happens for simulations, just copy
            tmpNeu=apriori.pntCrd;
         otherwise
            error('unknown pntCrdType option')        
      end
      
      % Interpolate to the points at hand (only in space, we assume the
      % epochs are the same (to be fixed!)
      
      % The following code failts when point are not inside the convex hull
      % or on it's border (ti and bc contain Inf entries). We replace it by
      % Matlabs scatteredInterpolant, but this does not allow error
      % propagation.
      
%       % Create a delaunay Triangulation of the scattered points 
%       DT = delaunayTriangulation(tmpNeu(:,1:2));
%       % Find the triangle that encloses each query point using the pointLocation 
%       % method. In the code below, ti contains the IDs of the enclosing triangles 
%       % and bc contains the barycentric coordinates associated with each triangle.
%       [ti,bc] = pointLocation(DT,pntNeu(:,1:2));
%       % Find the point id's of the vertices surrounding the triangle
%       triIdx=DT(ti,:);
%       % Calculate the sum of the weighted values of V(x,y) using the dot product.
%       for l=1:numEpochs
%          for i=1:3
%             v=apriori.obsData(:,l,i);
%             aprioriDispl(:,l,i)=dot(bc,v(triIdx),2);
%          end
%       end
      
      F=scatteredInterpolant(tmpNeu(:,1:2),tmpNeu(:,3));
      for l=1:numEpochs
         for i=1:3
            F.Values=apriori.obsData(:,l,i);
            aprioriDispl(:,l,i)=F(pntNeu(:,1:2));

         end
      end
      
   else
      error(['Dataset ' opt.aprioriDispl ' with a-priori displacements does not exist.'])
   end
end

%% Support information (overview)
%
%  +----------+               +----------+     SC2 > 1   +----------+     ++
%  |          |+              |          |+              |          |+    ||+
%  |   SC     ||    ------>   |   SC2    ||    ------>   | hasdispl ||    |||
%  |          ||       |      |          ||              |          ||    |||
%  +----------+|       |      +----------+|              +----------+|    ++|
%   +----------+       |       +----------+               +----------+     ++
%     SC-code          |         SC-code                    logical         numDisplPerPoint
%        |             |
%        |  SC > 0     |                                 +----------+
%        v             |             SC-code [012347]    +----------++
%  +----------+        |                                  +----------+
%  |          |+       |             0: no-obs                 numDisplPerEpoch   
%  |  parmask ||    ---+             1: no-obs/needed                     
%  |          ||                     2: losAsc                    
%  +----------+|  <-  epochMask      3: losDsc
%   +----------+  <-  pntMask        4: nearUp/East
%     logical     <- ...             7: full
%
%
%  SC        : SC-code for full set
%  SC2       : SC-code for subset of parameters (parmask applied to SC)
%
%  epochMask : logical array, true for the epochs that will be solved 
%  pntMask   : logical array, true for the points that will be solved
%
%  parmask   : logical mask for the parameters in normal equations (includes SC2 == 1 no-obs/needed parameters)  *), **) 
%  hasdispl  : logical mask for the parameters that can be solved, subset of parmask, with SC2 > 1  **)
%
%  *)  parmask is set from epochMask, pntMask and a lot of other options using the SC-codes to make decisions
%  **) The parameters with SC2 == 1 are initially included when the normal matrix is filled, though they will not be used later,
%      and removed before the Choleski factorization. parmask reflects the situation during the fill, hasdispl what is estimated
%      (including the computation base)
%
%  numDisplPerEpoch : number of displacements per Epoch (can be different for each component)  
%  numDisplPerPoint : number of displacements per Point (can be different for each component)  
%
%  selDim           : logical (numNomParTypes x 1), true for dimensions (components) that can be solved
%  ndim             : number of dimensions (components) that can be solved
%
%% Prepare data masks (filter)
%
% The software uses several masks to filter data
%
% - pntMask        Point mask         numPoints x 1 logical array
% - epochMask      Epoch mask         1 x numEpochs logical array
% - parMask        Point mask         numPoints x numEpochs x ndim logical array
%
% - obsMask *)     Observation mask   embedded in stm.obsData, NaN is false 
%
% True indicates can use, false is do not use. In order to use an observation 
% all conditions needs to be true (and relation)
%
% *) the obsMask is embeded in stm.obsData (stored in datasets{1...numDatasets}).
%    Thus each dataset has its own implicit "obsMask". After implementation of
%    testing we may store an explicit obsMask in dslocal{1...numDatasets}.

% Find points in ROI and epochs in POI

if isempty(opt.ROI)
   pntMask=true(numPoints,1);
else
   roi=bbox2poly(opt.ROI);
   pntMask=inpolygon(pntCrd(:,1),pntCrd(:,2),roi(:,1),roi(:,2)); 
end

epochMask=( epochDyears >= opt.POI(1) & epochDyears <= opt.POI(2) );

% Exclude point and epochs
%
% opt.excludePointIds - Cell array with pointIds to be excluded from the integration
% opt.excludeEpochIds - Cell array with epochIds to be excluded from the integration

[~,ia]=intersect(pntIds,opt.excludePointIds);
if ~isempty(ia)
    if opt.verbose > 0
       fprintf('\nPoints excluded from the integration (%d):\n\n',numel(ia))
       for k=1:numel(ia)
          if pntMask(ia) 
             fprintf('%s\n',pntName{ia(k)})
          else
             fprintf('%s  (outside ROI)\n',pntName{ia(k)})
          end
       end
    end
    pntMask(ia)=false;
end

[~,ia]=intersect(epochIds,opt.excludeEpochIds);
if ~isempty(ia)
    if opt.verbose > 0
       fprintf('\nEpochs excluded from the integration (%d):\n\n',numel(ia))
       for k=1:numel(ia)
          if pntMask(ia) 
             fprintf('%s\n',epochIds{is(k)})
          else
             fprintf('%s  (outside POI)\n',epochIds{ia(k)})
          end
       end
    end
    epochMask(ia)=false;
end


% Set the parameter mask 
%
% - Based on SC code (parMask has same dimension as SC)
%
%   0 - Not observed
%   1 - Not observed, but necessary for the other components (must have a-priori value)
%   2 - losAsc parameter, depends on a-priori values for two other components
%   3 - losDsc parameter, depends on a-priori values for two other components
%   4 - nearUp or nearEast parameter, depends on a-priori values for one other component
%   7 - Fully estimable parameter
%
% - Based on pntMask and epochMask

% Initialize parMask and filter out any parameter that is not observed

parMask=true(numPoints,numEpochs,numel(nomParTypes));
parMask = parMask & SC > 0;

numDispl=sum(sum(parMask));

%  If opt.filterLos true, skip displacements which have only ascending, or only decending, InSAR observations (los qualifier)

if opt.filterLos
   % select only nearUp and fullUp for height
   parMask(:,:,3)=parMask(:,:,3) & SC(:,:,3) > 3;    
   % reset "not observed but necessary" condition for the other two components
   parMask(:,:,2)=parMask(:,:,2) & ~( SC(:,:,2) == 1 & SC(:,:,3) <= 3 );    
   parMask(:,:,1)=parMask(:,:,1) & ~( SC(:,:,1) == 1 & SC(:,:,3) <= 3 );    
end

% If opt.filterNear true, skip displacements which have only InSAR observations (near qualifier)

if opt.filterNear
   % select only fullUp and fullEast
   parMask(:,:,3)=parMask(:,:,3) & SC(:,:,3) > 4;    
   parMask(:,:,2)=parMask(:,:,2) & SC(:,:,2) > 4;    
   % reset "not observed but necessary" condition for the North
   parMask(:,:,1)=parMask(:,:,1) & ~( SC(:,:,1) == 1 & ( SC(:,:,2) <= 4 | SC(:,:,3) <= 4 ) );    
end

% If opt.filterPriori... true, skip displacements which depend on a-priori values
%
% TO DO: These test need to be adjusted, when an a-priori dataset is available,
% only skip when for that displacement no data is available (in combination
% with extrapolate flag).

if opt.filterPrioriEast
   parMask(:,:,3)=parMask(:,:,3) & SC(:,:,2) > 1 ;    
   parMask(:,:,2)=parMask(:,:,2) & SC(:,:,2) > 1 ;    
   parMask(:,:,1)=parMask(:,:,1) & SC(:,:,2) > 1 ;    
end

if opt.filterPrioriNorth
   parMask(:,:,3)=parMask(:,:,3) & SC(:,:,1) > 1 ;    
   parMask(:,:,2)=parMask(:,:,2) & SC(:,:,1) > 1 ;    
   parMask(:,:,1)=parMask(:,:,1) & SC(:,:,1) > 1 ;    
end

% Filter out epochs and points in fewer than opt.minDatasets[Epoch|Point] 

epochMask2=epochMask  & sum(tmpcounte) >= opt.minDatasetsEpoch;
if any(epochMask ~= epochMask2)
   fprintf('\nEpoch mask has changed because there were some epochs with only one dataset ...\n%s\n\n', ...
        strjoin(epochIds( epochMask & ~epochMask2)) )
   epochMask=epochMask2;
end

% Force some points back in

forcedPoints=false(numPoints,1);
if ~strcmp(opt.forceDatasetPoints,'')
   excempt=find(contains(datasetIds,opt.forceDatasetPoints));
   for iexcempt=1:length(excempt)
       fprintf('\n\nAdd points only observed in dataset %s to the processing:\n',datasetIds{excempt(iexcempt)})
       forcedPoints = forcedPoints | tmpcountp(excempt(iexcempt),:)' >= 1;
   end
   forcedPoints(~pntMask)=false;
end
   if any(forcedPoints)
    
       % Forced points not only contains new, but also already selected points
    
       upcomp=3;
    
       if opt.doplots > 0
          figure
          imagesc(SO(forcedPoints,epochMask,upcomp))
          colorbar()
          title('forcedPoints (initial)')
       end
    
       % find epochs with at least three new points
       numForcedPoints0=sum(forcedPoints);
       forcedEpochs = epochMask & sum(SO(forcedPoints,:,upcomp),1) > 2;
       numForcedEpochs=sum(forcedEpochs);
    
       fprintf('- Number of epochs at least three newly added points: %d\n',numForcedEpochs)
        
       %forcedPoints = forcedPoints & sum(SO(:,epochMask,upcomp),2) > 1;
       forcedPoints = forcedPoints & sum(SO(:,forcedEpochs,upcomp),2) > 1;
       numForcedPoints=sum(forcedPoints);
    
       fprintf('- Number of non-redundant points that are added: %d  (out of %d candidates)\n',numForcedPoints,numForcedPoints0)
        
       if opt.doplots > 0
          figure
          imagesc(SO(forcedPoints,epochMask,upcomp))
          colorbar()
          title('forcedPoints (final)')

          figure
          imagesc(SO(forcedPoints,forcedEpochs,upcomp))
          colorbar()
       end
    
    end

%pntMask2=pntMask  & sum(tmpcountp) >= opt.minDatasetsPoint;
pntMask2=pntMask  & ( sum(tmpcountp)' >= opt.minDatasetsPoint | forcedPoints  ) ;
if any(pntMask ~= pntMask2)
   np=sum(pntMask & ~pntMask2);
   if opt.verbose > 2 || ( opt.verbose > 1 && np <= 200 ) || ( opt.verbose > 0 && np <= 50 )   
      fprintf('\nPoint mask has changed because there were some points observed in only one dataset ...\n%s\n\n', ...
          strjoin(pntName( pntMask & ~pntMask2)) )
   else
      fprintf('\nPoint mask has changed because there were some points observed in only one dataset, %d points removed.\n\n', np)
   end
   pntMask=pntMask2;
end

% Filter out epochs and points with no connections

epochMask2=epochMask & any(redEpo,1);
if any(epochMask ~= epochMask2)
   if opt.onlyRedundantEpochs 
       % Keep only epochs with redundant data, ie. with connecting points between datasets  
       fprintf('\nEpoch mask has changed because there were some epochs without connecting(common) points ...\n%s\n\n', ...
           strjoin(epochIds( epochMask & ~epochMask2)) )
       epochMask=epochMask2;
   else
       fprintf('\nThere are some epochs without connecting(common) points, but we will keep them ...\n%s\n\n', ...
           strjoin(epochIds( epochMask & ~epochMask2)) )
   end
end

%pntMask2=pntMask  & any(redPnt,2);
pntMask2=pntMask  & ( any(redPnt,2) | forcedPoints );
if any(pntMask ~= pntMask2)
   np=sum(pntMask & ~pntMask2);
   if opt.onlyRedundantPoints 
      % Keep only points with redundant data, ie. with connecting points between datasets  
      if opt.verbose > 2 || ( opt.verbose > 1 && np <= 200 ) || ( opt.verbose > 0 && np <= 50 )   
         fprintf('\nPoint mask has changed because there were some points without redundancy ...\n%s\n\n', ...
              strjoin(pntName( pntMask & ~pntMask2)) )
      else
         fprintf('\nPoint mask has changed because there were some points without redundancy, %d points removed.\n\n', np)
      end
      pntMask=pntMask2;
   else
      if opt.verbose > 2 || ( opt.verbose > 1 && np <= 200 ) || ( opt.verbose > 0 && np <= 50 )   
         fprintf('\nThere are some points without redundancy, but we will keep them ...\n%s\n\n', ...
              strjoin(pntName( pntMask & ~pntMask2)) )
      else
         fprintf('\nThere are %d points without redundancy, but we will keep them.\n\n', np)
      end
   end
end

% Filter out whole points and epochs

parMask(~pntMask,:,:)=false;
parMask(:,~epochMask,:)=false;

% Update the point and epoch masks
%
% Point and epoch mask needs updating, because points and/or epochs may 
% not anymore have enough observations (criterion is an input option).
% We recompute the pntMask and epochMask based on parMask.

epochMask2=any(sum(parMask,1) >= opt.minDisplEpoch,3);
pntMask2=any(sum(parMask,2) >= opt.minDisplPoint,3);

if any(epochMask ~= epochMask2)
   fprintf('\nEpoch mask has changed on account of the parameter mask, newly removed epochs are ...\n%s\n\n', ...
        strjoin(epochIds( epochMask & ~epochMask2)) )
   epochMask=epochMask2;
end

if any(pntMask ~= pntMask2)
   np=sum(pntMask & ~pntMask2);
   if opt.verbose > 2 || ( opt.verbose > 1 && np <= 200 ) || ( opt.verbose > 0 && np <= 50 )   
      fprintf('\nPoint mask has changed on account of the parameter mask, newly removed points are ...\n%s\n\n', ...
           strjoin(pntName( pntMask & ~pntMask2)) )
   else
      fprintf('\nPoint mask has changed on account of the parameter mask, %d points removed.\n\n', np)
   end
   pntMask=pntMask2;
end

% Check that enough epochs and points are left for the combination

numPoints2=sum(pntMask(:));
numEpochs2=sum(epochMask(:));
numDispl2=sum(sum(parMask));

fprintf('\nFilter statistics:\n\n')
fprintf('%d points remain out of %d\n',numPoints2,numPoints)
fprintf('%d epochs remain out of %d\n\n',numEpochs2,numEpochs)

for l=1:numel(nomParTypes)
   fprintf('%d %s displacements remain out of %d (%d max)\n',numDispl2(l),nomParTypes{l},numDispl(l),numPoints*numEpochs);
end

if numEpochs2 < opt.minEpochs || numPoints2 < opt.minPoints
   error('Not enough points or epochs left for the combination.')
end

% Recompute the dimension ndim

selDim=numDispl2 >= opt.minDispl;
ndim=sum(selDim);

if all(~selDim)
   error('None of the components have enough observable displacements.')
end
if any(~selDim)
   fprintf('The %s component(s) cannot be estimated, only %s component(s) are estimated, the dimension is reduced to %d.\n',strjoin(nomParTypes(~selDim)),strjoin(nomParTypes(selDim)),ndim)
end


%% Print overview of datasets, number of points and epochs 

numObsTot2=zeros(1,numel(expObsTypes));
ncDs=zeros(numDatasets,4);
fprintf('\n                                                                       Observation count per Type\n\n')
fprintf(  'DatasetId                              TechId   Pnts Epochs  Dim    North   East     Up ascLos dscLos    total\n')
fprintf(  '-------------------------------------  ------ ------ ------ ----   ------ ------ ------ ------ ------   ------\n')
for k=1:numDatasets
    % Point and epoch mask for the dataset 
    idxPnt=dslocal{k}.idxPnt;
    idxEpo=dslocal{k}.idxEpoch;
    pntMaskDs=pntMask(idxPnt);
    epochMaskDs=epochMask(idxEpo);
    pntMaskDs2=pntMaskDs & ~all(all(isnan(datasets{k}.obsData(:,epochMaskDs,:)),2),3);
    epochMaskDs2=epochMaskDs & ~all(all(isnan(datasets{k}.obsData(pntMaskDs,:,:)),1),3);
    % Count number of points and epochs that are actually used
    ncDs(k,1)=sum(pntMaskDs2);
    ncDs(k,2)=sum(epochMaskDs2);
    ncDs(k,3)=sum(pntMaskDs)-ncDs(k,1);
    ncDs(k,4)=sum(epochMaskDs)-ncDs(k,2);
    dsnote='';
    if ncDs(k,3) > 0
       dsnote=sprintf('%s (removed %d point(s) w/o obs)',dsnote,ncDs(k,3));
    end
    if ncDs(k,4) > 0
       dsnote=sprintf('%s (removed %d epoch(s) w/o obs)',dsnote,ncDs(k,4));
    end
    % Observation type count matrix (handles missing observations)
    numObsDs=zeros(1,numel(expObsTypes));
    for ll=1:numel(dslocal{k}.obsTypes)
       l=find(ismember(expObsTypes,dslocal{k}.obsTypes{ll}));
       numObsDs(l)=numel(find(~isnan(datasets{k}.obsData(pntMaskDs2,epochMaskDs2,ll))));
    end     
    fprintf('%-38s %-5s  %6d %6d  %3d   %6d %6d %6d %6d %6d   %6d   %s\n',dslocal{k}.datasetId,dslocal{k}.techniqueId, ...
       ncDs(k,1),ncDs(k,2),numel(dslocal{k}.obsTypes),numObsDs,sum(numObsDs),dsnote);
    numObsTot2=numObsTot2+numObsDs;
end
fprintf(  '                                              ------ ------        ------ ------ ------ ------ ------   ------\n')
fprintf('Distinct / total                              %6d %6d        %6d %6d %6d %6d %6d   %6d\n\n',numPoints2,numEpochs2,numObsTot2,sum(numObsTot2))

% Check that any of the datasets has enouch epochs and points

if any( ncDs(:,2) < opt.minEpochs ) || any( ncDs(:,1) < opt.minPoints )  
   error('One or more datasets do not have enough epochs and/or points for the combination.')
end

%% Least squares solution (via observation equations)

% Initialize main normal equations N11tot * x1 = b1tot for the 1st part.
% Ignoring missing elements, this would be something like this
%
%    N11tot=zeros(numPoints*numEpochs*numNomParTypes,numPoints*numEpochs*numNomParTypes);
%    b1tot=zeros(numPoints*numEpochs*numNomParTypes,1);
%
% Where we initialize with numNomParTypes instead of ndim, to make admin easier as it
% corresponds to the dimension of the sensitivityMatrix.
%
% However, this matrix would be HUGE, and does not fit in memory, though many elements are zero. We
% will skip many of the zero elements, especially those that are easy to
% identify. For this we build two index arrays (from parMask)
%
%   idxParMask :  ind=idxParMask(idx) returns the linear position "ind" in the [numPoints numEpochs numNomParTypes] matrix
%                 for the corresponding positions "idx" in N11tot, b1tot and x1tot ,
%                 [i,j,l]=ind2sub([numPoints numEpochs 3],ind) would return the subscript values, 
%                 idxParMask has the same length as b1tot and x1tot
%
%   indParMask :  idx=indParMask(ind) returns the linear position "idx" in N11tot, b1tot and x1tot for 
%                 the corresponding positions "ind" (given by ind=sub2ind([numPoints numEpochs NumNomParTypes],i,j,l) )
%                 in the [numPoints numEpochs numNomParTupes] matrix, indParMask has the same length as parMask(:), 
%                 indParMask has NaN's for the elements that do not exist in N11tot, b1tot and x1tot. 
%
% The dimension of the normal equations is now equal to the length of idxParMask, which is assigned 
% to the value numDisplPar.

idxParMask=find(parMask(:));
indParMask=nan(size(parMask(:)));
indParMask(idxParMask)=1:size(idxParMask,1);

% Initialize main normal equations N11tot * x1 = b1tot for the 1st part

numDisplPar=size(idxParMask,1);
N11tot=zeros(numDisplPar,numDisplPar,'single');
b1tot=zeros(numDisplPar,1,'single');

% Initialize cell array to hold intermediate least squares results

idxlsq=0;
for k=1:numDatasets
   dslocal{k}.idxlsq=idxlsq;
   idxlsq=idxlsq+dslocal{k}.numObsTypes;  
end
lsqsave=cell(idxlsq,1);

% Step 1 - provisional solution for offsets and transformation (A2)

fprintf('\nStep 1 - Pre-elimination of offset and transformation parameters\n')

if opt.doplots > 0 || opt.plotCond
    hnorm=figure('name','Normal matrix / Cholesky factor diagonals');
    hn1=subplot(2,1,1);
    hn2=subplot(2,1,2);
    hcond=figure('name','Cholesky factor diagonals');
    hc=subplot(2,1,1);
end

fprintf('\nDataset                                 Tech  Pnts Epochs   Dim     obsType   nObs  nPar2  Rd  RCOND_Qyy  RCOND_U22\n\n')

nobs=0;
npar2=0;

for k=1:numDatasets

    % Number of points (np), epochs (ne), observations (m) and unknowns (n) for this dataset

    np0=dslocal{k}.numPoints;     % Number of points in dataset
    ne0=dslocal{k}.numEpochs;     % Number of epochs in dataset
    nc=dslocal{k}.numObsTypes;    % Number of elements in dataset
    m=np0*ne0;                    % Number of observations in dataset
    n=np0+ne0;                    % Number of unknowns

    % Indices to points and epochs, mask points and epochs, recount
   
    idxPnt=dslocal{k}.idxPnt;
    idxEpo=dslocal{k}.idxEpoch;
    pntMaskDs=pntMask(idxPnt);
    epochMaskDs=epochMask(idxEpo);
%     if opt.minObsDispl < 2  
%         pntMaskDs=pntMaskDs & ~all(all(isnan(datasets{k}.obsData(:,epochMaskDs,:)),2),3);
%         epochMaskDs=epochMaskDs & ~all(all(isnan(datasets{k}.obsData(pntMaskDs,:,:)),1),3);
%     else
%         pntMaskDs=pntMaskDs & any(sum(~isnan(datasets{k}.obsData(:,epochMaskDs,:)),2) >= opt.minObsDispl,3);
%         epochMaskDs=epochMaskDs & any(sum(~isnan(datasets{k}.obsData(pntMaskDs,:,:)),1) >= opt.minObsDispl,3);
%     end
    % update point masks   (need to iterate over this...)
    numDisplPoint = sum(~isnan(datasets{k}.obsData(:,epochMaskDs,:)),2);
    pntMaskDs = pntMaskDs & any(numDisplPoint >= opt.minObsDispl,3) ;
    % update epoch masks
    numDisplEpoch = sum(~isnan(datasets{k}.obsData(pntMaskDs,:,:)),1);
    epochMaskDs2 = epochMaskDs & any(numDisplEpoch >= opt.minObsDispl,3);
    if any(xor(epochMaskDs2,epochMaskDs))
       epochMaskDs=epochMaskDs2;
       % epoch mask has changed, repeat selection
       numDisplPoint = sum(~isnan(datasets{k}.obsData(:,epochMaskDs,:)),2);
       pntMaskDs = pntMaskDs & any(numDisplPoint >= opt.minObsDispl,3) ;
    end
    if any(~pntMaskDs) || any(~epochMaskDs)
       % recount
       idxPnt=idxPnt(pntMaskDs);
       idxEpo=idxEpo(epochMaskDs);
       np=length(idxPnt(:));
       ne=length(idxEpo(:));
       m=np*ne;
       n=np+ne;
       if opt.verbose > 1
           fprintf('recount: np %d -> %d, ne: %d -> %d\n',np0,np,ne0,ne)
       end
    else
       np=np0;
       ne=ne0;
    end

    fprintf('%-38s %5s %5d %6d  %3d\n',dslocal{k}.datasetId,dslocal{k}.techniqueId, ...
       np,ne,nc)

    % Indices to main normal equations N11tot and b1tot
    
    ii=repmat(idxPnt(:),[1,ne]);
    jj=repmat(idxEpo(:)',[np,1]);
    idx=ii(:)+(jj(:)-1)*numPoints;
    
    % Build the (sparse) design matrix A2 for this dataset

    ii=repmat([1:np]',[1,ne]);
    jj=repmat([1:ne],[np,1]);
    A2=sparse( [1:m 1:m]',[ ii(:); jj(:)+np],ones(m*2,1),m,n,m*2);

    % Get the stochastic model parameters 

    stochModel=dslocal{k}.stochModel;
    stochData=datasets{k}.stochData;
    
    % Index to cell array to hold intermediate results

    idxlsq=dslocal{k}.idxlsq;
    
    % Loop over the observed elements

    for l=1:dslocal{k}.numObsTypes

       % Get space-time matrix with observations (from datasets{k})

       y=datasets{k}.obsData(pntMaskDs,epochMaskDs,l);

       % Sensitivity matrix (numPoints,3,l)
             
       sensitivityMatrix=datasets{k}.sensitivityMatrix(pntMaskDs,:,l);
       zerocols=all(sensitivityMatrix == 0 );
       onescols=all(sensitivityMatrix == 1 );
       if  ( all( onescols | zerocols ) && sum(onescols) == 1 )
           sensitivityMatrixType='unitcol';
           sensitivityCols=find(onescols);
       elseif any(zerocols)
           sensitivityMatrixType='selcols';
           sensitivityCols=find(~zerocols);
       else
           sensitivityMatrixType='full';
           sensitivityCols=1:size(sensitivityMatrix,2);
       end
       
       % Observed minus a-priori and ymask
       
       ymask=~isnan(y);
       switch sensitivityMatrixType
           case 'unitcol'
               y=y-aprioriDispl(idxPnt,idxEpo,sensitivityCols);
               ymask=ymask & parMask(idxPnt,idxEpo,sensitivityCols);
           case {'selcols','full'}
               for i=1:numel(sensitivityCols)
                  iCol=sensitivityCols(i);
                  d=repmat(sensitivityMatrix(:,iCol),[1 ne]);
                  y=y-d.*aprioriDispl(idxPnt,idxEpo,iCol);   
                  ymask=ymask & parMask(idxPnt,idxEpo,iCol);
               end
           otherwise
               error('Hey, this should never happen, we have an illegal case for the sensitivity matrix.')
       end
       
       % Find rows (points) and columns (epochs) without observations, will be removed later from A2, N22 and b2
       
       x2mask=true(np+ne,1);
       removecols = [ find(sum(ymask,2) == 0) ; np+find(sum(ymask,1) == 0) ];
       x2mask(removecols)=false;
       
       % Vectorize observed minus a-priori and mask the non-observed (nan) entries

       y=y(:);    
       ymask=ymask(:);
       y=y(ymask);
       idxmask=idx(ymask);
       my=length(y);

       % Compose co-variance matrix of observations
       
       %Qy=stmstochmodel(stochModel{l},stochData,dslocal{k}.pntNeu,dslocal{k}.epochDyear);
       %Qy=stmstochmodel(stochModel{l},stochData,dslocal{k}.pntNeu(pntMaskDs,:),dslocal{k}.epochDyear(epochMaskDs));
       if any(strncmp(stochModel{l},'tudinsar4rd',11))
          % stochModel{l}='tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)';  % there is a problem with the tudinsar4rd model (not positive definite)
          pntAttrib=datasets{k}.pntAttrib;
          epochAttrib=datasets{k}.epochAttrib;
          Qy=stmstochmodel(stochModel{l},stochData,dslocal{k}.pntNeu,dslocal{k}.epochDyear,{pntMaskDs,epochMaskDs},pntAttrib,epochAttrib);
          % Qy=eye(m);
       else
          Qy=stmstochmodel(stochModel{l},stochData,dslocal{k}.pntNeu,dslocal{k}.epochDyear,{pntMaskDs,epochMaskDs});
       end
       Qy=Qy(ymask,ymask);
       Qyopt='full'; 

       % provisional solution dataset specific parameters (x2a)
       %
       % Ry=chol(Qy) ->  Qy = Ry'*Ry  => inv(Qy) = inv(Ry)*inv(Ry')
       %
       %    A2n=inv(Ry')*A2            -> A2n=Ry'\A2
       %    yn=inv(Ry')*y              -> yn=Ry'\y
       %    b2a=A2'*inv(Qy)*y         -> b2a=A2n'*yn
       %    N22=A2'*inv(Qy)*A2        -> N22=A2n'*A2n
       %    x2a=inv(N22)*b2a          -> x2a=N22\b2a
       %
       % U = chol(N22)  -> N22 = U'*U  => inv(N22)=inv(U)*inv(U')
       %
       %    x2a=inv(N22)*b2a          -> x2a=N22\b2a=U\(U'\b2a)
       %    Qx=inv(n22)               -> Qx=N22\I==U\(U'\I)

       % Convert design-matrix A2 and observations to y handle inv(Qy) ->  A2n, yn

       switch lower(Qyopt)
          case 'diag'
             % Diagonal co-variance matrix
             wy = 1./sqrt(Qy);
             A2n = wy.*A2(ymask,x2mask);    % A2n=Qy^(-1/2)*A2
             yn = wy.*y;               % yn=Qy^(-1/2)*y
             Ry=chol(diag(Qy));        % <- NEED TO OPTIMIZE THIS PART
          case 'full'
             % full co-variance matrix, compute Cholesky factor
             [Ry,p]=chol(Qy);
             if p ~= 0
                warning(['Co-variance matrix of dataset ' num2str(k) ' is not positive definite'])
             end
             A2n=Ry'\A2(ymask,x2mask);  % A2n=Qy^(-1/2)*A2
             yn=Ry'\y;                  % yn=Qy^(-1/2)*y
           otherwise
             error(['Unsupported covariance matrix option ' Qopt ]) 
       end
 
       sqrcond1=min(diag(Ry))/max(diag(Ry));
       rankdefect=1;
       
       % Compute provisional solution x2a for the second part of the equations

       doqr=false;
       if doqr
          % Use QR decomposition to obtain the solution x2a
          if issparse(A2n)
             perm = colamd(A2n(:,1:end-1));
             [z,U22] = qr(A2n(:,perm),yn,0);
          else
             [Q,U22,perm] = qr(A2n(:,1:end-1),0);
             z = Q'*yn;
          end
          x2a([ perm  end ])= [ U22 \ z ; 0 ];
          sqrcond2=min(diag(U22))/max(diag(U22));
       else
          % Use Cholesky factorization to obtain the solution x2a 
          b2a=A2n'*yn;           % b2a=A2'*Qy^(-1)*y
          N22=A2n'*A2n;          % N22=A2'*Qy^(-1)*A2
          [U22,p]=chol(N22(1:end-rankdefect,1:end-rankdefect));            % N22 = U22'*U22 is rank-defect
          if p ~= 0
             rankdefect=numel(b2a)-p+1;
             disp(['Rank-defect of N22 for dataset ' num2str(k) ' is not as expected, position=' num2str(rankdefect) ', rankdefect=' num2str(rankdefect)])
             warning(['Rank-defect of N22 for dataset ' num2str(k) ' is not as expected, position=' num2str(rankdefect) ', rankdefect=' num2str(rankdefect)])
          end
          sqrcond2=min(diag(U22))/max(diag(U22));
          % Actual solve (cath warnings)
          origState = warning;
          warning('off')
          lastwarn('','')
          x2a=[ U22 \ ( U22' \ b2a(1:end-rankdefect) ) ; zeros(rankdefect,1) ];     % x2a = N22\b2a = inv(N22)*b2a     
          [msg,warnID] = lastwarn;
          warning(origState)          
          if opt.verbose > 0 && ~isempty(msg)
             fprintf ('%2d %s %s %10.3e %s  %s\n',k,dslocal{k}.datasetId,dslocal{k}.obsTypes{l},sqrcond2, ...
                msg,warnID)
          end
          % Add diagonal of Cholesky factor to plot (optional)
          if opt.doplots > 0 || opt.plotCond
              x2labels=[ pntName(idxPnt) ; epochIds(idxEpo)' ];
              x2labels=x2labels(x2mask);
              figure(hnorm);
              subplot(2,1,1)
              h=semilogy(diag(N22),'-+','displayName',['N22 ' dslocal{k}.datasetId ' ' dslocal{k}.obsTypes{l}]);
              h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dataset',repmat({h.DisplayName},size(h.XData)));
              h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('name',x2labels);
              h.DataTipTemplate.Interpreter = 'none';
              hold on
              subplot(2,1,2)
              h=semilogy(diag(U22),'-x','displayName',['U22 ' dslocal{k}.datasetId ' ' dslocal{k}.obsTypes{l}]);
              h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dataset',repmat({h.DisplayName},size(h.XData)));
              h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('name',x2labels(1:end-rankdefect));
              h.DataTipTemplate.Interpreter = 'none';
              hold on
              figure(hcond);
              subplot(2,1,1)
              [du22,iu22]=sort(diag(U22),'descend');
              h=semilogy(du22,'-x','displayName',[dslocal{k}.datasetId ' ' dslocal{k}.obsTypes{l}]);
              h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dataset',repmat({h.DisplayName},size(h.XData)));
              h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('name',x2labels(iu22));
              h.DataTipTemplate.Interpreter = 'none';
              hold on
          end
       end
       
       % Reduced normal equation of common parameters (x1) update
       %
       % Ry=chol(Qy) ->  Qy = Ry'*Ry  => inv(Qy) = inv(Ry)*inv(Ry')
       %
       % A1n=inv(Ry')*A1            -> A1n=Ry'\A1
       % b1=A1'*inv(Qy)*(y-A2*x2)   -> b1=A1n'*(yn-A2n*x2)
       %
       % N11=A1'*inv(Qy)*A1 - A1'*inv(Qy)*A2*inv(N22)*A2'*inv(Qy)*A1
       %                           -> N11=A1n'*A1n - A1n'*A2n*inv(U22'*U22)*A2n'*A1n
       %                           -> N11=A1n'*A1n - A1n'*A2n*inv(U22)*inv(U22)'*A2n'*A1n
       %
       % Three scenarios are possible
       % A. Sensitivity matrix indicates sensitivity in one element only
       %    One column of the sensitivity matrix is all ones, the other columns 
       %    are all zero. In this case A1 can be replaced with an identify matrix and
       %    only specific part of the normal equations need to be updated
       % B. Sensitivity matrix has zero columns, but also one or more real
       %    valued columns (not equal to one).
       % C. Sensitivity matrix has all real values (no columns with all
       %    ones or zeros)

       switch sensitivityMatrixType
           case 'unitcol'
        
               % A1 is e.g. [ O ; O ; I ], with an identity matrix for the component
               % that is observed (in this case height displacement) -> we don't 
               % need A1n = inv(Ry')*A1 explicitly, we can simply solve with Ry ...
               %
               % b1=A1n'*(yn-A2n*x2)  -> b1=inv(Ry)*(yn-A2n*x2) ->  b1=Ry\(yn-A2n*x2)  
               % N11=inv(Qy)-inv(Qy)*A2*inv(N22)*A2'*inv(Qy) ->
               %       -> N11=inv(Ry)*inv(Ry')-inv(Ry)*A2n*inv(N22)*A2n'*inv(Ry')
               %       -> N11=inv(Ry)*(I-A2n*inv(N22)*A2n')*inv(Ry')
               %       -> N11=inv(Ry)*(I-A2n*inv(U22'*U22)*A2n')*inv(Ry')
               %       -> N11=inv(Ry)*(I-A2n*inv(U22)*inv(U22')*A2n')*inv(Ry')
               %       -> N11=inv(Ry)*(I-A2n*inv(U22)*inv(U22')*A2n')*inv(Ry')
               %       -> N11=inv(Ry)*(I-A2s*A2s')*inv(Ry')
               % A2s=A2n*inv(U22)=A2n/U22

               b1=Ry\(yn-A2n*x2a);
               %N11=inv(Qy)-inv(Qy)*A2(:,1:end-1)*inv(N22(1:end-1,1:end-1))*A2(:,1:ends-1)'*inv(Qy);      
               A2s= A2n(:,1:end-rankdefect) / U22 ;                  % A2s= A2n*inv(U22) 
               N11 = ( Ry \ (eye(my)-A2s*A2s') ) / Ry' ;
       
               % Add to reduced normal equation update to reduced normal equations
               %
               % b1tot(idx)=b1tot(idx)+b1   
               % N11tot(idx,idx)=N11tot(idx,idx)+N11
      
               ll=(sensitivityCols-1)*numPoints*numEpochs;
               %lll=ll+idxmask;
               lll=indParMask(ll+idxmask);
               b1tot(lll)=b1tot(lll)+b1;   
               N11tot(lll,lll)=N11tot(lll,lll)+N11;
              
           case {'selcols','full'}

               % A1 is e.g. [ D1 ; D2 ; D3 ], with D1, D2 and D3 diagonal matrices,
               % with Di=diag(repmat(si,[1,ne])), and si the i'th column of the
               % sensitivity matrix. The diagonal of Di is composed of si, repeated
               % as many times as there are epochs. With
               %
               %    A1n=inv(Ry)'*A1 ->  A1n= [ inv(Ry)'*D1 ; inv(Ry)'*D2 ; inv(Ry)'*D3 ]
               %
               % it is possible to rewrite this into a modified 'unitcol' case, with
               % b1 and N11 from the 'unitcol' case, we have
               %
               % b1c=A1n'*(yn-A2n*x2)  -> b1c= [D1 D2 D3]*inv(Ry)*(yn-A2n*x2) = [D1 D2 D3]*b1
               %
               % N11c=A1'*inv(Qy)*A1 - A1'*inv(Qy)*A2*inv(N22)*A2'*inv(Qy)*A1 =
               %        A1' * ( inv(Qy) - inv(Qy)*A2*inv(N22)*A2'*inv(Qy) ) * A1 =
               %        [D1 D2 D3] * N11 * [D1' ; D2' ; D3']

               b1=Ry\(yn-A2n*x2a);
               %N11=inv(Qy)-inv(Qy)*A2(:,1:end-1)*inv(N22(1:end-1,1:end-1))*A2(:,1:ends-1)'*inv(Qy);      
               A2s= A2n(:,1:end-rankdefect) / U22 ;                  % A2s= A2n*inv(U22) 
               N11 = ( Ry \ (eye(my)-A2s*A2s') ) / Ry' ;
       
               % Add to reduced normal equation update to reduced normal equations
               %
               % b1tot(idx)=b1tot(idx)+b1   
               % N11tot(idx,idx)=N11tot(idx,idx)+N11

               for i=1:numel(sensitivityCols)
                  iCol=sensitivityCols(i);
                  d=repmat(sensitivityMatrix(:,iCol),[ne 1]);
                  ll=(iCol-1)*numPoints*numEpochs;
                  lll=indParMask(ll+idxmask);
                  %lll=ll+idxmask;
                  b1tot(lll)=b1tot(lll)+d(ymask).*b1;   
                  N11tot(lll,lll)=N11tot(lll,lll)+(d(ymask)*d(ymask)').*N11;
                  %b1tot(lll)=b1tot(lll)+d.*b1;   
                  %N11tot(lll,lll)=N11tot(lll,lll)+(d*d').*N11;
                  for j=1:i-1
                     jCol=sensitivityCols(j);
                     dj=repmat(sensitivityMatrix(:,jCol),[ne 1]);
                     lj=(jCol-1)*numPoints*numEpochs;
                     llj=indParMask(lj+idxmask);
                     %llj=lj+idxmask;
                     N11tot(lll,llj)=N11tot(lll,llj)+(d(ymask)*dj(ymask)').*N11;
                     N11tot(llj,lll)=N11tot(llj,lll)+(dj(ymask)*d(ymask)').*N11;
                     %N11tot(lll,llj)=N11tot(lll,llj)+(d*dj').*N11;
                     %N11tot(llj,lll)=N11tot(llj,lll)+(dj*d').*N11;
                  end
               end
               
           otherwise
               error('Hey, this should never happen, we have an illegal case for the sensitivity matrix.')
       end
       
       % Print number of observations, parameters and rank-defect

       nobs=nobs+my;
       nx2=numel(x2a);
       npar2=npar2+nx2-rankdefect;

       fprintf('                                                                    %6s %6d %6d %3d  %10.4e %10.4e %s\n',dslocal{k}.obsTypes{l},my,nx2-rankdefect,rankdefect,sqrcond1,sqrcond2,msg)

       % Save intermediate results for step 3

       lsqsave{idxlsq+l}.np=np;
       lsqsave{idxlsq+l}.ne=ne;      
       lsqsave{idxlsq+l}.y=y;
       lsqsave{idxlsq+l}.Ry=Ry;
       lsqsave{idxlsq+l}.yn=yn;
       lsqsave{idxlsq+l}.b2a=b2a;
       lsqsave{idxlsq+l}.A2n=A2n;
       lsqsave{idxlsq+l}.U22=U22;
       lsqsave{idxlsq+l}.A2s=A2s;  
       lsqsave{idxlsq+l}.sensitivityMatrixType=sensitivityMatrixType;  
       lsqsave{idxlsq+l}.sensitivityCols=sensitivityCols;  
       lsqsave{idxlsq+l}.sensitivityMatrix=sensitivityMatrix;  
       lsqsave{idxlsq+l}.N11=N11;  
       lsqsave{idxlsq+l}.x2a=x2a;   % saved, though not needed later
       %lsqsave{idxlsq+l}.A1n=A1n;  % A1n=Ry\A1
       lsqsave{idxlsq+l}.idxmask=idxmask;
       lsqsave{idxlsq+l}.ymask=ymask;
       lsqsave{idxlsq+l}.x2mask=x2mask;
       lsqsave{idxlsq+l}.idxPnt=idxPnt;
       lsqsave{idxlsq+l}.idxEpo=idxEpo;
       lsqsave{idxlsq+l}.epochMaskDs=epochMaskDs;
       lsqsave{idxlsq+l}.pntMaskDs=pntMaskDs;
       lsqsave{idxlsq+l}.rankdefect=rankdefect;
       
    end
   
end

fprintf('                                             -----  -----  ---             ------ ------\n')
%fprintf('                                           %5d %6d  %3d             %6d %6d\n\n',numPoints,numEpochs,ndim,nobs,npar2)
fprintf('                                             %5d %6d  %3d             %6d %6d\n\n',sum(pntMask),sum(epochMask),ndim,nobs,npar2)

if opt.doplots > 0 || opt.plotCond
   legend('interpreter','none')
   title('Cholesky factor U22 diagonal (sorted)')
   ylabel('[-]')
end


%% Step 2 - Solution of displacement part 

fprintf('\n\nStep 2 - Solution of the displacement part\n')

% Redundancy in space-time format (diagonal of N11)
%
% - all non-observed parameters (displacements) are zero
% - when sum over rows is zero the corresponding point (component) is not observed  
% - when sum over columns is zero the epoch (component) is not observed  
%
% - badly observed parameters (displacements) have small numbers, e.g. North 
%   components for points with InSAR only, and may require special
%   treatment or regularization (not yet fully implemented -> criterion
%   epsObs should depend on sqrt(diag(Qy)) and sigma0)

%red11=reshape(diag(N11tot),[ numPoints numEpochs numNomParTypes]);

red11=nan([ numPoints numEpochs numNomParTypes],'single');
red11(idxParMask)=diag(N11tot);

if opt.doplots > 0

    figure('Name','Sqrt normal matrix diagonal','Position',[680 558 1260 420])
%   for k=1:ndim
    for k=1:numNomParTypes
       subplot(1,numNomParTypes,k)
       imagesc(sqrt(red11(:,:,k)))
       colorbar
       title(nomParTypes{k})
    end
    
    red11N=red11(:,:,1);
    red11E=red11(:,:,2);
    red11U=red11(:,:,3);

    figure('Name','Sqrt normal matrix histograms','Position',[680 558 1260 420])

    subplot(1,4,1)
    histogram(sqrt(red11N(SC(:,:,1) == 1)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    hold on
    histogram(sqrt(red11E(SC(:,:,2) == 1)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    histogram(sqrt(red11U(SC(:,:,3) == 1)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    legend('N','E','U')
    title('underObs')

    subplot(1,4,2)
    histogram(sqrt(red11N(SC(:,:,1) == 2 | SC(:,:,1) == 3)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    hold on
    histogram(sqrt(red11E(SC(:,:,2) == 2 | SC(:,:,2) == 3)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    histogram(sqrt(red11U(SC(:,:,3) == 2 | SC(:,:,3) == 3)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    legend('N','E','U')
    title('losObs')

    subplot(1,4,3)
    histogram(sqrt(red11N(SC(:,:,1) == 4)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    hold on
    histogram(sqrt(red11E(SC(:,:,2) == 4)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    histogram(sqrt(red11U(SC(:,:,3) == 4)),'BinWidth',.01,'DisplayStyle','stairs','LineWidth',2)
    legend('N','E','U')
    title('nearObs')

    subplot(1,4,4)
    histogram(sqrt(red11N(SC(:,:,1) == 7)),'BinWidth',.05,'DisplayStyle','stairs','LineWidth',2)
    hold on
    histogram(sqrt(red11E(SC(:,:,2) == 7)),'BinWidth',.05,'DisplayStyle','stairs','LineWidth',2)
    histogram(sqrt(red11U(SC(:,:,3) == 7)),'BinWidth',.05,'DisplayStyle','stairs','LineWidth',2)
    legend('N','E','U')
    title('fullObs')

end

% Booleans with observed parameters

hasdispl= ( SC > 1 );
numDisplPerEpoch=squeeze(sum(hasdispl,1));
numDisplPerPoint=squeeze(sum(hasdispl,2));

% Counts of the displacements parameters (full dataset)

fprintf('\nDisplacement parameters (full dataset):\n')
fprintf('                                                                Computing base parameters\n')
fprintf('           noObs underObs   losPar  nearPar  fullPar   usable        Transf  Offset\n')
fprintf('        -------- -------- -------- -------- -------- --------       ------- -------\n')
for l=1:numNomParTypes
   noobs=sum(sum( SC(:,:,l) == 0 ));
   underobs=sum(sum( SC(:,:,l) == 1));
   losobs=sum(sum( SC(:,:,l) == 2 | SC(:,:,l) == 3 ));
   nearobs=sum(sum( SC(:,:,l) == 4 ));
   fullobs=sum(sum( SC(:,:,l) > 4 ));
   useobs=losobs+nearobs+fullobs;
   numtransform=sum(numDisplPerEpoch(:,l) > 0);
   numoffset=sum(numDisplPerPoint(:,l) > 0);
   fprintf('%5s   %8d %8d %8d %8d %8d %8d       %7d %7d\n',nomParTypes{l},noobs,underobs,losobs,nearobs,fullobs,useobs,numtransform,numoffset)
end
noobs=sum( SC(:) == 0 );
underobs=sum( SC(:) == 0 );
losobs=sum( SC(:) == 2 | SC(:) == 3 );
nearobs=sum( SC(:) == 4 );
fullobs=sum( SC(:) > 4 );
useobs=losobs+nearobs+fullobs;
numtransform=sum(numDisplPerEpoch(:) > 0);
numoffset=sum(numDisplPerPoint(:) > 0);
fprintf('        -------- -------- -------- -------- -------- --------       ------- -------\n')
fprintf('        %8d %8d %8d %8d %8d %8d       %7d %7d\n', noobs,underobs,losobs,nearobs,fullobs,useobs,numtransform,numoffset)
    
%% Counts of the displacements parameters (after filtering)

SC2=SC;
SC2(~parMask)=0;

hasdispl= ( SC2 > 1 );

numDisplPerEpoch=squeeze(sum(hasdispl,1));
numDisplPerPoint=squeeze(sum(hasdispl,2));

fprintf('\nDisplacement parameters (after filtering):\n')
fprintf('                                                                Computing base parameters\n')
fprintf('           noObs underObs   losPar  nearPar  fullPar   usable        Transf  Offset\n')
fprintf('        -------- -------- -------- -------- -------- --------       ------- -------\n')
for l=1:numNomParTypes
   noobs=sum(sum( SC2(:,:,l) == 0 ));
   underobs=sum(sum( SC2(:,:,l) == 1));
   losobs=sum(sum( SC2(:,:,l) == 2 | SC(:,:,l) == 3 ));
   nearobs=sum(sum( SC2(:,:,l) == 4 ));
   fullobs=sum(sum( SC2(:,:,l) > 4 ));
   useobs=losobs+nearobs+fullobs;
   numtransform=sum(numDisplPerEpoch(:,l) > 0);
   numoffset=sum(numDisplPerPoint(:,l) > 0);
   fprintf('%5s   %8d %8d %8d %8d %8d %8d       %7d %7d\n',nomParTypes{l},noobs,underobs,losobs,nearobs,fullobs,useobs,numtransform,numoffset)
end
noobs=sum( SC2(:) == 0 );
underobs=sum( SC2(:) == 0 );
losobs=sum( SC2(:) == 2 | SC2(:) == 3 );
nearobs=sum( SC2(:) == 4 );
fullobs=sum( SC(:) > 4 );
useobs=losobs+nearobs+fullobs;
numtransform=sum(numDisplPerEpoch(:) > 0);
numoffset=sum(numDisplPerPoint(:) > 0);
fprintf('        -------- -------- -------- -------- -------- --------       ------- -------\n')
fprintf('        %8d %8d %8d %8d %8d %8d       %7d %7d\n', noobs,underobs,losobs,nearobs,fullobs,useobs,numtransform,numoffset)
    
%% Choose the computing base
%
% Choose a computing base
%
%   opt.cbaseOffset    = { 'firstEpoch' , 'lastEpoch', 'firstCommonEpoch' , 'lastCommonEpoch' , <decimalYear> }  
%   opt.cbaseTransform = { 'firstPoint' , 'lastPoint', 'firstCommonPoint' , 'lastCommonPoint' , <pointId> }
%
% In case of '[first|last][Epoch|Point]' always the first or last epoch (point) is chosen with all components observed. 
% This is not necessarily the same point for each epoch, or vice versa.
%
% In case of ''[first|last]Common[Epoch|Point]' the first or last epoch (point) is chosen with all 
% components and points (epochs) observed. When this is not possible, the epoch (point) is selected with
% most observed points (epochs).
%
% When a specific epoch (point) is selected it must be an epoch with all components available and
% all points (epochs) observed for that epoch (point).

% opt.cbaseOffset='firstCommonEpoch';                   % { 'firstEpoch' , 'lastEpoch', 'firstCommonEpoch' , 'lastCommonEpoch' , <decimalYear> }
% opt.cbaseTransform='firstCommonPoint';                % { 'firstPoint' , 'lastPoint', 'firstCommonPoint' , 'lastCommonPoint' , <pointId> }

% First select the available dimensions
%
% - selDim is a logical array (1..numNomParTypes)
% - ndim is the number of available dimensions

selDim=any(numDisplPerEpoch > 0);
ndim=sum(selDim);

% Compute the first and last point (epoch) for each epoch (point)
%
% - first|lastEpochPerPoint is a numPoints x numNomParTypes array with index to first|last epoch
% - first|lastPointPerEpoch is a numEpoch x numNomParTypes array with index to first|last point

firstEpochPerPoint=nan(size(hasdispl,1),numNomParTypes);
lastEpochPerPoint=nan(size(hasdispl,1),numNomParTypes);
firstPointPerEpoch=nan(size(hasdispl,2),numNomParTypes);
lastPointPerEpoch=nan(size(hasdispl,2),numNomParTypes);
for l=1:numNomParTypes
   if ~selDim(l), continue; end
   for i=1:size(hasdispl,1)
     iii=find(hasdispl(i,:,l));
     if isempty(iii),continue; end
     firstEpochPerPoint(i,l)=iii(1);
     lastEpochPerPoint(i,l)=iii(end);
   end
   for j=1:size(hasdispl,2)
     jjj=find(hasdispl(:,j,l));
     if isempty(jjj),continue; end
     firstPointPerEpoch(j,l)=jjj(1);
     lastPointPerEpoch(j,l)=jjj(end);
   end
end

% Possible replacement code for the section below 
%
% maxEpochsPerPoint=numDisplPerPoint(:,selDim) == max(numDisplPerPoint(:,selDim))
% maxPointsPerEpoch=numDisplPerEpoch(:,selDim) == max(numDisplPerEpoch(:,selDim))
% needoffset=numDisplPerPoint(:,selDim) > 0
% needtransform=numDisplPerEpoch(:,selDim) > 0

% Find candidate points which have all components observed
ibaseCandidatePoints=find(all(numDisplPerPoint(:,selDim) > 0 ,2)); 
% Keep candidate points which have all epochs observed
ibaseCandidatePoints=ibaseCandidatePoints(all( all( hasdispl(ibaseCandidatePoints, any(numDisplPerEpoch(:,selDim) > 0, 2), selDim) ,2) ,3));
% Check if candidates remain, if not, generate warning and fall back to first/last options
if isempty(ibaseCandidatePoints) && contains(opt.cbaseTransform,'Common')
   ctmp=strrep(opt.cbaseTransform,'Common','');       
   fprintf('\nOption cbaseTransform changed from %s to %s\n',opt.cbaseTransform,ctmp) 
   opt.cbaseTransform=ctmp;
end

% Find candidate epochs which have all components observed
ibaseCandidateEpochs=find(all(numDisplPerEpoch(:,selDim) > 0, 2));
% Keep candidate epochs which have all points observed
ibaseCandidateEpochs=ibaseCandidateEpochs(all( all( hasdispl( any(numDisplPerPoint(:,selDim) > 0,2), ibaseCandidateEpochs, selDim) ,1) ,3));
% Check if candidates remain, if not, generate warning and fall back to first/last options
if isempty(ibaseCandidateEpochs) && contains(opt.cbaseOffset,'Common')
   ctmp=strrep(opt.cbaseOffset,'Common','');       
   fprintf('\nOption cbaseOffset changed from %s to %s\n',opt.cbaseOffset,ctmp) 
   opt.cbaseOffset=ctmp;       
end

% Set the reference time series (transform parameters)
%
%   opt.cbaseTransform = { 'firstPoint' , 'lastPoint', 'firstCommonPoint' , 'lastCommonPoint' , <pointId> }

basePointPerEpoch=nan(numEpochs,numNomParTypes);
switch lower(opt.cbaseTransform)
    case lower('firstPoint')
       basePointPerEpoch=firstPointPerEpoch;
    case lower('lastPoint')
       basePointPerEpoch=lastPointPerEpoch;
    case lower('firstCommonPoint')
       basePointPerEpoch(:,selDim)=repmat(ibaseCandidatePoints(1),[numEpochs sum(selDim)]);
       basePointPerEpoch(numDisplPerEpoch == 0)=NaN;
    case lower('lastCommonPoint')
       basePointPerEpoch(:,selDim)=repmat(ibaseCandidatePoints(end),[numEpochs sum(selDim)]);
       basePointPerEpoch(numDisplPerEpoch == 0)=NaN;
    otherwise
       % Use selected point name
       ibpoint=ismember(opt.cbaseTransform,pntIds);
       basePointPerEpoch(:,selDim)=repmat(ibpoint,[numPoints sum(selDim)]);      
       if ~all(all(hasdispl(ibpoint,any(numDisplPerEpoch > 0,2),selDim),2),3)
          warning(['Inappropriate refPoint selected ' opt.cbaseTransform ' ' num2str(ibpoint) ' ' pntName{ibpoint}  ' -> fall back to first point per epoch' ]) 
          basePointPerEpoch=firstPointPerEpoch;
       end
end

% Set the computation base for the timeseries (offset parameters)
%
%   opt.cbaseOffset    = { 'firstEpoch' , 'lastEpoch', 'firstCommonEpoch' , 'lastCommonEpoch' , <decimalYear> }  

baseEpochPerPoint=nan(numPoints,numNomParTypes);
switch lower(opt.cbaseOffset)
    case lower('firstEpoch')
       baseEpochPerPoint=firstEpochPerPoint;
    case lower('lastEpoch')
       baseEpochPerPoint=lastEpochPerPoint;
    case lower('firstCommonEpoch')
       baseEpochPerPoint(:,selDim)=repmat(ibaseCandidateEpochs(1),[numPoints sum(selDim)]);
       baseEpochPerPoint(numDisplPerPoint == 0)=NaN;
    case lower('lastCommonEpoch')
       baseEpochPerPoint(:,selDim)=repmat(ibaseCandidateEpochs(end),[numPoints sum(selDim)]);
       baseEpochPerPoint(numDisplPerPoint == 0)=NaN;
    otherwise
       % Use selected epoch 
       refepoch=regexp(opt.cbaseOffset,'[\d\.]*','match','once');
       [~,ibepoch]=min(abs(epochDyears-str2double(refepoch)));
       baseEpochPerPoint(:,selDim)=repmat(ibepoch,[numEpochs sum(selDim)]);      
       if ~all(all(hasdispl(any(numDisplPerPoint > 0,2),ibepoch,selDim),1),3)
           warning(['Inappropriate refEpoch selected ' opt.cbaseOffset ' ' num2str(ibepoch) ' ' num2str(epochDyears(ibepoch)) ' -> fall back to first epoch per point' ])
           baseEpochPerPoint=firstEpochPerPoint;
       end
end

% One, but exactly one, of the computation base parameters, must also be a transformation base parameters

for l=1:numNomParTypes
   if ~selDim(l), continue; end
   ii=find(~isnan(baseEpochPerPoint(:,l)));
   jj=find(~isnan(basePointPerEpoch(:,l)));
   %iscommonbase=(baseEpochPerPoint(basePointPerEpoch(jj,l),l) ==  find(epochMask)');
   iscommonbase=(baseEpochPerPoint(basePointPerEpoch(jj,l),l) ==  jj);
   if sum(iscommonbase) == 0
      % one of the offset and transformation parameters must be the same
      lpoint=basePointPerEpoch(jj(1),l);
      baseEpochPerPoint(lpoint)=1;
   elseif sum(iscommonbase) > 1
      % too many offset and transformation parameters are the same
      ll=find(iscommonbase);
      for lll=2:numel(ll)
          lepoch=jj(ll(lll));
          lpoint=basePointPerEpoch(lepoch,l);
          % take nearest point instead
          lcandidates=find(hasdispl(:,lepoch,l));
          lwcandidates=find(baseEpochPerPoint(lcandidates,l)' ~= lepoch);
          lw=find(lcandidates == lpoint);
          if lw < numel(lcandidates)/2
             basePointPerEpoch(lepoch,l)=lcandidates(lwcandidates(1));
          else
             basePointPerEpoch(lepoch,l)=lcandidates(lwcandidates(end));
          end        
      end
   end
end

% Print the selected reference points and epochs

fprintf('\n\nSelected computing base point(s) and epoch(s):\n\n')
fprintf('Component  pntName                                  pntNum      #Epochs\n')
for l=1:numNomParTypes
   if ~selDim(l), continue; end
   [ibpoint,ibpointcount]=grpstats(basePointPerEpoch(:,l),basePointPerEpoch(:,l),{'mean','numel'});
   for i=1:numel(ibpoint)
      fprintf('%-9s  %-38s %6d       %6d\n',nomParTypes{l},pntName{ibpoint(i)},ibpoint(i),ibpointcount(i))
   end
end
fprintf('\n')
fprintf('Component  epochId       dYear  epochNum  #Points\n')
for l=1:numNomParTypes
   if ~selDim(l), continue; end
   [ibepoch,ibepochcount]=grpstats(baseEpochPerPoint(:,l),baseEpochPerPoint(:,l),{'mean','numel'});
   for j=1:numel(ibepoch)
      fprintf('%-9s  %-12s%4.2f  %6d   %6d \n',nomParTypes{l},epochIds{ibepoch(j)},epochDyears(ibepoch(j)),ibepoch(j),ibepochcount(j))
   end
end   

% Compute logical space time matrix with computation base

compbasis=false(numPoints,numEpochs,numNomParTypes);
for l=1:numNomParTypes
   if ~selDim(l), continue; end
   ii=find(~isnan(baseEpochPerPoint(:,l)));
   jj=find(~isnan(basePointPerEpoch(:,l)));
   iii=sub2ind([ numPoints, numEpochs, numNomParTypes], ii , baseEpochPerPoint(ii,l) , repmat(l,[numel(ii) 1]) );
   jjj=sub2ind([ numPoints, numEpochs, numNomParTypes], basePointPerEpoch(jj,l) , jj , repmat(l,[numel(jj) 1]) );
   compbasis(iii)=true;
   compbasis(jjj)=true;
end

% Select the displacements for inversion

selpar=( hasdispl & ~compbasis );

fprintf('\nNumber of parameters in computation base:  (max per component %d + %d - 1 = %d)\n\n',numPoints,numEpochs,numPoints+numEpochs-1)
fprintf('        ndispl  npar1  nbase\n')
fprintf('        ------ ------ ------\n')
for l=1:numNomParTypes
    fprintf('%5s   %6d %6d %6d\n',nomParTypes{l},sum(sum(hasdispl(:,:,l))),sum(sum(selpar(:,:,l))),sum(sum(compbasis(:,:,l))))
end
fprintf('        ------ ------ ------\n')
fprintf('        %6d %6d %6d\n',sum(hasdispl(:)),sum(selpar(:)),sum(compbasis(:)))

selpar=selpar(:);
npar1=sum(selpar);

fprintf('\nNumber of observations:                           %6d\n',nobs)
fprintf(  'Number of parameters:                             %6d\n',npar1+npar2)
fprintf(  '   + transformations and offsets           %6d\n',npar2)
fprintf(  '   + displacements               %6d -> %6d\n',sum(hasdispl(:)),npar1)
fprintf(  '   - computation base            %6d\n',sum(compbasis(:)))
fprintf(  'Redundancy (dof):                                 %6d\n',nobs-npar1-npar2)

%% Compute solution for the displacement part (x1)
%
% Displacements that have not been observed will be NaN
% The computing base is set to zero, the same is done for the displacements
% which are badly observed (these values are needed later), but these
% should be set to NaN finallly
%
% The following administrative arrays play a role
%
% - parmask(numPoints,numEpochs,numNomParTypes)     logical mask for the parameters in normal equations  
%
%     x1tot, b1tot and N11tot all have dimension numDisplPar=sum(parmask(:))
%
%     idxParMask :  ind=idxParMask(idx) returns the linear position "ind" in the [numPoints numEpochs numNomParTypes] matrix
%                   for the corresponding positions "idx" in N11tot, b1tot and x1tot ,
%                   [i,j,l]=ind2sub([numPoints numEpochs numNomParTypes],ind) would return the subscript values, 
%                   idxParMask has the same length as b1tot and x1tot
%
%     indParMask :  idx=indParMask(ind) returns the linear position "idx" in N11tot, b1tot and x1tot for 
%                   the corresponding positions "ind" (given by ind=sub2ind([numPoints numEpochs NumNomParTypes],i,j,l) )
%                   in the [numPoints numEpochs numNomParTupes] matrix, indParMask has the same length as parMask(:), 
%                   indParMask has NaN's for the elements that do not exist in N11tot, b1tot and x1tot. 
%
% - hasdispl(numPoints,numEpochs,numNomParTypes)    logical mask for the parameters that can be solved
% - compbasis(numPoints,numEpochs,numNomParTypes)   logical mask for the computation base
%
% - selpar(numPoints*numEpochs*numNomParTypes,1)    logical mask for the parameters that can be solved (part of U11), excluding computation base
%
%     selpar=( hasdispl & ~compbasis );
%     selpar=selpar(:);
%
% - selpar2(numDisplPar,1)                          logical mask for the parameters that can be solved (part of U11), excluding computation base,
%                                                   applies to N11tot, b1tot and x1tot
%
%     idxParMask(selpar2)returns the linear position "ind" in the [numPoints numEpochs numNomParTypes] matrix
%     for the corresponding positions "idx" in U11 , 
%     [i,j,l]=ind2sub([numPoints numEpochs numNomParTypes],idxParMask(selpar2)) returns the corresponding subscript values     

x1tot=zeros(numDisplPar,1);
selpar2=hasdispl(idxParMask);
selpar2(compbasis(idxParMask))=false;
[U11,p]=chol(N11tot(selpar2,selpar2));                % N11 = U11'*U11 
if p ~= 0
   size(U11)
   size(N11tot(selpar2,selpar2))
   rankdefect=sum(selpar2)-p+1;
   disp(['Rank-defect of N11 is not as expected, excess rankdefect=' num2str(rankdefect)])
   warning(['Rank-defect of N11 is not as expected, excess rankdefect=' num2str(rankdefect)])
end
x1tot(selpar2)= U11 \ ( U11' \ b1tot(selpar2) );      % x1 = N11\b1 = inv(N11)*b1                

x1=zeros(numPoints*numEpochs*numNomParTypes,1);
x1(idxParMask) =x1tot;

% Condition number and optional plot of diagonal choleski factor

sqrcond=min(diag(U11))/max(diag(U11));
fprintf('\nCondition number displacement part: RCOND = %10.4e\n',sqrcond)

if opt.doplots > 0 || opt.plotCond
   figure(hcond)
   subplot(2,1,2);
   h=semilogy(sort(diag(U11),'descend'),'-x','displayName','Combined displacements');
   h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dataset',repmat({h.DisplayName},size(h.XData)));
   h.DataTipTemplate.Interpreter = 'none';
   title('Cholesky factor U11 diagonal (sorted)')
   ylabel('[-]')
end

% TO DO:
% - covariance matrix computation (we can do this also on demand for
%   selected elements...) -> normal matrix is more practical

%% Step 3 - Back subsitution part to compute final solutions for x2

fprintf('\n\nStep 3 - Solution of the transformation and offset parameters\n')

mse=0;
normy=0;
normrhs2a=0;
normrhs2=0;
nobs_=0;
npar2_=0;
dof=0;

lsqstat2=cell(size(lsqsave));

for k=1:numDatasets
    
    % Residuals and transformation parameters will be saved to dslocal{k}.lsq

    lsq=cell(dslocal{k}.numObsTypes,1);
    
    % Loop over the observation types
    
    for l=1:dslocal{k}.numObsTypes
       
       % reload/recompute idx, b2a A2n A1n N22 from provisional solution (x2a)  save A2n'*A1n as N21 ...
       
       idxlsq=dslocal{k}.idxlsq;

       np=lsqsave{idxlsq+l}.np;
       ne=lsqsave{idxlsq+l}.ne;      
       yn=lsqsave{idxlsq+l}.yn;
       Ry=lsqsave{idxlsq+l}.Ry;
       b2a=lsqsave{idxlsq+l}.b2a;
       A2n=lsqsave{idxlsq+l}.A2n;
       %A1n=lsqsave{idxlsq+l}.A1n;    % A1n=Ry\A1
       U22=lsqsave{idxlsq+l}.U22;
       idxmask=lsqsave{idxlsq+l}.idxmask;
       ymask=lsqsave{idxlsq+l}.ymask;
       sensitivityMatrixType=lsqsave{idxlsq+l}.sensitivityMatrixType;
       sensitivityCols=lsqsave{idxlsq+l}.sensitivityCols;  
       sensitivityMatrix=lsqsave{idxlsq+l}.sensitivityMatrix;

       switch sensitivityMatrixType
           case 'unitcol'
        
               % update the provisional solution -> final solution for x2
   
               % b2 = A2'*inv(Qy)*(y-A1*x1) = b2a-N21*x1 =  -> b2a-A2n'*A1n*x1(idx) 
               % b2 = b2a - A2'*inv(Ry)*inv(Ry)'*A1*x1 = 
               %      b2a - A2'*inv(Ry)*inv(Ry)'*I*x1 = 
               %      b2a - A2n'*inv(Ry)'*x1 = 
               %      b2a - A2n'*x1n
               %
               % x1n= inv(Ry)'*x1 = Ry'\x1

               ll=(sensitivityCols-1)*numPoints*numEpochs;
               x1n=Ry'\x1(ll+idxmask);
               b2=b2a-A2n'*x1n;                               
               %x2=[ U22 \ ( U22' \ b2(1:end-1) ) ; 0 ];     % x2 = N22\b2 = inv(N22)*b2                
      
           case {'selcols','full'}
           
               % b2 = A2'*inv(Qy)*(y-A1*x1) = b2a-N21*x1 =  -> b2a-A2n'*A1n*x1(idx) 
               % b2 = b2a - A2'*inv(Ry)*inv(Ry)'*A1*x1 = 
               %      b2a - A2'*inv(Ry)*inv(Ry)'*D*x1 = 
               %      b2a - A2n'*inv(Ry)'*D*x1 = 
               %      b2a - A2n'*x1n
               %
               % x1n= inv(Ry)'*D*x1 = Ry'\x1t   x1t=D*x1
               
               x1t=zeros(size(idxmask));              
               for i=1:numel(sensitivityCols)
                  iCol=sensitivityCols(i);
                  d=repmat(sensitivityMatrix(:,iCol),[ne 1]);
                  ll=(iCol-1)*numPoints*numEpochs;
                  x1t=x1t+d(ymask).*x1(ll+idxmask);
               end
               x1n=Ry'\x1t;
               b2=b2a-A2n'*x1n;                               
               %x2=[ U22 \ ( U22' \ b2(1:end-1) ) ; 0 ];     % x2 = N22\b2 = inv(N22)*b2        

           otherwise
               error('Hey, this should never happen, we have an illegal case for the sensitivity matrix.')
       end       

       origState = warning;
       warning('off')
       lastwarn('','')
       x2=[ U22 \ ( U22' \ b2(1:end-1) ) ; 0 ];     % x2 = N22\b2 = inv(N22)*b2        
       [msg,warnID] = lastwarn;
       warning(origState)          
       sqrcond2=min(diag(U22))/max(diag(U22));
       if opt.verbose > 1 && ~isempty(msg)
           fprintf ('%2d %s %s %10.3e %s  %s\n',k,dslocal{k}.datasetId,dslocal{k}.obsTypes{l},sqrcond2, ...
                 msg,warnID)
       end

       % compute the least squares residuals and update the sums of the squared resididuals

       % e = y - A1*x1 - A2*x2;
       % en = inv(Ry)'*e = Ry'\e   -> 
       % en = Ry'\y - Ry'\A1*x1 - Ry'\A2*x2 = 
       %      yn - A1n*x1 - A2n*x2 = 
       %      yn - Ry'\x1 - A2n*x2
       %      yn - x1n - A2n*x2

       en = yn - x1n - A2n*x2;
       msek=en'*en;

       % covariance matrix for x2 (inversion is not necessary, can be done also on demand, normal
       % matrix is enough)
       %
       % inv( N22-N21*inv(N11)*N12 )

       % save the results to lsqsave

       lsqsave{idxlsq+l}.b2=b2;
       lsqsave{idxlsq+l}.x2=x2;
       lsqsave{idxlsq+l}.en=en;

       % compute redundancy number taking into account offsets and transformation parameters (npar2) 
       tmp=zeros(numPoints,numEpochs);
       tmp(idxmask)=1;
       tmp=tmp - 1./repmat(sum(tmp,1),[size(tmp,1),1]) - 1./repmat(sum(tmp,2),[1,size(tmp,2)]) + 1/sum(sum(tmp));
       % adjust for displacement parameters (npar1) 
       rn=tmp(idxmask)*(nobs-npar1-npar2)/(nobs-npar2);

       % save the residuals and transformation parameters to dslocal{k}.lsq{l}
       %
       % transformation parameters (x2), the order is given by 
       % - off-sets
       % - transformation parameters
       idxPnt=lsqsave{k}.idxPnt;
       idxEpo=lsqsave{k}.idxEpo;

       lsq{l}.np=numel(idxPnt); 
       lsq{l}.ne=numel(idxEpo); 
       lsq{l}.idxPnt=idxPnt;
       lsq{l}.idxEpo=idxEpo;      
       lsq{l}.x2=x2;
       
       % lsq residuals, normalized residuals and redundancy number, the
       % corresponding pntIds and epochIds is given by idxmask: 
       %
       %    epochNumber=floor(idxmask./numPoints)  -> epochIds(epochNumber)
       %    pntNumber=idxmask-epochNumber.*numPoints  -> pntIds(pntNumber)
       
       lsq{l}.idxmask=idxmask;
       lsq{l}.en=en;
       lsq{l}.e=Ry*en;
       lsq{l}.rn=rn;

       % statistics per observation group (part 2)

       lsqstat2{idxlsq+l}.nobs=numel(yn);
       lsqstat2{idxlsq+l}.npar2=numel(x2)-1;
       lsqstat2{idxlsq+l}.msek=msek;
       lsqstat2{idxlsq+l}.normy=yn'*yn;
       lsqstat2{idxlsq+l}.normrhs2=x2'*b2;

       % update overall model statistics

       mse=mse+msek;
       normy=normy+yn'*yn;
       normrhs2a=normrhs2a+x2(1:end-1)'*b2a(1:end-1);
       normrhs2=normrhs2+x2'*b2;
       nobs_=nobs_+numel(yn);
       npar2_=npar2_+numel(x2)-1;
       dof=dof+numel(yn)-numel(x2)+1;    
       
    end

    dslocal{k}.lsq=lsq;

end

npar1_=numel(x1(selpar));
if npar1_ ~= npar1 || npar2_ ~= npar2 || nobs_ ~= nobs
   warning('mismatch in observation and parameter counts')
end

dof=dof-npar1;
normrhs1=x1(selpar)'*b1tot(selpar2);       

fprintf('\nNumber of observations:                   %6d\n',nobs)
fprintf(  'Number of parameters:                     %6d\n',npar1+npar2)
fprintf(  '   - transformations and offsets   %6d\n',npar2)
fprintf(  '   - displacements (w/o base)      %6d\n',npar1)
fprintf(  'Redundancy (dof):                         %6d\n',nobs-npar1-npar2)

% Overall model test

%omt1 = (e'*inv(Qy)*e)/(m-(n-r));
%omt2 = ( y' * inv(Qy) * y -  x'* A' * inv(Qy) *y ) / (m-(n-r));

omt1=(mse/dof);
omt2=((normy-normrhs2-normrhs1)/dof);

fprintf('\nOverall model test value:                 %.4g\n',omt1)

if (opt.verbose > 1)
   fprintf('\nOverall model test value (alt):           %.4g\n',omt2)
   fprintf('Overall model test value (normy/dof):       %.4g\n',normy/dof)
   fprintf('Overall model test value (normrhs1/dof):    %.4g\n',normrhs1/dof)
   fprintf('Overall model test value (normrhs2/dof):    %.4g\n',normrhs2/dof)
   fprintf('Overall model test value (normrhs2a/dof):    %.4g\n',normrhs2a/dof)
end
fprintf('\n')

% overall lsq statistics (part 1 and part 2)

lsqstat1.nobs=nobs;
lsqstat1.npar2=npar1;
lsqstat1.npar2=npar2;
lsqstat1.omt=omt1;
lsqstat1.dof=dof;
lsqstat1.mse=mse;
lsqstat1.normy=normy;
lsqstat1.normrhs1=normrhs1;
lsqstat1.normrhs2=normrhs2;


%% Save the residuals and transformation parameters (optional), and compute test statistics

[filepath,outputId] = fileparts(outputfilename);
if opt.saveResiduals
   save(fullfile(filepath,[outputId '_res.mat']), 'outputId', ...
       'pntName','pntIds', 'pntCrd', 'pntNeu', 'epochIds', ...
       'lsqstat1', 'lsqstat2', 'dslocal', ...
       'opt', ...
       '-v7.3');
end

%stmresiduals(pntIds, pntNeu, epochIds, dslocal, lsqstat1, opt, outputId)
[omtPoints2,dofPoints2,omtEpochs2,dofEpochs2,omtDatasets,dofDatasets] = ...
       stmresiduals(pntName, pntNeu, epochIds, dslocal, lsqstat1, opt, outputId);

%% Save intermediate results (optional)

if opt.saveIntermediate
   [filepath,outputId] = fileparts(outputfilename);
   save(fullfile(filepath,[outputId '_lsq.mat']), 'outputId', ...
       'pntName','pntIds', 'pntCrd', 'pntNeu', 'epochIds', ...
       'lsqstat1', 'lsqstat2', 'dslocal', ...
       'selpar','selpar2','hasdispl','x1','b1tot','U11', ...
       'lsqsave', ...
       'opt', ...
       '-v7.3');
end

%% Build the output space time matrix dataset
%
% - add estimated displacement corrections to a-priori displacement
% - create and write output stm data structure
% - save only the subset in pntMask and epochMask
% - transformation and offset parameters are saved in the residual files,
%   together with the n(normalized) least squares residua;s
% - the results of statistical testing and other important results are
%   saved in the space time dataset
 
% Reshape the solution (corrections) into space time matrix
%
% - missing elements are NaN
% - unsolved components (e.g. North components) are also NaN, but are 
%   set to zero if a-priori displacements are to be added
% - the computing base is zero
% - optional a-priori displacements are added
% 
% If no a-priori data is available, the dimension is reduced to actually
% estimated components

stmout=x1;
if isempty(opt.aprioriDispl) 
   stmout(~hasdispl)=nan;
   stmObsTypes=parTypes;
   stmObsTypeMask=selDim;
else
   stmout(~parMask)=nan;
   stmObsTypes=nomParTypes;
   stmObsTypeMask=true(size(nomParTypes));
end
stmout=reshape(stmout,[numPoints,numEpochs,numNomParTypes]);

% Add estimated displacement corrrections to a-priori displacement

stmout = aprioriDispl + stmout;

% Analyze the estimated displacement corrections (to the a-priori values)

if ~isempty(opt.aprioriDispl) 
   fprintf('\nRms difference with the a-priori dataset:\n')
   stmdiff(stmout)
   fprintf('\n')
end

% Prepare parameter type/class for output in stm (based on SC)

% Parameter code (absolute value) - to be stored in aux data
%   0 - Not observed
%   1 - Not observed, but necessary for the other components (must have a-priori value)
%   2 - losAsc parameter, depends on a-priori values for two other components
%   3 - losDsc parameter, depends on a-priori values for two other components
%   4 - nearUp or nearEast parameter, depends on a-priori values for one other component
%   7 - Fully estimable parameter
% If negative, not estimated as displacement (hasdispl)
% If > 8 then used for computing base

SC3=single(SC);

pntDisplType=squeeze(max(SC3(pntMask,epochMask,:),[],2));
epochDisplType=squeeze(max(SC3(pntMask,epochMask,:),[],1))';

SC3(compbasis)=SC(compbasis)+8;
SC3(~hasdispl)=-SC(~hasdispl);

% Create stm structure 

st = stm(projectId,'integrated','displ');    

% Point and epoch information (write only points and epochs that took part in the integration)

st.numPoints=numPoints2;
st.numEpochs=numEpochs2;
st.pntName=pntName(pntMask);
st.pntCrd=pntCrd(pntMask,:);
st.epochDyear=epochDyears(epochMask);

% Point and Epoch attributes (including harmonized pntId and epochId)
% - add point class as attributes
% - insert proper obsTypes

pntAttrib=[];
pntAttrib.pntId=pntIds(pntMask);
pntAttrib.pntNeu=pntNeu(pntMask,:);     
pntAttrib.numDisplPerPoint=numDisplPerPoint(pntMask,:);
pntAttrib.omt=omtPoints2;
pntAttrib.dof=dofPoints2;
st.pntAttrib=pntAttrib;

epochAttrib=[];
epochAttrib.epochId=epochIds(epochMask);
epochAttrib.numDisplPerEpoch=numDisplPerEpoch(epochMask,:)';
epochAttrib.omt=omtEpochs2;
epochAttrib.dof=dofEpochs2;
st.epochAttrib=epochAttrib;

% Dataset attributes, set projectFile name and creationDate(projectId is set upon initialization).
  
datasetAttrib=st.datasetAttrib;
datasetAttrib.softwareOptions=opt;
datasetAttrib.projectFile=projectFile;
datasetAttrib.projectFileDate=projectFileDate;
st.datasetAttrib=datasetAttrib;

% Integrated displacements and sensitivity matrix
%
% - the sensitivity matrix could(should?) be changed in case of nearEast and nearUp (and losUp) parameter
%   types, but this is not implemented
% - dimensions are reduced, if applicable, when no a-priori data is available (stmObsTypeMask)

st.obsTypes = stmObsTypes;
st.obsData = stmout(pntMask,epochMask,stmObsTypeMask);
st.sensitivityMatrix=nan(numPoints2,3,3);
st.sensitivityMatrix(:,:,1) = [ ones(numPoints2,1) zeros(numPoints2,1) zeros(numPoints2,1) ];
st.sensitivityMatrix(:,:,2) = [ zeros(numPoints2,1) ones(numPoints2,1) zeros(numPoints2,1) ];
st.sensitivityMatrix(:,:,3) = [ zeros(numPoints2,1) zeros(numPoints2,1) ones(numPoints2,1) ];     
st.sensitivityMatrix=st.sensitivityMatrix(:,:,stmObsTypeMask);

% Stochastic model
%
% - the Choleski factor is only available for the subset of parameters that are actually solved, 
%   this does NOT include the computation base and non-estimable (e.g. North) parameters
% - thus we provide an additional index array...
% 
%     idxParMask(selpar2) returns the linear position "ind" in the [numPoints numEpochs numNomParTypes] matrix
%     for the corresponding positions "idx" in U11 , [i,j,l]=ind2sub([numPoints numEpochs numNomParTypes],idxParMask(selpar2)) 
%     returns the corresponding subscript values     
%
%     subscripts and indices refer to the original arrays, so we still have to apply the pntMask and epochMask
%
% - the Choleski factor can be saved instead of the full inverse (optional), the full matrix is saved,
%   we have only to save the upper triangle of U11 but this is not implemented
%
% - the covariance matrix is computed using chol2stmcov, options for chol2stmcov are 
%
%      opt.byobstype=true;           % split output stochastic model in observation types (default true)
%      opt.addcb=true;               % add zero rows/columns for the computing base (default true)
%      opt.optimize=true;            % use optimized computation in case opt.byobstype=true;
%      opt.sreg = 0.0;               % regularization parameter
%      opt.plot=0;                   % plot intermediate results: 0 (none), 1 or 2 (default is 0)
%      opt.check=false;              % 
%
%   These options are set from stmintegrate options
%      opt.covmatrix                 % {'chol', 'full', 'byobstype'}
%      opt.sreg                      % regularization parameter for cov-matrices
%      opt.verbose
%      opt.doplots

tmp=false([numPoints numEpochs numNomParTypes]);
tmp(idxParMask(selpar2))=true;
tmp=tmp(pntMask,epochMask,stmObsTypeMask);
stochIndex=find(tmp(:));

switch opt.covmatrix
    case 'chol'
       %st.stochModel{1}={'cholfactor(format=index)'};  
       st.stochModel{1}=sprintf('cholfactor(format=index,size=[%d %d %d],lidx=%d,oidx=%d,odata=%d)', ...
           size(st.stochData),numel(stochIndex),0,numel(stochIndex));  
       st.stochData=[ stochIndex(:) ; U11(:) ];
    case 'byobstype'
       [st.stochModel,st.stochData]=chol2stmcov(st.obsTypes, st.obsData, U11, stochIndex, ...
         'byobstype',true,'sreg',opt.sreg,'addcb',true,'optimize',true, ...
         'plot',opt.doplots & opt.verbose>=2,'check',opt.verbose>=3);
    case 'full'
       [st.stochModel,st.stochData]=chol2stmcov(st.obsTypes, st.obsData, U11, stochIndex, ...
         'byobstype',false,'sreg',opt.sreg,'addcb',true, ...
         'plot',opt.doplots & opt.verbose>=2,'check',opt.verbose>=3);
    otherwise
       error('this should not never happen')
end

% Auxiliary data

st.auxTypes = {'SC-North' 'SC-East' 'SC-Up' } ;
numAuxTypes=numel(st.auxTypes);
st.auxData = SC3(pntMask,epochMask,:);

% Input dataset history

inputDatasets=st.inputDatasets;
for k=1:numDatasets
   inputDatasets(k).datasetId=datasets{k}.datasetId;
   inputDatasets(k).techniqueId=datasets{k}.techniqueId;
   inputDatasets(k).datasetAttrib=datasets{k}.datasetAttrib;
   inputDatasets(k).numPoints=datasets{k}.numPoints;
   inputDatasets(k).numEpochs=datasets{k}.numEpochs;
   inputDatasets(k).inputDatasets=datasets{k}.inputDatasets;
end
st.inputDatasets=inputDatasets;

% Save the integrated dataset

stmwrite(st,outputfilename);

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


end

%% Local functions (have their own workspace)

function poly=bbox2poly(bbox)

if size(bbox,1) == 2 && size(bbox,2) == 2 
   % make a polygon out of bounding box
   poly=[ bbox(1,1) bbox(1,2) ; ...    
          bbox(2,1) bbox(1,2) ; ...
          bbox(2,1) bbox(2,2) ; ...
          bbox(1,1) bbox(2,2) ; ...
          bbox(1,1) bbox(1,2) ];    
elseif size(bbox,1) > 2 && size(bbox,2) == 2
   % bbox is already a polygon
   poly=bbox;    
else
   error('This function expects a 2x2 bounding box  [latmin lonmin; latmax lonmax] or polygon.')    
end

% Close the polygon
if poly(1,1) ~= poly(end,1) || poly(1,2) ~= poly(end,2)
   poly = [poly ; poly(1,:) ];
end

end

