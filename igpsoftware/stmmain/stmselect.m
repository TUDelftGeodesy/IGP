function varargout=stmselect(inputfilenames,outputfilename,options,flag)
%stmselect   Select evaluation epochs and points for integrated processing.
%   STMSELECT(INPUTFILENAMES,OUTPUTFILENAME,OPTIONS) selects the evaluation
%   epochs and points for the integrated processing. INPUTFILENAMES must 
%   be the name of an file containing the input stm filenames, or a cell 
%   array with the input stm filenames. OUTPUTFILENAME is a character 
%   string with the name of the projectFile that contains the evaluation 
%   epochs and points for the reduction and integration to follow. 
%   OPTIONS is a cell array or structure with the processing options, if 
%   empty, or when fields are missing, default values are used.
%
%   STMSELECT(...,UPDATEFLAG) affects the processing in case OUTPUTFILENAME
%   exists. The following values of UPDATEFLAG are accepted
%
%     'create'     abort processing if the output dataset exists (default)
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
%
%   STAT=STMSELECT(...) returns a status code STAT. A status code of 0
%   indicates success.
% 
%   The evaluation points are typically based on existing levelling, gravity 
%   and GNSS points, with additional evaluation points selected from the
%   remaining datasets to fill the gaps. The evaluation epochs are based on 
%   typically existing levelling and campaign GNSS epochs, with additional 
%   epochs to fill gaps.
%
%   OPTIONS is a cell array or structure with the processing options, 
%   if empty, or when fields are missing, default values are used. Valid 
%   OPTIONS are 
%
%     projectId=''                Project Id (if empty, then the outputfilename is used as projectId) 
%
%     selectTechnique= { ... }    Cell array with techniques to include (default is {'lev' 'gnss' 'insar'})
%     selectDataset= ''           Cell array with datasetIds to select, you can give a pattern to select, default is all
%
%     selectEpochs={ ... }        Cell array with techniques or modes used for selecting epochs (remaining techniques 
%                                 will be used for desification), default is {'lev' 'grav' 'Campaign') 
%     selectPoints={ ... }        Cell array with techniques or modes used for selecting points (remaining techniques
%                                 will be used for densification), default is {'lev' 'gnss' 'grav'}
%
%     temporalTolerance = 0.1     Time tolerance [decimal years] for selection
%     numEpochsYear = 1;          Number of epochs per year for densification, if negative, no densification
%
%     spatialTolerance = 0.1      Spatial tolerance [km] for merging points in selection
%     numPointsKm2 = 0.5          Density for densification [km2], if negative, no densification
%     spatialMethod = 'grid'      Densification method, only 'grid' for now
%
%     opt.epochFmt = 'yyyymmdd'   Format for epochId, datestr format (e.g. 'yyyymmdd' or 'yyyymm') for dates or 
%                                 numerical (e.g. '%7.2f' or '%8.3f') for dYear. EpochIds always start with 'E'.
%
%     verbose=0                   Verbosity level, higher is more output, 0 is almost nothing
%     doplots=1                   Plot level, higher is more plots, 0 is no plots
%
%     ROI=[]                      Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, 
%                                 lat/lon polygon, or kml/shape file (Default all) *)
%     POI=[-Inf +Inf]             Period of interest [ dYearStart dYearEnd ] or 
%                                 [ dYearStart dYearEnd ; dYearStart dYearEnd ;...  ]   (Default all) *)
%
%     globalAttrib=[]             Structure with global attributes 
%
%   *) Options that are not yet implemented.
%
%   This module is the first of three steps:
%
%     1. STMSELECT     Select evaluation points and epochs
%     2. STMREDUCE, STMREDUCEGNSS and STMREDUCEINSAR
%                      Reduce the space time matrix to the evaluation
%                      epochs and points of the previous steps. Different
%                      functions are needed for GNSS CORS stations, InSAR
%                      and the other datasets (Levelling, GNSS Campaign,
%                      Gravimetry). 
%     3. STMINTEGRATE  Integrate the reduced datasets on the epochs and
%                      points selected by STMSELECT.
%
%   The reduce step for all datasets other than GNSS CORS and InSAR is 
%   simply an update of existing space time datasets with the common epoch 
%   id's and/or point id's, and project id. The reduce step for GNSS CORS is
%   more complicated as it involves a reduction in epochs based on a 
%   decomposition of the time series. The GNSS CORS decomposition and 
%   reduction are two different functions (decomposition is only once, the
%   reduction is often done multiple times). The reduce step for InSAR is
%   the most complicated as both a reduction in time, and space, is
%   necessary. The InSAR reduction and decomposition is done by a single
%   function, which needs to be called for multiple InSAR datasets.
% 
%   The main output of this function is a so-called projectFile. Multiple
%   projectFile can be generated from the same input files, using different
%   options, such that a multitude of integration runs are possible using
%   different evaluation (tie) epochs and points.
%
%   Examples:
%      stmselect('simtrue1_filenames.txt',simtrue1_combined.mat')      
%      stmselect('simtrue1_filenames.txt',simtrue1_combined.mat',options,'update')      
%      stmselect('simtrue1_filenames.txt',simtrue1_combined.mat')      
%
%   See also STMREDUCECAMPAIGN, STMDECOMPOSEGNSS, STMREDUCEGNSS, STMREDUCEINSAR and STMINTEGRATE.
%
%   (c) Freek van Leijen, Hans van der Marel, Delft University of Technology, 2020, 2021.

% Created:   19 October 2021 by Hans van der Marel 
% Modified:  18 September 2020 by Freek van Leijen
%               - original version called sel_evl
%            19 October 2021 by Hans van der Marel
%               - successor to sel_evl (replaces original stmreduce and sel_evl completely)
%               - copy paste from sel_evl, sel_epochs, create_epochs, sel_points
%                 and create_epochs, and edited code to prevent writing, and read
%                 meta data from datasets only once
%               - introduced projectFile (stm alike) to replace both
%                 evaluation*.mat files
%               - updating of actual datasets is deferred to stmreducecampaign 
%                 and is NOT done anymore in this function
%            20 October 2021 by Freek van Leijen
%               - Changed the clustering approach to cluster points within
%                 a dataset and earlier clustered points as well.
%            21 October 2021 by Hans van der Marel
%               - new function to select datasets, selection criteria are now
%                 options
%               - clean up the code and add missing documentation
%            28 October 2021 by Hans van der Marel
%               - ground-up rewrite of sel_epochs
%               - create_epochs changed
%               - epochIds are now done outside sel_epochs and
%                 create_epochs, format is an option
%               - plotting is improved and is now a seperate function, 
%                 using epochTable returned by sel_epochs and updated after 
%                 create_epochs  
%             9 April 2024 by Hans van der Marel
%               - If opt.numEpochsYear is negative, do no densification
%               - If opt.numPointskm2 is negative, do no densification
%              

%% Check the input arguments and options

if nargin < 3
     error('This function expects at least two input arguments.')
end

progname='stmselect';

% Default options 

opt.verbose=0;                                  % Default verbosity level, higher is more output, 0 is almost nothing
opt.doplots=1;                                  % Default plot level, higher is more plots, 0 is no plots

opt.projectId='';                               % Project Id (default is the outputfilename) 

opt.selectTechnique= { 'lev' 'gnss' 'insar' };  % Default is to include all supported techniques
opt.selectDataset= '';                          % Default is to not select on datasetIds, but you can give a pattern to select datasets

% Options that need to be implemented (see stmintegrate for example, pntMask, epochMask)
%
% opt.ROI=[];                     % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon, kml/shape file (Default none)
% opt.POI=[-Inf +Inf];            % Period of interest [ dYearStart dYearEnd ] or [ dYearStart dYearEnd ; dYearStart dYearEnd ; ... ] (Default none)

opt.selectEpochs={'lev' 'grav' 'Campaign'};     % Techniques or modes used for selecting epochs, remaining techniques will be 
                                                % for desification, default is {'lev' 'grav' 'Campaign') 
opt.temporalTolerance = 0.1;                    % [decimal years], time tolerance, for selection
opt.numEpochsYear = 1;                          % [decimal years], number of epochs per year, for densification, if negative, no densification

opt.selectPoints={'lev' 'gnss' 'grav'};         % techniques or modes used for selecting points, remaining techniques will be
                                                % used for densification, default is {'lev' 'gnss' 'grav'}
opt.spatialTolerance = 0.1;                     % spatial tolerance [km] for merging points in selection
opt.numPointsKm2 = 0.5;                         % density for densification [km2]
opt.spatialMethod = 'grid';                     % densification method, only 'grid' for now
opt.densificationROI = [];                      % ROI for the densification points, can be a techniqueID, datasetID, proper ROI or []
opt.shrink=0.5;                                 % shrink factor 0 < shrink < 1, larger is rougher boundary with less points

opt.clusterRadius=0.5;                          % cluster radius in [m]

opt.epochFmt = 'yyyymmdd';                      % Format for epochId, e.g. datestr format 'yyyymmdd' for dates or 
                                                % numerical '%7.2f' for dYear
%opt.epochFmt = '%7.2f';
%opt.epochFmt = 'mmmyyyy';
%opt.epochFmt = 'yyyyQQ'; 

opt.pntFmt = 'olc';                             % Format for pntId (only supported option now is Open Location Code (OLC)
opt.olcLsd = 3;                                 % number of least significant digits for OLC code (see olcencode) (Default is 3)

opt.globalAttrib=struct([]);                    % Global Attributes 

% Duplicate output to file, and catch errors, and start timing

try

[~,outputfileroot]=fileparts(outputfilename);
diary([ outputfileroot '_' datestr(now,30) '.log' ])

fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Check the options and if necessary overwrite the default values

[inputfilenames,outputfilename,opt]= ...
    stmcheckarguments(inputfilenames,outputfilename,opt,options,flag);
if isempty(outputfilename)
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    return;
end

%% Load space time matrix datasets 

% Load meta data from space time datasets 

datasets=[];
for k=1:numel(inputfilenames)
   %datasets{k}=stmread(inputfilenames{k},'NODATA');
   %datasets{k}=stmread(inputfilenames{k},'METADATA');  % <- new option for stmread, as below
   datasets{k}=load(inputfilenames{k}, ...
       'datasetId','datasetAttrib','techniqueId','techniqueAttrib','globalAttrib', ...
       'numPoints','numEpochs',...
       'pntName','pntCrd','pntAttrib', ...
       'epochDyear','epochAttrib',...
       'inputDatasets');
end

% Select the techniques to process

techniqueIds=cellfun(@(x) x.techniqueId,datasets,'UniformOutput',false);
datasetIds=cellfun(@(x) x.datasetId,datasets,'UniformOutput',false);
selected=ismember(techniqueIds,opt.selectTechnique) & contains(datasetIds,opt.selectDataset);

datasets(~selected)=[];

numDatasets=numel(datasets);

techniqueIds=cellfun(@(x) x.techniqueId,datasets,'UniformOutput',false);
datasetIds=cellfun(@(x) x.datasetId,datasets,'UniformOutput',false);

%% Print overview of datasets, with datasetId, TechniqueId and projectInfo
    
fprintf('\nDatasetId                              TechId   Pnts Epochs  createdBy       creationDate    softwareName    projectId       projectFile (projectFileDate)\n')
fprintf(  '-------------------------------------  ------ ------ ------  --------------- --------------- --------------- --------------- -------------------------------\n')
for k=1:numDatasets
   datasetId=datasets{k}.datasetId;
   if length(datasetId) > 38
      datasetId=[datasetId(1:35) '...'];
   end
   datasetAttrib=datasets{k}.datasetAttrib;
   fprintf('%-38s %-5s  %6d %6d  %-15s %-15s %-15s %-15s %s (%s)\n',datasetId,datasets{k}.techniqueId,datasets{k}.numPoints,datasets{k}.numEpochs, ...
          datasetAttrib.createdBy,datasetAttrib.creationDate,datasetAttrib.softwareName,datasetAttrib.projectId,datasetAttrib.projectFile,datasetAttrib.projectFileDate);
end
fprintf('\n')

%% Select the datasets for sel_epochs and sel_points

datasetMaskEpochs = sel_datasets(datasets,opt.selectEpochs); % Typically {'lev' 'grav' 'Campaign'} 
datasetMaskPoints = sel_datasets(datasets,opt.selectPoints); % Typically {'lev' 'gnss' 'grav'}


%% Get existing epochs from campaign datasets (datasetMaskEpochs)

fprintf('\nData reduction: selecting existing epochs from ... \n');
for k=1:numDatasets
   if datasetMaskEpochs(k)
      fprintf('- %s  (%s)\n',datasetIds{k},techniqueIds{k})
   end
end
fprintf('\n')

%[epochDyearExist,epochTable] = sel_epochs(datasets(datasetMaskEpochs),opt.temporalTolerance);
[epochDyearExist,epochTable] = sel_epochs(datasets(datasetMaskEpochs),opt);

if opt.doplots > 1
   plot_epoch_table(datasets(datasetMaskEpochs),epochTable)
end

%% Add new epochs

if opt.numEpochsYear > 0

   fprintf('\nData reduction: adding newly created epochs ... \n');

   % POI of is determined from datasets(~datasetMaskEpochs)

   [epochDyearTotal,epochDyearNew] = create_epochs(datasets(~datasetMaskEpochs), ...
      epochDyearExist,opt.numEpochsYear);
   %[epochDyearTotal,epochDyearNew] = create_epochs(~,epochDyearExist,opt.numEpochsYear,POI);

else
   epochDyearTotal=epochDyearExist;
   epochDyearNew=[];
end

% Update epochTable

epochTable2 = sortrows([ epochTable ; [ epochDyearNew' zeros(numel(epochDyearNew),1) nan(numel(epochDyearNew),size(epochTable,2)-2) ]],1);
if opt.verbose > 0
   epochTable2
end


%% Print summary output and optionally plot

% Summary

fprintf('\nNumber of epochs:   %d\n',numel(epochDyearTotal))
fprintf('- original   %5d\n',sum((epochTable2(:,2)==1)))
fprintf('- merged     %5d\n',sum((epochTable2(:,2)>1)))
fprintf('- new        %5d\n',numel(epochDyearNew))
epochIntervals=diff(epochDyearTotal);
fprintf('Min/mean/median/max interval between epochs: %.2f / %.2f / %.2f / %.2f [yr]\n',min(epochIntervals),mean(epochIntervals),median(epochIntervals),max(epochIntervals))
fprintf('MAD/StDev of intervals:  %.2f / %.2f [yr]\n',1.4826*median(abs( epochIntervals - median(epochIntervals))),std(epochIntervals))

% Plot

if opt.doplots > 0
   plot_epoch_table(datasets(datasetMaskEpochs),epochTable2)
end


%% Create array with epochIds

% Create epochId using the format specfied in opt.epochFmt

numEpochTotal=numel(epochDyearTotal);
if opt.epochFmt(1) == '%'        % EpochId is in decimal year format 
   epochIdTotal = cellstr([repmat('E',numEpochTotal,1) num2str(epochDyearTotal',opt.epochFmt)])';
else
   epochIdTotal = cellstr([repmat('E',numEpochTotal,1) datestr(dyear2date(epochDyearTotal'),opt.epochFmt)])';
end

% Find duplicate epochIds and rename

[epochIdUnique,~,idx2] = unique(epochIdTotal);
multipleEpochId = epochIdUnique(histcounts(idx2,1:numel(epochIdUnique)+1)>1);
for k = 1:numel(multipleEpochId)
    idx = strmatch(multipleEpochId{k},epochIdTotal);
    for m = 1:numel(idx)
      epochIdTotal(idx(m)) = {[epochIdTotal{idx(m)} char(double('a')+m-1)]};
    end
end

if opt.verbose > 0
   epochIdTotal
end


%% Get existing points from campaign and CORS datasets (datasetMaskPoints)

fprintf('\nData reduction: selecting existing points from ... \n');
for k=1:numDatasets
   if datasetMaskPoints(k)
      fprintf('- %s  (%s)\n',datasetIds{k},techniqueIds{k})
   end
end
fprintf('\n')

%[pntCrdExist,pntTable,pntNeuExist,pntCrdRef] = sel_points_legacy(datasets(datasetMaskPoints),opt.spatialTolerance);
[pntCrdExist,pntTable,pntNeuExist,pntCrdRef] = sel_points(datasets(datasetMaskPoints),opt.spatialTolerance);

if opt.verbose > 2
   pntTable
end

%% Print overview of merged points

numMerged=sum(pntTable(:,3)>1);
if numMerged > 0 && ( opt.verbose > 1 || numMerged < 150 )
   fprintf('\nNumber of existing points:   %d\n',size(pntCrdExist,1))
   fprintf('- original   %5d\n',sum((pntTable(:,3)==1)))
   fprintf('- merged     %5d\n',sum((pntTable(:,3)>1)))
   fprintf('\nMerged points ...\n\n')
   fprintf('  Latitude  Longitude Num  pntName   distance [m]  ...\n\n')
   ldsidx=find(datasetMaskPoints);
   for k=1:size(pntTable,1)
       if pntTable(k,3) <=1, continue; end
       fprintf('%10.6f %10.6f  %2d  ',pntCrdExist(k,1:2),pntTable(k,3));
       for l=1:numel(ldsidx)
            idxl=pntTable(k,l*3+1);
            if ~isnan(idxl)
                distl=sqrt((pntTable(k,l*3+2)-pntTable(k,1))^2 + (pntTable(k,l*3+3)-pntTable(k,2))^2 )*1000;
                fprintf('%-15s%7.1f  ',datasets{ldsidx(l)}.pntName{idxl},distl)
            else
                fprintf('                          ');
            end
       end
       fprintf('\n')
   end
end

%% Add new points

if opt.numPointsKm2 > 0

    fprintf('\nData reduction: adding newly created points ... \n');
    
    % [pntCrdTotal,pntNeuTotal,pntIdTotal] = create_points(datasets(~datasetMaskPoints), ...
    %     pntCrdExist,pntNeuExist,pntIdExist,pntCrdRef,opt.numPointsKm2,opt.spatialMethod);
    
    % Region of interest for densification points
    
    if ischar(opt.densificationROI) || iscell(opt.densificationROI) 
       % Check if it matches a techniqueID or datasetID , the set ROI based on extend of that technique
       techniqueIds=cellfun(@(x) x.techniqueId,datasets,'UniformOutput',false);
       datasetIds=cellfun(@(x) x.datasetId,datasets,'UniformOutput',false);
       selected=find(ismember(techniqueIds,opt.densificationROI) | contains(datasetIds,opt.densificationROI));
       dummyNEU=[];
       for k=1:numel(selected)
           dummyNEU=[dummyNEU ; plh2neusp(datasets{selected(k)}.pntCrd,pntCrdRef) ];
       end
       if ~isempty(dummyNEU)
    %      opt.polyBuffer=
    %      roi = polybuffer(dummyNEU(:,1:2),'points',opt.polyBuffer) 
          kk = boundary(dummyNEU(:,1:2),opt.shrink);
          roi=dummyNEU(kk,1:2);
       else
          roi=[];
       end
    elseif ~isempty(opt.densificationROI)
       % Can be a file or ROI in latitude / longitude, convert to NEU coordinates...
       roi=roi2poly(opt.densificationROI);
       roi=plh2neusp([ roi  zeros(size(roi,1),1) ],pntCrdRef);
       roi=roi(:,1:2);
    else
       roi=[];
    end

    % Compute new points
    
    switch opt.spatialMethod
       case 'grid'
          pntNeuNew = create_points_grid(pntNeuExist,opt.numPointsKm2,roi);
        case 'none'
          pntNeuNew = [];
       case {'triangulation','clustering'}
          error('Unimplemented point creation method.');
       otherwise
          error('You defined a non-existent point creation method.');
    end
    pntCrdNew=neu2plhsp(pntNeuNew,pntCrdRef);

    % Update pntCrdTotal and pntTable
    
    pntCrdTotal=[pntCrdExist;pntCrdNew];
    pntNeuTotal=[pntNeuExist;pntNeuNew];
    
    pntTable2 = [ pntTable ; [ pntNeuNew(:,1:2) ./ 1000 zeros(size(pntNeuNew,1),1) nan(size(pntNeuNew,1),size(pntTable,2)-3) ]];
    if opt.verbose > 1
       pntTable2
    end

else
    pntCrdTotal=pntCrdExist;
    pntNeuTotal=pntNeuExist;
    pntTable2=pntTable;
end

%pntCrdTotal=pntCrdExist;
%pntNeuTotal=pntNeuExist;
%pntTable2=pntTable;

%% Print summary output and optionally plot

% Summary

fprintf('\nNumber of points:   %d\n',size(pntCrdTotal,1))
fprintf('- original   %5d\n',sum((pntTable2(:,3)==1)))
fprintf('- merged     %5d\n',sum((pntTable2(:,3)>1)))
fprintf('- new        %5d\n',sum((pntTable2(:,3)==0)))
%epochIntervals=diff(epochDyearTotal);
%fprintf('Min/mean/median/max interval between epochs: %.2f / %.2f / %.2f / %.2f [yr]\n',min(epochIntervals),mean(epochIntervals),median(epochIntervals),max(epochIntervals))
%fprintf('MAD/StDev of intervals:  %.2f / %.2f [yr]\n',1.4826*median(abs( epochIntervals - median(epochIntervals))),std(epochIntervals))

% Plot

if opt.doplots > 0
   plot_pnt_table(datasets(datasetMaskPoints),pntTable2)
end


%% Create array with pntIDs

if strcmpi(opt.pntFmt,'olc')
   % Open Location Code with opt.lsd least significant digits (duplicates are handled by olcencode)
   pntIdTotal=cellstr(olcencode(pntCrdTotal(:,1),pntCrdTotal(:,2),opt.olcLsd)); % e.g. opt.olcLsd=3
else
   error('unsupported pntId format.')
end

% Sort on pntId

[pntIdTotal,idxTotal] = sort(pntIdTotal);
pntCrdTotal=pntCrdTotal(idxTotal,:);
pntNeuTotal=pntNeuTotal(idxTotal,:);
pntTable2=pntTable2(idxTotal,:);

%% Find point clusters (points that within a certain radius)
%
% opt.clusterRadius [m] is input
% T is an array with cluster indices
%

T=clusterdata(pntNeuTotal(:,1:2),'Criterion','distance','Cutoff',opt.clusterRadius);

[clusterCrd,clusterCount]=grpstats(pntCrdTotal,T,{'mean', 'numel'});
clusterCount=clusterCount(:,1);

clusterPntId(T)=pntIdTotal;              % clusterId  (same Id as one of the points in the cluster)
clusterPntIndex(T)=1:numel(pntIdTotal);  % index to the point in pntCrd
clusterPntCrd(T,:)=pntCrdTotal;          % If you want to use the same coordinates as the point with clusterId (instead of the mean)

pntClusterCount=clusterCount(T);         % Number of clustered points (same length as everything else in pntAttrib) 
pntClusterIndex=T;                       % Index to the point cluster
pntIsPrimary=false(numel(pntIdTotal),1);
pntIsPrimary(clusterPntIndex)=true;      % True if the point is the primary point of a cluster



%% Save projectId, epoch and point list to integration project definition file  
%
% The structure of the project definition file resembles that of the space
% time matrix, but does not contain any data
%
%   FieldName       Description                     Type                       Values/Examples/Notes
%   --------------  ------------------------------- ------------------------   -----------------------
%   datasetId       Solution Identifier             char string (free)         same as projectId !!
%   techniqueId     Technique Identifier            char string (reserved)     'projectFile' 
%   techniqueAttrib Technique Attributes            struct   
%   datasetAttrib   Dataset Attributes              struct                     
%   inputDatasets() Input dataset history           struct array 
%
%   numPoints       Number of Points                int scalar
%   numEpochs       Number of Epochs                int scalar
%
%   pntName         Point name                      cell array [numPoints]     
%   pntCrd          Point Cooordinates (lat,lon,h)  double [numPoints,3] matrix
%   pntAttrib       Point Attributes ....           table  [numPoints,* ]      
%      .pntId          Harmonized Point Identifier     cell array [numPoints]                 
%
%   epochDyear    Decimal year                      double [numEpochs] matrix 
%   epochAttrib   Epoch Attributes                  table [numEpochs,*]        
%      .epochId      Harmonized Epoch identifier       cell array [numEpochs]   
%
%  The above fields may change !!!! have to contain the same information
%  that is written to the old evaluation epochs and points

%save(fullfile(projectId,'evaluation_epochs'),'epochDyearTotal','epochIdTotal');
%save(fullfile(projectId, 'evaluation_points'),'pntCrdTotal','pntNeuTotal','pntIdTotal','pntCrdRef');
  
% Initialize projectID

if isempty(opt.projectId)
    [~,projectId]=fileparts(outputfilename);
else
    projectId=opt.projectId;
end

% datasetAttrib and techniqueAttrib structures
%
%   datasetAttrib   Dataset Attributes              struct
%      .createdBy      Username of the account running the software
%      .creationDate   File creation date in ISO format
%      .createdAs      File name at the time the file was created, with full path
%      .softwareName   Name of the software that created the dataset/file
%      .softwareOptionsOptions used to produce this dataset/file
%      .fileFormat     File format
%      .fileName       File name at the time the file was read
%      .projectId      Project Id or name  
%      .projectFile    Project file name 
%      .projectFileDate  Creation/modification date of the project file    
%                        in ISO format

%createdBy=getenv('username');
createdBy={getenv('USER'),getenv('USER'),'unknown'}; 
createdBy=createdBy{~cellfun('isempty',createdBy)};
creationDate=datestr(now,30);
dbstackCall=dbstack(1);
if ~isempty(dbstackCall)
   softwareName=dbstackCall.name;
else
   softwareName='';
end

datasetAttrib = struct( ...
    'createdBy',createdBy,...
    'creationDate',creationDate,...
    'createdAs',outputfilename,...
    'softwareName',softwareName,...
    'softwareOptions',opt,...
    'fileFormat','mat',...
    'projectId',projectId,...
    'projectFile',outputfilename,...
    'projectFileDate',creationDate );

techniqueAttrib=[];
techniqueAttrib.campaignDatasetIds=datasetIds(datasetMaskEpochs);
techniqueAttrib.markerDatasetIds=datasetIds(datasetMaskPoints);
techniqueAttrib.pntCrdRef=pntCrdRef;
%techniqueAttrib.projCrs=projCrs;

% Point attributes

dsPntIndex=pntTable2(:,4:3:end);          % array with index to pntCrd in each campaignDataset
dsMarkerNames = repmat({''},size(dsPntIndex));
ldsidx=find(datasetMaskPoints);
for k=1:numel(ldsidx)
   idxDs= dsPntIndex(:,ldsidx(k));
   idx= idxDs > 0;
   idxDs(~idx)=[];
   dsMarkerNames(idx,k) = datasets{ldsidx(k)}.pntName(idxDs);
end

pntAttrib=[];
pntAttrib.pntId=pntIdTotal;
pntAttrib.hasMarker=pntTable2(:,3);    % number of markerDatasets that link to this point 
pntAttrib.dsPntIndex=dsPntIndex;       % array with index to pntCrd in each campaignDataset
pntAttrib.dsMarkerNames=dsMarkerNames; % array with the original pntName for each markerDataset
pntAttrib.pntNeu=pntNeuTotal; 
pntAttrib.clusterCount=pntClusterCount;% Number of clustered points (same length as everything else in pntAttrib) 
pntAttrib.clusterIndex=pntClusterIndex;% Index to the point cluster
pntAttrib.isPrimary=pntIsPrimary;      % Primary point (in a cluster)
%pntAttrib.pntCrdProj=

clusterAttrib=[];
clusterAttrib.pntId=clusterPntId;         % clusterPntId  (same Id as one of the points in the cluster)
clusterAttrib.pntIndex=clusterPntIndex;   % index to the point in pntCrd
clusterAttrib.pntCrd=clusterPntCrd;       % If you want to use the same coordinates as the point with clusterId (instead of the mean)


% Epoch attributes

epochAttrib=[];
epochAttrib.epochId=epochIdTotal;
epochAttrib.hasCampaign=epochTable2(:,2)'; % number of campaignDatasets that link to this epoch 
epochAttrib.dsEpochIndex=epochTable2(:,3:2:end)'; % array with index to epochDyear in each campaignDataset
epochAttrib.dsEpochDyear=epochTable2(:,4:2:end)'; % array with the original epochDyear for each campaignDataset
%epochAttrib.campaignNames= % array with the campaign name (if available) for each campaignDataset

% inputDatasets()  Structure array with for each input dataset  struct array 
%      .datasetId
%      .techniqueId
%      .datasetAttrib
%      .numPoints
%      .numEpochs
%      .inputDatasets

inputDatasets=[];
for k=1:numDatasets
   inputDatasets(k).datasetId=datasets{k}.datasetId;
   inputDatasets(k).techniqueId=datasets{k}.techniqueId;
   inputDatasets(k).datasetAttrib=datasets{k}.datasetAttrib;
   inputDatasets(k).numPoints=datasets{k}.numPoints;
   inputDatasets(k).numEpochs=datasets{k}.numEpochs;
   inputDatasets(k).inputDatasets=datasets{k}.inputDatasets;
end

% Save the project file

projectFile = struct(  ...
    'datasetId',projectId, ...
    'techniqueId','projectFile',...
    'techniqueAttrib',techniqueAttrib,...
    'datasetAttrib',datasetAttrib ,...
    'numPoints',size(pntCrdTotal,1),...
    'numEpochs',size(epochDyearTotal,2),...
    'pntName',{pntIdTotal},...
    'pntCrd',{pntCrdTotal},...
    'pntAttrib',pntAttrib,...
    'clusterAttrib',clusterAttrib,...
    'epochDyear', {epochDyearTotal},...
    'epochAttrib',epochAttrib,...
    'inputDatasets',{ inputDatasets } ,...
    'globalAttrib',opt.globalAttrib);


projectFile


%save(outputfilename,'-v7.3','-struct','projectFile')   % v7.3 stores as hdf5, lot of overhead (10 Mb instead of 0.4 Mb file)
save(outputfilename,'-struct','projectFile') 

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

function datasetMask=sel_datasets(datasets,selectTechnique)
%SEL_DATASETS   Selection of datasets for.
% DATASETMASK=SEL_DATASETS(DATASETS,SELECTTECHNIQUE) selects datasets 
% from DATASETS for use by the selection function sel_epochs, create_epochs 
% and sel_points. The DATASETS variable is a structure array with meta 
% data for the datasets for ALL datasets in the project. 
%
% The relevant datasets are selected based on the techniqueId and 
% techniqueAttribs.mode (if available) using an "or" operation.
%
% The ouput is boolean array DATASETMASK with the selected datasets.
%
% Examples:
%
%    datasetMask = sel_datasets(datasets,{'lev' 'grav' 'Campaign'}); % for sel_epochs (only gnss campaigns)
%    datasetMask = sel_datasets(datasets,{'insar' 'CORS'});          % for create_epochs, is opposite of previous
%    datasetMask = sel_datasets(datasets,{'lev' 'gnss' 'grav'});     % for sel_points (all gnss)
%
% (c) Freek van Leijen, Hans van der Marel, Delft University of Technology, 2021. 

% Created:  21 October 2021 by Hans van der Marel
% Modified: 21 October 2021 by Hans van der Marel
%              - original code, written by Freek, adapted for this function
%              - added documentation

numDatasets = numel(datasets);

datasetMask=false(numDatasets,1);
for k = 1:numDatasets
    techniqueId=datasets{k}.techniqueId;
    techniqueAttrib=datasets{k}.techniqueAttrib;
    
    techniqueMode = {''}; % To extract the technique mode 
    if ~isempty(techniqueAttrib)
      if isfield(techniqueAttrib,'mode')
        techniqueMode = techniqueAttrib.mode;
      end
    end

    if ismember(techniqueId,selectTechnique) || ismember(techniqueMode,selectTechnique)
       datasetMask(k)=true;      
    end
    
end

end

function [epochDyearExist,epochTable] = sel_epochs(datasets,temporalTolerance) 
%SEL_EPOCHS   Selection of evaluation epochs.
% [EPOCHDYEAREXIST,EPOCHTABLE]=SEL_EPOCHS(DATASETS,TEMPORALTOLERANCE) selects 
% existing epochs based on levelling, gravity and GNSS campaign datasets 
% passed in DATASETS to this function. Epochs of different datasets within 
% the TEMPORALTOLERANCE are merged. Epochs within a dataset are not merged, 
% but can be merged with the same epochs in other datasets.
%
% The DATASETS variable is a structure array with the metadata of the datasets
% selected by the function sel_datasets prior to this function.
%
% The TEMPORALTOLERANCE is the maximum time span allowed between epochs
% (epochs with smaller time span should be merged). The value should be 
% given in decimal years, e.g., 0.1 .
%
% [EPOCHDYEAREXIST,EPOCHTABLE]=SEL_EPOCHS(DATASETS,OPT) does the same, but
% accepts as second argument a structure with options. Two fields are
% supported, 'verbose' and 'temporalTolerance'.
%
% The output consist of EPOCHDYEAREXIST with the decimal years of the
% unique existing and merged epochs, together with a matrix EPOCHTABLE
% linking the existing/merged epochs to the epochs in the respective
% datasets. The columns in EPOCHTABLE are respectively EPOCHDYEAREXIST,
% number of datasets that have this epoch, the epoch number in the first 
% dataset, the dYear of the first dataset, the epoch number in the second 
% dataset ... and so on. If epochs are not present in a dataset NaN's
% are inserted.
%
% Examples:
%
%    [epochDyearExist,epochTable] = sel_epochs(datasets,temporaltolerance) 
%    [epochDyearExist,epochTable] = sel_epochs(datasets(datasetMaskEpochs),opt.temporalTolerance);
%    [epochDyearExist,epochTable] = sel_epochs(datasets,opt.temporalTolerance)
%    [epochDyearExist,epochTable] = sel_epochs(datasets,opt)
%
% (c) Freek van Leijen, Hans van der Marel, Delft University of Technology, 2020, 2021. 

% Created:  06 September 2020 by Freek van Leijen
% Modified: 10 October 2020 by Freek van Leijen
%              - added header and annotation
%           21 October 2021 by Hans van der Marel
%              - became internal function to stmselect
%              - datasets are read outside this function (once) and passed
%                as structure array (relevant fields only)
%              - datasets are selected outside this function by sel_datasets 
%                and only relevant datasets are passed in DATASETS
%              - epochIds are now done outside this function, the second
%                output argument is an array linking the output epochs to
%                the relevant dataset epochs
%              - return EPOCHTABLE as 2nd output argument (for plotting)
%              - moved plotting outside of this function
%              - updating of actual datasets is deferred to stmreduce and
%                is NOT done anymore in this function
%           24 October 2021 by Hans van der Marel
%              - implemented new algorithm based on Matlab pdist2, this
%                solved problem linking nearby epochs and handles multiple
%                datasets properly
%              - new algorithm tested in test_sel_epochs_new.m
%              - second argument can be a structure with options

% Set the options 

opt.verbose=0;
if isstruct(temporalTolerance)
   opt=temporalTolerance;
   if ~isfield(opt,'temporalTolerance')
       error('Mandatory option temporalTolerance missing.')
   end
else
   opt.temporalTolerance = temporalTolerance;
end

% The new method of finding matching epochs works on pairs of datasets 
% and used Matlab's pdist2 function to find matching pairs with the 
% shortest time interval.
% 
% The distance and index between the values in X and Y are given by
%   
%   [D,I] = pdist2(X,Y,DISTANCE,'Smallest',1) returns a 1-by-MY matrix D
%   containing the smallest pairwise distances to observations in X for
%   each observation in Y, and returns a 1-by-MY matrix I containing 
%   indices of the observations in X corresponding to the smallest pairwise 
%   distances in D. To find out which observations are paired use
%
%      paired = ( dist' < opt.temporalTolerance) 
%   
%   Y(paired) and X(idx(paired)) should then be paired (Y should be the smaller of the 
%   two arrays for the algorithm to return only the minimum number of links)
% 
% An example implementation is given below:
%
%   epochDyear_i = datasets{i}.epochDyear;
%   epochDyear_j = datasets{j}.epochDyear;
%
%   if numel(epochDyear_j) <= numel(epochDyear_i) 
%      [dist,idx] = pdist2(epochDyear_i',epochDyear_j','euclidean','Smallest',1);
%      paired =( dist' < opt.temporalTolerance);
%      epochDyearMerged = ( epochDyear_i(idx(paired)) + epochDyear_j(paired) ) ./ 2 ;
%      epochDyear_i_notpaired=epochDyear_i;
%      epochDyear_i_notpaired(idx(paired))=[];
%      epochDyear_j_notpaired=epochDyear_j(~paired);
%   else
%      [dist,idx] = pdist2(epochDyear_j',epochDyear_i','euclidean','Smallest',1);
%      paired =( dist' < opt.temporalTolerance);
%      epochDyearMerged = ( epochDyear_i(paired) + epochDyear_j(idx(paired)) ) ./ 2 ;
%      epochDyear_i_notpaired=epochDyear_i(~paired);
%      epochDyear_j_notpaired=epochDyear_j;
%      epochDyear_j_notpaired(idx(paired))=[];
%   end
% 
%   epochDyearExist = sort([ epochDyear_i_notpaired epochDyearMerged epochDyear_j_notpaired ]) 
%
% For more than two datasets the above procedure should be iterated. We do
% this by forming every possible pair and updating two tables (per dataset)
%
%   for k=1:numDatasets
%      for l=k+1:numDatasets
%         % update for the datasets k and l tables with pointes and Dyears         
%      end
%   end
%
% The tables are for each dataset are
%
%   dsadmin(k).epochTable1[numEpochsDs,numDatasets]  with pointers to other datasets or zero if not a matching pair
%   dsadmin(k).epochTable2[numEpochsDs,numDatasets]  with nearby Dyears for other datasets or NaN
%
% The charm is that is you take the mean (ignoring NaN's) over epochTable2 
% you have the common epoch in Dyear.

% Initialize the temporary administrative structure array dsadmin

numDatasets=numel(datasets);

dsadmin=[];
for k=1:numDatasets
   % Copy epochDyears to dsadmin (store as column arrays for pdist2)
   numepochs=numel(datasets{k}.epochDyear);
   dsadmin(k).numepochs=numepochs;
   dsadmin(k).epochDyear=datasets{k}.epochDyear(:);     
   dsadmin(k).epochTable1=zeros(numepochs,numDatasets);
   dsadmin(k).epochTable2=nan(numepochs,numDatasets);
end

% Fill the tables in the adminstrative structure array dsadmin

for k=1:numDatasets
   for l=k+1:numDatasets
      % Make sure the 2nd dataset for pdist2 is the smaller dataset
      if dsadmin(l).numepochs <= dsadmin(k).numepochs
         i=k; j=l;
      else
         i=l; j=k;
      end
      [dist,idx] = pdist2(dsadmin(i).epochDyear,dsadmin(j).epochDyear,'euclidean','Smallest',1);
      paired =( dist' < opt.temporalTolerance);
      dsadmin(i).epochTable1(idx(paired),j)=find(paired);
      dsadmin(j).epochTable1(paired,i)=idx(paired);
      dsadmin(i).epochTable2(idx(paired),j)=dsadmin(j).epochDyear(paired);
      dsadmin(j).epochTable2(paired,i)=dsadmin(i).epochDyear(idx(paired));
   end
   dsadmin(k).epochTable2(:,k)=dsadmin(k).epochDyear;
end

if opt.verbose > 0
   for k=1:numDatasets
      disp(['epochTable1 and epochTable2 for dataset ' num2str(k)])
      dsadmin(k).epochTable1
      dsadmin(k).epochTable2
   end
end

% Create array with common epochs (this is a column vector).

epochDyear=[];
for k=1:numDatasets
   epochDyear = [ epochDyear ; nanmean(dsadmin(k).epochTable2,2) ];
end
epochDyear=unique(epochDyear);

% Create the epochTable

% The existing/merged epochs are linked to the epochs in the original 
% datasets through the matrix epochTable.
%
% The number of columns in epochTable is numDatasets*2 + 1, with the 
%   1. epochDyearExist 
%   2. number of dataset that have this epoch
%   3. the epoch number in the first dataset
%   4. the dYear of the first dataset
%   5. the epoch number in the second dataset 
%      ... and so on
%
% The number of rows in epochTable is initially the same as the number of
% rows in epochArray, but eventually becomes the same size as epochDyearExist

epochTable = nan(numel(epochDyear),numDatasets*2+2);
epochTable(:,1) = epochDyear';
epochTable(:,2) = zeros(numel(epochDyear,1));
for k=1:numDatasets
   % Distance and index to epochDyear
   dsEpochDyear=dsadmin(k).epochDyear;
   [dist,idx] = pdist2(epochDyear,dsEpochDyear,'euclidean','Smallest',1);
   % Update the table
   epochTable(idx,2)=epochTable(idx,2)+1;
   epochTable(idx,k*2+1)=[1:numel(dsEpochDyear)]';
   epochTable(idx,k*2+2)=dsEpochDyear;
   % Check for out of bounds
   outOfBounds =( dist' > opt.temporalTolerance);  
   fprintf('Number of matching ids  %d (out of %d) for dataset %d, max distance is %.4f years\n',sum(~outOfBounds),numel(dsEpochDyear),k,max(dist));
end

% Return epochDyearExist as row vector

epochDyearExist=epochDyear';

if opt.verbose > 0
   disp('epochDyearExist')
   epochDyearExist
   disp('epochTable')
   epochTable
end

end


function [epochDyearTotal,epochDyearNew] = create_epochs(datasets,epochDyearExist,  ...
              numEpochsYear,yearRange)
%CREATE_EPOCHS Create epochs to sub-sample the existing epochs.
%   [EPOCHDYEARTOTAL,EPOCHDYEARNEW] = CREATE_EPOCHS(DATASETS, ...
%      EPOCHDYEAREXIST,NUMEPOCHSYEAR,YEARRANGE)
%   creates additional evaluation epochs by subsampling the earlier
%   selected existing epochs EPOCHDYEAREXIST for the time span set by
%   YEARRANGE. The sampling frequency is given by NUMEPOCHSYEAR.
%
%   YEARRANGE is optional. If YEARRANGE is missing, or empty, the date
%   range is determined from the datasets passed in DATASETS (e.g. InSAR and 
%   CORS GNSS datasets selected by sel_datasets from the full range of 
%   input DATASETS).
%
%   The created epochs are given as EPOCHDYEARTOTAL and EPOCHDYEARNEW, the 
%   latter only includes newly created epochs.
%
%  (c) Freek van Leijen, Delft University of Technology, 2020. 

% Created:  18 September 2020 by Freek van Leijen
% Modified: 21 October 2021 by Hans van der Marel
%              - became internal function to stmselect
%              - datasets are read outside this function (once) and passed
%                as structure array (relevant fields only)
%              - datasets are selected outside this function by sel_datasets 
%                and only relevant datasets are passed in DATASETS
%              - epochId is now created outside this function, EPOCHIDEXIST
%                and EPOCHIDTOTAL are removed fron the function call
%              - added YEARRANGE to the input arguments
%              - added EPOCHDYEARNEW to output arguments

if nargin < 4
   yearRange=[]; 
end

numDatasets=numel(datasets);
deltaTime = 1/numEpochsYear;

%% Get relevant time span (typically only gnns cors and insar datasets are passed to this function)

if isempty(yearRange)
   minDate=inf;
   maxDate=-inf;
   for k=1:numDatasets
      epochDyear=datasets{k}.epochDyear;      
      minDate = min([minDate epochDyear]);
      maxDate = max([maxDate epochDyear]);
   end
else
   minDate=min(yearRange(:));
   maxDate=max(yearRange(:));
end

%% Check whether minDate and maxDate should be added
%
% x----|--x--x-x----x-x-xx-x-xx-x----x---x..........|
% e   min e  e e    e e ee e ee e    e   e         max
%
% minDate and or maxDate should be added if time to
% closest existing epoch is > deltaTime. Maybe not best
% solution, since data reduction will not be optimal.

epochDyearNew=[];
if min(abs(epochDyearExist-minDate))>deltaTime
  epochDyearExist=[epochDyearExist minDate];
  epochDyearNew=[epochDyearNew minDate];
end
if min(abs(epochDyearExist-maxDate))>deltaTime
  epochDyearExist=[epochDyearExist maxDate];
  epochDyearNew=[epochDyearNew maxDate];
end
epochDyearExist=sort(epochDyearExist);

 
%% Adapt to Period of Interest
% todo

%% Temporary remove epochs before and after minDate and maxDate
% (not useful to generate epochs before and after these dates.)
minIdx=find(epochDyearExist<minDate);
maxIdx=find(epochDyearExist>maxDate);
epochDyearOutside=[epochDyearExist(minIdx) epochDyearExist(maxIdx)];
epochDyearExist([minIdx maxIdx])=[];

%% Determine new epochs
maxDiff=inf;
while maxDiff>deltaTime
  epochDyearExist=sort(epochDyearExist);
  [maxDiff,maxIdx]= max(abs(diff(epochDyearExist)));
  if maxDiff>deltaTime
     epochDyearExist=[epochDyearExist mean([epochDyearExist(maxIdx) ...
         epochDyearExist(maxIdx+1)])];
     epochDyearNew=[epochDyearNew mean([epochDyearExist(maxIdx) ...
         epochDyearExist(maxIdx+1)])];
  end
end

%% Create output

epochDyearTotal=sort([epochDyearExist epochDyearOutside]);


end 



function [pntCrdExist,pntTable,pntNeuExist,pntCrdRef] = sel_points(datasets,spatialTolerance)
%SEL_POINTS  Selection of evaluation points. 
%  [PNTCRDEXIST,PNTTABLE,PNTNEUEXIST,PNTCRDREF]=SEL_POINTS(DATASETS,SPATIALTOLERANCE).
%  Selection of evaluation points based on typically the existing levelling, 
%  GNSS and gravity points. The selection is based on the space time matrix 
%  DATASETS passed to this function (a subset of all possible datasets 
%  selected by sel_datasets). Points of different datasets within the 
%  SPATIALTOLERANCE [km] are merged. Points within a dataset are not
%  merged, but can be merged with nearby epochs in other datasets.
%
%  The DATASETS variable is a structure array with the metadata of the datasets
%  selected by the function sel_datasets prior to this function.
%
%  The SPATIALTOLERANCE is the maximum distance allowed between points
%  (points with smaller distance should be merged). The value should be 
%  given in kilometers, e.g., 0.1 .
%
%  [...]=SEL_POINTS(DATASETS,OPT) does the same, but accepts as second 
%  argument a structure with options. Two fields are supported, 'verbose' 
%  and 'spatialTolerance'.
%
%  The output consist of PNTCRDEXIST with the latitude and longitude of the
%  unique existing and merged epochs, together with a matrix PNTTABLE
%  linking the existing/merged points to the points in the respective
%  datasets. The columns in PNTTABLE are respectively lat and lon in
%  PNTDYEAREXIST, number of datasets that have this points, the point number
%  in the first dataset, the lat and lon of the first dataset, the point
%  number in the second dataset ... and so on. If points are not present 
%  in a dataset NaN's are inserted.
%
%  Examples:
%
%     [pntCrdExist,pntTable,pntNeuExist,pntCrdRef] = sel_points(datasets,spatialtolerance) 
%     [pntCrdExist,pntTable,pntNeuExist,pntCrdRef] = sel_points(datasets(datasetMaskPoints),opt.spatialTolerance);
%     [pntCrdExist,pntTable,pntNeuExist,pntCrdRef] = sel_points(datasets,opt.temporalTolerance)
%     [pntCrdExist,pntTable,pntNeuExist,pntCrdRef] = sel_pointss(datasets,opt)
%
%  (c) Freek van Leijen, Hans van der Marel, Delft University of Technology, 2020, 2021. 

% Created:  18 September 2020 by Freek van Leijen
% Modified: 19 Oct 2021 by Freek van Leijen
%              - Changed the clustering approach to cluster points within
%                a dataset and earlier clustered points as well.
%           21 October 2021 by Hans van der Marel
%              - became internal function to stmselect
%              - datasets are read outside this function (once) and passed
%                as structure array (relevant fields only)
%              - datasets are selected outside this function by sel_datasets 
%                and only relevant datasets are passed in DATASETS
%              - changed units of SPATIALTOLERANCE to km
%              - updating of actual datasets is deferred to stmreduce and
%                is NOT done anymore in this function
%            1 November 2021 by Hans van der Marel
%              - pntIds are now done outside this function, the third
%                output argument replaced by an array linking the output points to
%                the relevant dataset points
%              - return PNTTABLE as 2nd output argument (for plotting)
%              - moved plotting outside of this function
%              - implemented new algorithm based on Matlab pdist2, this
%                solved problem linking nearby points and handles multiple
%                datasets properly
%              - new algorithm tested in test_sel_points_new.m
%              - second argument can be a structure with options

% Set the options 

opt.verbose=0;
if isstruct(spatialTolerance)
   opt=spatialTolerance;
   if ~isfield(opt,'spatialTolerance')
       error('Mandatory option spatialTolerance missing.')
   end
else
   opt.spatialTolerance = spatialTolerance;
end

% The new method of finding matching points works on pairs of datasets 
% and used Matlab's pdist2 function to find matching pairs with the 
% shortest distance.
%
% The algorithm is an adaptation from sel_epochs. For a more extensive
% description of the method see sel_epochs.

% Get reference position for pntNeu computation

pntCrd=[];
numDatasets=numel(datasets);
for k=1:numDatasets
   pntCrd = [pntCrd ; datasets{k}.pntCrd ];
end
pntCrdRef=mean(pntCrd,1); 

% Initialize the temporary administrative structure array dsadmin

dsadmin=[];
for k=1:numDatasets
   % Copy pntCrd and pntNeu to dsadmin 
   numpoints=size(datasets{k}.pntCrd,1);
   dsadmin(k).numpoints=numpoints;
   dsadmin(k).pntCrd=datasets{k}.pntCrd;
   dsadmin(k).pntNeu=plh2neusp(datasets{k}.pntCrd,pntCrdRef);   
   dsadmin(k).pntTable0=zeros(numpoints,numDatasets);
   dsadmin(k).pntTable1=nan(numpoints,numDatasets);
   dsadmin(k).pntTable2=nan(numpoints,numDatasets);
   dsadmin(k).pntTable3=nan(numpoints,numDatasets);
end

% Fill the tables in the adminstrative structure array dsadmin

for k=1:numDatasets
   for l=k+1:numDatasets
      % Make sure the 2nd dataset for pdist2 is the smaller dataset
      if dsadmin(l).numpoints <= dsadmin(k).numpoints
         i=k; j=l;
      else
         i=l; j=k;
      end
      [dist,idx] = pdist2(dsadmin(i).pntNeu(:,1:2),dsadmin(j).pntNeu(:,1:2),'euclidean','Smallest',1);
      paired =( dist' < opt.spatialTolerance * 1000);
      dsadmin(i).pntTable0(idx(paired),j)=find(paired);
      dsadmin(j).pntTable0(paired,i)=idx(paired);
      dsadmin(i).pntTable1(idx(paired),j)=dsadmin(j).pntCrd(paired,1);
      dsadmin(j).pntTable1(paired,i)=dsadmin(i).pntCrd(idx(paired),1);
      dsadmin(i).pntTable2(idx(paired),j)=dsadmin(j).pntCrd(paired,2);
      dsadmin(j).pntTable2(paired,i)=dsadmin(i).pntCrd(idx(paired),2);
      dsadmin(i).pntTable3(idx(paired),j)=dsadmin(j).pntCrd(paired,3);
      dsadmin(j).pntTable3(paired,i)=dsadmin(i).pntCrd(idx(paired),3);
   end
   dsadmin(k).pntTable1(:,k)=dsadmin(k).pntCrd(:,1);
   dsadmin(k).pntTable2(:,k)=dsadmin(k).pntCrd(:,2);
   dsadmin(k).pntTable3(:,k)=dsadmin(k).pntCrd(:,3);
end

if opt.verbose > 1
   for k=1:numDatasets
      disp(['pntTable0, pntTable1 and pntTable2 for dataset ' num2str(k)])
      dsadmin(k).pntTable0
      dsadmin(k).pntTable1
      dsadmin(k).pntTable2
      dsadmin(k).pntTable3
   end
end

% Create array pntCrd with common points and pntAdmin with index to original points

pntCrd=[];
pntAdmin=[];
for k=1:numDatasets
   pntCrd = [ pntCrd ; nanmean(dsadmin(k).pntTable1,2) nanmean(dsadmin(k).pntTable2,2) nanmean(dsadmin(k).pntTable3,2) ];
   dsadmink = dsadmin(k).pntTable0;
   dsadmink(:,k) = 1:size(dsadmink,1);
   pntAdmin = [ pntAdmin ; dsadmink ];
end
pntMerged=sum(pntAdmin > 0,2) > 1;
pntCrdMerged=[ pntCrd(pntMerged,:) pntAdmin(pntMerged,:)];
pntCrdMerged = unique(pntCrdMerged,'rows');
pntCrd= [ pntCrd(~pntMerged,:) ; pntCrdMerged(:,1:3) ];
pntAdmin= [ pntAdmin(~pntMerged,:) ; pntCrdMerged(:,4:end) ];

pntNeu=plh2neusp(pntCrd,pntCrdRef);   

if opt.verbose > 1
   size(pntCrd)
   size(pntNeu)
   size(pntAdmin)
   pntAdmin
end


% Create the pntTable

% The existing/merged points are linked to the points in the original 
% datasets through the matrix pntTable.
%
% The number of columns in pntTable is numDatasets*3 + 3, with the 
%   1. latitude of pntCrdExist 
%   2. longitude of pntCrdExist 
%   3. number of dataset that have this point
%   4. the point number in the first dataset
%   5. the latitude of the first dataset
%   6. the longitude of the first dataset
%   7. the point number in the second dataset 
%      ... and so on

numPoints=size(pntCrd,1);
pntTable = nan(numPoints,numDatasets*3+3);
pntTable(:,1:2) =pntNeu(:,1:2) ./ 1000;
pntTable(:,3) = zeros(numPoints,1);
for k=1:numDatasets
   % Original dataset coordinates
   dsPntNeu=dsadmin(k).pntNeu;
   % Distance and index to pntCrd
   %
   % Similar to 
   %
   %    [dist,idx] = pdist2(pntNeu(:,1:2),dsPntNeu(:,1:2),'euclidean','Smallest',1)
   %    idxDS = 1:size(dsPntNeu,1)
   %   
   % but handles points with same coordinates well
   idxDs= pntAdmin(:,k);
   idx= idxDs > 0;
   idxDs(~idx)=[];
   dist=sqrt( (pntNeu(idx,1)-dsPntNeu(idxDs,1)).^2 + (pntNeu(idx,2)-dsPntNeu(idxDs,2)).^2 );
   % Update the table
   pntTable(idx,3)=pntTable(idx,3)+1;
   pntTable(idx,k*3+1)=idxDs;
   pntTable(idx,k*3+2)=dsPntNeu(idxDs,1) ./ 1000;
   pntTable(idx,k*3+3)=dsPntNeu(idxDs,2) ./ 1000;
   % Check for out of bounds
   outOfBounds =( dist > opt.spatialTolerance * 1000);  
   fprintf('Number of matching ids  %d (out of %d) for dataset %d, max distance is %.4f m\n',sum(~outOfBounds),size(dsPntNeu,1),k,max(dist));
   %fprintf('Number of matching ids  %d (out of %d) for dataset %d, max distance is %.4f m\n',sum(~outOfBounds),size(dsPntNeu,1),k,max(dist(~outOfBounds)));
end

% Return epochDyearExist as row vector

pntCrdExist=pntCrd;
pntNeuExist=pntNeu;

if opt.verbose > 0
   disp('pntCrdExist')
   pntCrdExist
   disp('pntTable')
   pntTable
end

end


function pntNeuNew = create_points_grid(pntNeuExist,numPointsKm2,roi)
%CREATE_POINTS_GRID Creation of evaluation points using the grid approach.
%  PNTNEUNEW = CREATE_POINTS_GRID(PNTNEUEXIST,NUMPOINTSKM2,ROI) creates
%  evaluation points using the grid approach for the region of interest
%  ROI. The existing PNTNEUEXIST points are subsampled based on the 
%  NUMPOINTSKM2 density. The newly created points PNTNEUNEW are passed as 
%  output.
%
%  PNTNEUNEW = CREATE_POINTS_GRID(PNTNEUEXIST,NUMPOINTSKM2) uses a 
%  bounding box around PNTNEUEXIST as ROI. Same result as with ROI=[]; 
%
%  (c) Freek van Leijen, Delft University of Technology, 2020. 

% Created:  18 September 2020 by Freek van Leijen
% Modified: 21 October 2021 by Hans van der Marel
%              - became internal function to stmselect
%              - moved plotting outside this function
%           21 October 2021 by Hans van der Marel
%              - added ROI as input

% The creation of points is based a density indicated by the
% NUMPOINTSKM2 parameter, using one of the SPATIALMETHODs
% - grid
% - triangulation (NOT IMPLEMENTED YET)
% - clustering (NOT IMPLEMENTED YET) 
% This function provides the grid based method
%
% Grid based
%
% 1) create grid based on nPointsKm2. E.g, 4, then create grid of 500 m
% 2) determine gridcells NOT CONTAINING an existing benchmark

% Check input arguments

if nargin < 3
   roi=[];
end

%% Determine maximum extend for the densification grid

if isempty(roi)
   xmin = min(pntNeuExist(:,2));
   xmax = max(pntNeuExist(:,2));
   ymin = min(pntNeuExist(:,1));
   ymax = max(pntNeuExist(:,1));
   roi = [xmin ymin ; xmax ymax ];
else
   roi = roi2poly(roi);
   xmin = min(roi(:,2));
   xmax = max(roi(:,2));
   ymin = min(roi(:,1));
   ymax = max(roi(:,1));
end

%% Determine grid size 
% (125*2^x m)
% e.g. nPointsKm2 = 4
% 4/1000^2 = 1/x^2
% x = sqrt(1000^2/4)
% xx = floor(log2(x/125)) (or maybe round?)
% gridSize = 125*2^xx

gridSize = 125*2^floor(log2(sqrt(1000^2/numPointsKm2)/125));


%% Setup grid

xGridMin = floor(xmin/gridSize)*gridSize;
xGridMax = ceil(xmax/gridSize)*gridSize;
yGridMin = floor(ymin/gridSize)*gridSize;
yGridMax = ceil(ymax/gridSize)*gridSize;

Nx = round((xGridMax-xGridMin)/gridSize);
Ny = round((yGridMax-yGridMin)/gridSize);

gridIdx = [kron((1:Nx)',ones(Ny,1)) kron(ones(Nx,1),(1:Ny)')];


%% Determine empty grid cells

pointsIdx = [ceil((pntNeuExist(:,2)-xGridMin)/gridSize) ceil((pntNeuExist(:,1)-yGridMin)/gridSize)];

newIdx = setdiff(gridIdx,pointsIdx,'rows');

newX = newIdx(:,1)*gridSize + xGridMin - 0.5*gridSize;
newY = newIdx(:,2)*gridSize + yGridMin - 0.5*gridSize;
pntNeuNew=[newY newX zeros(size(newX))];
 
%% Select grid cells only within ROI

pntMask=getpntmask(pntNeuNew,roi);
pntNeuNew(~pntMask,:)=[];

end

function plot_epoch_table(datasets,epochTable)

figure;hold on;
numDatasets=numel(datasets);
numEpoch=size(epochTable,1);
epochDyearTotal=epochTable(:,1);
epochDyearNew=epochDyearTotal(epochTable(:,2)==0);
epochDyearSingle=epochDyearTotal(epochTable(:,2)==1);
epochDyearMerged=epochDyearTotal(epochTable(:,2)>1);
xx=[];
for k=1:numDatasets
   x = epochTable(:,k*2+2)';
   y = (numDatasets-k+1)*ones(1,numEpoch);
   plot(x,y,'x','Markersize',5,'DisplayName',datasets{k}.datasetId);
   xx = [xx x];
end   
plot(xx,zeros(1,numel(xx)),'k.','markersize',5,'DisplayName','Original');
plot(epochDyearSingle,zeros(1,numel(epochDyearSingle)),'rs','markersize',7,'DisplayName','Single');
plot(epochDyearMerged,zeros(1,numel(epochDyearMerged)),'ro','markersize',7,'DisplayName','Merged');
plot(epochDyearNew,zeros(1,numel(epochDyearNew)),'b+','markersize',7,'DisplayName','New');
plot(epochDyearTotal,-1*ones(1,numel(epochDyearTotal)),'k*','markersize',5,'DisplayName','Total');
ylim([-2 numDatasets+3])
legend('Interpreter','None');
xlabel('dyear')
title('Selected epochs (closeby epochs merged)');

end


function plot_pnt_table(datasets,pntTable)

figure;hold on;
numDatasets=numel(datasets);
pntCrdTotal=pntTable(:,1:2);
pntCrdNew=pntCrdTotal(pntTable(:,3)==0,:);
pntCrdSingle=pntCrdTotal(pntTable(:,3)==1,:);
pntCrdMerged=pntCrdTotal(pntTable(:,3)>1,:);
for k=1:numDatasets
   x = pntTable(:,k*3+3);
   y = pntTable(:,k*3+2);
   plot(x,y,'.','Markersize',5,'DisplayName',datasets{k}.datasetId);
end   
%plot(pntCrdSingle(:,2),pntCrdSingle(:,1),'rs','markersize',5,'DisplayName','Single');
plot(pntCrdMerged(:,2),pntCrdMerged(:,1),'ro','markersize',7,'DisplayName','Merged');
plot(pntCrdNew(:,2),pntCrdNew(:,1),'g+','markersize',3,'DisplayName','New');
%plot(pntCrdTotal(:,2),pntCrdTotal(:,1),'k*','markersize',5,'DisplayName','Total');
axis equal;
xlabel('East [km]')
ylabel('North [km]')
legend('Interpreter','None','Location','best');
title('Selected points (closeby points merged)');

end





