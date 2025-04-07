function varargout=stmreducegnss(inputfilename,outputfilename,varargin)
%stmreducegnss   Reduce GNSS time series.
%   STMREDUCEGNSS(INPUTFILENAME,OUTPUTFILENAME,OPTIONS) does a time series
%   reduction of the GNSS space time matrix dataset INPUTFILENAME and
%   creates an output space time matrix OUTPUTFILENAME with the reduction
%   to predefined epochs given in OPTIONS. OPTIONS is a cell array or 
%   structure with the processing options, if empty, or when fields are 
%   missing, default values are used.
%
%   STMREDUCEGNSS(...,UPDATEFLAG) affects the processing in case OUTPUTFILENAME
%   exists. The following values of UPDATEFLAG are accepted
%
%     'create'     abort processing if the output dataset exists
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
%
%   STAT=STMREDUCEGNSS(...) returns a status code STAT. A status code of 0
%   indicates success.
%
%   Options for the decomposition can be given in the structure, or cell 
%   array with option/value pairs, OPTIONS. Supported options are
%
%      verbose           Verbosity level, higher is more output, 0 is almost nothing (default 0)
%      doplots           Plot level, 0 is no plots, higher is more detailed plotting (default 0)
%
%      projectID         The projectID (if not specified, it is the output directory)
%      projectFile       Name of the projectFile if not based on projectId (Default <projectId>.mat) n 
%
%      tolEpoch          Tolerance for matching evaluation epochs in years (default 2 weeks = 14/365)
%
%      ROI               Region of interest, as [latmin lonmin ; latmax lonmax] bounding 
%                        box, or lat/lon polygon, or kml/shape file (Default all)
%      POI               Period of interest [ dYearStart dYearEnd ] or 
%                        [ dYearStart dYearEnd ; dYearStart dYearEnd ;...  ]  (Default all)
%
%      rmafilt           Number of days to be used for robust moving average filter in output space
%                        time matrix (default 0), if less or equal to one, the moving average filter 
%                        is not applied
%      datasetId         DatasetId for the output data space time matrix, if empty, "_reduced" is added 
%                        to the datasetId of the input space time matrix
%
%      globalAttrib      Struct with global attribute updates (default from input file)
% 
%   Examples:
%      stmreducegnss('groningen_GNSS_decomposed.mat','groningen_GNSS_reduced.mat')      
%      stmreducegnss('groningen_GNSS_decomposed.mat','groningen_GNSS_reduced.mat',options,'update')      
%
%   See also GNSS2STM, STMDECOMPOSEGNSS, STM, STMCHECKARGUMENTS, STMREAD and STMDIFF.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:   5 November 2020 by Hans van der Marel
% Modified: 10 March 2021 by Hans van der Marel 
%              - first alpha release combined with decomposition
%           29 August 2021 by Hans van der Marel
%              - split in decomposition and reduction functions
%            5 November 2021 by Hans van der Marel
%              - use the new projectFile

%% Check the input arguments and options

progname='stmreducegnss';

if nargin < 2
   error('This function expects at least two input arguments.')
end

% Default options

opt.verbose=0;                                  % Default verbosity level, higher is more output, 0 is almost nothing
opt.doplots=0;                                  % Default plot level, 0 is no plots, higher is more detailed plotting

opt.projectId='';                               % Project Id (if empty, projectId is determined from the outputfilename directory)
opt.projectFile='';                             % Project File (if empty, it is the same as the projectId with extension .mat),
                                                % if projectFile is used, projectID is set from the projectFile

opt.evaluationMethod='index';                   % Evaluation method, 'index' or 'dist' (default index), 'dist' requires also
opt.tolEpoch=14/365;                            % Tolerance for matching evaluation epochs in years (default 2 weeks)
opt.tolPointDist=.001;                          % Tolerance in distance for matching evaluation points in km (default 1 m)

opt.inputDir='';                                % Directory with the inputfiles

opt.ROI=[];                                     % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default none)
opt.POI=[-Inf +Inf];                            % Period of interest [ dYearStart dYearEnd ] (Default none)

opt.rmafilt=0;                                  % Number of days for final robust moving average filter (default 0), if less or equal to one, no filtering
opt.datasetId='';                               % DatasetId for the output data set, if empty, "_reduced" is added to the input datasetId.
opt.outputTarget='create';                      % Default mode for output STM {'create','overwrite','update'}, can be changed with flag

opt.globalAttrib=struct([]);                    % Struct with global attribute updates (default from input file)


%% Duplicate output to file, and catch errors, and start timing

try

[~,outputfileroot]=fileparts(outputfilename);
diary([ outputfileroot '_' datestr(now,30) '.log' ])

fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Check the options and if necessary overwrite the default values
%
% - on output opt contains the merged default an changed option values
% - on output outputfilename is empty in case no processing has to be done
%   (based on updateflag)
% - in case you have defined multiple input arguments, merge them into
%   a cell array for validation (make it a cell anyway for the checking
%   to work, if you pass inputfilenames as string, the function will
%   think it is a file with input names,

[inputfilename,outputfilename,opt]= ...
    stmcheckarguments({inputfilename},outputfilename,opt,varargin{:});
if isempty(outputfilename)
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    error('A name for the outputfile is expected.');
end
if numel(inputfilename) ~= 1
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    error('Only one inputfile is expected.');
end
inputfilename=char(inputfilename);

%% Get projectId and projectFile

[dataDir,datasetId,~] = fileparts(outputfilename);
if isempty(opt.projectId) 
   dataDirParts = regexp(dataDir, filesep, 'split');
   opt.projectId = dataDirParts{end};
end
projectId=opt.projectId;

if isempty(opt.projectFile)
   opt.projectFile = [opt.projectId '.mat'];
end


%% Get projectId, evaluation epochs and points from the projectFile

if isempty(opt.projectFile)
   error('The name of the project file could not be determined, specify projectId or projectFile as option.')    
end

fprintf('Reading project file %s ... \n',opt.projectFile);

projectFile=stmread(fullfile(opt.inputDir,opt.projectFile));
if ( projectFile.techniqueId ~= 'projectFile' )
   error('The specified projectFile is actually not a project file.')     
end
projectId=projectFile.datasetId;

numEvaluationEpochs=projectFile.numEpochs;
numEvaluationPoints=projectFile.numPoints;

techniqueAttribTotal=projectFile.techniqueAttrib;
campaignDatasetIds=techniqueAttribTotal.campaignDatasetIds;  % {'gron_levelling_flaggedOutliers_v2.3.02_20191030_commonIDs'}
markerDatasetIds=techniqueAttribTotal.markerDatasetIds;      % {'gron_levelling_flaggedOutliers_v2.3.02_20191030_commonIDs'  '06gps_nam_reduced_202107_reduced'}
pntCrdRef=techniqueAttribTotal.pntCrdRef;


pntCrdTotal=projectFile.pntCrd;
pntAttribTotal=projectFile.pntAttrib;
pntNeuTotal=pntAttribTotal.pntNeu;
pntIdTotal=pntAttribTotal.pntId;
             
epochDyearTotal=projectFile.epochDyear;
epochAttribTotal=projectFile.epochAttrib;
epochIdTotal=epochAttribTotal.epochId;


%% Read the GNSS STM (with decomposed time series)

fprintf('Reading input GNSS STM %s ... \n',inputfilename);

st = stmread(inputfilename);

%% Set some frequently used variables and optionally select subset of points and epochs

datasetId=st.datasetId;

numPoints=st.numPoints;
numEpochs=st.numEpochs;

pntCrd=st.pntCrd;
pntAttrib=st.pntAttrib;

epochDyears=st.epochDyear;
epochDate=dyear2date(epochDyears);
epochAttrib=st.epochAttrib;

%% Set Harmonized points ids

fprintf('Reducing GNSS STM ... \n');

% Select regions of interest  (ROI)

fprintf('Number of points in input space time dataset %d\n',numPoints);
if isempty(opt.ROI)
   pntMask=true(numPoints,1);
else
   pntMask=getpntmask(pntCrd,opt.ROI);
   fprintf('Number of points inside ROI %d (out of %d)\n',sum(pntMask),numPoints);
end

if strcmpi(opt.evaluationMethod,'index')
   idxDs=find(strcmpi(datasetId,markerDatasetIds));
   if isempty(idxDs)
      error('No matching datasetId for markers found.')
   end
   %dsPntIndex=pntAttribTotal.dsPntIndex(logical(pntAttribTotal.hasMarker(:,idxDs)),idxDs);
   dsPntIndex=pntAttribTotal.dsPntIndex(:,idxDs);
   hasMarker=~isnan(dsPntIndex);
   pntAttrib.pntId(dsPntIndex(hasMarker),1)=pntIdTotal(hasMarker);    
elseif strcmpi(opt.evaluationMethod,'dist')
   % Compute north, east coordinates using the same reference point as for evaluation points
   pntNeu = plh2neusp(pntCrd,pntCrdRef);  % units are [m] (but stm units are [km]!)
   % Distance and index to total evaluation points (units pf pntNeuTotal are [m], not [km])
   [dist,idx] = pdist2(pntNeuTotal(:,1:2),pntNeu(:,1:2),'euclidean','Smallest',1);
   % get pointIDs and update pointmask
   pntId = pntIdTotal(idx);
   outOfBounds =( dist' > opt.tolPointDist * 1000);  % tolerance in [km]
   pntMask = pntMask & ~outOfBounds;
   pntId(outOfBounds) = num2cell(repmat(' ',[sum(outOfBounds),1]));
   pntAttrib.pntId=pntId;
   pntAttrib.evalPntDist=dist';
   fprintf('Number of points with matching pntIds %d (out of %d)\n',sum(~outOfBounds),numPoints);
   fprintf('Number of points inside ROI and with matching pntIds %d (out of %d)\n',sum(pntMask),numPoints);
else
   error(['Unknown evaluation method option' opt.evaluationMethod])
end

% Check that we have pntIds in the attributes

if ~isfield(pntAttrib,'pntId')
   warning('The input space time dataset contains no pntIds in the attributes and no file with harmonized pntIds was specified we process the data, but later integration runs will not be possible.')    
end

%% Reduce the dataset for integration

% The structure evaluationEpochs contain two fields 
%
%   epochDyearTotal 
%   epochIdTotal

epochIds = epochIdTotal;

% Select period of interest (POI)

fprintf('Number of epochs in input space time dataset %d\n',numEpochs);
if isempty(opt.POI)
   epochMask=true(1,numEpochs);
else
   epochMask=getepochmask(epochDyearTotal,opt.POI);
   fprintf('Number of epochs inside POI %d (out of %d)\n',sum(epochMask),numEpochs);
end

epochDyearTotal=epochDyearTotal(epochMask);
epochIds=epochIds(epochMask);

% Find nearest entry on space time matrix to epochDyearTotal 

[d,epochIdx]=pdist2(st.epochDyear(:),epochDyearTotal(:),'euclidean','smallest',1);

% Remove evaluation epochs that are not close to actual data 

keepers= ( abs(d) < opt.tolEpoch );

epochIds=epochIds(keepers);
epochIdx=epochIdx(keepers); 

fprintf('Number of evaluation epochs %d (out of %d)\n',numel(epochIdx),numEpochs);

% Optionally smooth the residuals using a moving average filter (opt.rmafiltDays > 1)
%      opt.rmafilt=21  number of days in moving average filter
%      opt.rmacrit=5   ....
%      opt.rmastepcrit=0  step criterion (0 is ignore steps)
% Unlike stmdecomposegnss the robust moving average filter operates over
% the whole of the signal (thus including harmonics when included)

obsData=st.obsData;

if opt.rmafilt > 1
   fprintf('Robust moving average filtering (%d days)\n',opt.rmafilt);
   opt.rmacrit=5;    % hardwired option for rmafilt
   opt.rmastepcrit=0;
   for k=1:numPoints
      if pntMask(k)
         tmp=obsData(k,:,:);
         tmp(:,1)=rmafilt(tmp(:,1),opt.rmafilt,'crit',opt.rmacrit,'stepcrit',opt.rmastepcrit);
         tmp(:,2)=rmafilt(tmp(:,2),opt.rmafilt,'crit',opt.rmacrit,'stepcrit',opt.rmastepcrit);
         tmp(:,3)=rmafilt(tmp(:,3),opt.rmafilt,'crit',opt.rmacrit,'stepcrit',opt.rmastepcrit);
         obsData(k,:,:)=tmp;
      end
   end
end

%% Create output space time matrix

fprintf('Creating reduced GNSS STM %s ... \n',outputfilename);

if ~isempty(opt.datasetId)
   datasetIdOut = opt.datasetId;
else
   datasetIdOut = datasetId;
end
if strcmpi(datasetIdOut,st.datasetId)
   datasetIdOut = [ st.datasetId '_reduced' ];
end
   
stout = stm(projectId,datasetIdOut,'gnss');  

techniqueAttrib=stout.techniqueAttrib;
techniqueAttrib.status='Reduced';
stout.techniqueAttrib=techniqueAttrib;

datasetAttrib=stout.datasetAttrib;
datasetAttrib.softwareOptions=opt;
stout.datasetAttrib=datasetAttrib;

stout.numPoints=sum(pntMask);
stout.numEpochs=numel(epochIdx);
stout.pntName=st.pntName(pntMask);
stout.pntCrd=st.pntCrd(pntMask,:);
stout.epochDyear=epochDyears(epochIdx);
 
if ~isempty(pntAttrib)
   pntAttribFields = fieldnames(pntAttrib);
   for k=1:numel(pntAttribFields)
      pntAttribField=pntAttribFields{k};
      tmp=pntAttrib.(pntAttribField);
      pntAttrib.(pntAttribField)=tmp(pntMask,:);
   end
end
stout.pntAttrib=pntAttrib;

epochAttribFields = fieldnames(epochAttrib);
for k=1:numel(epochAttribFields)
   epochAttribField=epochAttribFields{k};
   tmp=epochAttrib.(epochAttribField);
   epochAttrib.(epochAttribField)=tmp(:,epochIdx);
end
epochAttrib.epochId=epochIds;
stout.epochAttrib=epochAttrib;

stout.obsTypes = st.obsTypes;
stout.obsData = obsData(pntMask,epochIdx,:);
stout.sensitivityMatrix = st.sensitivityMatrix(pntMask,:,:);
stout.stochModel = st.stochModel;
stout.stochData = st.stochData;

stout.auxTypes = {};
stout.auxData = [];

% Structure array with data on each input dataset

inputDatasets=struct([]);
inputDatasets(1).datasetId=st.datasetId;
inputDatasets(1).techniqueId=st.techniqueId;
datasetAttrib=st.datasetAttrib;
datasetAttrib.fileName=inputfilename;
inputDatasets(1).datasetAttrib=datasetAttrib;
inputDatasets(1).numPoints=st.numPoints;
inputDatasets(1).numEpochs=st.numEpochs;
inputDatasets(1).inputDatasets=st.inputDatasets;

stout.inputDatasets=inputDatasets;

% Global attributes (copy from input dataset and optionally change).

globalAttrib=st.globalAttrib;
allfields=fieldnames(struct(opt.globalAttrib));
for k=1:numel(allfields)
  globalAttrib.(allfields{k})=opt.globalAttrib.(allfields{k});
end
stout.globalAttrib=globalAttrib;


%% Save the reduced dataset

stmwrite(stout,outputfilename);


%% Finish the function

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

% [End of main]

