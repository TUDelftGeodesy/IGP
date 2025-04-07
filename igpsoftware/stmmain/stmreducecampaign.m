function varargout=stmreducecampaign(inputfilename,outputfilename,varargin)
%stmreducecampaign   Reduce campaign datasets.
%   STMREDUCECAMPAIGN(INPUTFILENAME,OUTPUTFILENAME,OPTIONS) does a reduction 
%   of Levelling, GNSS or Gravity campaign space time matrix dataset 
%   INPUTFILENAME and creates an output space time matrix OUTPUTFILENAME 
%   with the proper point and campaign ids. The output space time matrix
%   is basically a copy of the input matrix, the only values that are
%   added are the pntAttrib.pntId and epochAttrid.epochId. OPTIONS is a cell 
%   array or structure with the processing options, if empty, or when 
%   fields are missing, default values are used.
%
%   STMREDUCECAMPAIGN(...,UPDATEFLAG) affects the processing in case OUTPUTFILENAME
%   exists. The following values of UPDATEFLAG are accepted
%
%     'create'     abort processing if the output dataset exists
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
%
%   STAT=STMREDUCECAMPAIGN(...) returns a status code STAT. A status code of 0
%   indicates success.
%
%   Options for the reduction can be given in the structure, or cell 
%   array with option/value pairs, OPTIONS. Supported options are
%
%      verbose           Verbosity level, higher is more output, 0 is almost nothing (default 0)
%      doplots           Plot level, 0 is no plots, higher is more detailed plotting (default 0)
%
%      projectID         The projectID (if not specified, it is the output directory)
%      projectFile       Name of the projectFile if not based on projectId (Default <projectId>.mat) n 
%
%      evaluationMethod  Evaluation method, 'index' or 'dist' (default index), 'dist' requires 
%                        also the two parameters below
%      tolEpoch          Tolerance for matching evaluation epochs in years (default 2 weeks = 14/365)
%      tolPointDist      Tolerance in distance for matching evaluation points in km (default 0.001 km)
%
%      ROI               Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, 
%                        or lat/lon polygon (Default none)
%      POI               Period of interest [ dYearStart dYearEnd ] (Default [-Inf +Inf])
%
%      datasetId         DatasetId for the output data space time matrix, if empty, "_reduced" is added 
%                        to the datasetId of the input space time matrix
%
%      globalAttrib      Struct with global attribute updates (default from input file)
% 
%   Examples:
%      stmreducecampaign('groningen_levelling.mat','groningen_levelling_reduced.mat')      
%      stmreducecampaign('groningen_levelling.mat','groningen_levelling_reduced.mat',options,'update')      
%
%   See also STMSELECT, STMREDUCEGNSS and STMREDUCEINSAR.
%
%  (c) Hans van der Marel, Delft University of Technology, 2021.

% Created:   4 November 2021 by Hans van der Marel
% Modified: 

%% Check the input arguments and options

progname='stmreducecampaign';

%if nargin < 2
%   error('This function expects at least two input arguments.')
%end

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


%% Read the campaign STM 

fprintf('Reading input GNSS STM %s ... \n',inputfilename);

st = stmread(inputfilename);

%% Set some frequently used variables and optionally select subset of points and epochs

datasetId=st.datasetId;

numPoints=st.numPoints;
numEpochs=st.numEpochs;

pntCrd=st.pntCrd;
pntAttrib=st.pntAttrib;

epochDyear=st.epochDyear;
epochAttrib=st.epochAttrib;

% Select regions of interest  (ROI)

fprintf('Number of points in input space time dataset %d\n',numPoints);
if isempty(opt.ROI)
   pntMask=true(numPoints,1);
else
   pntMask=getpntmask(pntCrd,opt.ROI);
   fprintf('Number of points inside ROI %d (out of %d)\n',sum(pntMask),numPoints);
end

% Select period of interest (POI)

fprintf('Number of epochs in input space time dataset %d\n',numEpochs);
if isempty(opt.POI)
   epochMask=true(1,numEpochs);
else
   epochMask=getepochmask(epochDyear,opt.POI);
   fprintf('Number of epochs inside POI %d (out of %d)\n',sum(epochMask),numEpochs);
end

%% Set Harmonized points ids

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


%% Set harmonized epoch ids

if strcmpi(opt.evaluationMethod,'index')
   idxDs=find(strcmpi(datasetId,campaignDatasetIds));
   if isempty(idxDs)
      error('No matching datasetId for epochs found')
   end
   %dsEpochIndex=epochAttribTotal.dsEpochIndex(idxDs,logical(epochAttribTotal.hasCampaign(idxDs,:)));
   dsEpochIndex=epochAttribTotal.dsEpochIndex(idxDs,:);
   hasCampaign=~isnan(dsEpochIndex);
   epochAttrib.epochId(dsEpochIndex(hasCampaign))=epochIdTotal(hasCampaign);    
elseif strcmpi(opt.evaluationMethod,'dist')
   % Find nearest entry on space time matrix to epochDyearTotal 
   [dist,idx]=pdist2(epochDyear(:),epochDyearTotal(:),'euclidean','smallest',1);
   % get epochIDs and update epochMask
   epochId = epochIdTotal(idx);
   outOfBounds =( dist' > opt.tolEpoch);  % tolerance in [dyears]
   epochMask = epochMask & ~outOfBounds;
   epochId(outOfBounds) = num2cell(repmat(' ',[sum(outOfBounds),1]));
   epochAttrib.epochId=epochId;
   epochAttrib.evalEpochDist=dist;
   fprintf('Number of epochs with matching epochIds %d (out of %d)\n',sum(~outOfBounds),numEpochs);
   fprintf('Number of epochs inside POI and with matching epochIds %d (out of %d)\n',sum(epochMask),numEpochs);
else
   error(['Unknown evaluation method option' opt.evaluationMethod])
end


%% Create output space time matrix

fprintf('Creating reduced campaign STM %s ... \n',outputfilename);

if ~isempty(opt.datasetId)
   datasetIdOut = opt.datasetId;
else
   datasetIdOut = [ datasetId '_reduced' ];
end
   
stout = stm(projectId,datasetIdOut,st.techniqueId);  

techniqueAttrib=stout.techniqueAttrib;
techniqueAttrib.status='Reduced';
stout.techniqueAttrib=techniqueAttrib;

datasetAttrib=stout.datasetAttrib;
datasetAttrib.softwareOptions=opt;
stout.datasetAttrib=datasetAttrib;

stout.numPoints=sum(pntMask);
stout.numEpochs=numel(epochMask);
stout.pntName=st.pntName(pntMask);
stout.pntCrd=st.pntCrd(pntMask,:);
stout.epochDyear=epochDyear(epochMask);
 
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
   epochAttrib.(epochAttribField)=tmp(:,epochMask);
end
stout.epochAttrib=epochAttrib;

stout.obsTypes = st.obsTypes;
stout.obsData = st.obsData(pntMask,epochMask,:);
stout.sensitivityMatrix = st.sensitivityMatrix(pntMask,:,:);
stout.stochModel = st.stochModel;
stout.stochData = st.stochData;
if ~isempty(st.stochData) && ( ~all(pntMask) || ~all(epochMask) )
   error('POI and ROI not implemented for stochastic data.')    
end
stout.auxTypes = stout.auxTypes;
if ~isempty(st.auxTypes)
   stout.auxData = st.auxData(pntMask,epochMask,:);
else
   stout.auxData = [];
end

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

