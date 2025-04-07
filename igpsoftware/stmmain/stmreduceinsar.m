function varargout = stmreduceinsar(inputfilename,outputfilename,varargin)
%stmreduceinsar  Data reduction InSAR (including decomposition).
%   STMREDUCEINSAR(INPUTFILENAME,OUTPUTFILENAME,OPTIONS) Function to
%   reduce an InSAR dataset, specified with INPUTFILENAME, to a 
%   pre-determined list of evaluation points and epochs. These
%   points and epochs are loaded via the [projectId]/evaluation_points
%   and [projectId]/evaluation_epochs .mat files. The resulting
%   space-time dataset is saved as the OUTPUTFILENAME file.
%
%   INPUTFILENAME is a cell array with the input filename. 
%   OUTPUTFILENAME is a character string with the name of the output 
%   Space Time Matrix dataset. OPTIONS is a cell array or
%   structure with the processing options, if empty, or when fields are 
%   missing, default values are used.
%
%   STMREDUCEINSAR(...,UPDATEFLAG) affects the processing in case OUTPUTFILENAME
%   exists. The following values of UPDATEFLAG are accepted
%
%     'create'     abort processing if the output dataset exists
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
% 
%   STAT=STMREDUCEINSAR(...) returns a status code STAT. A status code of 0
%   indicates success.
%
%   Options for the data reduction can be given in the structure, or cell 
%   array with option/value pairs, OPTIONS. Supported options are
%
%      verbose              Verbosity level, higher is more output, 0 is almost nothing (default 0)
%      doplots              Plot level, 0 is no plots, higher is more detailed plotting (default 0)
%
%      projectId            The projectId (if not specified, it is the output directory)
%      projectFile          Name of the projectFile if not based on projectId (Default <projectId>.mat) 
%
%      inputDir             Directory with the inputfiles
%      datasetId            DatasetId for the output data space time matrix, if empty, "_reduced" is added 
%                           to the datasetId of the input space time matrix
%      outputTarget         Mode for output STM {'create','overwrite','update'}, can be changed with flag,
%                           default is 'create';
%
%      ROI                  Region of interest, as [latmin lonmin ; latmax lonmax] bounding 
%                           box, or lat/lon polygon, or kml/shape file (Default all)
%      POI                  Period of interest [ dYearStart dYearEnd ] or 
%                           [ dYearStart dYearEnd ; dYearStart dYearEnd ;...  ]  (Default all)
%
%      insarRadius          Spatial search radius for spatial reduction [m], e.g. 500 (default empty)
%      diffEpochsMax        Maximum absolute epoch difference for temporal reduction [year], e.g. 0.1.
%                           If empty, this criterion is not used. (default empty)
%      numEpochsMax         Maximum number of epochs for temporal reduction, e.g. 5. If empty, this
%                           criterion is not used. (default empty)
%
%      splitDefoRegimes     Apply split in deformation regimes (deep/shallow), true (default) or false.
%
%      InSARModels          InSAR models for selection of best model, e.g.
%                           {...
%                           {'constant','linear'};...
%                           {'constant','linear','periodic'};...
%                           {'constant','linear','breakpoint'};...
%                           {'constant','linear','heavyside'};...
%                           {'constant','linear','periodic','breakpoint'};...
%                           {'constant','linear','periodic','heavyside'};...
%                           };
%
%                           or
%
%                           {...
%                           {'constant','linear'};...
%                           {'constant','linear','periodic'};...
%                           };
%
%      InSARModelsMinEpochs Minimum number of epochs for partial model (breakpoint, heavyside), e.g. 5.
%
%      stochModelParameterFile  Name of the file with stochastic model parameters
%                               (Default 'insarStochasticModelParameters.txt') 
%
%      globalAttrib         Struct with global attribute updates (default from input file)
%
%   Examples:
%      stmreduceinsar({'input_filename.mat'},'output_filename.mat')      
%      stmreduceinsar({'input_filename.mat'},'output_filename.mat',options,'update')      
%
%   The following 6 steps can be distinguished in the module:
%   1) Reading of input space-time dataset, and evaluation points and epochs.
%   2) Reduction to evaluation points. Hereby, the parameter opt.insarRadius
%      is used to find InSAR points in the surroundings of the evaluation points.
%      For now, the representative time series is obtained by an unweighted
%      average. This may be improved in the future (e.g., weighting, outlier
%      removal). If a separation in deep and total deformation regimes is applied
%      while creating the InSAR space-time dataset, a representative time series is
%      calculated for both deep and shallow (deep-total) processes.
%   3) Estimation of the best fitting funtional model on the time series per 
%      evaluation point by multi-hypothesis testing, using the mht toolbox. This
%      can be seen as a decomposition. Subsequently, the residuals are calculated.
%   4) Reduction to the evaluation epochs based on an unweighted mean of the 
%      time series residuals. The residuals to be used in the averaging are
%      selected based on 1) a binning of all original epochs around the evaluation
%      epochs to avoid overlap (hence, obtaining uncorrelated estimates), and 
%      optionally 2) a maximum number of epochs to use (set by the opt.numEpochsMax
%      parameter) and/or a maximum absolute difference between the observation 
%      and the evaluation epochs. Hereby, unnecessary averaging over long time 
%      spans is avoided. Note that both, one of the, or none of the parameters can
%      be used.
%   5) Restore of the estimated deformation model to obtain the final time series
%      at the evalution points.
%   6) Writing of the reduced space-time dataset. This includes the stochastic
%      model for reduced InSAR datasets, together with the parameters needed by
%      the model.
%
%   (c) Freek van Leijen, Delft University of Technology, 2020. 

% Created:  12 Nov 2020 by Freek van Leijen
% Modified: 08 Mar 2021 by Freek van Leijen
%           - Changed dataset.obsType 'losCompaction' to 'losShallow'
%           19 Oct 2021 by Freek van Leijen
%           - Added the opt.diffEpochsMax parameter, the max abs epoch difference 
%             for temporal reduction [optional]
%           - Made the opt.numEpochsMax optional, just like the opt.diffEpochsMax
%             parameter.
%           21 Oct 2021 by Freek van Leijen
%           - Renamed from rd_insar.m to stmreduceinsar.m
%           - Aligned setup with the stmreducegnss.m function
%           - Added the decomp_insar_mht and eval_insar_model function as
%             subfunctions.
%           5 Nov 2021 by Hans van der Marel
%           - read and use new project file
%           - use getpointmask and getepochmask for ROI and POI
%           - global attributes
%           7 Jan 2022 by Hans van der Marel
%           - stochModelParameterFile is input option
%           4 Feb 2023 by Freek van Leijen
%           - fixed bug regarding position of loop for shallow signal
%           - changed check on defoRegime (because now always initialized with zeros)
%          27 oct 2023 by Hans van der Marel
%           - fixed minor bug in output datasetId
%           8 jan 2024 by FvL?
%           - added option for splitDefoRegimes
%           7 march 2024 by Hans van der Marel
%           - respect doplots option
%

%% Check the input arguments and options

progname = 'stmreduceinsar';

if nargin < 2
    error('This function expects at least two input arguments.')
end

% Default options (none)

opt.verbose = 0;                 % Default verbosity level, higher is more output, 0 is almost nothing
opt.doplots = 0;                 % Default plot level, 0 is no plots, higher is more detailed plotting

opt.projectId = '';              % Project Id (if empty, projectId is determined from the outputfilename directory)
opt.projectFile = '';            % Project File (if empty, it is the same as the projectId with extension .mat),
                                 % if projectFile is used, projectID is set from the projectFile

opt.inputDir = '';               % Directory with the inputfiles
opt.datasetId = '';              % DatasetId for the output data set, if empty, "_reduced" is added to the input datasetId.
opt.outputTarget = 'create';     % Default mode for output STM {'create','overwrite','update'}, can be changed with flag

opt.ROI = [];                    % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default none)
opt.POI = [-Inf +Inf];           % Period of interest [ dYearStart dYearEnd ] (Default none)

opt.insarRadius = [];            % [m], spatial search radius for spatial reduction
opt.diffEpochsMax = [];          % [year], max abs epoch difference for temporal reduction [optional]
opt.numEpochsMax = [];           % maximum number of epochs for temporal reduction [optional]

opt.splitDefoRegimes = true;     % Apply split in deformation regimes (deep/shallow), default true.

opt.InSARModels = {};            % InSAR models for selection of best model
opt.InSARModelsMinEpochs = [];   % Minimum number of epochs for partial model (breakpoint, heavyside)

opt.stochModelParameterFile='insarStochasticModelParameters.txt';  % Name of the file with stochastic model parameters

opt.globalAttrib = struct([]);   % Struct with global attribute updates (default from input file)


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

techniqueAttribTotal=projectFile.techniqueAttrib;
campaignDatasetIds=techniqueAttribTotal.campaignDatasetIds;  % {'gron_levelling_flaggedOutliers_v2.3.02_20191030_commonIDs'}
markerDatasetIds=techniqueAttribTotal.markerDatasetIds;      % {'gron_levelling_flaggedOutliers_v2.3.02_20191030_commonIDs'  '06gps_nam_reduced_202107_reduced'}
pntCrdRef=techniqueAttribTotal.pntCrdRef;

numEvaluationEpochs=projectFile.numEpochs;
numEvaluationPoints=projectFile.numPoints;

evaluationPoints.pntCrdTotal=projectFile.pntCrd;       % emulate evaluationPoints stucture;
pntAttribTotal=projectFile.pntAttrib;  
evaluationPoints.pntNeuTotal=pntAttribTotal.pntNeu;
evaluationPoints.pntIdTotal=pntAttribTotal.pntId;
             
evaluationEpochs.epochDyearTotal=projectFile.epochDyear;  %  replaces evaluationEpochs.epochDyearTotal
epochAttribTotal=projectFile.epochAttrib;
evaluationEpochs.epochIdTotal=epochAttribTotal.epochId;

%% Read InSAR STM

fprintf('Reading input InSAR STM %s ... \n',inputfilename);

st = stmread(inputfilename);


%% Optionally select subset of points and epochs

% Select regions of interest  (ROI)

fprintf('Number of evaluation points %d\n',numEvaluationPoints);
if isempty(opt.ROI)
   pntMask = true(numEvaluationPoints,1);
else
   pntMask = getpointmask(evaluationPoints.pntCrdTotal,opt.ROI);
   fprintf('Number of evaluation points inside ROI %d (out of %d)\n',sum(pntMask),numEvaluationPoints);

   % Update evaluation points
   evaluationPoints.pntCrdTotal = evaluationPoints.pntCrdTotal(pntMask,:);
   evaluationPoints.pntNeuTotal = evaluationPoints.pntNeuTotal(pntMask,:);
   evaluationPoints.pntIdTotal = evaluationPoints.pntIdTotal(pntMask);
   numEvaluationPoints = sum(pntMask);
end

% Select periods of interest (POI)

fprintf('Number of evaluation epochs %d\n',numEvaluationEpochs);
if isempty(opt.POI)
   epochMask = true(numEvaluationEpochs,1);
else
   epochMask = getepochmask(evaluationEpochs.epochDyearTotal,opt.POI);
   fprintf('Number of evaluation epochs inside POI %d (out of %d)\n',sum(epochMask),numEvaluationEpochs);

   % Update evaluation epochs
   evaluationEpochs.epochDyearTotal = evaluationEpochs.epochDyearTotal(epochMask);
   evaluationEpochs.epochIdTotal = evaluationEpochs.epochIdTotal(epochMask);
   numEvaluationEpochs = sum(epochMask);
end


%% Reduction to evaluation points
% RS = Reduced Spatial, RST = Reduced Spatial Temporal

fprintf('\nReduction to evaluation points ... \n');

% Compute north, east coordinates using the same reference point as for evaluation points
pntNeu = plh2neu(st.pntCrd,pntCrdRef,'deg');

if opt.splitDefoRegimes
  numObsTypes = 2;
else
  numObsTypes = 1;
end
obsDataRS = NaN(numEvaluationPoints,st.numEpochs,numObsTypes,'single');
cellCountPoints = NaN(numEvaluationPoints,numObsTypes,'single');
cellAvgDist = NaN(numEvaluationPoints,numObsTypes,'single');
sensitivityMatrixRS = NaN(numEvaluationPoints,3,'single'); % generic, based on deep only
incAngle = NaN(numEvaluationPoints,1,'single'); % generic, based on deep only
azAngle = NaN(numEvaluationPoints,1,'single'); % generic, based on deep only

for k=1:numEvaluationPoints

  % Find points within bounding box radius
  idx1 = find(pntNeu(:,1)>evaluationPoints.pntNeuTotal(k,1)-opt.insarRadius & ...
              pntNeu(:,1)<evaluationPoints.pntNeuTotal(k,1)+opt.insarRadius & ...
              pntNeu(:,2)>evaluationPoints.pntNeuTotal(k,2)-opt.insarRadius & ...
              pntNeu(:,2)<evaluationPoints.pntNeuTotal(k,2)+opt.insarRadius);

  if ~isempty(idx1)
    % Find points within radius
    pntDist = hypot(pntNeu(idx1,1)-evaluationPoints.pntNeuTotal(k,1),pntNeu(idx1,2)-evaluationPoints.pntNeuTotal(k,2));
    idx2 = find(pntDist<opt.insarRadius);

    if ~isempty(idx2)

      % Check for IGRS
      if isfield(st.pntAttrib,'pntClass') 
         idx3 = find(st.pntAttrib.pntClass(idx1(idx2))==5);
      else
         idx3 = [];
      end

      %if numel(idx3)==1 & pntDist(idx2(idx3))<10 % CR, AR, or IGRS found
      if numel(idx3)==1 % CR, AR, or IGRS found

        obsDataRS(k,:,1) = st.obsData(idx1(idx2(idx3)),1:st.numEpochs); % Deep only
        cellCountPoints(k,1) = 1; % Deep only, so numAveragePoints(k,2) remains NaN
        cellAvgDist(k,1) = Inf; % Set to Inf, for proper stochastic model evaluation
        sensitivityMatrixRS(k,:) = st.sensitivityMatrix(idx1(idx2(idx3)),:);
        incAngle(k) = st.pntAttrib.incAngle(idx1(idx2(idx3)));
        azAngle(k) = st.pntAttrib.azAngle(idx1(idx2(idx3)));
        %TODO still get total (and thereby shallow) deformation from surroundings?

      elseif numel(idx3)>1

        error('More than one CR/AR/IGRS is found in vincinity.');
        %TODO There could be close-by CR, e.g. in Wassenaar ...

      else

        if ~isempty(find(st.pntAttrib.defoRegime>0)) & opt.splitDefoRegimes % deformation regime classification
          idx4a = find(st.pntAttrib.defoRegime(idx1(idx2))==1); % Deep
          idx4b = find(st.pntAttrib.defoRegime(idx1(idx2))==2); % Total
        else % all points are unclassified
          idx4a = (1:numel(idx2))'; 
          idx4b = [];
        end  

        if ~isempty(idx4a) % Points with Deep signal exist

          %obsData1 = NaN(numel(idx4a),st.numEpochs);
          %for v = 1:numel(idx4a)
            %v
            %tic
            %%obsData1(v,:) = stmread(st,'OBSDATA',idx1(idx2(idx4a(v))),1:st.numEpochs,1);
            %obsData1(v,:) = st.obsData(idx1(idx2(idx4a(v))),1:st.numEpochs);
            %toc
          %end
          %obsDataRS(k,:,1) = mean(obsData1,1);
          %sensitivityMatrixRS(k,:) = mean(sensitivity,1);
            
          obsDataRS(k,:,1) = mean(st.obsData(idx1(idx2(idx4a)),:),1);
          cellCountPoints(k,1) = numel(idx4a);

          if cellCountPoints(k,1)==1
            cellAvgDist(k,1) = Inf; % Set to Inf, for proper stochastic model evaluation
          else
            cellAvgDist(k,1) = mean(pdist(pntNeu(idx1(idx2(idx4a)),1:2)));
          end

          sensitivityMatrixRS(k,:) = mean(st.sensitivityMatrix(idx1(idx2(idx4a)),:),1);

          incAngle(k) = nanmean(st.pntAttrib.incAngle(idx1(idx2(idx4a))));
          azAngle(k) = nanmean(st.pntAttrib.azAngle(idx1(idx2(idx4a))));

          if ~isempty(idx4b) % no need if idx4a is empty, no deep points for integration
            % Shallow = total - deep;
            obsDataRS(k,:,2) = mean(st.obsData(idx1(idx2(idx4b)),:),1) - obsDataRS(k,:,1);
            cellCountPoints(k,2) = numel(idx4b); 

            if cellCountPoints(k,2)==1
              cellAvgDist(k,2) = Inf; % Set to Inf, for proper stochastic model evaluation
            else
              cellAvgDist(k,2) = mean(pdist(pntNeu(idx1(idx2(idx4b)),1:2)));
            end

          end

        end

      end
    end
  end
end

% Remove points without value
nanIdx = find(~isnan(cellCountPoints(:,1)));
obsDataRS = obsDataRS(nanIdx,:,:);
sensitivityMatrixRS = sensitivityMatrixRS(nanIdx,:);
cellCountPoints = cellCountPoints(nanIdx,:);
cellAvgDist = cellAvgDist(nanIdx,:);
incAngle = incAngle(nanIdx,:);
azAngle = azAngle(nanIdx,:);
pntNeuRS = evaluationPoints.pntNeuTotal(nanIdx,:);
pntCrdRS = evaluationPoints.pntCrdTotal(nanIdx,:);
pntIdRS = evaluationPoints.pntIdTotal(nanIdx);
numPointsRS = size(obsDataRS,1);

if opt.doplots > 1
    h = NaN(3,1);
    figure;hold on;
    h(1) = plot(pntNeu(:,2),pntNeu(:,1),'r.');
    h(2) = plot(evaluationPoints.pntNeuTotal(:,2),evaluationPoints.pntNeuTotal(:,1),'b.');
    h(3) = plot(pntNeuRS(:,2),pntNeuRS(:,1),'go');
    axis equal;
    legend(h,'InSAR points','Original+created evaluation points','Available evaluation points');
    xlabel('East [m]')
    ylabel('North [m]')
    title('InSAR spatial reduction');
end

fprintf('\nEstimation optimal displacement model (decomposition) ... \n');
%% Estimation optimal displacement model

% Multi-hypothesis testing to find optimal model
[displModelLib,displModel] = decomp_insar_mht(obsDataRS,st.epochDyear,st.stochModel,st.stochData,opt);

% Get residuals wrt displacement model
[modelDataRS,modelDataRST] = eval_insar_model(obsDataRS,st.epochDyear,evaluationEpochs.epochDyearTotal,st.stochModel,st.stochData,displModelLib,displModel,opt);
resDataRS = obsDataRS-modelDataRS;


%% Evaluation at epochs
%dataEpochs = evl_insar_time(dataSep,epochs);

fprintf('\nReduction to evaluation epochs ... \n');

% Determine number of observations within bins around evaluation epochs
edges = [-inf 0.5*(evaluationEpochs.epochDyearTotal(1:end-1)+evaluationEpochs.epochDyearTotal(2:end)) inf];
[~,~,binIdx] = histcounts(st.epochDyear,edges);

resDataRST = NaN(numPointsRS,numEvaluationEpochs,numObsTypes,'single');
cellCountEpochs = NaN(numPointsRS,numEvaluationEpochs,numObsTypes,'single');
%cellAvgPeriod = NaN(numEvaluationEpochs,2,'single'); % Equal for all points
cellAvgPeriod = NaN(numPointsRS,numEvaluationEpochs,numObsTypes,'single');
for k = 1:numEvaluationEpochs

  idx1 = find(binIdx==k);
  if ~isempty(idx1)

    % Select epochs
    epochDiff = abs(st.epochDyear(idx1)-evaluationEpochs.epochDyearTotal(k));
    [~,idx2] = sort(epochDiff);

    % Optionally apply a maximum epoch time span
    if ~isempty(opt.diffEpochsMax)
      idx3 = find(epochDiff(idx2)<opt.diffEpochsMax);
      idx2 = idx2(idx3); % update idx2 if diffEpochsMax
    end

    % Optionally apply a maximum number of epochs
    if ~isempty(opt.numEpochsMax)
      numSamp = min(numel(idx2),opt.numEpochsMax); % Number of samples
    else
      numSamp = numel(idx2);
    end

    resData = resDataRS(:,idx1(idx2(1:numSamp)),:);
    resDataRST(:,k,:) = nanmean(resData,2); % Consider appeared/disappeared scatterers, e.g., IGRS
    nanMaskOrig = ~isnan(resData);
    cellCountEpochs(:,k,:) = sum(nanMaskOrig,2); % Again, appeared/disappeared scatterers

    %%for all points
    %[ta,tb] = meshgrid(st.epochDyear(idx1(idx2(1:numSamp))),st.epochDyear(idx1(idx2(1:numSamp))));
    %cellPeriod = abs(ta-tb)+triu(NaN(numSamp)); % NaN for diagonal and upper triangular part
    %cellAvgPeriod(k,1) = nanmean(cellPeriod(:));

    %unique per point, bin, deep and shallow ....
    [ta,tb] = meshgrid(st.epochDyear(idx1(idx2(1:numSamp))),st.epochDyear(idx1(idx2(1:numSamp))));
    cellPeriodOrig = abs(ta-tb)+triu(NaN(numSamp)); % NaN for diagonal and upper triangular part
    cellPeriodOrig = repmat(cellPeriodOrig,1,1,numPointsRS); % Setup array for each point

    for l=1:numObsTypes
      nanMask = nanMaskOrig(:,:,l);
      cellPeriod = cellPeriodOrig;
      nanMaska = kron(nanMask,ones(1,numSamp))'; % setup mask in column direction
      nanMaska = reshape(nanMaska(:),numSamp,numSamp,numPointsRS);
      nanMaskb = kron(ones(1,numSamp),nanMask)'; % setup mask in row direction
      nanMaskb = reshape(nanMaskb(:),numSamp,numSamp,numPointsRS);
      cellPeriod(~nanMaska) = NaN;
      cellPeriod(~nanMaskb) = NaN;
      cellPeriod = reshape(permute(cellPeriod,[3 1 2]),numPointsRS,numSamp*numSamp);
      cellAvgPeriod(:,k,l) = nanmean(cellPeriod,2); %TODO check for NaN if number samples is 1, should be Inf 
    end
  end
  
end

% Add model back
obsDataRST = modelDataRST+resDataRST;

% Remove epochs without value
nanIdx = find(sum(~isnan(obsDataRST(:,:,1)),1)~=0); % Based on Deep
obsDataRST = obsDataRST(:,nanIdx,:);
cellCountEpochs = cellCountEpochs(:,nanIdx,:);
cellAvgPeriod = cellAvgPeriod(:,nanIdx,:);
epochDyearFinal = evaluationEpochs.epochDyearTotal(nanIdx);
epochIdFinal = evaluationEpochs.epochIdTotal(nanIdx);
numEpochsFinal = size(obsDataRST,2);

% Remove points without epochs with value (to be sure...)
nanIdx = find(sum(~isnan(obsDataRST(:,:,1)),2)~=0); % Based on Deep
obsDataFinal = obsDataRST(nanIdx,:,:);
sensitivityMatrixFinal = sensitivityMatrixRS(nanIdx,:);
cellCountPoints = cellCountPoints(nanIdx,:);
cellAvgDist = cellAvgDist(nanIdx,:);
cellCountEpochs = cellCountEpochs(nanIdx,:,:);
cellAvgPeriod = cellAvgPeriod(nanIdx,:,:);
incAngle = incAngle(nanIdx,:);
azAngle = azAngle(nanIdx,:);
pntNeuFinal = pntNeuRS(nanIdx,:);
pntCrdFinal = pntCrdRS(nanIdx,:);
pntIdFinal = pntIdRS(nanIdx,:);
numPointsFinal = size(obsDataFinal,1);

%for k = 1:numPointsFinal
%    h = NaN(4,1);
%    figure;hold on
%    h(1) = plot(st.epochDyear,obsDataRS(k,:,1),'b-o');
%    h(2) = plot(epochDyearFinal,obsDataFinal(k,:,1),'r-o');
%    h(3) = plot(st.epochDyear,obsDataRS(k,:,2),'b-.');
%    h(4) = plot(epochDyearFinal,obsDataFinal(k,:,2),'r-.');
%    grid on
%    ylabel('Displacement [mm]');
%    legend(h,'Original deep','Reduced deep','Original shallow','Reduced shallow'); 
%    display('pause')
%    pause
%end


%% Get InSAR stochastic model parameters

if ~exist(opt.stochModelParameterFile,'file')
    error('File with stochastic model parameters does not exist.');
end 

fid = fopen(opt.stochModelParameterFile,'r');
modelParam = textscan(fid,['%s%s%s%f32%f32%f32%f32%f32'], ...
     'Delimiter',',','Headerlines',1);
fclose(fid);

idx1 = strmatch(st.techniqueAttrib.system,modelParam{1});
idx2 = strmatch(st.techniqueAttrib.systemMode,modelParam{2}(idx1));
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



%% Create output space time matrix

fprintf('Creating reduced InSAR STM %s ... \n',outputfilename);

if ~isempty(opt.datasetId)
   datasetIdOut = opt.datasetId;
else
   datasetIdOut = [ st.datasetId '_reduced' ];
end
   
stout = stm(projectId,datasetIdOut,'insar');  

% general attributes
techniqueAttrib = stout.techniqueAttrib;
techniqueAttrib.status = 'Reduced';
stout.techniqueAttrib = techniqueAttrib;

datasetAttrib = stout.datasetAttrib;
datasetAttrib.softwareOptions = opt;
stout.datasetAttrib = datasetAttrib;

% epoch attributes
stout.numEpochs = numEpochsFinal;
stout.epochDyear = epochDyearFinal;

% Store pntAttrib first in temporary array, before moving it to dataset (objects only support one level of indexing)
epochAttrib = [];
epochAttrib.epochId = epochIdFinal;
%epochAttrib.cellCount = cellCountEpochs; %TODO change to stochData.cellCountEpochs ???
%epochAttrib.cellSize = cellAvgPeriod; %TODO change to stochData.cellAvgPeriod ???
epochAttrib.cellCount = mode(cellCountEpochs(:,:,1),1); %now generic for all points, deep only, TODO change to stochData.cellCountEpochs ???
epochAttrib.cellSize = mode(cellAvgPeriod(:,:,1),1); %now generic for all points, deep only, TODO change to stochData.cellAvgPeriod ???
stout.epochAttrib = epochAttrib;

% point attributes
stout.numPoints = numPointsFinal;
stout.pntName = pntIdFinal; % Same as pntId, unique identifier for each measurement point.
stout.pntCrd = pntCrdFinal; % Lat [Deg], Lon [Deg], h [m]

% Store pntAttrib first in temporary array, before moving it to dataset (objects only support one level of indexing)
pntAttrib = [];
pntAttrib.pntId = pntIdFinal; % Unique identifier for each measurement point.
pntAttrib.incAngle = incAngle;
pntAttrib.azAngle = azAngle;
pntAttrib.cellCount = cellCountPoints(:,1); %TODO shallow???, change to stochData.cellCountPoints ???
pntAttrib.cellSize = cellAvgDist(:,1); %TODO shallow???, change to stochData.cellAvgDist ???
pntAttrib.displModel = displModel;
stout.pntAttrib = pntAttrib;

stout.displModelLib = displModelLib;

% observation attributes
stout.parTypes = {'North' 'East' 'Up'};

% check whether any
if opt.splitDefoRegimes
  if any(any(obsDataFinal(:,:,2))) % If only NaNs
    numObsTypes = 1;
    obsDataFinal = obsDataFinal(:,:,1);
  end
end

if numObsTypes == 1
  stout.obsTypes = {'los'};
elseif numObsTypes == 2  
  stout.obsTypes = {'los','losShallow'};
  sensitivityMatrixFinal = cat(3,sensitivityMatrixFinal,sensitivityMatrixFinal);
end  
stout.obsData = obsDataFinal;
stout.sensitivityMatrix = sensitivityMatrixFinal;

stout.auxTypes = [];
stout.auxData = [];


%% Set stochasic model

%stout.stochModel={'tudinsar4rd(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Nt=data1,Dt=data2,Rs=1.11,Ns=data3,Ds=data4)'};
%stout.stochData= ..  % [M-by-N-4] matrix

%stout.stochModel = {[[stochModelParam.model 'rd'] ...
%                    '(s20=' num2str(stochModelParam.s20) ...
%                    ',s2t=' num2str(stochModelParam.s2t) ...
%                    ',s2s=' num2str(stochModelParam.s2s) ...
%                    ',Rt=' num2str(stochModelParam.Rt) ...
%                    ',Nt=data1' ... %cellCount
%                    ',Dt=data2' ... %cellSize
%                    ',Rs=' num2str(stochModelParam.Rs) ...
%                    ',Ns=data3' ... %cellCount
%                    ',Ds=data4)']}; %cellSize
%
%stochData = NaN(stout.numpoints,stout.numepochs,4);
%stochData(:,:,1) = cellCountEpochs;
%stochData(:,:,2) = cellAvgPeriod;
%stochData(:,:,3) = cellCountPoints;
%stochData(:,:,4) = cellAvgDist;
%stout.Stochdata = stochData;

stochModel = {[[stochModelParam.model 'rd'] ...
                '(s20=' num2str(stochModelParam.s20) ...
                ',s2t=' num2str(stochModelParam.s2t) ...
                ',s2s=' num2str(stochModelParam.s2s) ...
                ',Rt=' num2str(stochModelParam.Rt) ...
                ',Nt=cellCount' ... %cellCountEpochs?
                ',Dt=cellSize' ... %cellAvgPeriod?
                ',Rs=' num2str(stochModelParam.Rs) ...
                ',Ns=cellCount' ... %cellCountPoints?
                ',Ds=cellSize)']}; %cellAvgDist?

if numObsTypes == 1
  stout.stochModel = {stochModel};
elseif numObsTypes == 2
  stout.stochModel = {stochModel,stochModel}; %TODO create different models for deep and shallow
end

stout.stochData = [];


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


%% Write the reduced dataset

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


function [H0,pntModel] = decomp_insar_mht(obsData,epochDyear,stochModel,stochData,opt)
%DECOMP_INSAR_MHT Decompose InSAR time series based on multi-hypothesis testing.
%  [H0,PNTMODEL] = DECOMP_INSAR_MHT(OBSDATA,EPOCHDYEAR,STOCHMODEL,STOCHDATA,
%  OPT) decomposes InSAR time series in OBSDATA in an optimal set of 
%  displacement models based on multi-hypothesis testing. The time series
%  epochs are given by EPOCHDYEAR. The stochastic model defined by STOCHMODEL
%  and STOCHDATA is used. The following parameters are parsed via the OPT
%  variable:
%  - opt.InSARModels              Array with displacement models to evaluate
%  - opt.InSARModelsMinEpochs     Minimal number of epochs at the beginning
%                                 and end of the displacement time series,
%                                 applicable for breakpoint and heavyside
%                                 models only.
%
%  The outputs are an array H0 with the specifications of the models used, and
%  PNTMODEL with an index in the H0 model array for each point.
%
%  This function is called by the RD_INSAR.m function, and uses the MHT toolbox.
%
%
%  (c) Freek van Leijen, Delft University of Technology, 2020. 

% Created:  18 September 2020 by Freek van Leijen
% Modified: 21 Februari 2021 by Freek van Leijen
%           - options parsing via opt variable
%

% TODO 
% - add functionality for temperature
% - check the stochastic model used
%

%% Set parameters
[numPoints,numEpochs,numLayers] = size(obsData);

alpha0 = 1/(2*numEpochs); % see Chang and Hanssen, hardcoded for now
gamma0 = 0.5; % Best choice for MHT

bufferSize = 100; % Number of points per buffer, trade of between number of loops
                  % and array sizes. Tests show that (for current hardware) 100 is 
                  % better than 10 or 1000.


%% Setup model
Btemp = epochDyear-epochDyear(1);
BT = []; %No temperature model yet

% Setup covariance matrix
Qyy = stmstochmodel(stochModel,stochData,[0 0;1 1],epochDyear);
Qyy = Qyy(1:numEpochs,1:numEpochs);
invQyy = inv(Qyy);

% Setup models
[H0,numParMax] = mht_setup_models(opt.InSARModels,Btemp,BT,opt.InSARModelsMinEpochs);
numH0 = size(H0,2);

% Get critical values
kOMT = NaN(numH0,1);
alphaOMT = NaN(numH0,1);
  
for v = 1:numH0
    % lam0 : non-centre parameter of chi-square distribution
    % k1: critical value for 1d test
    % kb: critical value of b-dimensional test
    % level of significance for b-dimensional test
   [lam0,k1,kOMT(v),alphaOMT(v)] = pretest(numEpochs-H0(v).q,alpha0,gamma0);
end


%% Test H0
OMT = NaN(numPoints,numH0,numLayers,'single');
TOMT = NaN(numPoints,numH0,numLayers,'single');
varfac = NaN(numPoints,numH0,numLayers,'single');

for v = 1:numLayers

  validMask = ~isnan(obsData(:,:,v)); % Check for temporary scatterers or missing data (layer >1)
  validIdx = find(sum(validMask,2)==numEpochs);
  nanIdx = setdiff(1:numPoints,validIdx);
  numValid = numel(validIdx);
  numNan = numel(nanIdx);

  % Analyze full time series
  if numValid > 0

    for w = 1:numH0 %Loop over different hypothesis
      A0 = H0(w).Cj;
      Qxx = inv(A0'*invQyy*A0);
      xhat = Qxx*A0'*invQyy*obsData(validIdx,:,v)';
      ehat = obsData(validIdx,:,v)' - A0*xhat;
      Qee = Qyy-A0*Qxx*A0';

      Nbuffer = floor(numValid/bufferSize);
      remBufferSize = rem(numValid,bufferSize);
      if Nbuffer==0
        z = 0;
      else
        for z = 1:Nbuffer
          OMT(validIdx(1+(z-1)*bufferSize:z*bufferSize),w,v) = ... 
          diag(ehat(:,1+(z-1)*bufferSize:z*bufferSize)'*invQyy* ...
          ehat(:,1+(z-1)*bufferSize:z*bufferSize));
        end
      end
      if remBufferSize~=0
        OMT(validIdx(z*bufferSize+1:z*bufferSize+remBufferSize),w,v) = ...
         diag(ehat(:,z*bufferSize+1:z*bufferSize+remBufferSize)'*invQyy* ...
         ehat(:,z*bufferSize+1:z*bufferSize+remBufferSize));
      end

      TOMT(validIdx,w,v) = OMT(validIdx,w,v)/kOMT(w);

      varfac(validIdx,w,v) = OMT(validIdx,w,v)/(numEpochs-H0(w).q);

    end
  end

  % Analyze temporary scatterers
  for k = 1:numNan  % Skipped if numNan==0
    if length(find(validMask(nanIdx(k),:))) >= 2*numParMax % should be something to estimate
      Qyy1 = Qyy(validMask(nanIdx(k),:),validMask(nanIdx(k),:));
      invQyy1 = inv(Qyy1);
      numEpochs1 = size(Qyy1,1);

      for w = 1:2 %numH0 %Loop over different hypothesis
        %TODO breakpoint and heavyside only work within time span
        A1 = H0(w).Cj(validMask(nanIdx(k),:)',:);
        Qxx1 = inv(A1'*invQyy1*A1);
        xhat1 = Qxx1*A1'*invQyy1*obsData(nanIdx(k),validMask(nanIdx(k),:),v)';
        ehat1 = obsData(nanIdx(k),validMask(nanIdx(k),:),v)' - A1*xhat1;
        Qee1 = Qyy1-A1*Qxx1*A1';

        OMT(nanIdx(k),w,v) = ehat1'*invQyy1*ehat1;

        [~,~,kOMT1,~] = pretest(numEpochs1-H0(w).q,alpha0,gamma0); %TODO base this on q only
        TOMT(nanIdx(k),w,v) = OMT(nanIdx(k),w,v)/kOMT1;
        varfac(nanIdx(k),w,v) = OMT(nanIdx(k),w,v)/(numEpochs1-H0(w).q);
      end
    end
  end
end

[minOMT,pntModel] = min(TOMT,[],2);
minOMT = reshape(minOMT,numPoints,numLayers);
pntModel = reshape(pntModel,numPoints,numLayers);
pntModel(isnan(minOMT)) = NaN; % Avoid index 1 for NaN rows

end %end decomp_insar_mht


function [modelData,modelDataRST] = eval_insar_model(obsData,epochDyear,epochDyearRST,stochModel,stochData,displModel,pntModel,opt)
%EVAL_INSAR_MODEL Evaluates the optimal displacement models at the evaluation epochs.
%  [MODELDATA,MODELDATARST] = EVAL_INSAR_MODEL(OBSDATA,EPOCHDYEAR,EPOCHDYEARRST,
%  STOCHMODEL,STOCHDATA,DISPLMODEL,PNTMODEL,OPT) evaluates the displacement
%  model for each point, based on the optimal model selected in the
%  DECOMP_INSAR_MHT.m function. These optimal models are indicted with an index
%  PNTMODEL in the DISPLMODEL array. The evaluation is done both for the original 
%  epochs EPOCHDYEAR and the evalutation EPOCHDYEARRST. The stochastic model defined
%  by STOCHMODEL and STOCHDATA is used. The following parameters are parsed via the
%  OPT variable:
%  - opt.InSARModels              Array with displacement models to evaluate
%  - opt.InSARModelsMinEpochs     Minimal number of epochs at the beginning
%                                 and end of the displacement time series,
%                                 applicable for breakpoint and heavyside
%                                 models only.
%
%  The outputs the evalated models per point, both for the original epochs MODELDATA
%  and the evaluation epochs MODELDATARST.
%
%  This function is called by the STMREDUCEINSAR.m function.
%
%
%  (c) Freek van Leijen, Delft University of Technology, 2020. 

% Created:  18 September 2020 by Freek van Leijen
% Modified: 21 Februari 2021 by Freek van Leijen
%           - options parsing via opt variable
%

% TODO 
% - check the stochastic model used
%

%% Set parameters
[numPoints,numEpochs,numLayers] = size(obsData);
numEpochsRST = numel(epochDyearRST);


% Setup covariance matrix
Qyy = stmstochmodel(stochModel,stochData,[0 0;1 1],epochDyear);
Qyy = Qyy(1:numEpochs,1:numEpochs);
invQyy = inv(Qyy);


% Setup models for reduced dataset
BtempRST = epochDyearRST-epochDyear(1); %Btemp with same reference as original dataset
BTRST = []; %No temperature model yet

[displModelRST,~] = mht_setup_models(opt.InSARModels,BtempRST,BTRST,opt.InSARModelsMinEpochs);


%% Estimate model data
modelData = NaN(numPoints,numEpochs,numLayers);
modelDataRST = NaN(numPoints,numEpochsRST,numLayers);

for v = 1:numLayers

  validMask = ~isnan(obsData(:,:,v)); % Check for temporary scatterers or missing data (layer >1)
  validIdx = find(sum(validMask,2)==numEpochs);
  nanIdx = setdiff(1:numPoints,validIdx);
  numValid = numel(validIdx);
  numNan = numel(nanIdx);

  % Analyze full time series
  if numValid > 0
    uniqueModels = unique(pntModel(validIdx,v));
    numUniqueModels = numel(uniqueModels);

    for w = 1:numUniqueModels

      % Model original epochs
      pntIdx = find(pntModel(validIdx,v)==uniqueModels(w));
      A0 = displModel(uniqueModels(w)).Cj;
      Qxx = inv(A0'*invQyy*A0);
      xhat = Qxx*A0'*invQyy*obsData(validIdx(pntIdx),:,v)';
      modelData(validIdx(pntIdx),:,v) = (A0*xhat)';

      % Model reduced epochs
      B0 = displModelRST(uniqueModels(w)).Cj;
      modelDataRST(validIdx(pntIdx),:,v) = (B0*xhat)';

    end
  end

  % Analyze temporary scatterers
  for k = 1:numNan  % Skipped if numNan==0
    if ~isnan(pntModel(nanIdx(k),v))
      Qyy1 = Qyy(validMask(nanIdx(k),:),validMask(nanIdx(k),:));
      invQyy1 = inv(Qyy1);
      numEpochs1 = size(Qyy1,1);

      % Model original epochs
      A1 = displModel(pntModel(nanIdx(k),v)).Cj(validMask(nanIdx(k),:)',:);
      Qxx1 = inv(A1'*invQyy1*A1);
      xhat1 = Qxx1*A1'*invQyy1*obsData(nanIdx(k),validMask(nanIdx(k),:),v)';
      modelData(nanIdx(k),validMask(nanIdx(k),:),v) = (A1*xhat1)';

      % Model reduced epochs
      B1 = displModelRST(pntModel(nanIdx(k),v)).Cj;
      modelDataRST(nanIdx(k),:,v) = (B1*xhat1)';

    end
  end

end

end % end eval_insar_model

