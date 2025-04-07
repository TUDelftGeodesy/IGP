function [mask,commonPoints] = stmclustering(stm,varargin)
%STMCLUSTERING Clustering of nearby observation points.
%  MASK = STMCLUSTERING(STM) performs a clustering of points in the
%  Space-Time Matrix STM and creates a MASK with the 'best' point
%  per cluster. The clustering is based on a maximum distance
%  (spatialTolerance) betweens points, or on benchmark ID. See options 
%  below.The selection of the 'best' point of a cluster is based on 
%  subsequently
%
%  1) most epochs
%  2) longest time span
%  3) random choice (first in list)
%
%  The clustering is based on a maximum distance (spatialTolerance) 
%  betweens points, or on benchmark ID. See options below.
%
%  [MASK,COMMONPOINTS] = STMCLUSTERING(STM) also outputs a structure
%  with information on the points clustered.
%
%  [...]=STMCLUSTERING(STM,'option',value,...) allows to specify options
%  for the clustering
%
%    'clusteringMethod'   clustering method ('distance' (default) or 'benchmark')
%    'spatialTolerance'   spatial tolerance [m] for clustering based on 'distance'
%                         (default 1 m), can also be 0 m.
%    'doplots'            create plots of results (0 or 1, default 0)
%
%  Examples:
%
%    mask = stmclustering('simtest0b_combined.mat')
%
%  (c) Freek van Leijen, Delft University of Technology, 2024.

% Created:  30 Aug 2024 by Freek van Leijen
% Modified: 
%

% Check input arguments and process options

if nargin < 1
   error('This function expects at least one input argument.')
end

opt.clusteringMethod = 'distance';   % clustering method, {'benchmark','distance'}
                                     % (default 'distance')
opt.spatialTolerance = 1;            % [m], maximum distance between benchmarks for 
                                     % 'distance' method, default 0.1 m
opt.doplots = 0;                     % make plots, default 0

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end

% Initialize
numObs = numel(stm.obsTypes);  % Number of elements in dataset

if isfield(stm.pntAttrib,'mapCrd')
  xy=stm.pntAttrib.mapCrd;
elseif isfield(stm.pntAttrib,'pntNeu')
  xy = stm.pntAttrib.pntNeu(:,[2 1])*1000; % convert to meters, original in km
else
  if exist('plh0','var') % use existing reference point
    pntNeu = plh2neusp(stm.pntCrd,plh0);  
  else
    [pntNeu,plh0] = plh2neusp(stm.pntCrd);
  end
  xy = pntNeu(:,[2 1]);
end

if opt.doplots > 0
  figure;
  scatter(xy(:,1),xy(:,2));
end


%% Detect common/clustered points

switch lower(opt.clusteringMethod)

  case 'benchmark'   % Detection based on benchmark id
  
    if isfield(stm.pntAttrib,'pntId')

      % Find unique point Ids
      pntId = char(stm.pntAttrib.pntId);
      pntId2 = [pntId(:,1:8) pntId(:,10:14)];
      [pntIdUnique,idxUnique] = unique(pntId2,'rows','stable'); %stable for original order of pntIds
      numUnique = numel(idxUnique);

      % Create structure of common points
      commonPoints = [];
      count = 0;
      for v = 1:numUnique
        idxCommon = strmatch(pntIdUnique(v,:),pntId2);

        if numel(idxCommon)>1
          count = count + 1;
          commonPoints(count).idx = idxCommon;
          commonPoints(count).numPoints = numel(idxCommon);
        end
    
      end
  
    else
  
      error('No pntId available, cannot find common points.');
  
    end

  
  case 'distance'   % Detection based on distance
  
    if opt.spatialTolerance == 0
  
      % Find unique points (based on exact same coordinates)
      [xyUnique,idxUnique] = unique(xy,'rows','stable'); %stable for original order of xy
      numUnique = numel(idxUnique);

      % Create structure of common points
      commonPoints = [];
      count = 0;
      for v = 1:numUnique
        idxCommon = find(xy(:,1) == xy(idxUnique(v),1) & ...
                         xy(:,2) == xy(idxUnique(v),2)  ...
                        );
    
        if numel(idxCommon)>1
          count = count + 1;
          commonPoints(count).idx = idxCommon;
          commonPoints(count).numPoints = numel(idxCommon);
        end
    
      end

    else % Detection of point clusters based on distance
  
      %%  Pairwise distance between observations
      %
      % Returns a vector d containing the Euclidean distances between each pair 
      % of observations. Use 'squareform' to make this a matrix.

      d = pdist(xy);

      %% Cluster data using distance criterion

      % Create a hierarchical cluster tree using the 'average' method and the 'euclidean' metric.
      %
      % Single linkage (default method), also called nearest neighbor, uses the smallest 
      % distance between objects in the two clusters. Average linkage uses the average 
      % distance between all pairs of objects in any two clusters.Value	
      % 
      % The metric is what is used by 'pdist' for evaluating distance. We use the 
      % 'pdist' default 'euclidean'. ('pdist' is used by 'linkage')
      %
      % Z is a matrix that encodes a tree containing hierarchical clusters of the rows 
      % of the input data matrix.

      Z = linkage(xy,'average','euclidean');

      % Cluster the data using a threshold of 1.5 for the 'distance' criterion.

      T = cluster(Z,'cutoff',opt.spatialTolerance,'Criterion','distance');

      % T contains numbers that correspond to the cluster assignments. Find the number of classes that cluster identifies.


      %% Create centroids for each cluster and count using grpstats

      %%[xycentroid,count]=grpstats(xy,T,{"mean","numel"});
      [xycentroid,count]=grpstats(xy,T,{'mean','numel'});
      count=count(:,1);

      %%% Plot the centroids with count

      numgroups=length(unique(count));

      figure;
      clr = hsv(numgroups);
      gscatter(xycentroid(:,1),xycentroid(:,2),count,clr,".",[1:numgroups]*3);

      TUnique = unique(T);
      numUnique = numel(TUnique);

      % Create structure of common points    
      commonPoints = [];
      numPoints = NaN(stm.numPoints,1);
      count = 0;
      for v = 1:numUnique
        idxCommon = find(T==TUnique(v));
        numPoints(idxCommon) = numel(idxCommon);
        if numel(idxCommon)>1
          count = count + 1;
          commonPoints(count).idx = idxCommon;
          commonPoints(count).numPoints = numel(idxCommon);
        end
    
      end
    
    end

  otherwise

    error('You specified an unsupported unique point clustering method.');
    
end

numCommonPoints = numel(commonPoints);

mask = true(stm.numPoints,1);

if numCommonPoints > 0

  %% Plot common time series (optional)

  if opt.doplots > 0
    for v = 1:numCommonPoints
  
      figure;
      h1 = NaN(commonPoints(v).numPoints,1);
      h2 = NaN(commonPoints(v).numPoints,1);
      h3 = NaN(commonPoints(v).numPoints,1);
      for w = 1:commonPoints(v).numPoints
        hh1 = subplot(3,1,1); h1(w) = plot(stm.epochDyear,stm.obsData(commonPoints(v).idx(w),:,1)); hold on;
        hh2 = subplot(3,1,2); h2(w) = plot(stm.epochDyear,stm.obsData(commonPoints(v).idx(w),:,2)); hold on;
        hh3 = subplot(3,1,3); h3(w) = plot(stm.epochDyear,stm.obsData(commonPoints(v).idx(w),:,3)); hold on;     
      end
      legend(h1,stm.pntName(commonPoints(v).idx));
      legend(h2,stm.pntName(commonPoints(v).idx));
      legend(h3,stm.pntName(commonPoints(v).idx));   
      ylabel(hh1,[stm.obsTypes{1} ' [mm]'],'fontweight','bold');
      ylabel(hh2,[stm.obsTypes{2} ' [mm]'],'fontweight','bold');
      ylabel(hh3,[stm.obsTypes{3} ' [mm]'],'fontweight','bold');
         
    end
  end


  %% Analyze common points

  % Select 'best' point based on
  % 1) most epochs
  % 2) longest time span
  % 3) random choice (first in list)
  % and create point mask

  for v = 1:numCommonPoints
    commonPointsIdx = commonPoints(v).idx;
    numEpochs = NaN(commonPoints(v).numPoints,numObs);
    timeSpan = NaN(commonPoints(v).numPoints,numObs);
    for w = 1:commonPoints(v).numPoints
      for z = 1:numObs
        pointMask = ~isnan(stm.obsData(commonPointsIdx(w),:,z));
        numEpochs(w,z) = nnz(pointMask);
        if numEpochs(w,z)~=0
          timeSpan(w,z) = max(stm.epochDyear(pointMask)) - min(stm.epochDyear(pointMask));
        end
      end
    end
  
    % Select 'best' point: 1) most epochs
    numEpochsTot = sum(numEpochs,2,'omitnan');
    maxEpochs = max(numEpochsTot);
    maxEpochsIdx = find(numEpochsTot==maxEpochs);
  
    if numel(maxEpochsIdx) > 1
      % Select 'best' point: 2) longest time span
      timeSpanTot = sum(timeSpan(maxEpochsIdx,:),2,'omitnan');
      maxTimeSpan = max(timeSpanTot);
      maxTimeSpanIdx = find(timeSpanTot==maxTimeSpan);
    
      if numel(maxTimeSpanIdx) > 1
        % Select 'best' point: 3) random choice (first in list)
        commonPoints(v).idxPrimary = commonPointsIdx(maxEpochsIdx(maxTimeSpanIdx(1))); %Store index primary/best point
        commonPointsIdx(maxEpochsIdx(maxTimeSpanIdx(1))) = []; % remove selected point ( 3) random choice) 
      else
        commonPoints(v).idxPrimary = commonPointsIdx(maxEpochsIdx(maxTimeSpanIdx)); %Store index primary/best point
        commonPointsIdx(maxEpochsIdx(maxTimeSpanIdx)) = []; % remove selected point ( 2) longest time span)
      end
    else
      commonPoints(v).idxPrimary = commonPointsIdx(maxEpochsIdx); %Store index primary/best point
      commonPointsIdx(maxEpochsIdx) = []; % remove selected point ( 1) most epochs) 
    end
  
    mask(commonPointsIdx) = false; % mask removed points
    commonPoints(v).idxSecondary = commonPointsIdx; %Store index secondary point(s)
      
  end

end

numRemoved = numel(find(mask==false));
display(['Number of point clusters: ' num2str(numCommonPoints)]);
display(['Number of removed points: ' num2str(numRemoved)]);

if numRemoved < numCommonPoints
  error('The number of removed points is smaller than the number of point clusters. This should not be possible.');
end

%% Plot the centroids with count (optional)

if opt.doplots > 0

  numgroups=length(unique(count));

  figure;
  clr = hsv(numgroups);
  gscatter(xy(mask,1),xy(mask,2),count,clr,".",[1:numgroups]*3);

end


