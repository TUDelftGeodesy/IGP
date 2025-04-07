function [predSignal,predSignalErrorVar,QzyQyyinvAFull] = stmprediction(stmIn,stmOut,emat,Qy,Atrend,varargin)
%STMPREDICTION   Spatio-temporal prediction based on a space time matrix.
%  [PREDSIGNAL,PREDSIGNALERRORVAR,QZYQYYINVAFULL] = STMPREDICTION(STMIN,STMOUT,EMAT,QY, ...
%  ATREND,OPTIONS) Function to predict displacements at a list of desired points and epochs 
%  based on an input space-time matrix STMIN and an output space-time matrix STMOUT. STMIN
%  contains the input data. STMOUT is assumed to contain at least the prediction points,
%  prediction epochs, prediction layers, temporal reference and the stochastic signal
%  model. The prediction is based on the residuals EMAT obtained after removal of a trend,
%  estimated using the ATREND design matrix and the QY covariance matrix. ATREND is used
%  here for error propagation.
%
%  The prediction is based on an iterative scheme, where predictions at the observation
%  points and epochs are used to converge to a solution (see OPTIONS). A 
%  The output of the function consists of the predicted signal PREDSIGNAL, the error
%  variance of the predicted signal PREDSIGNALERRORVAR, and the matrix QZYQYYINVAFULL,
%  which contains information needed for error propagation later on.
%
%  [...] = STMPREDICTION(STMIN,STMOUT,EMAT,QY,ATREND,'option',value,...) 
%  allows to specify options for the prediction
%
%    'convThres'             % Convergence threshold for iterative prediction [mm]
%    'maxIter'               % Maximum number of iterations for iterative prediction
%
%  (c) Freek van Leijen, Delft University of Technology, 2023.

% Created:   07 Nov 2023 by Freek van Leijen
% Modified: 
%            08 Jan 2024 by Freek van Leijen
%            - fixed bug preventing xMax, yMax, tMax to be included in the results
%            21 Nov 2024 by Freek van Leijen
%            - removed unnecessary call to parsestmstochmodel.m

% Check input arguments and process options

if nargin < 5
   error('This function expects at least three input arguments.')
end

fprintf('\n');
fprintf('Start prediction ...\n');

% Default options
opt.convThres = 10000;             % Convergence threshold for iterative prediction (sum of
                                   % absolute values) [mm].
opt.maxIter = 10;                  % Maximum number of iterations for iterative prediction

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end

%% Initialize
obsNeu = stmIn.pntAttrib.pntNeu;   % Observation points (Neu coordinates)
tIn = stmIn.epochDyear;            % Observation epochs

pntNeu = stmOut.pntAttrib.pntNeu;  % Prediction points (Neu coordinates)
tOut = stmOut.epochDyear;          % Prediction epochs
layers = stmOut.techniqueAttrib.obsIdx;  % Prediction layers
t0 = stmOut.temporalRef;           % Temporal reference
signalModel = stmOut.techniqueAttrib.signalModel; % Stochastic signal model for prediction

numPointsIn = stmIn.numPoints;     % Number of points in dataset 
numEpochsIn = stmIn.numEpochs;     % Number of epochs in dataset
numObsIn = numel(stmIn.obsTypes);  % Number of elements in dataset

m = numPointsIn * numEpochsIn;     % Number of observations per layer
mTot = m * numObsIn;               % Total number of observations in dataset
n = numPointsIn*2 + numEpochsIn;   % Number of unknowns per layer (rang defect is 1)

numPointsOut = size(pntNeu,1);     % Number of prediction points
numEpochsOut = numel(tOut);        % Number of prediction epochs
numLayers = numel(layers);         % Number of prediciton layers
numObsOut = max(layers);           % Maximum number of layers (for output array)

mp = numPointsOut * numEpochsOut;  % Number of predictions per layer

               
%% Get stochastic model signal and find maximum ranges to determine buffer size
if ~isempty(signalModel)
  maxRs = 0;
  maxRt = 0;
  for k = 1:numel(signalModel)
    [signalModelType{k},signalModelPar{k}] = parsestochmodel(signalModel{k}{1});
    if isfield(signalModelPar{k},'Rs')
      maxRs = max(maxRs,signalModelPar{k}.Rs);
    elseif isfield(signalModelPar{k},'L')
      maxRs = max(maxRs,signalModelPar{k}.L);
    end
    if isfield(signalModelPar{k},'Rt')
      maxRt = max(maxRt,signalModelPar{k}.Rt);
    end
  end
  if maxRs==0
    maxRs = Inf; %No spatial buffering
  end
  if maxRt==0
    maxRt = Inf; %No temporal buffering
  end
else
  error('You have to provide a covariance function for the signal prediction.');
end

% Set buffer of influence space and time
spaceBuffer = 3*maxRs;
timeBufer = 3*maxRt;
%spaceBuffer = maxRs;
%timeBufer = maxRt;
%spaceBuffer = 1;
%timeBuffer = 1;

%% Setup evaluation grid
xMin = min(pntNeu(:,2));
xMax = max(pntNeu(:,2));
yMin = min(pntNeu(:,1));
yMax = max(pntNeu(:,1));
tMin = min(tOut);
tMax = max(tOut);

xRange = [xMin:maxRs:xMax xMax+0.1]; % Grid of maxRs km, results in [xMin xMax] if maxRS==Inf, +0.1 to include xMax in loop
yRange = [yMin:maxRs:yMax yMax+0.1]; % Grid of maxRs km, +0.1 to include yMax in loop
tRange = [tMin:maxRt:tMax tMax+0.1]; % Grid of maxRt year, +0.1 to include tMax in loop
Nx = numel(xRange)-1;
Ny = numel(yRange)-1;
Nt = numel(tRange)-1;
      
%% Setup output arrays
predSignal = zeros(numPointsOut,numEpochsOut,numObsOut); %zeros because of iteration
predSignalErrorVar = NaN(numPointsOut,numEpochsOut,numObsOut);
predObs = zeros(numPointsIn,numEpochsIn,numObsIn);

numPointsPred = zeros(Nx,Ny,numObsOut);
numEpochsPred = zeros(Nx,Ny,numObsOut);
numPred = zeros(Nx,Ny,numObsOut);

%% Loop over layers
for k = layers

  % Get residuals for layer
  res = emat(:,:,k);
  
  % Create masks
  obsMask = ~isnan(res);
  pntMask = any(obsMask,2);
  epochMask = any(obsMask,1); 
  
  % Apply mask to residuals  
  res = res(obsMask(:));
  numObs = length(res);

  
  if numObs>0 %skip layer if no observations
    
    %%% Prediction
    fprintf('\n');
    fprintf('Predicting layer %d ...\n',k);

    %% Setup prediction system
    QzyQyyinvA = NaN(numPointsOut*numEpochsOut,numPointsIn); %For error propagation later on

    % Get covariance matrices
    QzyObs = stmstochmodel(signalModel{k},[],obsNeu(:,1:2),tIn,obsMask,t0);
    %QyyObs = Qy{k} + stmstochmodel(signalModel{count},[],obsNeu(:,1:2),dtIn,obsMask,t0);
    QyyObs = Qy{k} + QzyObs; % For observations, QzyObs == Qss
        
    % Calculate prediction weights for observations
    invQyyObs = inv(QyyObs);
    Wobs = QzyObs*invQyyObs;

    dObs = inf;
    numIter = 0;

    % Iterate until maxIter or sum(abs(increments) below convThres
    while dObs(end) > opt.convThres & numIter <= opt.maxIter
      numIter = numIter+1;

      dpredObs = zeros(numPointsIn,numEpochsIn);
      dpredSignal = NaN(numPointsOut,numEpochsOut);

      % Get observations
      predObsLayer = predObs(:,:,k);
      obs = res - predObsLayer(obsMask(:));
            
      % Least-squares prediction observations
      dpredObs(obsMask) = Wobs*obs;

      dObs = [dObs; sum(abs(dpredObs(:)))];
      predObs(:,:,k) = predObs(:,:,k) + dpredObs;


      %% Prediction

      for v = 1:Nx
  
        fprintf('Predicting buffer %3d of %3d for layer %3d, iteration %3d ....\n',v,Nx,k,numIter);

        idx1a = find(pntNeu(:,2)>=xRange(v) & pntNeu(:,2)<xRange(v+1));

        if ~isempty(idx1a)

          for w = 1:Ny

            idx1b = find(pntNeu(idx1a,1)>=yRange(w) & pntNeu(idx1a,1)<yRange(w+1));

            if ~isempty(idx1b)

              % Get number of points, epochs
              numPointsPred(v,w,k) = numel(idx1b);
              numEpochsPred(v,w,k) = numEpochsOut;
              numPred(v,w,k) = numPointsPred(v,w,k)*numEpochsPred(v,w,k);

              % Create prediction mask
              pntIdxPred = idx1a(idx1b);
              pntMaskPred = logical(zeros(numPointsOut,1));
              pntMaskPred(pntIdxPred) = 1;
              epochIdxPred = 1:numEpochsOut;
              epochMaskPred = logical(ones(1,numEpochsOut));
              maskPred = {pntMaskPred epochMaskPred};
    
              % Get prediction matrix
              Qzy = stmstochmodel(signalModel{k},[],obsNeu(:,1:2),tIn,obsMask,t0, ...
                                        pntNeu(:,1:2),tOut,maskPred);


              % Least-squares prediction
              W = Qzy*invQyyObs;
              dpredSignal(pntMaskPred,epochMaskPred) = reshape(W*obs, ...
                  numPointsPred(v,w,k),numEpochsPred(v,w,k));

              if numIter == 1 % Only needed once
              
                % Construct full prediction matrix for later use (error prop).
                maskPredFull = false(numPointsOut,numEpochsOut);
                maskPredFull(pntMaskPred,epochMaskPred) = true;

                QzyQyyinvA(maskPredFull(:),pntMask) = W*Atrend{k}(:,1:numel(find(pntMask))); %only velocity part of Atrend
                
                % Prediction error variance                    
                Qzz = stmstochmodel(signalModel{k},[],pntNeu(:,1:2),tOut,maskPred,t0);
                Pzhat = Qzz - W*Qzy';

                predSignalErrorVar(pntMaskPred,epochMaskPred,k) = ...
                  reshape(diag(Pzhat), ...
                  numPointsPred(v,w,k),numEpochsPred(v,w,k));
                    
              end

            end
          end % for Ny
        end
      end % for Nx

      predSignal(:,:,k) = predSignal(:,:,k) + dpredSignal;

    end %end while
    
    % Store for later use in prediction error calculation (see stmrestore.m)
    QzyQyyinvAFull{k} = QzyQyyinvA;

  end % end if observations
  
end % end loop layers

end


% Sub-functions

function [name,params]=parsestochmodel(str)
%PARSESTOCHMODEL  Parse string with stochastic model specification.
%   [name,pars]=PARSESTOCHMODEL(STR) parses the string STR with the stochastic
%   model specification and returns the model name NAME and a structure
%   PARAMS with the stochastic model parameters.
%
%   Example
%
%      [stochModel,stochParameters]=parsestochmodel('tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)')
%
%   (c) Hans van der Marel, Delft University of Technology, 2020

% Created:  10 September 2020 by Hans van der Marel
% Modified:

str=strrep(str,' ','');

ii=regexp(str,'[()]');
name=str(1:ii(1)-1);
pargs=strsplit(str(ii(1)+1:ii(end)-1),',');

params=[];
for k=1:numel(pargs)
   if ~isempty(pargs{k})
      p=strsplit(pargs{k},'=');
      value=str2num(p{2});
      if isempty(value)
         value=p{2};
      end
      params.(p{1})=value;
   end
end


end

