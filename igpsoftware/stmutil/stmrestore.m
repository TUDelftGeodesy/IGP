function [predTrend,predTrendErrorVar,predTotal,predTotalErrorVar,vel_int] = stmrestore(stmIn,stmOut,vel,Qvel,predSignal,predSignalErrorVar,QzyQyyinvAFull,varargin)
%stmrestore   Restore the trend after prediction.
%  [PREDTREND,PREDTRENDERRORVAR,PREDTOTAL,PREDTOTALERRORVAR,VEL_INT] =
%  STMRESTORE(STMIN,STMOUT,VEL,QVEL,PREDSIGNAL,PREDSIGNALERRORVAR,QZYQYYINVAFULL,OPTIONS)
%  Function to restore the trend after prediction. Using the metadata in STMIN and
%  STMOUT, the velocities in VEL are interpolated to the prediction points to obtain
%  VEL_INT. The interpolation method can be specified via the OPTIONS (see below). Based
%  on these, the trend at the prediction points PREDTREND is evaluated. Using
%  the covariance matrix of the velocities QVEL and QZYQYYINVAFULL (calculated in the 
%  prediction.m function), the variances of the interpolated trend PREDTRENDERRORVAR are
%  calculated (since the number of prediction points can in principle be very large, it
%  is not feasible to calculate the full covariance matrix). 
%
%  Finally, the total prediction can be calculated as PREDTOTAL = PREDTREND + PREDSIGNAL.
%  Similarly, for the variance holds PREDTOTALERRORVAR = PREDTRENDERRORVAR + 
%  PREDSIGNALERRORVAR.
%
%  [...] = STMRESTORE(...,'option',value,...) allows to specify options
%  for the restore step
%
%  'interpMethod'               interpolation method, 'linear' for triangulation-based 
%                               linear interpolation or 'idw' for inverse distance weighting
%                               (default 'linear')
%  'idwPower'                   power for inverse distance weighting (if used) (default 2)
%  't0'                         reference epoch {'mean' , 'best', <dyear>} (default 'best')
%
%  The 'linear' and 'idw' interpolation functions are provided as subfunctions.
%  In case the 'linear' interpolation method is used, a Delaunay triangulation is used.
%  This means that no interpolation (extrapolation) is possible outside the convex hull of
%  the observation points (resulting in NaN for prediction points outside the convex hull).
%  With the 'idw', also values outside the convex hull are obtained.
%
%  (c) Freek van Leijen, Delft University of Technology, 2023.

% Created:   16 Nov 2023 by Freek van Leijen
% Modified: 
%            08 Jan 2024 by Freek van Leijen
%            - fixed bug resulting in zerors instead of NaNs in restoring trend
%

% Check input arguments and process options

if nargin < 7
   error('This function expects at least six input arguments.')
end

fprintf('\n');
fprintf('Start restore ...\n');

opt.interpMethod = 'linear';       % interpolation method, 'linear' for triangulation-based 
                                   % linear interpolation or 'idw' for inverse distance weighting
opt.idwPower = 2;                  % power for inverse distance weighting (if used)
opt.t0 = 'best';                   % reference epoch {'mean' , 'best', <dyear>} (default 'best')

for k = 1:2:numel(varargin)
    opt.(varargin{k}) = varargin{k+1};
end


%% Initialize
obsNeu = stmIn.pntAttrib.pntNeu;   % Observation points (NEU coordinates)
pntNeu = stmOut.pntAttrib.pntNeu;  % Prediction points (NEU coordinates)
tIn = stmIn.epochDyear;            % Observation epochs
tOut = stmOut.epochDyear;          % Prediction epochs

layers = stmOut.techniqueAttrib.obsIdx;  % Prediction layers
numLayers = numel(layers);

numPointsIn = size(obsNeu,1);
numEpochsIn = numel(tIn);
numPointsOut = size(pntNeu,1);
numEpochsOut = numel(tOut);
numOut = numPointsOut*numEpochsOut;

pntMaskAll = true(numPointsIn,1);
epochMaskAll = true(1,numEpochsIn);

% Reference epoch and delta times
if isnumeric(opt.t0)
  t0 = opt.t0;
else

  switch opt.t0
    case 'mean'
       t0 = mean(tIn(epochMaskAll));
    case 'best'
       epochCount = sum(sum(~isnan(stmIn.obsData(pntMaskAll,epochMaskAll,:)),3),1);
       t0 = sum(tIn(epochMaskAll) .* epochCount) / sum(epochCount);
    otherwise
       error('You specified an invalid t0.');
  end
  
end
dtOut = tOut(:) - t0;


%% Interpolation of velocities
vel_int = NaN(numPointsOut,numLayers);
coverageMask = true(numPointsOut,3);

switch opt.interpMethod
  case 'idw'

    % Inverse distance interpolation
    for l = 1:numLayers
      pntMask{l} = any(vel(:,l),2); % Mask can be different for every layer
      weights{l} = idw(obsNeu(pntMask{l},1:2),pntNeu(:,1:2),opt.idwPower); %inverse distance weighting
      vel_int(:,l) = weights{l}*vel(pntMask{l},l);
    end

  case 'linear'

    % Triangulation-based linear interpolation 
    for l = 1:numLayers
      pntMask{l} = any(vel(:,l),2); % Mask can be different for every layer
      [weights{l},coverageMask(:,l)] = triangulation_linear(obsNeu(pntMask{l},1:2),pntNeu(:,1:2));
      vel_int(coverageMask(:,l),l) = weights{l}(coverageMask(:,l),:)*vel(pntMask{l},l);
    end

  otherwise
    error('You specified an unsupported interpolation technique.');
end


%% Construct velocity-based time series
i0p = repmat((1:numPointsOut)',[numEpochsOut 1]);
a0p = reshape(repmat(dtOut',[numPointsOut 1]), numOut, 1);
Ap = sparse((1:numOut)' , i0p, a0p, numOut, numPointsOut, numOut);

predTrend = NaN(numPointsOut,numEpochsOut,numLayers);
predTrendErrorVar = NaN(numPointsOut,numEpochsOut,numLayers);
for l = 1:numLayers

  % Calculate time series
  predTrend(:,:,l) = reshape(Ap*vel_int(:,l),numPointsOut,numEpochsOut);
  predTrend(isnan(vel_int(:,l)),:,l) = NaN; % Needed because sparse 0 * NaN gives 0 instead of NaN
  
  % Calculate Trend error variance
  R = Ap*weights{l} - QzyQyyinvAFull{l}(:,pntMask{l});
  
  % original implementation, takes too much memory
  %predTrendErrorVar(:,:,l) = ...
  %                reshape(diag(R*Qvel{l}(1:numel(find(pntMask{l})),1:numel(find(pntMask{l})))*R'), ...
  %                numPointsOut,numEpochsOut);

  R2 = R*Qvel{l}(1:numel(find(pntMask{l})),1:numel(find(pntMask{l})));
  Rdiag = sum(R2.*R,2);
  predTrendErrorVar(:,:,l) = reshape(Rdiag,numPointsOut,numEpochsOut);
  predTrendErrorVar(~coverageMask(:,l),:,l) = NaN; % Apply coverageMask for same coverage as PredTrend

end

%% Restore trend
predTotal = predTrend + predSignal;
predTotalErrorVar = predTrendErrorVar + predSignalErrorVar;

end

function weights = idw(xy,xyp,p)
% IDW Inverse distance weights.
%
% Input:  - xy       spatial coordinates observation points
%         - xyp      spatial coordinates prediction points
%         - p        power of the distance. E.g.,
%                    p = 0: mean
%                    p = 1: inverse distance weighting
%                    p = 2: inverse square distance weighting
%                    p can take any value >= 0.
%
% Output: - weights  calculated weight matrix
% 
%   (c) Freek van Leijen, Delft University of Technology, 2023

% Created:  03 January 2023 by Freek van Leijen
% Modified:

d = pdist2(xyp,xy);
w = d.^-p;
sumw = sum(w,2);
weights = w.*repmat(1./sumw,1,size(xy,1));

end

function [weights,coverageMask] = triangulation_linear(xy,xyp)
% Triangulation-based linear interpolation.
%
%
% Input:  - xy       spatial coordinates observation points
%         - xyp      spatial coordinates prediction points
%
% Output: - weights  calculated weight matrix
% 
%   (c) Freek van Leijen, Delft University of Technology, 2023

% Created:  29 October 2023 by Freek van Leijen
% Modified:

N = size(xy,1);
Np = size(xyp,1);

dt = delaunayTriangulation(xy);
[ti,bc] = pointLocation(dt,xyp); % triangle index, weights

coverageMask = ~isnan(ti); % points within convex hull
Nc = sum(coverageMask);
tp = dt(ti(coverageMask),:); % triangulation points
bc = bc(coverageMask,:); % triangulation weights

weights = sparse(Np,N);
weights(sub2ind([Np N],repmat(find(coverageMask),3,1),tp(:))) = bc(:);

end

