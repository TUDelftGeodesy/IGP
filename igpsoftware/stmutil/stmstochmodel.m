function R=stmstochmodel(stochModel,stochData,xy,t,mask,varargin)
%STMSTOCHMODEL  Compute (square-root) covariance/precision matrix. 
%   R=STMSTOCHMODEL(STOCHMODEL,STOCHDATA) builds the covariance matrix 
%   R from STOCHMODEL and STOCHDATA. STOCHMODEL is a character string (or 
%   cell array) with the model specification and STOCHDATA is a numerical 
%   array/matrix with the stochastic data. The output R is a covariance 
%   matrix matching the order of observations in the (vectorized) space 
%   time matrix. This form is only possible when the co-variance matrix is 
%   stored explicitly, for parameterized model computations use the form
%   below.
%
%   R=STMSTOCHMODEL(STOCHMODEL,STOCHDATA,XY,T,MASK) computes the 
%   covariance matrix R from the stochastic model data in STOCHMODEL, for 
%   points with coordinates in XY and times in T, masking the entries in 
%   MASK. STOCHMODEL is a character string or cell array which specifies 
%   the stochastic model, stochastic model parameters and/or how the 
%   stochastic model data is stored in the numerical array/matrix STOCHDATA.
%   STOCHDATA can be an empty array when the co-variance matrix is generated 
%   from parameterized model data only. MASK is a logical matrix of the same
%   dimensions as the space time matrix, with true for the elements to
%   keep, or a cell array MASK={pntMask,epochMask[,mask]} with logical vectors
%   pntMask and epochMask with true for points and epochs to keep. If MASK 
%   is empty, or missing, no entries are excluded. The output R is a 
%   covariance matrix matching the order of observations in the (vectorized) 
%   space time matrix, as determined by the size of XY, T and the contents 
%   of MASK.
%
%   R=STMSTOCHMODEL(...,PNTATTRIB,EPOCHATTRIB) allows to specify extra
%   attributes required by some models. Details are provided with each model.
%
%   R=STMSTOCHMODEL(...,TREFORIG) allows to specify extra attribute trefOrig
%   required by sure models in case tlag is used as parameter. Details are
%   provided with each model.
%
%   R=STMSTOCHMODEL(...,XY0,T0,MASK0) enables the calculation of the
%   prediction matrix for prediction, with XY0 and T0 the prediction
%   points and epochs, and MASK0 a cell array MASK={pntMask0,epochMask0} 
%   with logical vectors pntMask0 and epochMask0 with true for points and 
%   epochs to keep.  The prediction matrix with is returned in R. 
%   This syntax is only supported for the tudpred1 and sure models. 
%
%   R=STMSTOCHMODEL(...,OUT) with OUT = [ 'COV' | 'PREC' ] returns in R the 
%   covariance matrix (Default or 'COV') or the precision matrix ('PREC').
%   OUT must be the last parameter in the argument list. The precision matrix 
%   is the inverse of the co-variance matrix and sometimes can be computed 
%   directly without numerical inversion. For some models, the precision 
%   matrix can be sparse, whereas the co-variance matrix is full. This 
%   option is not yet inplemented, output is always a covariance matrix.
%
%   XY a is a matrix with the X and Y coordinates for the points 1...M, 
%   the first column contains the X coordinate, the second column must contain 
%   the Y coordinate. The units are [km]. The number of rows is M. 
%
%   T is an array with the decimal year for the epochs 1...N. The length of
%   T is N. 
%
%   MASK is a M-by-N logical matrix (or array) with false for 
%   missing observations. If multiple observed components are involved 
%   (e.g. for North, East, and Up components) MASK is a three dimensional 
%   M-by-N-by-P logical matrix, with P the number of components.
%   If MASK is empty, or missing, no entries are excluded. 
%
%   The output R is a r-by-r matrix, with r=M*N or r=M*N*P, less the number 
%   of masked elements in the MASK logical matrix. For some models the
%   output R can be a sparse matrix (this is not yet implemented, output
%   is always a full matrix).
%
%   STOCHMODEL is a character string or cell array of character strings
%   with the model specification. The model specification contains
%   the stochastic model name and the parameters for each model. When 
%   STOCHMODEL is a cell array the covariance matrix is constructed using
%   several layers, each with it own model specification. The order of
%   the stochastic layers is important. As a general rule of thumb start
%   with 1) the temporal stochastic models, 2) spatial models and 3) combined
%   spatio-temporal models. 
%
%   Supported stochastic models are:
%
%   stochModel (temporal)                         Attributes
%   --------------------------------------------- ---------------
%   WN(sigma=#[,fs=#,ndays=#])                    t[,epochAttrib]
%   FOGM(rho=#,sigma=#[,fs=#,ndays=#])            t[,epochAttrib]
%   PL(alpha=#,sigma=#[,fs=#,ndays=#])            t[,epochAttrib]
%   iGGM(alphaR=#,alphaL=#,gamma=#,sigma=#        t[,epochAttrib]
%                            [,fs=#,ndays=#])
%   houtenbos(s2t=#,pt=#)                         t
%
%   sigma acts as a scale factor (amplitude), its unit is linked to the 
%   square root of the covariance matrix, e.g. mm. 
%
%   rho is the exponential decay in [1/days] for the FOGM, or AR(1), process.
%
%   alpha, alphaR and alphaL are spectral indices for the PL and iGGM 
%   respectively, these can be real valued number between -3 and 1. For
%   the range -3 < alpha < -1 we have non-stationary fractional Brownian 
%   motion, for -1 < alpha < 1 stationary fractional "white" noise. For
%   integer alpha=0 we have white noise, alpha=-1 flicker noise and alpha=-2 
%   random walk. The factor gamma (0 -> 1) for iGGM defines the cross-over 
%   between alphaR (high frequency) and alphaL (low frequency) part. 
%
%   fs (default 1 sample/day) and ndays (default 1 day) are optional
%   parameters. The value of ndays is numeric or a string with the
%   fieldname in epochAttrib.
%
%   s2t is the temporal variance in the houtenbos model.
%
%   pt is the exponential decay in the houtenbos model.
%
%   stochModel (spatial)                          Attributes
%   --------------------------------------------- ---------------
%   spatialcov(rho=#[,q=#])
%
%   rho is the exponential decay in [1/km] and q and optional exponent,
%   with default q=1 (exponential co-variance function), q=2 returns the
%   co-variance matrix for a Gaussian co-variance function.   
%
%   stochModel (spatio-temporal)                  Attributes
%   --------------------------------------------- ---------------
%   tudinsar4(s20=#,s2t=#,s2s=#,Rt=#,Rs=#)        xy,t
%   tudinsar4rd(s20=#,s2t=#,s2s=#,Rt=#,Rs=#)      xy,t,pntAttrib,epochAttrib
%   tudpred1(s20=#,s2t=#,s2s=#,Rt=#,Rs=#)         xy,t[,xy0,t0,mask0]
%   sure(s2z=#,L=#,q=#,tref=#)                    xy,t[,xy0,t0,mask0]
%   sure(s2z=#,L=#,q=#,tlag=#)                    xy,t,trefOrig[,xy0,t0,mask0]
%
%   pntAttrib and epochAttrib are required for tudinsar4rd model,
%   xy0 and t0 are optional attributes in case the tudpred1 and sure
%   models are used for prediction, without xy0 and t0 these models
%   only provide the co-variance matrix. For the sure model you can
%   choose between tref or tlag, if tlag is used, and extra attribute
%   trefOrig is required. 
%
%   stochModel (data)                             Attributes
%   --------------------------------------------- ---------------
%   stdev([sigma=#,size=#,odata=#])               stochData        
%   covmatrix(format=[full|lower|upper])          stochData
%   covmatrix(format=index,lidx=#,oidx=#,odata=#) stochData
%   covmatrix(format=symdiags,[nd=#])             stochData
%   covmatrix(format=blkdiag,[nd=#],[sreg=#])     stochData
%   precmatrix(format=...)                        stochData
%
%   The value of format can be 'full', 'index', 'lower', 'blkdiag' or 'symdiags'. 
%
%   All formats have optional parameters 'size=#' and 'odata=#'
%   - 'size' is a scalar with the number of columns and rows in the matrix,
%      it can also be given as array with number of points, epochs (and 
%      observation types) in which case size is the product of the
%      elements.
%   - 'odata'+1 is the position in stochData with the first element of 
%      the covariance matrix starts. 
%   The last element of the element in stochdata depends on the format. 
%   For a full matrix the last elements is in 'odata'+'size'^2.
%
%   For format 'index' two additional parameters are required, 'ldix' and
%   'odix', and 'odata' is mandatory. 'lidx' is the lenght of the index and
%   'odix' is where the index starts in stochData. The values of index
%   are the row/column numbers in the matrix. The co-variance matrix (with
%   'lidx rows/columns' can be found in stochData starting at 'odata'+1 and
%   ending at 'odata'+'ldix'^2.
%   
%   For 'blkdiag' and 'symdiags' an extra parameter nd=# with the number of 
%   diagonal blocks respectively diagonals is required (optional for
%   'symdiags'). Block diagonal matrix are assumed to have a rankdefect of
%   1 for each block and will be regularized using the parameter 'sreg'. 
%   Default value for 'sreg' is 0.3e-4 (sqrt(1e-9)). If they are full-rank, 
%   use 'sreg=0' to skip regularization, though the regularization won't 
%   hurt even in case of full rank.
%
%   The co-variance matrix is a symmetric r-by-r matrix. The full covariance
%   matrix can be stored in (r+1)/2 stm layers ( r-by-(r+1)/2 matrix ).
%   With symdiags we store the matrix in diagonals, with lower it is stored
%   as columnwise lower triangle, and with blkdiag only the diagonal block
%   for each epoch is stored.   
%
%   Examples:
%
%     C=stmstochmodel('stdev()',stochData);
%     C=stmstochmodel('stdev(sigma=.5)',stochData);
%     C=stmstochmodel('covmatrix(format=full)',stochData);    
%     C=stmstochmodel(['covmatrix(format=blkdiag,nd=' num2str(numEpochs) ')'],stochData);    
%
%     C=stmstochmodel('tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11'), ...
%           [],xy,t);                      % Radarsat-2
%     C=stmstochmodel('tudinsar4(s20=9.49,s2t=4.53,s2s=4.96,Rt=0.70,Rs=1.09'), ...
%           [],xy,t,mask);                 % Sentinel-1 with observation mask
%
%     C=stmstochmodel('tudpred1(s20=3,s2t=2,s2s=2,Rt=1.0,Rs=3.5'), ...
%           [],xy,t,mask);                 % Covariance matrix
%     C=stmstochmodel('tudpred1(s20=3,s2t=2,s2s=2,Rt=1.0,Rs=3.5'), ...
%           [],xy,t,mask,xy0,t0,mask0);    % Prediction matrix
%
%     C=stmstochmodel({ ...
%           'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=1)', ...  
%           'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=1)' ...
%           'spatialcov(rho=0.0887)'}, ...
%           [],xy,t,mask);                 % recommended NAM/06-GPS model continuous GNSS
%     C=stmstochmodel({ ...
%           'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=ndays)', ...  
%           'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=ndays)' ...
%           'spatialcov(rho=0.0887)'}, ...
%           'WN(ndays=setup)'}, ...
%           [],xy,t,mask,pntAttrib,epochAttrib); % recommended NAM/06-GPS model campaign GNSS
%
%   Known limitations:
%   - fieldnames for tudinsar4rd are fixed (better to have these as variables)
%   - units of rho in FOGM should be verified
%   - The output is always a full covariance matrix (support for precision
%     matrix and sparse matrices not yet implemented)
%   - mask does not work in situations where stochData is combined with models
%
%   References:
%
%   [1] F.J. van Leijen, S. Samiei-Esfahany, H. van der Marel and R.F. Hanssen,
%       Improving the Functional and Stochastic Model of InSAR, NAM INSAR FM SM
%       Project, Report for NAM, 28 April 2020.
%
%   [2] Simon Williams, Description of GPS uncertainties within the Long Term 
%       Study on Anomalous Time-Dependent Subsidence, Report for NAM, 2015.
%
%   [3] A.P.E.M. Houtenbos, 2004. Subsidence Residual Modeling, SURE user manual, 
%       A.P.E.M. Houtenbos Geodetic Consultancy.
%
%   [4] A.P.E.M. Houtenbos and F. Kenselaar, 2001. Peilmerk hoogte variaties,
%       Stochastische analyse van peilmerkbeweging in Nederland, Delft University
%       of Technology (in Dutch).
%
%
%   GNSS temporal stochastic models for Vertical Continuous data [2]
%   ------------------------------------------------------------
%
%   Model 1: General 06-GPS/NAM model for Vertical Continuous data
%
%      iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sig=0.0967,fs=24,ndays=1)  
%
%   Model 2: General 06-GPS nam model plus 0.5 mm/sqrt(year) of random walk
%
%      iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sig=0.0967,fs=24,ndays=1)  
%      PL(alpha=-2,sig=0.00534,fs=24,ndays=1)
%
%   with sig= 0.00534 = 0.5 / sqrt(24*3600*365.25/Ts) = 0.5 / sqrt(24*365.25) = 
%         = 0.5 / sqrt(fs*365.25) ( fs= samples/day ).
%
%   Model 3 is model 2 with 1mm/sqrt(year) random walk.
% 
%   Model 4 is model 1 with sig=0.0967*3.
%
%   Model 5: 06-GPS/NAM modified at low frequencies to mimic the amplitude of the 
%   regional noise. 
%
%      iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.997380,sig=0.0967,fs=24,ndays=1)  
%
%   The main change is the parameter gamma (in the report the value of gamma 
%   is slightly different, 0.99238), which shifts the cross over frequency to a lower 
%   value which then increases the power at low frequencies. 
%
%   Model 6: General 06-GPS nam model with an additional FIGGM to mimic the 
%   regional noise. 
%
%      iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sig=0.0967,fs=24,ndays=1)  
%      iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sig=0.04,fs=24,ndays=1)  
%
%   This is similar to the argument for additional random walk (of model 3)
%   so that only powers at low frequencies are affected except the spectrum 
%   will tend to 0.8525 instead of -2 at very low frequencies (in the report 
%   the value of gamma is slightly different, 0.99806)
%
%   Model 7: A model based on the daily solutions from the regional filtered data 
%   from NGL. 
%
%      PL(alpha=-0.8525,sig=1.5,fs=1,ndays=1)
%      WN(sig=2.13)
%
%   This model is a power law with spectral index = -0.8525 and amplitude of 
%   1.50 mm (5.26mm/yrK4) in CATS notation) and white noise of amplitude 2.13 mm. 
%
%
%   GNSS temporal stochastic models for Horizontal Continuous data [2]
%   --------------------------------------------------------------
%
%   Model 1: General 06-GPS/NAM model for Horizontal Continuous data
%
%      iGGM(alphaR=-1.91,alphaL=-0.70,gamma=0.9137,sig=0.1481,fs=24,ndays=1)  
%
%   Model 2: General 06-GPS/NAM model modified at low frequencies to mimic the 
%   amplitude of the regional noise (same reasoning as Model 5 for vertical)
%
%      iGGM(alphaR=-1.91,alphaL=-0.8281,gamma=0.9771,sig=0.1481,fs=24,,ndays=1)  
%
%   The main change is the parameter gamma, which shifts the cross over frequency 
%   to a lower value which then increases the power at low frequencies. 
%
%   Model 3: General  06-GPS/NAM model with an additional FIGGM to mimic the 
%   regional noise.  
%
%      iGGM(alphaR=-1.91,alphaL=-0.70,gamma=0.9137,sig=0.1481,fs=24,,ndays=1)  
%      iGGM(alphaR=-2,alphaL=-0.8281,gamma=0.992,sig=0.06,fs=24,ndays=1)  
%
%   Model 4: A model based on the daily solutions from the regional filtered data 
%   from NGL. 
%
%      PL(alpha=-0.8281,sig=0.76,fs=1,ndays=1)
%      WN(sig=0.63)
%
%   This model is a power law with spectral index = -0.8281 and amplitude of 
%   0.76 mm (2.58 mm/yrK4) in CATS notation) and white noise of amplitude 0.63 mm. 
%
%
%   The spatial co-variance for GNSS is generated from a unit spatial
%   co-variance matrix with an exponential decay (comparable to FOGM model
%   for the temporal model). The model specifications from the LTS2 study
%   are
%
%      spatialcov(rho=0.0827)  
%      spatialcov(rho=0.1291)  
%      spatialcov(rho=0.0887)
%
%   for respectively the North, East and Vertical (Up) components. The 
%   unit for the exponential decay rho is [km^-1].
%
%   (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:   3 Sep 2020 by Hans van der Marel
% Modified: 14 Sep 2020 by Hans van der Marel
%            - Release for internal review 
%           28 Oct 2020 by Hans van der Marel
%            - extra parameters for tudinsar4rd model
%            - added input option for mask as cell
%           02 Nov 2020 by Hans van der Marel
%            - patch for NaN's in tudinsar4rd model and single/double
%            - optimized pntMask and epochMask operation on data, now 
%              also implemented efficiently for blkdiagonal matrices
%           07 Nov 2020 by Freek van Leijen
%            - added tudpred1 model for prediction (full spatio-
%              temporal covariance needed for prediction and calculation
%              of prediction matrix).
%            - added local function predexp for the calculation of 
%              components needed for the prediction matrix
%           20 Nov 2020 by Freek van Leijen
%            - added sure model for prediction
%            - extension of the covexp and predexp functions for
%              Gaussian covariance function
%           22 Dec 2020 by Freek van Leijen and Hans van der Marel 
%            - imported undocumented changes from Freek
%            - correct implementation of tudinsar4rd model?
%           10 Feb 2021 by Freek van Leijen
%            - implemented houtenbos model and added documentation
%           18 Oct 2021 by Hans van der Marel
%            - added regularization of block diagonal matrices (default)
%            - updated description
%           14 Nov 2021 by Freek van Leijen
%            - enabled observation mask for prediction matrix
%           07 Nov 2023 by Freek van Leijen
%            - changed the SURE model implementation using a time lag
%              in combination with the reference time. 
%            - time reference inserted as extra input parameter.
%           21 Nov 2023 by Hans van der Marel
%            - support both tref and tlag with sure model
%            - made the function backwards compatible 
%            - added documentation
%            5 Dec 2023 by Hans van der Marel
%            - fixed standard fieldname logic in tudinsar4rd model
%            - major change in covmatrix and precmatrix, support indexed
%              models, offsets and size specification
%           16 Jan 2024 by Hans van der Marel
%            - minor changes to improve logic to covmatrix, index is now a 
%              format, more general support for point and epoch masks
%            - removed variance as model
%            - added size and odata parameters to stdev and covmatrix
%            - added documentation
%           30 Jan 2024 by Hans van der Marel
%            - changed specification of size parameter

% Check the number and type of input arguments

if nargin < 2
   error('This function expects at least two input arguments.')
end
if ischar(stochModel) || isstring(stochModel)
   stochModel=cellstr(stochModel);
elseif ~iscell(stochModel)
   error('stochModel must be a cell array or string.')
end
if ~isnumeric(stochData) && ~isempty(stochData)
   error('stochData must be numeric or an empty array.')
end

% Number of points and epochs (only when xy and t present)

if nargin > 3
  [m,k]=size(xy);
  if ~isempty(xy) && k ~= 2, error('xy must have columns with X and Y coordinates.'); end
  n=numel(t);
  t=t(:);
  pntMask=true(m,1);
  epochMask=true(n,1);
else
  pntMask=[]; 
  epochMask=[];
end

% Check MASK

if nargin < 5 || isempty(mask)
   mask=[];
elseif iscell(mask)
   pntMask=mask{1};
   epochMask=mask{2};
   if numel(mask) > 2
      mask=mask{3};
   else
      mask=[];
   end
   if isempty(xy) 
      m=numel(pntMask);
   elseif ( numel(pntMask) == m )
      xy=xy(pntMask,:);
      m=size(xy,1);
   else
      error('Dimensions of pntMask and XY do not match.')
   end
   if isempty(t)
      n=numel(epochMask);
   elseif ( numel(epochMask) == n )
      t=t(epochMask);
      n=numel(t);
   else
      error('Dimensions of epochMask and T do not match.')
   end
%    % For data model compute the mask
%    if isempty(xy) && isempty(t)
%       mask0=true(m,n);
%       mask0(~pntMask,:)=false;
%       mask0(:,~epochMask)=false;
%    else
%       mask0=[];
%    end
else
   if m ~= size(mask,1) || n ~= size(mask,2)
       error('Dimensions of MASK must match the dimensions of  XY and T.') 
   end
end

% Check output type

out='COV';
if nargin == 3 && ischar(xy)
   out=xy;
elseif nargin > 5 && ~isempty(varargin) && ischar(varargin{end})
   out=varargin{end};
end
if ~any(strcmpi(out,{'COV','PREC'}))
   error('Invalid value for OUT.') 
end


% Build the stochastic model

% We need to keep track of where we are
%  
%    stage = 1  - Build up of temporal model (temporal submatrix numEpochsx numEpochs )
%    stage = 2  - Build up of spatial correlation (temporal submatrix -> full spatio-temporal matrix) 
%    stage = 3  - Spatio-temporal components (full matrix)

stage=0;

for k=1:numel(stochModel)

  [modelType,modelParameters]=parsestochmodel(stochModel{k});

  % Evaluate the model

  switch modelType
       
    case {'PL','iGGM','FOGM','WN'}
        
      %   stochModel (temporal)                         Attributes
      %   --------------------------------------------- ---------------
      %   WN(sigma=#[,fs=#,ndays=#])                    t[,epochAttrib]
      %   FOGM(rho=#,sigma=#[,fs=#,ndays=#])            t[,epochAttrib]
      %   PL(alpha=#,sigma=#[,fs=#,ndays=#])            t[,epochAttrib]
      %   iGGM(alphaR=#,alphaL=#,gamma=#,sigma=#        t[,epochAttrib]
      %                            [,fs=#,ndays=#])
      %
      %   fs (default 1 sample/day) and ndays (default 1 day) are optional
      %   parameters. The value of ndays is numeric or a string with the
      %   fieldname in epochAttrib.

      % checkstochmodel(modelParameters,{'..'}) is done inside tempcov 
      
      if stage <=0               
         % Initizalize stage 1 processing
          C=tempcov(modelType,modelParameters,t,varargin{:});
          stage=1;
      elseif stage == 1
         % Add layer to temporal co-variance matrix
         C=C+tempcov(modelType,modelParameters,t,varargin{:});
      elseif stage > 1
         % Add layer to full co-variance matrix
         C=tempcov(modelType,modelParameters,t,varargin{:});
         R=R+kron(eye(size(xy,1)),C);
      end        

    case 'spatialcov'
      
      %   stochModel (spatial)                          Attributes
      %   --------------------------------------------- ---------------
      %   spatialcov(rho=#[,q=#])

      % Spatial unit co-variance

      checkstochmodel(modelParameters,{'rho'},{'q'})
      if isfield(modelParameters,'q')
         S=covexp(xy,modelParameters.rho,modelParameters.q);
      else
         S=covexp(xy,modelParameters.rho);
      end

      if stage == 0
         if numel(stochModel) == 1
            % Just return the spatial covariance
            R=S;
         else
            error('spatialcov cannot be used as first model.')
         end
      elseif stage == 1               
         % Combine temporal and spatial covariance matrix
         R=kron(C,S);
         stage=2;
      elseif stage == 2
         error('spatialcov cannot be called twice.')
      elseif stage == 3
         % apply spatial covariance to full matrix, implement as loop, does
         % not make much sense
         error('spatialcov cannot be used in the full spatio-temporal covariance stage.')
      end
       
    case {'tudinsar4','tudinsar4rd','stdev','variance','covmatrix','precmatrix','tudpred1','sure','houtenbos'}

      %   stochModel (temporal)                         Attributes
      %   --------------------------------------------- ---------------
      %   houtenbos(s2t=#,pt=#)                         t

      %   stochModel (spatio-temporal)                  Attributes
      %   --------------------------------------------- ---------------
      %   tudinsar4(s20=#,s2t=#,s2s=#,Rt=#,Rs=#)        xy,t
      %   tudinsar4rd(s20=#,s2t=#,s2s=#,Rt=#,Rs=#,Ns=<>,Nt=<>,Ds=<>,Dt=<>)     
      %                                                 xy,t,pntAttrib,epochAttrib
      %   tudpred1(s20=#,s2t=#,s2s=#,Rt=#,Rs=#)         xy,t[,xy0,t0,mask0]
      %   sure(s2z=#,L=#,q=#,tref=#)                    xy,t[,xy0,t0,mask0]
      %   sure(s2z=#,L=#,q=#,tlag=#)                    xy,t,trefOrig[,xy0,t0,mask0]

      %   stochModel (data)                             Attributes
      %   --------------------------------------------- ---------------
      %   stdev()                                       stochData        
      %   covmatrix(format=...)                         stochData
      %   precmatrix(format=...)                        stochData

      if stage == 0
         % Initialize the covariance matrix
         R=0;
      elseif stage == 1
         % Expand temporal covariance matrix
         R=kron(C,eye(m));
      end
      stage=3;    
      
      switch modelType
          
        case 'tudinsar4'
            
          % TUD InSAR model 4 for non-reduced datasets (Ref: InSAR_FM_SM_v1.1.pdf, Eq. 4.4, p. 16)  
       
          % Check stochastic model parameters
          checkstochmodel(modelParameters,{'s20','s2s','s2t','Rs','Rt'})

          % Generate the spatial and temporal components      
          S=covexp(xy,1./modelParameters.Rs);
          T=covexp(t,1./modelParameters.Rt);

          % Covariance matrix is the sum of three components
          R= R + modelParameters.s20 * eye(m*n) +  ...
                 modelParameters.s2s * kron(eye(n),S) + ...
                 modelParameters.s2t * kron(T,eye(m));

        case 'tudinsar4rd'
        
          % TUD InSAR model 4 for reduced (averaged) datasets (InSAR_FM_SM_v1.1.pdf, Table 4.4, p. 21)  

          % This model requires extra point and epoch attributes given in pntAttrib and epochAttrib:
          %
          %    Par  Default                Description
          %    ---  ---------------------- ----------------------------------------------
          %    Ns   pntAttrib.cellCount    number of points spatial averaging cell (mi in Table 4.4)
          %    Ds   pntAttrib.cellSize     averaged distances among the points of each spatial cell (hi in Table 4.4)
          %    Nt   epochAttrib.cellCount  number of epochs/acquisitions in the averaging temporal intervals (np in Table 4.4) 
          %    Dt   epochAttrib.cellSize   averaged time differences between dates of the acquisitions in each averaging temporal interval (Tp in Table 4.4).
          %
          % We use the modified method with correlations to guarantee the co-variance matrix is positive definite. Note that this
          % method can also be implemented on the spatial and temporal  componentes, but the result is not the same, but it is something to
          % investigate.
        
          % Check stochastic model parameters
          
          checkstochmodel(modelParameters,{'s20','s2s','s2t','Rs','Rt'},{'Ns','Nt','Ds','Dt'})
          if ~isfield(modelParameters,'Ns')
             modelParameters.Ns='cellCount';
          end
          if ~isfield(modelParameters,'Ds')
             modelParameters.Ds='cellSize';
          end
          if ~isfield(modelParameters,'Nt')
             modelParameters.Nt='cellCount';
          end
          if ~isfield(modelParameters,'Dt')
             modelParameters.Dt='cellSize';
          end

          % Check inputs
        
          if numel(varargin) < 2
              error('The tudinsar4rd model expects two additional input arguments: pntAttrib and epochAttrib.')
          end
          pntAttrib=varargin{1};
          epochAttrib=varargin{2};
          if ~isfield(pntAttrib,modelParameters.Ns) || ~isfield(pntAttrib,modelParameters.Ds)
              error('Missing field in pntAttrib for tudins4rd model.')
          end
          if ~isfield(epochAttrib,modelParameters.Nt) || ~isfield(epochAttrib,modelParameters.Dt')
              error('Missing field in epochAttrib for tudins4rd model.')
          end

          mi=pntAttrib.(modelParameters.Ns);
          hi=pntAttrib.(modelParameters.Ds);
          np=epochAttrib.(modelParameters.Nt);
          tp=epochAttrib.(modelParameters.Dt);
          
          mi=mi(pntMask);
          hi=hi(pntMask);
          np=np(epochMask);
          tp=tp(epochMask);
          
          hi(isnan(hi))=Inf;     % Patch for NaN ...
          tp(isnan(tp))=Inf;

          % Create diagonal components
    
          % Diagonal components according to report and applied to full
          % matrix. Could also create positive definite temporal and spatial
          % components. The sum is also guaranteed to be positive definite.
        
          f=repmat(1./mi(:),[1 n]).*repmat(1./np(:)',[m 1]);
          s2d0 = modelParameters.s20 .* f;
          s2ds = modelParameters.s2s .* f .* repmat( 1 + (mi(:)-1).*exp( -hi(:)./modelParameters.Rs ) , [1 n] );
          s2dt = modelParameters.s2t .* f .* repmat( 1 + (np(:)'-1).*exp( -tp(:)'./modelParameters.Rt ) , [m 1 ] );
        
          % Generate the spatial and temporal correlation components

          S=covexp(xy,1./modelParameters.Rs);
          T=covexp(t,1./modelParameters.Rt);

          % Covariance matrix

          %s2d=sqrt(s2d0(:) + s2ds(:) + s2dt(:));
          %R= R + (s2d*s2d') .* ( kron(eye(n),S) + kron(T,eye(m))  ) ;
          
          s2ds = sqrt(s2ds(:));
          s2dt = sqrt(s2dt(:));
          R = R + diag(s2d0(:)) + (s2ds*s2ds').*kron(eye(n),S) + (s2dt*s2dt').*kron(T,eye(m));

          R=double(R);    % Patch 
          
        case 'tudpred1'
            
          % TUD Prediction model 1 (similar as tudinsar4, but now with full spatio-temporal ...
          % correlation, and calculation of prediction matrix, to enable prediction)  
       
          % Check stochastic model parameters
          checkstochmodel(modelParameters,{'s20','s2s','s2t','Rs','Rt'})
 
          if numel(varargin) == 3 % create prediction matrix
            xy0 = varargin{1};
            t0 = varargin{2};
            mask0 = varargin{3};
            t0 = t0(:);

            xy0 = xy0(mask0{1},:);
            t0 = t0(mask0{2});
            m0 = size(xy0,1);
            n0 = numel(t0);

            xy = kron(ones(n,1),xy);
            xy0 = kron(ones(n0,1),xy0);
            t = kron(t,ones(m,1));
            t0 = kron(t0,ones(m0,1));

            % Generate the spatial and temporal components      
            S = predexp(xy0,xy,1./modelParameters.Rs);
            T = predexp(t0,t,1./modelParameters.Rt);

            % Prediction matrix is the sum of two components (no observation noise)
            R = modelParameters.s2s * S + ...
                modelParameters.s2t * T;

          else % create covariance matrix

            % Generate the spatial and temporal components      
            S=covexp(xy,1./modelParameters.Rs);
            T=covexp(t,1./modelParameters.Rt);

            % Covariance matrix is the sum of two components (no observation ...
            % noise, full covariance)
            R = R + modelParameters.s2s * kron(ones(n),S) + ...
                    modelParameters.s2t * kron(T,ones(m));
            %
            % Covariance matrix is the sum of three components (full covariance)
            %R = R + modelParameters.s20 * eye(m*n) +  ...
            %        modelParameters.s2s * kron(ones(n),S) + ...
            %        modelParameters.s2t * kron(T,ones(m));

          end

        case {'sure'}

          % SURE prediction model   
       
          % Check stochastic model parameters
          checkstochmodel(modelParameters,{'s2z','L','q'},{'tref','tlag'})
          if isfield(modelParameters,'tref')
              tref = modelParameters.tref;   
              hasOrig=0;
          elseif isfield(modelParameters,'tlag')
              if numel(varargin) < 1
                 error('You need to specify the reference time in case of the SURE model.')
              end
              trefOrig=varargin{1};
              hasOrig=1;
              tref = trefOrig - modelParameters.tlag;
          else
              %tref = min(t);
              error('You need to specify tlag or tref in case of the SURE model.')
          end
 
          if numel(varargin) == 3 || numel(varargin) == 4     % create prediction matrix
            xy0 = varargin{hasOrig+1};
            t0 = varargin{hasOrig+2};
            mask0 = varargin{hasOrig+3};
            t0 = t0(:);

            xy0 = xy0(mask0{1},:);
            t0 = t0(mask0{2});
            m0 = size(xy0,1);
            n0 = numel(t0);

            xy = kron(ones(n,1),xy);
            xy0 = kron(ones(n0,1),xy0);
            t = kron(t,ones(m,1));
            t0 = kron(t0,ones(m0,1));

            % Generate the spatial and temporal components
            S = predexp(xy0,xy,1./modelParameters.L,2);
            T1 = kron(abs(t0-tref).^(2*modelParameters.q),ones(1,n*m));
            T2 = kron(abs(t'-tref).^(2*modelParameters.q),ones(n0*m0,1));
            T = T1 - pdist2(t0,t).^(2*modelParameters.q) + T2;

            % Prediction matrix
            R = 0.5*modelParameters.s2z * T .* S;

          else % create covariance matrix

            % Generate the spatial and temporal components
            S = covexp(xy,1./modelParameters.L,2);
            T1 = kron(abs(t-tref).^(2*modelParameters.q),ones(1,n));
            T = T1 - squareform(pdist(t).^(2*modelParameters.q)) + T1';

            % Covariance matrix is the sum of two components (no observation ...
            % noise, full covariance)
            R = 0.5*modelParameters.s2z * kron(T,ones(m)) .* kron(ones(n),S);
            
          end

        case {'houtenbos'}

          % Houtenbos idealization model (Houtenbos and Kenselaar, 2001)
       
          % Check stochastic model parameters
          checkstochmodel(modelParameters,{'s2t','pt'})

          T = 2*modelParameters.s2t*squareform(pdist(t).^(2*modelParameters.pt));
          R = kron(T,eye(m));

        case {'stdev'}

          % Check stochastic model parameters
          
          checkstochmodel(modelParameters,{},{'sigma' 'size' 'odata'})
          if ~isfield(modelParameters,'odata')
             modelParameters.odata=0;
          end
          if ~isfield(modelParameters,'size')
             if nargin > 3
                modelParameters.size=[m n];
             else
                error('could not determine missing size in stochmodel description ')
             end
          end
          modelParameters.size_orig=modelParameters.size;
          modelParameters.size=prod(modelParameters.size);
   
          % Diagonal co-variance matrix (standard deviation)

          if  ( isempty(pntMask) || isempty(epochMask) || (  all(pntMask) && all(epochMask) ) ) 
             if ~isfield(modelParameters,'sigma')
                R = R + diag(stochData(modelParameters.odata+1:modelParameters.odata+modelParameters.size).^2);
             else
                R = R + modelParameters.sigma^2 * eye(modelParameters.size);
             end
          else
             if ~isfield(modelParameters,'sigma')
                R1 = stochData(modelParameters.odata+1:modelParameters.odata+modelParameters.size).^2;
                R1 = reshape(R1,[m n]);
                R1 = R1(pntMask,epochMask);
                R = R + diag(R1);
             else
                R = R + modelParameters.sigma^2 * eye(sum(pntMask)*sum(epochMask));
             end
          end
        
        case {'covmatrix','precmatrix'}

          % Check stochastic model parameters
          if strcmpi(modelParameters.format,'index')
              % Indexed -> check stochastic model parameters  (lidx, oidx, odata are mandatory) 
              checkstochmodel(modelParameters,{'format' 'lidx' 'oidx' 'odata'},{'size' 'sreg'})
          elseif strcmpi(modelParameters.format,'blkdiag') 
              % Blkdiag -> nd is required now
              checkstochmodel(modelParameters,{'format' 'nd'},{'sreg' })
              if ~isfield(modelParameters,'sreg')
                  % Default value for regularization parameter
                  modelParameters.sreg=sqrt(1e-9);
              end
          else
              % Non idexed ->  Check stochastic model parameters  (size, odata are optional) 
              checkstochmodel(modelParameters,{'format'},{'size' 'odata' 'nd' 'sreg'})
              if ~isfield(modelParameters,'odata')
                  modelParameters.odata=0;
              end
          end
          if ~isfield(modelParameters,'size')
             if nargin > 3
                modelParameters.size=[m n];
             else
                error('could not determine missing size in stochmodel description ')
             end
          end
          modelParameters.size_orig=modelParameters.size;
          modelParameters.size=prod(modelParameters.size_orig);

          switch modelParameters.format

            case {'full'}

              if modelParameters.odata+modelParameters.size^2 > length(stochData)
                  error('stochData has incorrect length')
              end
              R1 = stochData(modelParameters.odata+1:modelParameters.odata+modelParameters.size^2);
              R1 = reshape(R1,[modelParameters.size modelParameters.size]);

            case {'index'}

              R1=zeros(modelParameters.size,modelParameters.size);

              % Insert indexed covariance data 
              if modelParameters.odata+modelParameters.lidx^2 > length(stochData)
                  error('stochData has incorrect length')
              end
              idx=stochData(modelParameters.oidx+1:modelParameters.oidx+modelParameters.lidx);
              R1(idx,idx)=reshape(stochData(modelParameters.odata+1:modelParameters.odata+modelParameters.lidx^2), ...
                     [modelParameters.lidx modelParameters.lidx]);

            case {'upper', 'lower'}

              % Upper (or lower) triangular matrix stored in a single vector

              % stochData is a 1-by-(M*(M+1)/2) vector corresponding to the M*(M+1)/2 
              % unique pairs of covariances in the M-by-M co-variance matrix.
              %
              % stochData is arranged in the order of ((1,1),(2,1),(3,1),..., (M,1),
              % (2,2),...(M,2),.....(M,M)), i.e. the lower left triangle of the full
              % M-by-M co-variance matrix in column order.  
              %
              % To get the co-variance between the Ith and Jth observations (I < J), 
              % use the formula D((I-1)*(M-I/2)+J-I), or use the code below to
              % return a square and symmetric co-variance matrix.
    
              len2=length(stochData(:));
              len1=floor(sqrt(2*len2));
              if len1*(len1+1)/2 ~= len2
                 error('Hey, this should never happen!!');
              end
              R1 = zeros(len1);
              R1(tril(true(len1),0)) = stochData(:);
              R1 = R1 + triu(R1',1);
              
            case {'blkdiag'}

              % Block diagonal matrix

              % stochData is a M-by-M-by-N matrix, with N diagonal blocks,
              % of each a M-by-M matrix.
                            
              % Find out the length and reshape stochData
              nd=modelParameters.nd;
              len2=length(stochData(:));
              len1=floor(sqrt(len2/nd));
              if nd*len1^2 ~= len2
                 error('Hey, this should never happen!!');
              end
              stochData=reshape(stochData,[len1 len1 nd]);

              % prepare epochMask and pntMask
              if isempty(epochMask), epochMask=true(1,nd); end
              if isempty(pntMask), pntMask=true(len1,1); end
              len1mask=sum(pntMask);
              ndmask=sum(epochMask);
              
              % Copy the diagonal blocks into the covariance matrix 
              R1 = zeros(len1mask*ndmask);
              k2=0;
              for kk=1:nd
                 if ~epochMask(kk), continue; end
                 k1=k2+1;
                 k2=k2+len1mask;
                 if  modelParameters.sreg <= 0
                     % Blocks are full rank, no regularization required
                     R1(k1:k2,k1:k2)=stochData(pntMask,pntMask,kk);
                 else
                     % Co-variance matrix for each block is rank defect, regularize 
                     % it (doesn't harm the estimate, even if full rank)
                     R1(k1:k2,k1:k2)=stochData(pntMask,pntMask,kk)+ones(len1mask,len1mask)*modelParameters.sreg.^2;
                 end
              end

            case {'symdiags'}

              % Matrix with cyclic diagonals [B,d]
              
              %   A= SYMDIAGS(B) reconstucts the symmetric m-by-m matrix A from the 
              %   m-by-m/2+1 matrix of cyclic diagonals B.  The output matrix is
              %   a full matrix with only the upper triangle.
              %        
              %   A= SYMDIAGS(B,d) reconstucts the symmetric m-by-m matrix A from the p 
              %   cyclic diagonals in the m-by-p matrix B, with the cyclic diagonals 
              %   specified in d. The values of d must be in the range [0:m/2].
              %   The output is sparse matrix with only the upper triangle with
              %   diagonals d(1:p) and m-d(1:p).
        
              if isfield(modelParameters,'nd')              
                 % Return as sparse matrix?
                 R1 = symdiags(stochData,modelParameters.nd);
              else
                 R1 = symdiags(stochData);
              end
                  
            otherwise
                  
              error('Illegal format for covariance matrix.')
              
          end

          % Blank out pntMask and epochMask data (optionally)
          if  ( isempty(pntMask) || isempty(epochMask) || (  all(pntMask) && all(epochMask) ) || strcmpi(modelParameters.format,'blkdiag') )
              R = R + R1;
          else
              mask0=true(numel(pntMask),numel(epochMask));
              mask0(~pntMask(:),:)=false;
              mask0(:,~epochMask(:))=false;
              R = R + R1(mask0(:),mask0(:));
          end
          
          % Logic for precision matrix is not yet implemented
          if strcmp(modelType,'precmatrix')
            error('Precision matrix not yet implemented.')
          end
          
      end
      
    otherwise
          
      error('Unknown stochastic model type.')
  
  end

end

% If stage is one, return only temporal component

if stage == 1
   R=C;
end

% Apply observation mask

if ~isempty(mask) && stage > 1
  if size(R,1)==size(R,2) % Covariance matrix
    R=R(mask(:),mask(:));
  else                          % Prediction matrix
    R=R(:,mask(:));
  end
end

end


%% Internal functions

% % Code examples to test parsing of stochastic model ...
% 
% [a,b]=parsestochmodel('tudinsar4(Rt=34.4, s2s=8, s2t=2, s20=2, Rs=0)')
% checkstochmodel(b,{'s20','s2s','s2t','Rs','Rt'},{'p'})
% [a,b]=parsestochmodel('tudinsar4(Rt=34.4, s2s=8, s2t=2, s20=2, Rs=0, p=Hans)')
% checkstochmodel(b,{'s20','s2s','s2t','Rs','Rt'},{'p'})
% [a,b]=parsestochmodel('tudinsar4(Rt=34.4, s2s=8, s2t=2, s20=2, Rs=0, a=Hans)')
% checkstochmodel(b,{'s20','s2s','s2t','Rs','Rt'},{'p'})
% [a,b]=parsestochmodel('tudinsar4(Rt=34.4, s20=2, Rs=0)')
% checkstochmodel(b,{'s20','s2s','s2t','Rs','Rt'},{'p'})

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
% Modified: 30 January 2024 by Hans van der Marel
%              - do not remove blanks from parameter values

%str=strrep(str,' ','');

ii=regexp(str,'[()]');
name=str(1:ii(1)-1);
pargs=strsplit(str(ii(1)+1:ii(end)-1),',');

params=[];
for k=1:numel(pargs)
   if ~isempty(pargs{k})
      p=strsplit(pargs{k},'=');
      p{1}=strrep(p{1},' ','');
      value=str2num(p{2});
      if isempty(value)
         value=p{2};
      end
      params.(p{1})=value;
   end
end


end

function checkstochmodel(modelParameters,expected,optional)
%CHECKSTOCHMODEL  Check stochastic model parameters.
%   CHECKSTOCHMODEL(MODELPARAMETERS,EXPECTED,OPTIONAL) checks the stochastic
%   model parameters in the stucture MODELPARAMETERS against the expected
%   parameter fields in the cell array EXPECTED and optional fields in
%   the cell array OPTIONAL.
%
%   Examples
%
%      [stochModel,stochParameters]=parsestochmodel('tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)')
%      checkstochmodel(stochParameters,{'s20','s2s','s2t','Rs','Rt'})
%
%      [stochModel,stochParameters]=parsestochmodel('PL(alpha=-0.8281,sigma=0.76,fs=1,ndays=1)')
%      checkstochmodel(stochParameters,{'alpha','sigma'},{'fs','ndays'})
%
%   (c) Hans van der Marel, Delft University of Technology, 2020

% Created:  3 September 2020 by Hans van der Marel
% Modified:

if nargin < 2
    error('At least two arguments expected.')
end
if nargin < 3
   optional=[];
end

parameterNames= fieldnames(modelParameters);
[~,~,ib]=setxor(parameterNames,expected);
if ~isempty(ib)
   error(['Not all required parameters (' strjoin(expected(ib)) ') have been specified in the stochastic model, quiting.'])
end
[~,ia]=setxor(parameterNames,[ expected optional ]);
if ~isempty(ia)
   warning([ 'The stochastic model contains unknown parameters (' strjoin(parameterNames(ia)) '), these will be ignored.'])
end

end

function C=tempcov(modelType,modelParameters,t,varargin)
%TEMPCOV   Covariance matrix for the temporal co-variance.
%  C=TEMPCOV(MODELTYPE,MODELPARAMETERS,T) computes the co-variance matrix C 
%  for a powerlaw (PL), Fractionally Integrated Generalised Gauss Markov 
%  (iGGM), First-Order Gauss Markov (FOGM) exponential decay model, or
%  white noise (WN) model. MODELTYPE is either PL, iGGM, FOGM or WN, 
%  MODELPARAMETERS is a structure with the model parameters, and T is a 
%  vector with decimal years.
%
%  C=TEMPCOV(...,EPOCHATTRIB) provides the optional epoch attributes, in
%  case the model uses the 'ndays' parameter to specify a field in the
%  epoch attributes.
%
%  MODELTYPE and MODELPARAMETERS are output of parsestochmodel, and the
%  parameters should be checked by checkstochmodel before the call to
%  this function.
%
%  See also PARSESTOCHMODEL, CHECKSTOCHMODEL, iGGM, HOSKING and COVEXP.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:  10 Sep 2020 by Hans van der Marel
% Modified:

% Check input arguments

if nargin < 3
    error('This function needs at least three input arguments.')
end

% Get the sample rate (number of samples per day)

if isfield(modelParameters,'fs')
    fs=modelParameters.fs;
else
    fs=1;
end

% Get the (moving) averaging window length (in days) 

if isfield(modelParameters,'ndays')
   ndays=modelParameters.ndays;
   if ~isnumeric(ndays)
       % read ndays from attribute field
       if numel(varargin) < 1
          error('Epoch attribute required, but missing in argument list.')
       end
       epochAttrib=varargin{1};
       if ~isfield(epochAttrib,ndays)
          error(['Expecting field ' ndays ' in epoch attributes, but missing.'])
       end
       ndays=epochAttrib.(ndays);
   end
else
   ndays=1;
end

% Set logical noaveraging to true when it is straight sampling without averaging

noaveraging = isscalar(ndays) && ( abs(ndays*fs - 1) < 1e-6 );


% Get the optional lead-up time (in days)  

if isfield(modelParameters,'leadup')
   leadup=modelParameters.leadup;
else
   leadup=0;
end

% Convert time from decimal years to MJD

mjd=51544.5+(t(:)-2000)*365.25;

% If ndays scalar, convert into array

nepochs=size(t(:),1);
if isscalar(ndays)
    ndays=ndays*ones(nepochs,1);
end

% Create a full continous integer sequence [0:nfull-1] from leadup days 
% before first mjd to last mjd, with sample period dt=1/fs. Only required
% for PL, iGGM and FOGM models (not used for WN, but computed anyhow)

mjd_start = floor(min(mjd)-max(ndays)/2-leadup); 
mjd_end   = ceil(max(mjd)+max(ndays)/2);       

nfull=(mjd_end-mjd_start)*fs+1;

% Convert time into sample number in the full continuous integer sequence (mc),
% and determine for each epoch the first (m1) and last (m2) sample number.
% Only required for PL, iGGM and FOGM models (not used for WN, but computed anyhow)

mc = round( (mjd-mjd_start)*fs ) + 1;
m1 = mc - floor((ndays*fs-1)/2);
m2 = mc + floor(ndays*fs/2);

% Compute the co-variance matrix

switch modelType

    case {'iGGM','PL'}
        
       %   PL(alpha=#,sigma=#[,fs=#,ndays=#])            t[,epochAttrib]
       %   iGGM(alphaR=#,alphaL=#,gamma=#,sigma=#        t[,epochAttrib]
       %                            [,fs=#,ndays=#])

       % Generate differencing/integrating vector
       
       switch modelType
          case 'iGGM'
             checkstochmodel(modelParameters,{'alphaR','alphaL','gamma','sigma'},{'fs','ndays','leadup'});
             h = iGGM(nfull,modelParameters.alphaR,modelParameters.alphaL,modelParameters.gamma) * modelParameters.sigma;
          case 'PL'
             checkstochmodel(modelParameters,{'alpha','sigma'},{'fs','ndays','leadup'});
	         h = hosking(nfull,modelParameters.alpha) * modelParameters.sigma;
       end
       
       if noaveraging
           % No averaging needed, use Toeplitz matrix in computation
           S = tril(toeplitz(h));
  		   C = S(mc,:)*S(mc,:)';
           % Other method would be to multiply the columns of the Toeplitz
           % matrix with the t(i)-t(i-1), using nepochs instead of nfull.
       else
           % Averaging over ndays*fs samples
           np = 2^nextpow2(nfull);
           c = [h(1);zeros(nfull-1,1);zeros(2*np-(2*nfull-1),1);flipud(h(2:end))];
           fftc=fft(c);
           S = zeros(2*np,nepochs);
           for j = 1:nepochs
              ndaysj=m2(j)-m1(j)+1;
              T = zeros(2*np,1);
              T(m1(j):m2(j)) = ones(ndaysj,1) ./ ndaysj;
              S(:,j) = ifft(fftc.*fft(T));
           end
           S = S(1:nfull,:)';
           C = S*S';
       end
      
    case 'FOGM'
        
       %   FOGM(rho=#,sigma=#[,fs=#,ndays=#])            t[,epochAttrib]
       
       checkstochmodel(modelParameters,{'rho','sigma'},{'fs','ndays'});

       if noaveraging
          C=covexp(mc,modelParameters.rho)*modelParameters.sigma^2; 
       else
          C=covexp((0:nfull-1)',modelParameters.rho) * modelParameters.sigma^2;
          T=zeros(nepochs,nfull);
          for j = 1:nepochs
             ndaysj=m2(j)-m1(j)+1;
             T(j,m1(j):m2(j)) = ones(ndaysj,1) ./ ndaysj;
          end
          C=T*C*T';
       end
        
    case 'WN'
        
       %   WN(sigma=#[,fs=#,ndays=#])                    t[,epochAttrib]
       
       checkstochmodel(modelParameters,{'sigma'},{'fs','ndays'});
       ndays(ndays==0)=+Inf;  % If ndays is zero, no white noise is added! In this way we can play with the "setup" noise
       C= diag(modelParameters.sigma^2 ./ (fs*ndays) );
         
    otherwise
        
       error('Invalid modelType.')

end

end


function C=covexp(x,rho,q)
%COVEXP  Covariance matrix for with exponential co-variance function.
%  C=COVEXP(X) computes the co-variance matrix C for the exponential 
%  co-variance function with co-variance C(i,j)=exp(-d(i,j)). The distance
%  d(i,j) is |x(i)-x(j)| when X is a vector or the euclidian distance
%  sqrt((x(i)-x(j))^2+(y(i)-y(j))^2) when X is a matrix with two columns.
%
%  C=COVEXP(X,RHO) computes the co-variance matrix C for the exponential 
%  co-variance function with co-variance C(i,j)=exp(-d(i,j)*RHO). 
%
%  C=COVEXP(X,RHO,Q) computes the co-variance matrix C for the 
%  co-variance function with co-variance C(i,j)=exp((-d(i,j)*RHO)^Q). If 
%  Q=2, the co-variance function is Gaussian
%
%  See also PRECEXP and PREDEXP.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:   3 Sep 2020 by Hans van der Marel
% Modified: 20 Nov 2020 by Freek van Leijden and Hans van der Marel
%              - added third argument q to enable also Gaussian co-variance
%                function (q=2)

% Check input arguments

if nargin < 1, error('This function expects at least one input argument.'); end
if size(x,1)==1
   C = 1;
   return
end

% Compute vector d with distance between each pair of observations in X.

% We use pdist. This is an efficient implementation using a mex function.
% The ouput d is a 1-by-(M*(M-1)/2) row vector, corresponding to the M*(M-1)/2 pairs 
% of observations in X.

d=pdist(x);
if nargin > 1
  d=d.*rho;
end
if nargin > 2
  d=d.^q;
end

% Co-variance matrix C  (similar to C=exp(-squareform(d)) );

[~,n]=size(d);
m = ceil(sqrt(2*n));
if m*(m-1)/2 ~= n
   error('Hey, this should never happen!!');
end

C = zeros(m);
C(tril(true(m),-1)) = exp(-d);
C = eye(m) + C + C';

end

function C=predexp(x1,x2,rho,q)
%PREDEXP  Prediction matrix for with exponential co-variance function.
%  C=PREDEXP(X1,X2) computes the prediction matrix C for the exponential 
%  co-variance function with co-variance C(i,j)=exp(-d(i,j)). The distance
%  d(i,j) is |x2(i)-x1(j)| when X1 and X2 are vectors, or, the euclidian 
%  distance sqrt((x2(i)-x1(j))^2+(y2(i)-y1(j))^2) when X1 and X2 are
%  matrices with two columns. X1 and X2 can have different lengths.
%
%  C=PREDEXP(X1,X2,RHO) computes the prediction matrix C for the exponential 
%  co-variance function with co-variance C(i,j)=exp(-d(i,j)*RHO). 
%
%  C=PREDEXP(X1,X2,RHO,Q) computes the co-variance matrix C for the 
%  co-variance function with co-variance C(i,j)=exp((-d(i,j)*RHO)^Q). If 
%  Q=2, the co-variance function is Gaussian
%
%  See also COVEXP and PRECEXP.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:   7 Nov 2020 by Freek van Leijen, based on COVEXP
% Modified: 22 Dec 2020 by Hans van der Marel
%              - created the help text 

% Check input arguments

if nargin < 2, error('This function expects at least two input argument.'); end

% Compute vector d with distance between each pair of observations in X. 

% We use pdist2. This is an efficient implementation using a mex function.
% The ouput d is a 1-by-(M*(M-1)/2) row vector, corresponding to the M*(M-1)/2 pairs 
% of observations in X.

d=pdist2(x1,x2);

if nargin > 2
  d=d.*rho;
end

if nargin > 3
  d=d.^q;
end

% Co-variance matrix C
C = exp(-d);

end

function [p0,p1]=precexp(x)
%PRECEXP   Tridiagonal precision matrix for exponential co-variance function.
%  [P0,P1]=PRECEXP(X) computes the main P0 and first diagonal P1 of the 
%  precision matrix for the exponential co-variance function with
%  co-variance C(i,j)=exp(-|x(i)-x(j)|)). The precision matrix P is the 
%  inverse of the co-variance matrix C.  The input coordinates must be
%  sorted with x(i-1) <= x(i) and scaled with correlation length (X are
%  normalized unitless coordinates). X does not have to be equidistant.
%
%  The precision matrix P is tri-diagonal and symmetric matrix with 
%  P=diag(P0)+diag(P1,1)+diag(P1,-1). P0 is an array of length
%  N, P1 an array of length N-1, with N the length of X. To create a
%  sparse precision matrix P use spdiags, with
%
%     P = spdiags([ [p1;0] p0 [0;p1]],[-1 0 1],n,n)
%
%  Note that when the distance X(i)-X(i+1) between two coordinates becomes  
%  very small (large correlation), the corresponding elements in P0 and P1 
%  will become very large, indicating a strong dependency. Since the computation
%  in this function is analytic, the result is numerically more stable 
%  than obtained by inverting the covariance matrix numerically.
%
%  See also spdiags.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:  29 Aug 2020 by Hans van der Marel
% Modified:

x=x(:);
dx=diff(x);

pdx=exp(dx);
qdx=exp(-dx);

pqdx=pdx-qdx;
p2dx=exp(x(3:end)-x(1:end-2));
q2dx=exp(x(1:end-2)-x(3:end));

p0 = [ pdx(1)./pqdx(1) ; ...
      (p2dx-q2dx)./(pqdx(1:end-1).*pqdx(2:end)) ; ...
       pdx(end)./pqdx(end) ];
p1 = -1./pqdx;

end

function [h]=hosking(length,spectral_index)
%HOSKING   Differencing/integrating vector using Hosking algorithm
%   H=HOSKING(LENGTH,SPECTRAL_INDEX) Generate a differencing/integrating 
%   vector H of length LENGTH and with spectral_index SPECTRAL_INDEX.
%   The SPECTRAL_INDEX can be of any order including fractional.  Thise are 
%   the discrete-time versions of fractional differentiation and integration. 
%
%   See Hosking, R. M., 1981, Biometrika, 68, pp. 165-176

% Created:  28 May 1998 by Hadley Johnson, hjohnson@ucsd.edu.
% Modified: 13 Jul 2007 by Simon Williams, sdwil@noc.ac.uk
%            8 Sep 2020 by Hans van der Marel
%              - changed help text
%              - minor syntax change (precedence by () instead of [])

degree = spectral_index / 2;
%h=cumprod([1;[0:1:length-2]'-degree]./[1;[1:1:length-1]']);
h=cumprod([1;(0:1:length-2)'-degree]./[1;(1:1:length-1)']);

end

function [h] = iGGM(N,alphaR,alphaL,gamma)
%iGGM   Fractionally Integrated Generalised Gauss Markov (FIGGM) model.
%   H=iGGM(iGGM(N,alphaR,alphaL,gamma) generates a differencing/integrating
%   vector H of length N for the FIGGM model with spectral index ALPHAR
%   to right of cross-over frequency and spectral index ALPHAL to the left
%   of the cross-over frequency. The factor GAMMA (0 -> 1) defines the 
%   cross-over for AR(1).  As GAMMA tend to 1 the cross-over frequency 
%   decreases in frequency. 
%
%   A Fractionally Integrated Generalized Gauss Markov (FIGGM) model has 
%   two regions separated at some frequency, fc, which is related to the 
%   AR(1) parameter, GAMMA. To the left of fc, the spectrum has a slope in 
%   log-log space equal to ALPHAL and a slope of ALPHAR to the right of fc. 
%   All of the other models can be considered special cases of FIGGM. 
%
%   Model           ALPHAL     ALPHAR     GAMMA
%   -------------   -------    -------    -------
%   AR(1)/FOGM        0          -2       GAMMA
%   ARFIMA/ARFI     ALPHAL     ALPHAL-2   GAMMA
%   POWER-LAW       ALPHA      ALPHA        0       (ALPHA->ALPHAL)
%   POWER-LAW       ALPHA      ALPHA        1       (ALPHA->APLHAR)
%   GGM               0        ALPHAR     GAMMA
%   FIGGM           ALPHAL     ALPHAR     GAMMA
%
%   For instance AR(1) (also called First-Order Gauss Markov (FOGM) noise) 
%   simply has one parameter, GAMMA, which describes where the break point 
%   is; ALPHAL is 0 and ALPHAR is -2. Power-law noise which simply has one 
%   slope can arise from a number of combinations of the three parameters. 
%
%   Ref: https://www.psmsl.org/products/trends/methods.php
%
%   The autocovariance function for AR(1) decays with a decay time (also 
%   called time constant) TAU=-1/LOG(GAMMA). The cross-over frequency is
%   associated to 1/TAU = -LOG(GAMMA). [THIS IS TO BE CONFIRMED]
%
%   See also hosking.

% Created:  June 2015 by Simon Williams
% Modified: 8 Sep 2020 by Hans van der Marel
%              - added help text with all the details on cross-over

hl = hosking(N,alphaL);
alpha1 = alphaR - alphaL;
hr = hosking(N,alpha1);

g = gamma.^(0:N-1)';

h = fftfilt(hl,g.*hr);

end

