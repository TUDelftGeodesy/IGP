function [v,t0,offsets,refseries,omt,emat,dmat,qxhat,qehat,qyobs,Amatrices]=stmvelocity(stm,varargin)
%STMVELOCITY  Velocity estimation for space time matrix.
%  V=STMVELOCITY(STM) estimates the velocity V of the points in the space-time
%  matrix structure STM. The function return a plain matrix with dimensions
%  [numPoints x numObsTypes] with points and observation type in the same
%  order as the input space-time matrix STM. 
% 
%  The output V is filled with NaN's when no velocity could be determined 
%  because of insufficient data.
%
%  [V,T0,OFFSET,REFSERIES,OMT]=STMVELOCITY(STM) also returns the reference 
%  epoch TO (decimal years) with the estimated OFFSETS, the estimated
%  reference time series REFSERIES, and overall model test value OMT. 
%  OFFSET is a [numPoints x numObsTypes] matrix (same dimensions as V), 
%  REFSERIES is a matrix with dimensions [numEpochs x numObsTypes] and
%  OMT a vector with dimension [numObsTypes].
%
%  The reference time series REFSERIES is only estimated if information
%  content in the input dataset corresponds to 'relative displacements'.
%  This is typical for InSAR. If the information content in the input 
%  dataset corresponds to 'displacements' in a consistent reference frame 
%  over time, then no reference time-series needs to be estimated. This is
%  typical for some GNSS datasets. In this case REFSERIES will be filled 
%  with zeros. The particular mode is set with the boolean option 'relative').
% 
%  The matrices OFFSETS and REFSERIES are NaN filled when there was not 
%  enough data to estimate the parameters. REFSERIES is also filled with 
%  zeros when no reference time series was estimated (option 'relative' is
%  false).
% 
%  [V,T0,OFFSET,REFSERIES,OMT,EMAT,DMAT]=STMVELOCITY(STM) returns also the
%  [numPoints x numEpochs x numObsTypes] matrices EMAT and DMAT. EMAT is
%  a matrix with the residuals and NaNs in places without observations and
%  unused observations. DMAT contains the model displacements computed from
%  the estimated velocities.  
%
%  [V,T0,OFFSET,REFSERIES,OMT,EMAT,DMAT,QXHAT,QEHAT,QY]=STMVELOCITY(STM) gives 
%  also cell arrays QXHAT, QEHAT and QY. QXHAT contains the covariance matrix 
%  of the estimated parameters for each of the observation types. The order of 
%  the parameters in each covariance matrix is the same as in the outputs 
%  V, OFFSETS, REFSERIES (if applicable), and in this order 
%  [ numPoints numPoints [numEpochs]]. The matrix is not NaN filled, so in 
%  all likelyhood, it will be smaller dimensions than the corresponding 
%  parameter estimates. QEHAT and QY are respectively the covariance matrix of the
%  residuals (EMAT) and the underlying observations (only those that have been 
%  actually used in the estimation. These are again cells arrays, with in
%  each cell a covariance matrix. 
%
%  [V,T0,OFFSET,REFSERIES,OMT,EMAT,DMAT,QXHAT,QEHAT,QY,AMATRICES]=STMVELOCITY(STM)
%  provides additioal output of the design matrices AMATRICES. 
%
%  [...]=STMVELOCITY(STM,'option',value,...) allows to specify options
%  for the velocity estimation
%
%    'trendmodel'         trend model (default and only option 'linear')
%    't0'                 reference epoch {'mean' , 'best', <dyear>} (default 'best')
%
%    'refsystem'          Choice of reference system (default 'min-refseries'):
%                         - 'inherit'        
%                              Same as 'relative'=false 
%                         - 'minimum-norm'   
%                              Minimum norm solution (use only when offsets are close to zero)
%                         - 'min-refseries'  
%                              Minimize the reference series corrections (default)
%                         - pntRefMask, pntRefNumbers, refROI 
%                              Implicitly sets 'refsystem' to 'min-velocities'. This sets
%                              the mean velocity of the points specified by prnRefMask, 
%                              pntRefNumbers or refROI to the value of 'cvelocity'.
%
%    'relative'           true if displacements in STM are relative (default true),
%                         when 'relative' is set to false the 'refsystem', 'cvelocity'
%                         and 'cweight' options become meaningless. Same as
%                         'refsystem'='inherit'.
%
%    'cvelocity'          Velocity for point velocity constraints (Default 0)
%    'cweight'            Weight for pseudo-observation constraints (Default 1)
%
%    'ROI'                Region of interest, as [latmin lonmin ; latmax lonmax] 
%                         bounding box, or lat/lon polygon (Default [], is all)
%    'POI'                Period of interest [ dYearStart dYearEnd ] (Default [-Inf +Inf])
%    'pntMask'            Point mask (default [], is all), will be combined with ROI
%    'epochMask'          Epoch mask (default [], is all), will be combined with POI
%
%    'minObsDisplPoint'   Minimum number observed displacements per point, 
%                         if below, point is removed (default 2)
%    'minObsDisplEpoch'   Minimum number observed displacements per epoch, 
%                         if below, epoch is removed (default is 2 for relative, 
%                         1 otherwise)
%
%    'verbose'            Verbose level (default 0)
%    'check'              If true, check some of the computions (default is false)
%    'ignoreStochModel'   If true, ignore stochastic model from stm, use 
%                         unit matrix instead (default false) 
%    'pntCrdType'         Point coordinate type ['deg/m','km/m'] (default 'deg/m')
%
%  In case of relative displacements (option 'relative' is true) the 
%  underlying model is rank-defect. The rank defect is 2. The rank defects
%  are between the
%
%    - offsets and reference series: the corresponding null-space is
%
%         N1 = [ zeros(1,numPoints) ones(1,numPoints) ones(1,numEpochs) ] 
%
%    - velocities and reference series: the corresponding null-space is
%         
%         N2 = [ ones(1,numPoints) zeros(1,numPoints) t-t0 ] 
% 
%      with t-t0 a [1 x numEpoch] vector.
%
%  The rank defect is solved by adding two conditions to the observation
%  equations
%  
%       AC1 * [v ; offsets; refseries ] = 0
%       AC2 * [v ; offsets; refseries ] = 0
%
%  depending on the option 'refsystem':
%
%     - 'min-refseries' (default):   
%             AC1= [ zeros(1,2*numPoints) ones(1,numEpochs) ]
%             AC2= [ zeros(1,2*numPoints) t-t0 ]
%     - 'min-velocities
%             AC1= [ zeros(1,2*numPoints) ones(1,numEpochs) ]
%             AC2= [ ones(refPntMask)/nref zeros(1,numPoints+numepochs)  ]
%     - 'minimum-norm':            
%             AC1=N1, AC2=N2
%     - 'inherit':            
%             sets relative=false, full-rank, no constraints necessary
%
%  The rank defects are solved with pseudo-observations (instead of hard
%  constraints). 
%
%  In case of the model without reference time series, (option 'relative' 
%  is false) the model has full rank. This is because the input space time
%  mattix uses a consistent reference frame over time.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:  28 Sep 2023 by Hans van der Marel
% Modified: 29 Sep 2023 by Hans van der Marel
%            - working version without full options
%            - alpha release 
%            1 Oct 2023 by Hans van der Marel
%            - extended options and minor changes
%            - min-velocities refsystem option added
%            - 2nd alpha release
%            4 Oct 2023 by Hans van der Marel
%            - Compute Qe_hat, output Qe_hat and Qy, useful for prediction 
%            - 3rd alpha release
%           09 Jan 2024 by Freek van Leijen
%            - Additional output of design matrices, for later reuse
%           30 Jan 2024 by Hans van der Marel
%            - hardcoded tolerances moved to (undocumented) options
%            - freeze of arguments and options
%            - first public release

% Consider:  xx xxx xxxx by Hans van der Marel
%            - Remove dmat from this function, implement in own function
%              with options to gap-fill dmat and extend/extrapolate dmat,
%              extrapolation range to be specified in dyears
%            - Use Lagrange multipliers to solve for the constraints
%              instead of, or as option, using pseudo-observations,
%              check correct scaling of cweight when using
%              pseudo-observations

% Check input arguments and process options

if nargin < 1
   error('This function expects at least one input argument.')
end

opt.trendmodel='linear';      % trend model (default 'linear')
opt.t0='best';                % reference epoch {'mean' , 'best', <dyear>} (default 'best')

opt.relative=true;            % true if displacements in STM are relative (default true)
opt.refsystem='min-refseries';% Solution to rankdefect {'minimum-norm','min-refseries'} (default 'min-refseries')
opt.cvelocity=0.;             % Velocity for point velocity constraints (Default 0)
opt.cweight=1;                % Weight for the pseudo-observation constraints (Default 1)

opt.ROI=[];                   % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default none)
opt.POI=[-Inf +Inf];          % Period of interest [ dYearStart dYearEnd ] (Default none)
opt.pntMask=[];               % Point mask (default [], is all), will be combined with ROI
opt.epochMask=[];             % Epoch mask (default [], is all), will be combined with POI

opt.minObsDisplPoint=2;       % Minimum number observed displacements per point, if below, point is removed
opt.minObsDisplEpoch=1;       % Minimum number observed displacements per epoch, if below, epoch is removed (default is 2 for relative, 1 otherwise)
opt.verbose=0;                % Verbose level (default 0)
opt.check=false;              % If true, check some of the computions (default is false)

opt.pntCrdType='deg/m';       % Coordinate type in <stm>.pntCrd {'deg/m', 'km/m'}, default 'deg/m' (Latitude, Longitude, Height)
opt.ignoreStochModel=false;   % If true, ignore stochasticmodel from stm, use unit matrix instead 

opt.mincond=1e-3;             % Warning and corrective action if rcond < opt.mincond
opt.inversenorm=1e-5;         % Warning if inverse norm exceeds this criterion

opt.sregularize=1e0;         % regularization factor if covariance matrix turns out singular

opt.doqr=false;               % If true, do QR factorization (not yet supported)

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end

% Number of points (np), epochs (ne) and epoch times for this dataset

numPoints=stm.numPoints;     % Number of points in dataset
numEpochs=stm.numEpochs;     % Number of epochs in dataset
pntCrd=stm.pntCrd;           % Point coordinates
epochDyear=stm.epochDyear;   % Epoch times in dyear

% Set point and epochs masks from ROI and epochs in POI

if isempty(opt.ROI)
   pntMask=true(numPoints,1);
else
   roi=bbox2poly(opt.ROI);
   pntMask=inpolygon(pntCrd(:,1),pntCrd(:,2),roi(:,1),roi(:,2)); 
end
epochMask=( epochDyear >= opt.POI(1) & epochDyear <= opt.POI(2) );

% Set/update pntMask and/or epochMask from pntMask/epochMask option

if ~isempty(opt.pntMask)
    pntMask = pntMask & opt.pntMask;
end
if ~isempty(opt.epochMask)
    epochMask = epochMask & opt.epochMask;
end

% Prepare reference system option

if islogical(opt.refsystem) || isnumeric(opt.refsystem) 
   % set pntRefMask for the refsystem='min-velocities' 
   if islogical(opt.refsystem)
     % mask with reference points
     pntRefMask = opt.refsystem;
   elseif any(size(opt.refsystem)==1) 
     % point number or array with reference point numbers
     pntRefMask = false(numPoints,1);
     pntRefMask(opt.refsystem)=true;
   else 
     % bounding box or polygon for the reference points
     refroi=bbox2poly(opt.refsystem);
     pntRefMask=inpolygon(pntCrd(:,1),pntCrd(:,2),refroi(:,1),refroi(:,2)); 
   end
   fprintf('Number of reference point candidates (velocity constraints): %d\n',sum(pntRefMask))
   pntRefMask = pntRefMask & pntMask;
   fprintf('Number of reference point remaining after applying points mask: %d\n',sum(pntRefMask))
   opt.refsystem='min-velocities';
elseif strcmpi(opt.refsystem,'inherit')
    % refsystem='inherit' is the same as relative=false
   opt.relative=false;
end

if opt.relative && opt.minObsDisplEpoch == 1
   opt.minObsDisplEpoch=2;     % Must have at least two displacements per epoch in case or relative displacement, if below, epoch is removed 
end

% Convert latitude and longitude to km for use with stochastic model (if necessary)

pntAttrib=stm.pntAttrib;
if isfield(pntAttrib,'pntNeu')
    pntNeu=pntAttrib.pntNeu;
else
    switch lower(opt.pntCrdType)
       case 'deg/m'
          % convert latitude/longitude into local topocentric coordinates
          [pntNeu,~] = plh2neusp(pntCrd);     % deg/deg/m -> m/m/m
          pntNeu(:,1:2)=pntNeu(:,1:2)./1000;     % m/m/m -> km/km/m
       case 'km/m'
          % coordinates are already in the right units, this is exceptional,
          % and only happens for simulations, just copy
          pntNeu=pntCrd;
    otherwise
       error('unknown pntCrdType option')        
    end
end
pntNeu=pntNeu(:,1:2);

% Reference epoch and delta times

switch opt.t0
    case 'mean'
       t0 = mean(epochDyear(epochMask));
    case 'best'
       epochCount = sum(sum(~isnan(stm.obsData(pntMask,epochMask,:)),3),1);
       t0 = sum(epochDyear(epochMask) .* epochCount) / sum(epochCount);
    otherwise
       t0 = opt.t0;
end
dt = epochDyear(:) - t0;

% Initialize output variables

numObsTypes=length(stm.obsTypes);

v=nan(numPoints,numObsTypes);
offsets=nan(numPoints,numObsTypes);
refseries=nan(numEpochs,numObsTypes);

qxhat=cell(numObsTypes,1);
qehat=cell(numObsTypes,1);
qyobs=cell(numObsTypes,1);

omt=nan(1,numObsTypes);
emat=nan(numPoints,numEpochs,numObsTypes);
dmat=nan(numPoints,numEpochs,numObsTypes);

% Loop over observation types

for l=1:numObsTypes

  % Set up observation equations

  switch opt.trendmodel

    case 'linear'
  
       % Unknowns are                   
       % - velocities (numPoints)
       % - offsets at t0 (numPoints)
       % - reference time series (numEpochs) if relative==true
       % Total is 
       % - numPoints*2 + numEpochs, with a rank defect of one, if relative is true 
       % - numPoints*2  if relative is false

       % update point masks   (need to iterate over this...)

       numDisplPoint = sum(~isnan(stm.obsData(:,epochMask,l)),2);
       pntMaskDs = pntMask & numDisplPoint >= opt.minObsDisplPoint ;

       numDisplEpoch = sum(~isnan(stm.obsData(pntMaskDs,:,l)),1);
       epochMaskDs = epochMask & numDisplEpoch >= opt.minObsDisplEpoch;
       if any(xor(epochMaskDs,epochMask))
           % epoch mask has changed, repeat selection
           numDisplPoint = sum(~isnan(stm.obsData(:,epochMaskDs,l)),2);
           pntMaskDs = pntMask & numDisplPoint >= opt.minObsDisplPoint ;
       end

       if opt.verbose > 1
           disp('numDisplPoints:')
           disp(numDisplPoint')
           disp('pntMaskDs:')
           disp(pntMaskDs')
           disp('numDisplEpochs:')
           disp(numDisplEpoch)
           disp('epochMaskDs:')
           disp(epochMaskDs)
       end

       % Build the nominal sparse design matrix A for this dataset

       np=sum(pntMaskDs);    % number of valid points
       ne=sum(epochMaskDs);  % number of valid epochs
       m = np * ne;          % max number of possible observations

       if np == 0 || ne == 0
           warning([ stm.obsTypes{l} ' has not sufficient points or epochs... '])
           continue; 
       end 
       
       ia0 = repmat(1:np,[ne,1])'; ia0=ia0(:);
       a0 = repmat(dt(epochMaskDs),[1,np])'; a0=a0(:);
       ia1 = np +ia0;
       a1 = repmat(ones(ne,1),[1,np]); a1=a1(:); 
       ia2 = np*2 + repmat(1:ne,[np,1]);ia2=ia2(:);
       a2 = -a1;

       if opt.relative 
          n = np * 2 + ne;   % Max unknowns
          rankdefect=2;
          A = sparse( [1:m 1:m 1:m ]' , [ia0 ; ia1 ; ia2 ], [a0 ; a1; a2], m, n, m*3);
       else
          n = np * 2;        % Max unknowns
          rankdefect=0;
          A = sparse( [1:m 1:m ]' , [ia0 ; ia1 ], [a0 ; a1 ], m, n, m*2);
       end

       % Get space-time matrix with observations 

       y=stm.obsData(pntMaskDs,epochMaskDs,l);
             
       % Vectorize observed minus a-priori and mask the non-observed (nan) entries

       y=y(:);    
       ymask=~isnan(y);
       y=y(ymask);
       my=length(y);

       Amatrices{l} = A(ymask,:); %store for later use in error propagation

       % Compose co-variance matrix of observations

       lstoch=min(l,numel(stm.stochModel));   % PATCH FOR INSAR DATASETS
       if opt.ignoreStochModel
          Qy=eye(my,my);
       elseif any(strncmp(stm.stochModel{lstoch},'tudinsar4rd',11))
          % stochModel{l}='tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)';  % there is a problem with the tudinsar4rd model (not positive definite)
          pntAttrib=stm.pntAttrib;
          epochAttrib=stm.epochAttrib;
          Qy=stmstochmodel(stm.stochModel{lstoch},stm.stochData,pntNeu,epochDyear,{pntMaskDs,epochMaskDs},pntAttrib,epochAttrib);
          Qy=Qy(ymask,ymask);
       else
          Qy=stmstochmodel(stm.stochModel{lstoch},stm.stochData,pntNeu,epochDyear,{pntMaskDs,epochMaskDs});
          Qy=Qy(ymask,ymask);
       end
       Qyopt='full'; 

       % set parameter mask

       % x2mask=true(n,1);

       % Least squares solution (BLUE) 
       %
       % Ry=chol(Qy) ->  Qy = Ry'*Ry  => inv(Qy) = inv(Ry)*inv(Ry')
       %
       %    A2n=inv(Ry')*A2           -> A2n=Ry'\A2
       %    yn=inv(Ry')*y             -> yn=Ry'\y
       %    b2a=A2'*inv(Qy)*y         -> b2a=A2n'*yn
       %    N22=A2'*inv(Qy)*A2        -> N22=A2n'*A2n
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
             A2n = wy.*A(ymask,:);     % A2n=Qy^(-1/2)'*A2
             yn = wy.*y;               % yn=Qy^(-1/2)'*y
             %Ry=chol(diag(Qy));       % Not needed now (otherwise could be done more efficiently)
          case 'full'
             % full co-variance matrix, compute Cholesky factor
             [Ry,flag]=chol(Qy);
             if flag ~= 0
                warning(['Co-variance matrix of dataset ' num2str(l) ' is not positive definite, do regularization and continue ...'])
                %Qy = Qy + opt.sregularize*ones(size(Qy));
                Qy = Qy + opt.sregularize*eye(size(Qy));
                [Ry,flag]=chol(Qy);
                if flag ~= 0
                   error(['Co-variance matrix of dataset ' num2str(l) ' is not positive definite even after regularization, quiting.'])
                end
             end
             A2n=Ry'\A(ymask,:);        % A2n=Qy^(-1/2)'*A2=inv(Ry')*A
             yn=Ry'\y;                  % yn=Qy^(-1/2)'*y=inv(Ry')*y
           otherwise
             error(['Unsupported covariance matrix option ' Qopt ]) 
       end
       
       % Compute solution

       x2a=nan(n,1);

       if opt.doqr
          % Use QR decomposition to obtain the solution x2a
          if issparse(A2n)
             perm = colamd(A2n(:,1:end-rankdefect));
             [z,U22] = qr(A2n(:,perm),yn,0);
          else
             [Q,U22,perm] = qr(A2n(:,1:end-rankdefect),0);
             z = Q'*yn;
          end
          if rankdefect == 0
             x2a(perm)= U22 \ z ;
          else
             x2a([ perm  end-rankdefect:end ])= [ U22 \ z ; zeros(rankdefect,1) ];
          end
       else
          % Normal equations 
          b2a=A2n'*yn;             % b2a=A2'*Qy^(-1)*y=A2*inv(Ry)*inv(Ry')*y=A2n'*yn
          N22=A2n'*A2n;            % N22=A2'*Qy^(-1)*A2=A2*inv(Ry)*inv(Ry')*A=A2n'*A2n
          % Solve rankdefect by adding constraints (pseudo-observations)
          if opt.relative
              % Null space
              null1=[ zeros(np,1) ; ones(np,1) ; ones(ne,1) ];
              null2=[ ones(np,1) ; zeros(np,1) ; dt(epochMaskDs) ];
              % Check Null-space
              if any(A(ymask,:)*null1 ~= 0) || any(A(ymask,:)*null2 ~= 0) 
                  warning(['Null-space is not as expected in component ' num2str(l) ])
                  disp(null1'*A(ymask,:)')
                  disp(null2'*A(ymask,:)')
              end
              % Solve rank defect by adding pseudo observations
              switch opt.refsystem
                  case 'minimum-norm'
                     % minimum norm solution (makes only sense if expected
                     % value of offsets is zero!)
                     ac1=null1';
                     ac2=null2';
                     cvel=0;
                     cweight=opt.cweight;
                  case 'min-refseries'
                     % minimize refseries corrections 
                     ac1=[ zeros(np,1) ; zeros(np,1) ; ones(ne,1) ]';
                     ac2=[ zeros(np,1) ; zeros(np,1) ; dt(epochMaskDs) ]';
                     cvel=0;
                     cweight=opt.cweight;
                  case 'min-velocities'
                     % minimize corrections to
                     % - refseries
                     % - velocities of selected points (by ROI or logical array)
                     pntRefMaskDs = pntMaskDs & pntRefMask;
                     fprintf('Number of velocity constraints for %s: %d (out of %d)\n',stm.obsTypes{l},sum(pntRefMaskDs),sum(pntMaskDs))
                     ac1=[ zeros(np,1) ; zeros(np,1) ; ones(ne,1) ]';
                     if sum(pntRefMaskDs) >= 1
                        ac21=zeros(numPoints,1);
                        ac21(pntRefMaskDs)=1/sum(pntRefMaskDs);
                        ac2=[ ac21(pntMaskDs) ; zeros(np,1) ; zeros(ne,1) ]';
                        cweight=opt.cweight*sum(pntRefMaskDs)^2;
                     else
                        warning('at least one velocity constraint must be within the pntMask, continue with min-refseries...')
                        ac2=[ zeros(np,1) ; zeros(np,1) ; dt(epochMaskDs) ]';
                        cweight=opt.cweight;
                     end
                     cvel=opt.cvelocity;
                  otherwise
                     error(['unknown refsystem option' opt.refsystem ])
              end
              N22=N22 + ac1'*cweight*ac1 + ac2'*cweight*ac2;
              b2a=b2a + ac2'*cweight*cvel;
              rankdefect=0;
          end
          % Use Cholesky factorization to obtain the solution x2a
          %
          % U22 = chol(N22)  -> N22 = U22'*U22  => inv(N22)=inv(U22)*inv(U22')
          %
          %    x2a=inv(N22)*b2a          -> x2a=N22\b2a=U22\(U22'\b2a)
          %    Qx=inv(N22)               -> Qx=N22\I=U22\(U22'\I)
          %   
          % Determine order for Cholesky factorization for maximum stability
          [~,perm]=sort(diag(N22),'descend');
          [U22P,flag]=chol(N22(perm,perm)); 
          rcond=diag(U22P);
          if opt.verbose > 0
              disp('Square root of Cholesky factor diagonal:')
              disp(rcond')
              disp('Sequential decrease in square root of Cholesky factor diagonal:')
              disp(rcond(1:end-1)'./rcond(2:end)')
          end
          if flag ~=0 && n-flag+1 ~= rankdefect
              warning(['Rank-defect (check 1) is not as expected in component ' num2str(l) ', actual ' num2str(n-flag+1) ', expected ' num2str(rankdefect), ', trying to continue...' ])
              rankdefect=n-flag+1;
          end
          flag2=find(rcond < opt.mincond,1); if isempty(flag2), flag2=0; end
          if flag2 ~=0 && n-flag2+1 ~= rankdefect
              warning(['Rank-defect (check 2) is not as expected in component ' num2str(l) ', actual ' num2str(n-flag2+1) ', expected ' num2str(rankdefect), ', min(rcond) ' num2str(min(rcond)) ', trying to continue...' ])
              rankdefect=n-flag2+1;
              U22P(:,flag2:end)=[];
              U22P(flag2:end,:)=[];
          end
          b2p=b2a(perm);
          x2p=[ U22P \ ( U22P' \ b2p(1:end-rankdefect) ) ; zeros(rankdefect,1) ];     % x2a = N22\b2a = inv(N22)*b2a                
          x2a(perm)=x2p;
       end

       % update velocity vector etc.

       v(pntMaskDs,l)=x2a(1:np);
       offsets(pntMaskDs,l)=x2a(np+1:np*2);
       if opt.relative 
           refseries(epochMaskDs,l)=x2a(np*2+1:end);
       else
           refseries(epochMaskDs,l)=zeros(ne,1);
       end

       % residuals and omt
       %
       %     e = y - A2*x2;
       %     en = inv(Ry)*e = Ry\e   -> 
       %     en = Ry\y - Ry\A2*x2 = yn - A2n*x2 

       en = yn - A2n*x2a;
       mse=en'*en;
       omt(1,l)=mse/(m-n+rankdefect);

       if opt.verbose > 0
           disp('Overall model test (OMT) values:')
           disp(omt)
       end

       e = nan(np,ne);
       e(ymask) = y - A(ymask,:)*x2a;
       emat(pntMaskDs,epochMaskDs,l) = e;

       % model values (offsets and reference series removed)
       %
       % To do:
       % - gap filling and optional extrapolation
       % - option to add reference series?
       %
       % Maybe we should do this outside this function

       s = nan(np,ne);
       s(ymask) = A(ymask,1:np)*x2a(1:np);
       dmat(pntMaskDs,epochMaskDs,l) = s;

       % compute co-variance matrix Qxhat (optional)

       if nargout > 7

           Q22=nan(n,n);
           if rankdefect > 0
              Q22(perm,perm) = [ U22P \ ( U22P' \ eye(n-rankdefect,n-rankdefect) )  zeros(n-rankdefect,rankdefect) ;
                              zeros(rankdefect,n ) ];
           else
              Q22(perm,perm) = U22P \ ( U22P' \ eye(n,n) );
           end
           inversenorm=norm(eye(n,n)-Q22*N22);
           if opt.verbose > 0
              disp(['norm(Q*N-I) for component ' num2str(l) ',is ' num2str(inversenorm) ])
           end
           if  inversenorm > opt.inversenorm
              warning(['Possible problem with inverse in component ' num2str(l) ', norm(Q*N-I) is ' num2str(inversenorm) ', continuing...' ])
           end
           qxhat{l} = Q22;
       end

       % output co-variance matrix Qehat and Qy (optional)

       if nargout > 8

           % Compute Qehat using straightforward method
           % 
           %    Qehat = Qy - A * Qx A'

           qehat{l} = Qy - A(ymask,:)*Q22*A(ymask,:)';
           
           if opt.check

               % Compute Qehat using more stable method (but also more resource intensive)
               %
               %    PAorthogonal = eye(m,m) - A(ymask,:)*Q22*A(ymask,:)'*inv(Qy)      inv(Qy)=inv(Ry)*inv(Ry');
               %    PAorthogonal * Ry'  = Ry' - A(ymask,:)*Q22*A(ymask,:)'*inv(Ry);
               %    PAorthogonal * Ry'  = Ry' - A(ymask,:)*Q22*A2n'    ->  PAon
               %
               %    Qehat = PAorthogonal * Qy * PAorthogonal' = PAon * PAon';

               PAn = Ry' - A(ymask,:)*Q22*A2n';
               qealt = PAn * PAn';

               dnorm=norm(qealt-qehat{l});
               disp(['Norm(Qehat1-Qehat2) component ' num2str(l) ': ' num2str(dnorm) ] )
               if dnorm > 1e-6
                  warning(['Possible problem with Qehat computation, component ' num2str(l) ', norm(Qehat2-Qehat1) is ' num2str(dnorm) ', continuing...' ])
               end
               qehat{l}=qealt;

           end

           % Save covariance matrix of observations (same dimensions as Qehat) for later use           

           qyobs{l} = Qy;

       end

  otherwise
    error (['unsupported trend model ' opt.trendmodel])
  end

end

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

end


