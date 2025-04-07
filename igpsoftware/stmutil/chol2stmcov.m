function [stochModel,stochData]=chol2stmcov(obsTypes, obsData, U, idxU, varargin)
%CHOL2STMCOV  Compute covariance matrix from Cholesky factor. 
%   [STOCHMODEL,STOCHDATA]=CHOL2STMCOV(OBSTYPES,OBSDATA,U,IDXU) 
%   computes the covariance matrix Q from the cholesky factor U, and return it in
%   the covariance matrix in STOCHMODEL and STOCHDATA. The input is the
%   Cholesky factor U and index IDXU to the observation data OBSDATA with
%   observation types OBSTYPES. The output STOCHMODEL is a character string 
%   with the model specification and output STOCHDATA is a numerical array/matrix 
%   with the co-variance matrix, possibly including an index that gives
%   for each row/column in the covariance matrix the corresponding 
%   element in OBSDATA. This index is similar to IDXU, but not necessarily 
%   identical, because the computation base (zero columns/rows) can be 
%   inserted in the co-variance matrix. 
%
%   [STOCHMODEL,STOCHDATA]=CHOL2STMCOV(...,'option',value) supports the
%   following option/value pairs
%
%     option       default description
%     ------------ ------- -----------------------------------------
%     byobstype    true    Return co-variance matrix for each observation
%                          type. If false, one big co-variance matrix is
%                          returned for all observation types together.
%     addcb        true    Include zero rows/columns for the computing base.
%     optimized    true    Use optimized computation method when byobstype==true
%     plot         0       Plot intermediate results: 0 (none), 1 or 2
%     sreg         0.0     regularization parameter
%
%   The option 'byobstype' is an important one. If true, it returns the
%   (indexed) full covariance matrix for each observation type, but not the
%   covariances between observation types. The stochModel is something like
% 
%      stochModel = { 
%         { 'covmatrix(format=index,lidx=#,oidx=#,odata=#,[size=#],[sreg=#])' } ,
%         { 'covmatrix(format=index,lidx=#,oidx=#,odata=#,[size=#],[sreg=#])' } ,
%         { 'covmatrix(format=index,lidx=#,oidx=#,odata=#,[size=#],[sreg=#])' } 
%      }
% 
%   with a full indexed co-variance matrix for each observation type. The parameter
%   odata refers to start position (odata+1) in stochData (default odata=0), 
%   size is the dimension of the covariance matrix and sreg is the regularisation 
%   parameter that was used to regularize each matrix. The index, which
%   contains the column numbers, is given in stochData(oidx+1:oidx+lidx),
%   the co-variance matrix is in stochData(odata+1:odata+lidx^2).
%
%   If you need the covariances between observation types, then the option
%   'byobstype' should be set to false. The output stochModel is then
%
%      stochModel = { { 'covmatrix(format=full,[sreg=#])' } }
%
%   However, this type of stochModel is not supported by modules (such as
%   stmvelocity or stmintegrate which process data observation type by
%   observation type).
%
%   Default value for 'sreg' is 0.3e-4 (sqrt(1e-9)). If they are full-rank, 
%   use 'sreg=0' to skip regularization, though the regularization won't 
%   hurt even in case of full rank.
%
%   Example:
%
%     [stm.stochModel, stm.stochData] = chol2stmcov(obsType, obsData, U, idxU)
%     [stm.stochModel, stm.stochData] = ... 
%          chol2stmcov(obsType, obsData, U, idxU, 'plot', true)
%  
%   (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:  14 Nov 2023 by Hans van der Marel
% Modified: 27 Nov 2023 by Hans van der Marel
%            - Development version 
%           20 Dec 2023 by Hans van der Marel
%            - Release for testing with stmStochModel
%           25 Jan 2024 by Hans van der Marel
%            - Final parameter names for use with redesigned stmStochModel
%           30 Jan 2024 by Hans van der Marel
%            - size parameter is now an array with [#pnts #epochs]


%% Check input arguments and process options

% check number of input parameters

if nargin < 4
   error('This function expects at least three input arguments.')
end

% process options

opt.byobstype=true;           % split output stochastic model in observation types (default true)
opt.addcb=true;               % add zero rows/columns for the computing base (default true)
opt.optimize=true;            % use optimized computation in case opt.byobstype=true;
opt.sreg = 0.0;               % regularization parameter
opt.plot=0;                   % plot intermediate results: 0 (none), 1 or 2 (default is 0)
opt.check=false;              % 

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end

%% Set up index arrays for the Cholesky factor and computing base
%
% Cholesky factor
%
% - the Cholesky factor is only available for the subset of parameters that are actually solved, 
%   this does NOT include the computation base and non-estimable (e.g. North) parameters
% - thus we provide an additional index array to refer to the position in the space time matrix
% - the computing base is not included in the Cholesky factor

% Optionally plot the Cholesky factor

if opt.plot > 0
    figure
    imagesc(U)
    colorbar
    cmax=quantile(abs(U(:)),.95);
    clim([-cmax cmax])
    text(0.1,0.1,['min: ' , num2str(min(U(:))), ', max: ' , num2str(max(U(:))), ', 95%: ' num2str(cmax) ],'Units','normalized')
    title('Cholesky factor')
end

% Get the dimension of the space time matrix and observation mask

[numPoints,numEpochs,numObsTypes]=size(obsData);
ymask=~isnan(obsData(:));

% get subscript numbers for the points, epochs and observation types

[idxPoints,idxEpochs,idxObsTypes]=ind2sub([numPoints,numEpochs,numObsTypes],idxU);

% the Cholseky factor is blocked, get the block boundaries in subObsTypes

subObsTypes=[0 ; find(diff(idxObsTypes)> 0) ; numel(idxObsTypes) ];
if numel(subObsTypes) ~= numObsTypes+1
   error('Error subObsTypes !!!')
end
for l=1:numObsTypes
   if ~all(idxObsTypes(subObsTypes(l)+1:subObsTypes(l+1))==l)
       error('Error idxObsTypes!!')
   end
end

% get indices for the computing base (these are zero in the Cholesky factor, bot not saved, and therefore not present in U)

baseMask=ymask;
baseMask(idxU)=false;
idxBase=find(baseMask);
[idxBasePoints,idxBaseEpochs,idxBaseObsTypes]=ind2sub([numPoints,numEpochs,numObsTypes],idxBase);

% get the block boundaries for the computing base in subBaseObsTypes

subBaseObsTypes=[0 ; find(diff(idxBaseObsTypes)> 0) ; numel(idxBaseObsTypes) ];
if numel(subBaseObsTypes) ~= numObsTypes+1
   error('Error subBaseObsTypes!!!')
end
for l=1:numObsTypes
   if ~all(idxBaseObsTypes(subBaseObsTypes(l)+1:subBaseObsTypes(l+1))==l)
       error('Error idxBaseObsTypes!!')
   end
end

% compute new set of indices with computing base included

if opt.addcb
   npar=numel(idxU);
   nbase=numel(idxBase);
    
   idxQ = [idxU ; idxBase ];
   [idxQ,iperm]=sort(idxQ);
    
   % get subscript numbers for the points, epochs and observation types
    
   [idxPoints2,idxEpochs2,idxObsTypes2]=ind2sub([numPoints,numEpochs,numObsTypes],idxQ);
    
   % the Covariance matrix is blocked, get the block boundaries in subObsTypes
    
   subObsTypes2=[0 ; find(diff(idxObsTypes2)> 0) ; numel(idxObsTypes2) ];
   if numel(subObsTypes2) ~= numObsTypes+1
      error('Error subObsTypes2 !!!')
   end
   for l=1:numObsTypes
      if ~all(idxObsTypes2(subObsTypes2(l)+1:subObsTypes2(l+1))==l)
          error('Error idxObsTypes2!!')
      end
   end
else
   subObsTypes2=subObsTypes;
   idxPoints2=idxPoints;
   idxEpochs2=idxEpochs;
   idxObsTypes2=idxObsTypes;
   idxQ=idxU;
end 

% Optionally plot the indices 

if opt.plot > 0
    for l=1:numObsTypes
       figure;
       plot(idxEpochs(idxObsTypes==l),idxPoints(idxObsTypes==l),'.')
       hold on;
       plot(idxBaseEpochs(idxBaseObsTypes==l),idxBasePoints(idxBaseObsTypes==l),'*r')
       xlabel('Epochs')
       ylabel('Points')
       title(['Space time matrix elements ' obsTypes{l}])
       legend('observed','computing base')
    end
end


%% Compute the full covariance matrix

if ~opt.optimize || ~opt.byobstype || opt.check

   %  Compute the full covariance matrix from the Cholesky factor

   tic
   Q= U \ ( U' \ eye( size(U ) ) );
   toc

   % Plotting of the full covariance matrix (optional)

   if opt.plot > 1
       plotfullcov(Q,numPoints,numEpochs,numObsTypes,obsTypes,subObsTypes,idxPoints,idxEpochs,subBaseObsTypes,idxBasePoints,idxBaseEpochs)
   end

   % Insert zeros for computing base

   if opt.addcb
       Q = [ Q zeros(npar,nbase) ; zeros(nbase,npar) zeros(nbase,nbase) ];
       Q = Q(iperm,iperm);
   end 

   % Plotting of the full covariance matrix (optional)

   if opt.plot > 0
       plotfullcov(Q,numPoints,numEpochs,numObsTypes,obsTypes,subObsTypes2,idxPoints2,idxEpochs2)
   end

   % Regularization

   if opt.sreg > 0  
       N = kron( eye(numObsTypes,numObsTypes) , ...
          [ repmat(eye(numPoints,numPoints),[numEpochs,1]) kron(eye(numEpochs,numEpochs),ones(numPoints,1)) ]);
       N = N(ymask,:);
       Q = Q + N*N'*opt.sreg^2;
   end

end


%% Compute the covariance matrix in parts

% Covariance matrix of the up component (see Van der Marel, 1988, pp 128-132,137-141,215):
%
%  Normal matrix
%
%    N = U'*U 
% 
%    N = | N11  N12 | = | U11'      | * | U11 U12 |  =  | U11'*U11      U11'*U12       |  
%        | N12' N22 |   | U12' U22' |   |     U22 |     | U12'*U11  U22'*U22+U12'*U12  |  
%
%    U11 = chol(N11) ,   U12 = inv(U11)'*N12 = U11'\N12 ,  
%                               U22 = chol(N22-U12'*U12) = chol(N22-N12'*inv(N11)*N12)
%
%  Cholesky factor         Inverse Cholesky factor (not needed explicitly)
%
%    U = | U11 U12 |          inv(U) = | inv(U11)  -inv(U11)*U12*inv(U22) |       
%        |     U22 |                   |                   inv(U22)       |
%
%                             inv(U) = U \ I ,  inv(U') = U' \ I
%  Inverse
%
%    Q = inv(N) = inv(U)*inv(U') = U \ ( U' \ I ) 
%
%       = | inv(U11) -inv(U11)*U12*inv(U22) |  *  |      inv(U11)'                       | 
%         |                 inv(U22)        |     | -inv(U22)'*U12'*inv(11)'   inv(U22)' |
%
%       = | inv(U11)*inv(U11)' + inv(U11)*U12*inv(U22)*inv(U22)'*U12'*inv(U11)'  -inv(U11)*U12*inv(U22)*inv(U22)' | 
%         |            -inv(U22)*inv(U22)'*U12'*inv(U11)'                            inv(U22)*inv(U22)'           |
%
%       = | inv(N11) + inv(N11)*N12*inv(N22)*N12'*inv(N11)   -inv(N11)*N12*inv(N22) | 
%         |            -inv(N22)*N12'*inv(N11)                   inv(N22)           |
%
%       = | inv(N11)*(N11 + N12*inv(N22)*N12')*inv(N11)   -inv(N11)*N12*inv(N22) | 
%         |            -inv(N22)*N12'*inv(N11)                inv(N22)           |
%
% So, to get the inverse of the Up component (which is at the end of the
% Cholesky factor), is easy:
%
%   Q22 = U22 \ ( U22' \ I )
%
% The inverse of the first part is
%
%   Q11 = inv(U11)*inv(U11)' + inv(U11)*U12*inv(U22)*inv(U22)'*U12'*inv(U11)' =
%         inv(U11)*(I+U12*inv(U22)*inv(U22)'*U12')*inv(U11)' =
%         inv(U11)*(I+V12*V12')*inv(U11)' =  U11 \ ( (I+V12*V12')*inv(U11)' ) =
%         U11 \ ( (I+V12*V12') / U11' )
%
%   Q12 = -inv(U11)*U12*inv(U22)*inv(U22)' =
%         -inv(U11)*V12*inv(22)' =
%         -inv(U11)*W12 = U11 \ W12
%
% with 
%   V12' = inv(U22)'*U12' = U22' \ U12'
%   W12' = inv(U22) * V12' = inv(U22)*inv(U22)'*U12' = U22 \ ( U22' \ U12' ) =  U22 \ V12'
%

if opt.byobstype && opt.optimize 

   % Split in two blocks

   l=numObsTypes;
   k1=subObsTypes(l);
   k2=subObsTypes(l+1);

   % Compute the second block (with Up)

   Q22 = U(k1+1:k2,k1+1:k2) \ ( U(k1+1:k2,k1+1:k2)' \ eye(k2-k1) );
   
   % Compute the first block (with North and East)
    
   V12T = ( U(k1+1:k2,k1+1:k2)' \  U(1:k1,k1+1:k2)' );
   Q11 = U(1:k1,1:k1) \ ( ( eye( k1,k1) + V12T'*V12T ) / U(1:k1,1:k1)' );

   % extract the sub-matrices into cell arrays

   Qsub=cell(numObsTypes,1);
   idxQsub=cell(numObsTypes,1);
   for l=1:numObsTypes-1
      kk1=subObsTypes(l)+1;
      kk2=subObsTypes(l+1);
      Qsub{l}=Q11(kk1:kk2,kk1:kk2);
      idxQsub{l}=idxU(kk1:kk2)-(l-1)*(numPoints*numEpochs);
      %(l-1)*(numPoints*numEpochs)
   end
   l=numObsTypes;
   kk1=subObsTypes(l)+1;
   kk2=subObsTypes(l+1);
   Qsub{l}=Q22;
   idxQsub{l}=idxU(kk1:kk2)-(l-1)*(numPoints*numEpochs);
   % (l-1)*(numPoints*numEpochs)

   % Insert zeros for computing base

   if opt.addcb
       for l=1:numObsTypes
          k1=subObsTypes(l)+1;
          k2=subObsTypes(l+1);
          kb1=subBaseObsTypes(l)+1;
          kb2=subBaseObsTypes(l+1);
          npar=k2-k1+1;
          nbase=kb2-kb1+1;
          idxQtmp = [idxU(k1:k2) ; idxBase(kb1:kb2) ];
          [idxQtmp,iperm]=sort(idxQtmp);
          tmp = [ Qsub{l} zeros(npar,nbase) ; zeros(nbase,npar) zeros(nbase,nbase) ];
          Qsub{l} = tmp(iperm,iperm);
          idxQsub{l}=idxQtmp-(l-1)*(numPoints*numEpochs);
       end
   end

   % Regularization

   if opt.sreg > 0  
       l2=0;
       for l=1:numObsTypes
           l1 = l2+1;
           l2 = l2 + numPoints*numEpochs;
           N = [ repmat(eye(numPoints,numPoints),[numEpochs,1]) kron(eye(numEpochs,numEpochs),ones(numPoints,1)) ];
           N = N(ymask(l1:l2),:);
           Qsub{l} = Qsub{l} + N*N'*opt.sreg^2;
       end
   end

   % Plot and compare
    

   if opt.plot > 0 && opt.check

      l=numObsTypes;
      k1=subObsTypes(l);
      k2=subObsTypes(l+1);
      icomp=numObsTypes;

      figure
      s=sqrt(diag(Q22));
      s2=sqrt(diag(Q(k1+1:k2,k1+1:k2)));
      plot([s s2 s-s2])
      smax=quantile(s,.995);
      ylim([-1 10])
      text(0.01,0.2,['min: ' , num2str(min(s)), ', max: ' , num2str(max(s)), ', 99.5%: ' num2str(smax) , ', rms-diff: ', num2str(std(s-s2))],'Units','normalized')
      legend('Q22','Q','Q22-Q')
      title([ 'Sqrt(diag) Covariance Matrix ' obsTypes{icomp}])

      figure
      s=sqrt(diag(Q11));
      s1=sqrt(diag(Q(1:k1,1:k1)));
      plot([s s1 s-s1])
      smax=quantile(s,.995);
      ylim([-1 10])
      text(0.01,0.2,['min: ' , num2str(min(s)), ', max: ' , num2str(max(s)), ', 99.5%: ' num2str(smax) , ', rms-diff: ', num2str(std(s-s1))],'Units','normalized')
      legend('Q11','Q','Q11-Q')
      title([ 'Sqrt(diag) Covariance Matrix ' strjoin(obsTypes(1:icomp-1), '+')])

      Q22diff=Q22-Q(k1+1:k2,k1+1:k2);

      figure
      imagesc(Q22diff)
      colorbar
      cmax=quantile(abs(Q22diff(:)),.999);
      clim([-cmax cmax])
      text(0.01,0.1,['min: ' , num2str(min(Q22diff(:))), ', max: ' , num2str(max(Q22diff(:))), ', 99.9%: ' num2str(cmax), ', rms-diff: ', num2str(std(Q22diff(:))) ],'Units','normalized')
      title(['Covariance matrix difference ' obsTypes{icomp}])
       
      Q11diff=Q11-Q(1:k1,1:k1);
    
      figure
      imagesc(Q11diff)
      colorbar
      cmax=quantile(abs(Q11diff(:)),.999);
      clim([-cmax cmax])
      text(0.01,0.1,['min: ' , num2str(min(Q11diff(:))), ', max: ' , num2str(max(Q11diff(:))), ', 99.9%: ' num2str(cmax), ', rms-diff: ', num2str(std(Q11diff(:))) ],'Units','normalized')
      title(['Covariance matrix difference ' strjoin(obsTypes(1:icomp-1), '+')])
    
   end


   if opt.plot > 0
    
      figure
      imagesc(Q22)
      colorbar
      cmax=quantile(abs(Q22(:)),.999);
      %clim([-1 20])
      clim([-cmax cmax])
      if opt.check
         text(0.01,0.1,['min: ' , num2str(min(Q22(:))), ', max: ' , num2str(max(Q22(:))), ', 99.9%: ', num2str(cmax), ', rms-diff: ', num2str(std(Q22diff(:))) ],'Units','normalized')
      else
         text(0.01,0.1,['min: ' , num2str(min(Q22(:))), ', max: ' , num2str(max(Q22(:))), ', 99.9%: ', num2str(cmax) ],'Units','normalized')
      end
      title(['Covariance matrix ' obsTypes{numObsTypes}])
    
      figure
      imagesc(V12T)
      colorbar
      %clim([-0.0001 0.0001])
      cmax=quantile(abs(V12T(:)),.99);
      clim([-cmax cmax])
      text(0.01,0.1,['min: ' , num2str(min(V12T(:))), ', max: ' , num2str(max(V12T(:))) ],'Units','normalized')
      title('V12T')

      figure
      imagesc(Q11)
      colorbar
      cmax=quantile(abs(Q11(:)),.999);
      %clim([-1 20])
      clim([-cmax cmax])
      if opt.check
         text(0.01,0.1,['min: ' , num2str(min(Q11(:))), ', max: ' , num2str(max(Q11(:))), ', 99.9%: ', num2str(cmax), ', rms-diff: ', num2str(std(Q11diff(:))) ],'Units','normalized')
      else
         text(0.01,0.1,['min: ' , num2str(min(Q11(:))), ', max: ' , num2str(max(Q11(:))), ', 99.9%: ', num2str(cmax) ],'Units','normalized')
      end
      title(['Covariance matrix ' strjoin(obsTypes(1:numObsTypes-1), '+')])

   end

end

%% Prepare stochModel string and stochData 
%
%  The format for stochModel is either ('byobstype'=true)
%
%      stochModel = { 
%         { 'covmatrix(format=index,[size=#],lidx=#,oidx=#,odata=#,[sreg=#])' } ,
%         { 'covmatrix(format=index,[size=#],lidx=#,oidx=#,odata=#,[sreg=#])' } ,
%         { 'covmatrix(format=index,[size=#],lidx=#,oidx=#,odata=#,[sreg=#])' } 
%      }
% 
%   with a full co-variance matrix for each observation type, or ('byobstype'=false)
%
%      stochModel = { { 'covmatrix(format=index,[size=#],lidx=#,oidx=#,odata=#,[sreg=#])' } }
%
%   with a single full co-variance matrix for all observation types.
%   If you need the covariances between observation types, then this
%   is the way to go, but this is not (yet) supported by the other modules 
%   such as stmvelocity or stmintegrate which process data observation type by
%   observation type.
%
%   The only supported value for format is 'index'. Other formats (e.g.
%   lower or upper) are not yet supported.
%
%   The meaning of the optional parameters in the model specification is 
%
%      parameter  default  type     description
%
%      size                integer   number of rows/columns in the co-variance matrix
%                                    or array with #pnts and #epochs
%      odata               integer   offset co-variance matrix in stochData
%      lidx                integer   length of index in stochData
%      oidx                integer   offset of index in stochData
%      sreg       0        float     regularisation parameter that was used to regularize
%                                    the co-variance matrix.
%
%   The part of stochData that is used are
%
%       stochData(oidx+1:oidx+lidx)          index
%       stochData(odata+1:odata+lidx^2)      co-variance matrix

stochModel={};
if opt.byobstype
   stochData=[];
   odata=0;
   if opt.optimize
      % extract the sub-matricesfrom the first part 
      for l=1:numObsTypes
         k1=subObsTypes2(l)+1;
         k2=subObsTypes2(l+1);
         Qtmp=Qsub{l};
         oidx=odata;
         odata=odata+k2-k1+1;
         stochModel{l}{1}=sprintf('covmatrix(format=index,size=[%d %d],lidx=%d,oidx=%d,odata=%d,sreg=%s)',numPoints,numEpochs,k2-k1+1,oidx,odata,num2str(opt.sreg));  
         odata=odata+(k2-k1+1)^2;
         stochData=[ stochData ; idxQsub{l} ; Qtmp(:) ];
      end
   else
      % extract the sub-matrices from the single-full covariance matrix
      for l=1:numObsTypes
         k1=subObsTypes2(l)+1;
         k2=subObsTypes2(l+1);
         Qsub=Q(k1:k2,k1:k2);
         oidx=odata;
         odata=odata+k2-k1+1;
         stochModel{l}{1}=sprintf('covmatrix(format=index,size=[%d %d],lidx=%d,oidx=%d,odata=%d,sreg=%s)',numPoints,numEpochs,k2-k1+1,oidx,odata,num2str(opt.sreg));  
         odata=odata+(k2-k1+1)^2;
         stochData=[ stochData ; idxQ(k1:k2)-(l-1)*(numPoints*numEpochs) ; Qsub(:) ];
      end
   end
else
   % Single covariance matrix with all the components 
   oidx=0;
   odata=numel(idxQ);
   stochModel{1}{1}=sprintf('covmatrix(format=index,size=[%d %d %d],lidx=%d,oidx=%d,odata=%d,sreg=%s)',numPoints,numEpochs,numObsTypes,size(Q,1),oidx,odata,num2str(opt.sreg));  
   stochData=[idxQ ; Q{:} ];
end


end

%%

function plotfullcov(Q,numPoints,numEpochs,numObsTypes,obsTypes,subObsTypes,idxPoints,idxEpochs,subBaseObsTypes,idxBasePoints,idxBaseEpochs)

extratitle='';
if nargin > 8
    extratitle=' w/o base';
end

% Plot the covariance matrix

figure
imagesc(Q)
colorbar
cmax=quantile(abs(Q(:)),.99);
%clim([-1 20])
clim([-cmax cmax])
text(0.1,0.1,['min: ' , num2str(min(Q(:))), ', max: ' , num2str(max(Q(:))), ', 99%: ' num2str(cmax) ],'Units','normalized')
title(['Covariance matrix' extratitle])

% Plot the diagonals

figure
s=sqrt(diag(Q));
plot(s)
smax=quantile(s,.995);
ylim([0 10])
text(0.1,0.1,['min: ' , num2str(min(s)), ', max: ' , num2str(max(s)), ', 99.5%: ' num2str(smax) ],'Units','normalized')
title(['Sqrt(diag) Covariance Matrix' extratitle])

% Plot the covariance matrix by component

for l=1:numObsTypes
    figure
    k1=subObsTypes(l)+1;
    k2=subObsTypes(l+1);
    Qsub=Q(k1:k2,k1:k2);
    imagesc(Qsub)
    colorbar
    cmax=quantile(abs(Qsub(:)),.999);
    %clim([-1 20])
    clim([-cmax cmax])
    text(0.1,0.1,['min: ' , num2str(min(Qsub(:))), ', max: ' , num2str(max(Qsub(:))), ', 99.9%: ' num2str(cmax) ],'Units','normalized')
    title(['Covariance matrix ' obsTypes{l} extratitle])
end

% Plot the square root of the diagonal

s=sqrt(diag(Q));

figure
for l=1:numObsTypes
    k1=subObsTypes(l)+1;
    k2=subObsTypes(l+1);
    plot(idxPoints(k1:k2),s(k1:k2),'.')
    hold on
end
ylim([0 10])
xlabel('Points')
ylabel('St.dev.')
legend(obsTypes)
title(['Sqrt(diag) Covariance Matrix' extratitle ])

figure
for l=1:numObsTypes
    k1=subObsTypes(l)+1;
    k2=subObsTypes(l+1);
    plot(idxEpochs(k1:k2),s(k1:k2),'.')
    hold on
end
ylim([0 10])
xlabel('Epochs')
ylabel('St.dev.')
legend(obsTypes)
title(['Sqrt(diag) Covariance Matrix' extratitle])
        
% Plot square root of diagonal in space time matrix format

s=sqrt(diag(Q));

for l=1:numObsTypes
    figure
    tmp=nan(numPoints,numEpochs);
    k1=subObsTypes(l)+1;
    k2=subObsTypes(l+1);
    tmp(sub2ind([numPoints numEpochs],idxPoints(k1:k2),idxEpochs(k1:k2)))=s(k1:k2);
    if nargin > 8
       k1=subBaseObsTypes(l)+1;
       k2=subBaseObsTypes(l+1);
       tmp(sub2ind([numPoints numEpochs],idxBasePoints(k1:k2),idxBaseEpochs(k1:k2)))=zeros(k2-k1+1,1);
    end
    imagesc(tmp)
    colorbar
    clim([-5 10])
    xlabel('Epochs')
    ylabel('Points')
    title(['Sqrt(diag) Covariance Matrix ' obsTypes{l} extratitle])
end

% Dito, but focus on the extreme values

for l=1:numObsTypes
    figure
    tmp=nan(numPoints,numEpochs);
    k1=subObsTypes(l)+1;
    k2=subObsTypes(l+1);
    tmp(sub2ind([numPoints numEpochs],idxPoints(k1:k2),idxEpochs(k1:k2)))=log10(s(k1:k2));
    imagesc(tmp)
    colorbar
    xlabel('Epochs')
    ylabel('Points')
    title(['Log10(Sqrt(diag)) Covariance Matrix ' obsTypes{l} extratitle])
end

end