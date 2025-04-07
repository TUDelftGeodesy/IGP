function varargout=stmplotcov(stm,varargin)
%STMPLOTCOV  Plot space time matrix (STM) covariance matrix.
%  H=STMPLOTCOV(STM) creates a pseudo-color image of the covariance
%  matrix for the space time matrix structure STM. STM can be a a 
%  space time matrix structure or the name of a space time dataset. The
%  function returns the plot handle(s) H.
%
%  H=STMPLOTCOV(...,'option',value,...) allows for extra options
%
%    'type'               Output type {'covmatrix','cholesky','precmatrix'} (default 'covmatrix')
%
%    'ROI'                Region of interest, as [latmin lonmin ; latmax lonmax] 
%                         bounding box, or lat/lon polygon (Default [], is all)
%    'POI'                Period of interest [ dYearStart dYearEnd ] (Default [-Inf +Inf])
%    'pntMask'            Point mask (default [], is all), will be combined with ROI
%    'epochMask'          Epoch mask (default [], is all), will be combined with POI
%
%    'verbose'            Verbose level (default 0)
%
%  Example:
%
%    stmplotcov('../2_reduce/gron_levelling_flaggedOutliers_v2.3.02_20191030_reduced.mat')
%
%  See also: stmstochmodel
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:   2 Oct 2023 by Hans van der Marel
% Modified: 

% Check input arguments

if nargin < 1
   error('This function expects at least one STM input argument.')
end

if ischar(stm)
   stm=stmread(stm);
elseif ~isstruct(stm)
   error('The first argument must be a character string with the filename or space time matrix structure.')
end

% Process options

opt.type='covmatrix';         % Output type {'covmatrix','cholesky','precmatrix'} (default 'covmatrix')
opt.chol=false;

opt.ROI=[];                   % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default none)
opt.POI=[-Inf +Inf];          % Period of interest [ dYearStart dYearEnd ] (Default none)
opt.pntMask=[];               % Point mask (default [], is all), will be combined with ROI
opt.epochMask=[];             % Epoch mask (default [], is all), will be combined with POI
opt.verbose=0;                % Verbose level (default 0)

opt.pntCrdType='deg/m';       % Coordinate type in <stm>.pntCrd {'deg/m', 'km/m'}, default 'deg/m' (Latitude, Longitude, Height)

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end

% Number of points, epochs, coordinates and epoch times for this dataset

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

% Prepatory work for stmstochmodel

numObsTypes=length(stm.obsTypes);

numStochModels=numel(stm.stochModel);
if numStochModels ~= numObsTypes
   fprintf('Number of stochatic models (%d) is smaller than number of observation types (%d)\n',numStochModels,numObsTypes)
end

if strncmp(stm.stochModel{1},'cholfactor',10)
   if numStochModels ~= 1 && numObsTypes > 1
       fprintf('Single Cholesky factor for multiple observation types (%d).\n',numObsTypes)
   end
   U = stm.stochData;
   %U = U + tril(nan(size(U)),-1);
   
   figure
   imagesc(U, 'AlphaData', 1-0.5*isnan(U))

   size(U)
   size(stm.obsData)
   ymasks=cell(1,3);
   for l=1:numObsTypes
       %y=stm.obsData(pntMask,epochMask,l);
       y=stm.obsData(pntMask,epochMask,l);
       y=y(:);    
       ymasks{l}=~isnan(y);
   end
   opt.chol=false
end

% Plot covariance matrices (or Cholesky factors) for the different observations types

h(1)=figure('Units','normalized','OuterPosition',[0.03 0.03 0.3*numObsTypes 0.95]);
t=tiledlayout(6,numObsTypes,'TileSpacing','compact','Padding','compact');  

for l=1:numObsTypes

    % Get space-time matrix with observations 
    
    y=stm.obsData(pntMask,epochMask,l);
         
    % Vectorize observed minus a-priori and mask the non-observed (nan) entries
    
    y=y(:);    
    ymask=~isnan(y);
    y=y(ymask);
    my=length(y);
    
    % Point and Epoch number matrices
    
    yp=repmat([1:numPoints]',[1,numEpochs]);
    ye=repmat(1:numEpochs,[numPoints,1]);
    yp=yp(pntMask,epochMask);
    yp=yp(ymask);
    ye=ye(pntMask,epochMask);
    ye=ye(ymask);
    
    
    % Compose co-variance matrix of observations
    
    lstoch=min(l,numel(stm.stochModel));   % PATCH FOR INSAR DATASETS
    stm.stochModel{lstoch}
    
    if any(strncmp(stm.stochModel{lstoch},'tudinsar4rd',11))
       % stochModel{l}='tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)';  % there is a problem with the tudinsar4rd model (not positive definite)
       pntAttrib=stm.pntAttrib;
       epochAttrib=stm.epochAttrib;
       Qy=stmstochmodel(stm.stochModel{lstoch},stm.stochData,pntNeu,epochDyear,{pntMask,epochMask},pntAttrib,epochAttrib);
       Qy=Qy(ymask,ymask);
    elseif strncmp(stm.stochModel{lstoch},'cholfactor',10)
       U=stm.stochData;
       opt.chol=false
       Qy=U + tril(nan(size(U)),-1);
       size(U)
    else 
       Qy=stmstochmodel(stm.stochModel{lstoch},stm.stochData,pntNeu,epochDyear,{pntMask,epochMask});
    end
    Qy=Qy(ymask,ymask);

    % Select between Covariance matrix and Cholesky factor

    if opt.chol
        Qy=chol(Qy) + tril(nan(size(Qy)),-1);
    end

    % Pseudo-color plot of covariance matrix
    
    nexttile(t,l,[3,1])
    imagesc(Qy,'AlphaData', 1-0.5*isnan(Qy));
    colorbar
    title(stm.obsTypes{l},'interpreter','none')
    
    % Plot point and epoch numbers

    nexttile(t,l+3*numObsTypes)
    plot(yp,'.')
    xlim([1,my])
    ylabel('Point number')
    
    nexttile(t,l+4*numObsTypes)
    plot(ye,'.')
    xlim([1,my])
    ylabel('Epoch number')

    % Plot square root of diagonal

    nexttile(t,l+5*numObsTypes)
    plot(sqrt(diag(Qy)),'-')
    xlim([1,my])
    ylabel('sqrt(diag)')
    
end

title(t,stm.datasetId ,'interpreter','none');

% Special analysis for block diagonal matrices

if any(strncmp(stm.stochModel{1},'covmatrix(format=blkdiag',24))
    h(2)=figure('Units','normalized','OuterPosition',[0.03 0.03 0.9 0.95]);

    t=tiledlayout('Flow','TileSpacing','compact','Padding','compact');  

    for l=1:numObsTypes

       fprintf('epo  #pnts  rd (cols..)    rcond1   rcond2    rank defect; all zero; diagonal only; sub-net\n')
       for i=1:stm.numEpochs
          C=covplot(stm,i);
       end
    end
end

if nargout > 0
    varargout{1}=h;
end

end

function C=covplot(stm,iepoch)

epochlabel= [ 'Epoch ' num2str(iepoch) ': ' stm.epochAttrib.prjName{iepoch}  '  (' num2str(stm.epochDyear(iepoch)) ')'];

ymask=~isnan(stm.obsData(:,iepoch,1));
C=stm.stochData(ymask,ymask,iepoch);
pntName=stm.pntName(ymask);

nexttile
imagesc(C)
colorbar
hold on
plot(size(C,1)-sqrt(diag(C))*size(C,1)/max(sqrt(diag(C))),'-w')
title(epochlabel,'interpreter','none')

nc=size(C,1);

IC = C == 0;

% Find columns/rows with all zero's
ir = find( all(IC,1) & all(IC,2)' );

% Find columns/rows with only diagonal elements
id = find( sum(IC,1) == nc-1 & sum(IC,2)' == nc-1 );

% Find columns/rows with many zero's 
ixcrit=floor(size(IC,1)/10);
lx= sum(IC,1) > ixcrit & sum(IC,2)' > ixcrit ;
ix=find(lx);
% Check that outside there are only zero elements
% any(any(IC(~lx,lx) ~=1)) & any(any(IC(~lx,lx) ~=1))
% remove all zero's and diagonals
lx(ir)=false;
lx(id)=false;
ix=find(lx);

if length(ix) > 2
   figure
   imagesc(C(ix,ix))
   title(['sub-net' num2str(iepoch)])
   colorbar
end

% determine rank-defects

p=Inf;
p1=0;
psave=[];
while p > 0
  [U,p]=chol(C(p1+1:end,p1+1:end));
  psave= [psave p];
  p1=p;  
end
psave=psave(1:end-1);
psave=cumsum(psave);
rd=length(psave);
rdcols=['(' num2str(psave) ')'];

% regularization

[U,p]=chol(C+ones(nc,nc)*1e-9);
if p == 0
    rcond=min(diag(U))/max(diag(U));
else
    fprintf('Epoch %d: Unexpected rankdefect after regularization: pivot %d\n',iepoch,p)
    rcond=0;
end

% idealization precision

[U,p]=chol(C+eye(nc,nc)*1e-9);
if p == 0
    rcondi=min(diag(U))/max(diag(U));
else
    fprintf('Epoch %d: Very unexpected rankdefect with diagonal regularization: pivot %d\n',iepoch,p)
    rcondi=0;
end

fprintf('%3d  %5d  %2d %-10s %7.5f  %7.5f    %s  ;  %s  ;  %s  ;  %s\n',iepoch,nc,rd,rdcols,rcond,rcondi,strjoin(pntName(psave),', '),strjoin(pntName(ir),', '),strjoin(pntName(id),', '),strjoin(pntName(ix),', '))

end
