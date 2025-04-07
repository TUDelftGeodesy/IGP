%% Test script for stmstochmodel - Compare with LTS2 NAM/06-GPS model
%
% Hans van der Marel
%
%% Simulate coordinates [km] and time [dyears]

numPoints=70;   
xy=[ -2+rand(numPoints,2)*4 ];

t=2001:2010;
t=t+rand(size(t))*.5+.3;

%% Generate observation matrix

npoints=size(xy,1);
obsindex=zeros(npoints*numel(t),2);
k2=0;
for k=1:numel(t)
   k1=k2+1;
   k2=k2+npoints;
   obsindex(k1:k2,:)=[ [1:npoints]' k*ones(npoints,1) ];
end

%% Test case with 3 day averages

% LTS2 implementation of recommended NAM/06-GPS model for the vertical

cov1=gpscov1(obsindex,xy*1000,t,3,'up');

% Recommended NAM/06-GPS model for the vertical in STM notation 
%
% - to be compatible we have a leadup time of 180 days
% - to be compatible we have to take into account that the LTS2 implementation 
%   uses one extra sample in the average

C=stmstochmodel({ ...
      'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=3+1/24,leadup=180)', ...  
      'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=3+1/24,leadup=180)' ...
      'spatialcov(rho=0.0887)'}, ...
      [],xy,t);   
figure;imagesc(C);colorbar
max(max(abs(C-cov1*1000*1000)))
figure;imagesc(C-cov1*1000*1000);colorbar

%% Test case with 1 day averages

% LTS2 implementation of recommended NAM/06-GPS model for the vertical

cov1=gpscov1(obsindex,xy*1000,t,1,'up');

% Recommended NAM/06-GPS model for the vertical in STM notation 
%
% - to be compatible we have a leadup time of 180 days
% - to be compatible we have to take into account that the LTS2 implementation 
%   uses one extra sample in the average

C=stmstochmodel({ ...
      'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=1+1/24,leadup=180)', ...  
      'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=1+1/24,leadup=180)' ...
      'spatialcov(rho=0.0887)'}, ...
      [],xy,t);   
figure;imagesc(C);colorbar
max(max(abs(C-cov1*1000*1000)))
figure;imagesc(C-cov1*1000*1000);colorbar

%% Functions from LTS project

function cov1=gpscov1(obsindex,pntcoord,tyear,ndays,comp)
%GPSCOV1   Make a GPS Covariance matrix.
%   COV1=GPSCOV1(OBSINDEX,PNTCOORD,TYEAR,NDAYS,OPT) creates a GPS 
%   for the observations identified by OBSINDEX. Each row in OBSINDEX 
%   represents an GPS observation: the first column contains the station
%   number and the second row contains the project number. The station 
%   numbers refer to the station coordinates PNTCOORD (in [m]), and the
%   project numbers to the project date in TYEAR (in decimal years).
%   NDAYS is the number of days a point is observed. OPT is either
%   a string with the component ('up','east','north') or a structure 
%   with the explicit parameters for the covariance model.
%
%   Example:
%       cov1=gpscov1(obsindex,pntcoord,tyear,1,'up'); 
%       cov1=gpscov1(obsindex,pntcoord,tyear,1,config.gpscov(3)); 

%% Set the default models parameters

SWmodel = [ 7 7 6 ];          % Default Simon Williams model for the North,East,Up component
rho= [ 0.0827 0.1291 0.0887]; % Default exponential decay for the North, East,Up component

doplots=true;                 % Make plots for debugging and reporting purposes

%%  Check the input arguments.

if nargin ~= 5
    error('The function requires five input arguments, see help.')
end

%% Check dimensions

numpnt=numel(pntcoord);
numprj=numel(tyear);
numobs=size(obsindex,1);

if size(obsindex,2) < 2
   error('OBSINDEX must contain at least two columns with point numbers and project numbers.')
end

if isscalar(ndays)
   ndays=ndays*ones(size(tyear));
end
if numel(ndays) ~= numel(tyear)
   error('NDAYS must be scalar or an array of the same length as TYEAR.')
end

%% Set the model parameters

if ischar(comp)
  switch comp
    case {'up','UP','u','U'}
       opt.SWmodel=SWmodel(3);
       opt.rho=rho(3);
       opt.comp='up';
       opt.doplots=doplots;
    case {'east','EAST','e','E'}
       opt.SWmodel=SWmodel(2);
       opt.rho=rho(2);
       opt.comp='east';
       opt.doplots=doplots;
    case {'north','NORTH','n','N'}
       opt.SWmodel=SWmodel(1);
       opt.rho=rho(1);
       opt.comp='north';
       opt.doplots=doplots;
    otherwise
       error('Illegal value for comp')
  end
elseif isstruct(comp)
   if ~all(ismember(fieldnames(comp),{'SWmodel','rho','comp','setupnoise','setuplevellingnoise','doplots'}))
      error('Illegal field name in option string.')
   end
   opt=comp;
else
   error('Invalid data type for opt argument.')
end

%% Temporal covariance matrix for a single GPS station
%
% Compute temporal GPS covariance matrix for a single station without
% setup noise. Setup noise will be added later as the setup noise
% will be spatially and temporally uncorrelated. 

Q=gpstemporalcov(tyear,'setup',0,'ndays',ndays,'model',opt.SWmodel);

if opt.doplots
 
   % Plot the evolution of the standard deviations

   v=sqrt(diag(Q));

   figure('Name','GPSStDev'); 
   plot(tyear,v*1e3,'-s')
   ylabel('[mm]')
   title('GPS st.dev. (single station)')

   figure('Name','GPSTempCov'); 
   imagesc(Q.*1e6); hc=colorbar;ylabel(hc,'[mm^2]')
   title('GPS Covariance matrix (single station)')

   figure('Name','GPSTempCor'); 
   imagesc(Q./(v*v')); hc=colorbar;ylabel(hc,'[-]')
   title('GPS Correlation matrix (single station)')

   % Temporal difference analysis

   % In this section all possible temporal single differences are computed 
   % and plotted

   figure('Name','GPS SD'); 
   plot(tyear,v*1e3,'-s','DisplayName','Undifferenced')
   hold on

   markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

   n=size(Q,1);
   for p=1:n
     sd=sdpivot(n,p);
     Qsd=sd*Q*sd';
     idx=1:n;idx(p)=[];
     vsd(idx)=sqrt(diag(Qsd))*1e3;
     vsd(p)=NaN;
     plot(tyear,vsd,[':' markers{mod(p,numel(markers))+1}],'DisplayName',['Pivot ',num2str(p)])
   end
   ylabel('[mm]')
   title('GPS Single-Difference st.dev. (single station)')
   legend('show')

end

%% Spatial correlations GPS stations

sCor= gpsspatialcov(pntcoord(:,1),pntcoord(:,2),ones(size(pntcoord,1),1),opt.rho);

if opt.doplots

   figure('Name','GPSSpatialCor'); 
   imagesc(sCor); hc=colorbar; ylabel(hc,'[-]')
   title('GPS Spatial Correlation matrix')

end

%% Combine spatial and temporal covariance matrices

cov1=kron(sCor,Q);

if opt.doplots

   figure('Name','GPSCampCov1a'); 
   imagesc(cov1.*1e6); hc=colorbar;ylabel(hc,'[mm^2]')
   title('GPS Covariance matrix (w/o setup noise) sorted by station')

   figure('Name','GPSCampCov1b'); 
   imagesc(kron(Q,sCor).*1e6); hc=colorbar;ylabel(hc,'[mm^2]')
   title('GPS Covariance matrix (w/o setup noise) sorted by epoch')

end

%% Remove unused observations

% The computed covariance matrix assumes numpnt*numprj observations,
% sorted in blocks of numprj, repeated numpnt times (sorted by station),
% ie. numprj is the fastest running index.
%
% The actual number of observations is numobs.

% Compute linear index for the observations vector

%obsmask=sub2ind([numpnt,numprj],obstable(:,1),obstable(:,2));
obsmask=sub2ind([numprj,numpnt],obsindex(:,2),obsindex(:,1));

cov1=cov1(obsmask,obsmask);

if opt.doplots

   figure('Name','GPSCampCov1'); 
   imagesc(cov1.*1e6); hc=colorbar;ylabel(hc,'[mm^2]')
   title('GPS Covariance matrix (w/o setup noise)')

end

end

function cov2=gpssetupcov(obsindex,cov1,opt,iscamp)
%GPSSETUPCOV   Add setup noise to GPS Covariance matrix.
%   COV2=GPSSETUPCOV(OBSDATA,COV1,SETUP,ISCAMP) adds the setup noise SETUP 
%   to the GPS covariance matrix COV1. Setup noise is only added to the
%   observations from GPS campaign stations identified by ISCAMP. OBSINDEX 
%   is the table with point id's and ISCAMP is a logical array with true 
%   for campaign stations and false for CORS stations. SETUP may also be a 
%   matlab structure with the co-variance parameters.
%
%   Example:
%       cov2=gpssetupcov(obsindex,cov1,0.0015,~iscors); 
%       cov2=gpssetupcov(obsindex,cov1,config.gpscov(3),~iscors); 

%% Set the defaults

doplots=true;            % Make plots for debugging and reporting purposes

%%  Check the input arguments.

if nargin ~= 4
   error('The function requires four input arguments, see help.')
end

%% Check dimensions

if size(obsindex,1) ~= size(cov1,1)
   error('OBSINDEX and COV1 dimensions do not match.')
end

%% Check the third parameter

if isstruct(opt)
   if ~all(ismember(fieldnames(opt),{'SWmodel','rho','comp','setupnoise','doplots'}))
      error('Illegal field name in option string.')
   end
   doplots=opt.doplots;
   setup=opt.setupnoise;
else
   setup=opt;
end


%% Add setup noise

if any(iscamp)
  % generate (diagonal) matrix with setup noise
  setup_cov1 = eye(size(cov1)) .* setup.^2;
  % blank out the CORS stations
  idx=find(~iscamp);
  numidx=numel(idx);
  setup_cov1(idx,idx)=zeros(numidx,numidx);
  % Add setup noise to input matrix
  cov2 = cov1 + setup_cov1;
else
  cov2=cov1;
end


%%  Optional plotting

if doplots

   figure('Name','GPSCampCov2'); 
   imagesc(cov2.*1e6); hc=colorbar;ylabel(hc,'[mm^2]')
   title('GPS Covariance matrix (w/ setup noise)')

   figure('Name','GPSCampCov2 (log)'); 
   imagesc(log10(abs(cov2.*1e6))); hc=colorbar;ylabel(hc,'[log mm^2]')
   title('GPS Covariance matrix (w/ setup noise)')

end

end

function cov = gpsspatialcov(x,y,dv,direction)
%GPSSPATIALCOV  Spatial covariance matrix for GPS
%  COV=GPSSPATIALCOV(X,Y,DIRECTION) calculates the covariance matrix for 
%  a set of sites, given their X and Y coordinates in meters, for the  
%  directional component DIRECTION ('north','east' or 'up'), using defaults
%  for the exponential decay.
%
%  COV=GPSSPATIALCOV(X,Y,RHO) calculates the covariance matrix for 
%  a set of sites, given their X and Y coordinates in meters, with an
%  explicit exponential decay RHO [1/km].
%
%  DV are the block diagonal matrices

if ischar(direction)
   validatestring(direction,{'north','east','up','vertical'},{'gpsspatialcov'});
   if strncmp(direction,'nor',3) 
	  rho = 0.0827;
   elseif strncmp(direction,'eas',3)
	  rho = 0.1291;
   else 
	  rho = 0.0887;
   end
elseif isfloat(direction)
   rho=direction;
else
   error('Either specify character string with direction or RHO explicitly.')
end

y = y(:);
x = x(:);
dv  = dv(:);

if length(y) ~= length(x)
	error('x and y vectors must be the same length');
end
if length(y) ~= length(dv)
	error('lat, lon and dv vectors must be the same length');
end

for j = 1:length(y)
	for k = 1:length(y)
		rng = sqrt((y(j)-y(k))^2+(x(j)-x(k))^2) * 1e-3;  % range in km
		cov(j,k) = exp(-rho * rng) * dv(j) * dv(k);
	end
end

end 

function Q =  gpstemporalcov(dyear,varargin)
%GPSTEMPORALCOV  Compute Covariance matrix for GPS station.
%   Q=GPSTEMPORALCOV(DYEAR,varargin) computes the covariance matrix C for
%   campaign GPS data based on the noise model derived for the 06-GPS/NAM 
%   datasets. At minimum it requires the time of the campaigns
%   in decimal years.
%
%   Additonal arguments are input as pairs with the first argument
%   indicating the parameter to set and the second the parameter itself.
%   The options are :
% 
%      model : specify which stochastic model to choose from between the range [1 7]
%              refer to the document "Description of GPS Uncertainties within the
%              Long Term Study on Anomalous Time-Dependent Subsidence". The default
%              model is #6
% 
%      ndays : The assumption is that every campaign is 5 days long. If you know the
%              measurement length for each campaign then you can input it as a vector
%              of the same length as the time and data vectors. The vector should be
%              a set of integers specifying the number of days.
%
%      setup : set the white noise setup error. By default this is assumed to be ~1.5 mm.
%
%      leadup: set the lead-up time (days) for the full continuous sequence. By default 
%              this is assumed to be 180 days.
%
%   Example:
%   
%      Q = gpstemporalcov(dyear,'model',1,'setup',0);
%
% Created:  June 2015 by Simon Williams
% Modified: 22 August 2016 by Hans van der Marel
%             - stripped estimation part, only generation of covariance matrix
%             - renamed to gpstemporalcov
%           31 October 2016 by Hans van der Marel
%             - Fixed bug in original Simon Williams code
%             - Modified period of full continuous sequence and accept
%               lead-up time as input argument
%             - Update documentation

nvar = nargin - 1;

if mod(nvar,2) == 1
	error('Incorrect number of optional arguments : arguments need to be in pairs.');
end

dyear = dyear(:);

npoints = length(dyear);

gotNdays  = 0;
gotModel  = 0;
gotSetup  = 0;
gotLeadup = 0; 
gotAnnual = 0; 

for j = 1:2:nvar
	switch lower(varargin{j})
		case 'ndays'
			ndays = varargin{j+1};
			ndays = ndays(:);
			validateattributes(ndays,{'numeric'},{'nonempty','size',[npoints 1],'integer'});
                        gotNdays = 1;
		case 'model'
			model = varargin{j+1};
			validateattributes(model,{'numeric'},{'nonempty','>=',1,'<=',7,'size',[1 1]});
                        gotModel = 1;
		case 'setup'
			setup = varargin{j+1};
			validateattributes(setup,{'numeric'},{'nonempty','size',[1 1]});
                        gotSetup = 1;
		case 'leadup'
			leadup = varargin{j+1};
			validateattributes(leadup,{'numeric'},{'nonempty','size',[1 1]});
                        gotLeadup = 1;
		case 'periodic'
			annual = varargin{j+1};
			validateattributes(annual,{'numeric'},{'nonempty','size',[1 1]});
                        gotAnnual = 1;
            disp('This option is not supported')
		otherwise
            error(sprintf('Unknown input argument : %s',varargin{j}))
        end
end


if gotNdays == 0
	% without specifying the number of campaign days per campaign then
	% we are assuming 5 days 
	ndays = ones(npoints,1)*5;
end

if gotModel == 0
	% the default model is #6 from the document
	% "Description of GPS uncertainties within
	% "the Long Term Study on Anomalous Time Dependent
	% Subsidence"
	model = 6;
end

if gotSetup == 0
	% the default amplitude of setup error is ~5 mm
	%setup = 4.4970;  % original sdwil
    setup = 1.5;      % modified by H. Baehr (NAM); see footnote in report, page 41
end

if gotLeadup == 0
	% if you do not want to take into account the effect of the annual signal then 
	% just dont set the annual in the arguments
	leadup = 180;
end

if gotAnnual == 0
	% if you do not want to take into account the effect of the annual signal then 
	% just dont set the annual in the arguments
	annual = 0;
end

% Convert decimal year to MJD

mjd=51544.5+(dyear-2000)*365.25;
%mjd=fix(mjd+1e-7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mjd_all = [];
for j = 1:npoints
	% for simplicity we will assume the hours will always be centred around the date given
	% and because ndays * 24 is always a even number we will add one more hour
	id = (0:ndays(j)*24)' - ndays(j)*12;
	mjd_all = [mjd_all ; mjd(j) + id ./ 24];
end

% we create a full continous sequence from leadup days before first mjd to 2 days after last mjd

mjd_start = min(mjd_all)-leadup;    % -180 in original code
mjd_end   = max(mjd_all)+1;         % +180 in original code;

m = (mjd_start:1/24:mjd_end)';

nfull = length(m);

T = zeros(npoints,nfull);
for j = 1:npoints
	id = (0:ndays(j)*24)' - ndays(j)*12;
	m1 = mjd(j) + id ./ 24;
	% ij = round(m1 - mjd_start) .* 24 + 1;           % bug corrected 31/10/2016
    ij = unique(round((m1 - mjd_start) .* 24) + 1);  
	T(j,ij) = ones(length(id),1) ./ length(id);
end

if model == 1
	h1 = iGGM(nfull,-1.91,-0.74,0.9861) * 0.0967;
	C = multiply(T,h1);
elseif model == 2
	h1 = iGGM(nfull,-1.91,-0.74,0.9861) * 0.0967;
	h2 = hosking(nfull,-2) * 0.005340;
	C = multiply(T,h1) + multiply(T,h2);
elseif model == 3
	h1 = iGGM(nfull,-1.91,-0.74,0.9861) * 0.0967;
	h3 = hosking(nfull,-2) * 0.005340*2;
	C = multiply(T,h1) + multiply(T,h3);
elseif model == 4
	h4 = iGGM(nfull,-1.91,-0.74,0.9861) * 0.0967 * 3;
	C = multiply(T,h4);
elseif model == 5
	h5 = iGGM(nfull,-1.91,-0.8525,0.997380) * 0.0967;
	C = multiply(T,h5);
elseif model == 6
	h1 = iGGM(nfull,-1.91,-0.74,0.9861) * 0.0967;
	h6 = iGGM(nfull,-1.91,-0.8525,0.999096) * 0.04;
	C = multiply(T,h1) + multiply(T,h6);
elseif model == 7                                      % Is Model 3 for horizontal data 
	h1 = iGGM(nfull,-1.91,-0.70,0.9137) * 0.1481;
	h3 = iGGM(nfull,-2,-0.8281,0.992) * 0.06;
	C = multiply(T,h1) + multiply(T,h3);
else
	error('Model out of acceptable range [1 6]');
end

Q = C ./ 1000 ./ 1000;
Q = Q + eye(size(Q)) .* (setup./1000).^2;

end

function [h]=hosking(length,spectral_index)

%    [h]=hosking(length,spectral_index)
%
%    Generate a differencing/integrating vector of any order including
%    fractional.  These are the discrete-time versions of fractional
%    differentiation and integration.
%
%    length    size of vector
%    spectral_index spectral index 
%
%    h    output vector

% See Hosking, R. M., 1981, Biometrika, 68, pp. 165-176

% Hadley Johnson, hjohnson@ucsd.edu, 28-May-1998 version
% Simon Williams  sdwil@noc.ac.uk, 13-July-2007

degree = spectral_index / 2;
h=cumprod([1;[0:1:length-2]'-degree]./[1;[1:1:length-1]']);

end

function [h] = iGGM(N,alphaR,alphaL,gamma)

%    [h]=iGGM(N,alphaR,alphaL,gamma)
%
%    N    size of vector
%    alphaR spectral index to right of cross-over 
%    alphaL spectral index to left of cross-over 
%    gamma gamma for ar(1) (0 -> 1) 
%
%    h    output vector
%
%    Simon Williams : June 2015

import sdwil.*;  % added by H. Baehr (NAM)

hl = hosking(N,alphaL);
alpha1 = alphaR - alphaL;
hr = hosking(N,alpha1);

g = gamma.^(0:N-1)';

h = fftfilt(hl,g.*hr);

end

function C = multiply(T,h)

N = length(h);
M = size(T,1);
nc = size(T,2);

if length(h) ~= size(T,2)
	error('Dimensions of h and T are incompatable')
end

np = 2^nextpow2(N);

%tic
c = [h(1);zeros(N-1,1);zeros(2*np-(2*N-1),1);flipud(h(2:end))];

S = zeros(2*np,M);
for j = 1:M
	A = zeros(2*np,1);
	A(1:N) = T(j,:)';
	S(:,j) = ifft(fft(c).*fft(A));
end

S = S(1:N,:)';

C = S*S';
%toc 

end 

function [ sd, id] = sdpivot(n,ip)
%SDPIVOT  Generate single difference transformation matrix.
%   SD=SDPIVOT(N,IP) conputes the single difference transformation
%   matrix SD with dimensions (N-1) x N , with IP the pivot number.
%
%   [SD,ID]=SDPIVOT(N,IP) also returns the array ID with the non-pivots.

sd=eye(n-1,n);
sd=[ sd(:,1:ip-1) -ones(n-1,1) sd(:,ip:n-1)];

id=ones(n,1);
id(ip)=[];

end



