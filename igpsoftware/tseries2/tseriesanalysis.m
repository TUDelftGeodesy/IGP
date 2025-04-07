function tseries=tseriesanalysis(year,neu,sneu,model,env,harmonic,tjump)
%TSERIESANALYSIS Analysis of timeseries of GPS station coordinates.
%   TSERIESANALYSIS(YEAR,NEU,SNEU,MODEL,ENV,HARMONIC,TJUMP) analyses a GPS 
%   timeseries NEU with North, East, Up components. YEAR is an array with 
%   decimal years. SNEU are the (optional) standard deviations of the time 
%   series. The function determines the rate of change of the coordinates, 
%   or a more complicated model as specified in MODEL, estimates the effect 
%   temperature and atmospheric loading when temperature and/or pressure
%   are specified in ENV, estimates the harmonic periods (in years ) as 
%   specified in the array HARMONIC and removes any  jumps specified in 
%   TJUMP. The arguments SNEU, MODEL, ENV, HARMONIC and TJUMP are optional. 
%   If an empty array is used the default is assumed.
%
%   TSERIES=TSERIESANALYSIS(...) saves the results in a structure.
%
%   Examples:
%      tseriesanalysis(year,neu);
%      tseriesanalysis(year,neu,sneu);
%      tseriesanalysis(year,neu,sneu,[],[],[],[2003.207]);
%      tseriesanalysis(year,neu,sneu,2,[],[1.0]);
%
%   (c) Hans van der Marel, Delft University of Technology, 2004-2016.

%   Created:    20 April 2004 by Hans van der Marel
%   Modified:   15 June 2005 by Hans van der Marel
%                  - adapted to JPL timeseries
%               30 April 2015 by Hans van der Marel
%                  - complete overhaul and rewrite as tseriesanalysis
%                  - added pressure and temperature effects, and
%                    corrections from residual stacks and common mode
%                  - file reading part removed
%               30 September 2016 by Hans van der Marel
%                  - store reference values for pressure and temperature
%                    in structure
%               23 March 2018 by Hans van der Marel
%                  - Add spline order and period as options to the trend
%                    model, e.g. spline3y2 means spline, 3 year intervals
%                    and of order 2 (piecewise linear). Default is cubic 
%                    (order 3) with 1y interval (spline1y3)

% Check mandatory arguments

if nargin < 2, error('too few arguments (specify at least YEAR and NEU)');end

if size(year,1) ~= size(neu,1), error('YEAR and NEU must have same column length'); end
if size(year,2) ~= 1, error('YEAR must be a vector'); end
if size(neu,2) ~= 3, error('NEU must have three columns'); end

% Check optional third argument

if nargin < 3 || isempty(sneu)
  sneu=ones(size(neu));
end
if size(neu,1) ~= size(sneu,1) || size(neu,2) ~= size(sneu,2), error('SNEU and NEU must have same size'); end

% Check optional fouth argument (default is veolocity model)

if nargin < 4 || isempty(model)
   model=1;
end

% Check other optional arguments (Environment, Harmonics and Jumps) 

if nargin < 5; env=[]; end
if nargin < 6; harmonic=[]; end
if nargin < 7; tjump=[]; end

% Latitude, longitude and height

lat=neu(:,1);siglat=sneu(:,1);
lon=neu(:,2);siglon=sneu(:,2);
rad=neu(:,3);sigrad=sneu(:,3);

% Let's find out what the time limits are, set mean time, and delta time

tfirst=min(year);
tlast=max(year);

t0=round(mean(year));
dt=year-t0;
m=length(dt);


% set up the trend model

if isscalar(model)
  method='polynomial';
  pdegree=model;
elseif ischar(model)
  smodel=regexp(model,'\d*\D*','match');
  switch lower(smodel{1})
     case {'auto'}
       if tlast-tfirst > 2
         method='spline';
         pdegree=0;
         kspline=2;
         pspline=2;
       else
         method='polynomial';
         pdegree=1;
       end
     case {'vel','velocity','linear'}
       method='polynomial';
       pdegree=1;
     case {'acc','acceleration','quadratic'}
       method='polynomial';
       pdegree=2;
     case {'spline'}
       method='spline';
       pdegree=0;
       pspline=1;
       kspline=3;
       if ~isempty(smodel{2}), pspline=str2double(regexpi(smodel{2},'\d+','match')); end
       if ~isempty(smodel{3}), kspline=str2double(regexpi(smodel{3},'\d+','match')); end
    otherwise
       error('Unknown trend model.')
  end
else
  error('MODEL must be scalar or character');
end
  
switch lower(method)
  case {'polynomial'}

    % polynomial model
  
    alat=[ones(m,1) dt ];
    alon=[ones(m,1) dt ];
    arad=[ones(m,1) dt ];
    parlabel{1}='p(0)'; 
    parlabel{2}='p(1)';  
    parunit{1}='mm';
    parunit{2}='mm/y';

    n=2;
    for k=2:pdegree
      alat=[ alat dt.^k ];
      alon=[ alon dt.^k ];
      arad=[ arad dt.^k ];
      parlabel{n+1}=sprintf('p(%d)',k);
      parunit{n+1}=sprintf('mm/y%d',k);
      n=n+1;
    end
  
  case {'spline'}
    t1=min(dt);t2=max(dt);
    ts=t1:(t2-t1)/round((t2-t1)/pspline):t2;
    knots = augknt(ts,kspline);
    %[~,ik]=min(abs(repmat(dt,[1 numel(knots)-2*k])-repmat(knots(k+1:end-k),[numel(dt) 1])));
    %newknots=[knots(1:k) dt(ik)' knots(end-k+1:end) ];
    %knots
    %knots-newknots
    alat = spcol(knots,kspline,dt);
    alon = alat;
    arad = alat;
    %save('sp.mat','kspline','knots','dt','alat')
    %figure;plot(dt,alat)
    %hold on
    %plot(knots,zeros(size(knots)),'*')
    %figure;
    %plot(diff(dt))
    n=size(alat,2);
    for i=1:n, parlabel{i}=sprintf('s(%d)',i);parunit{i}='mm'; end
  otherwise
    error('Unknown trend model.')
end
parindex(1)=n;

% Set up the temperature and atmospheric loading model

if ~isempty(env)
  P=env{1};
  T=env{2};
  P0=mean(P);
  T0=mean(T);
else
  T=[];
  P=[];
  P0=0;
  T0=0;
end
if ~isempty(P)
  dP=0.1*(P-P0);
  alat = [alat dP ];
  alon = [alon dP ];
  arad = [arad dP ];
  parlabel{n+1}='AtmLd';
  parunit{n+1}='mm/kPa';
  n=n+1;
end
if ~isempty(T)
  dT=0.1*(T-T0);
  alat = [alat dT ];
  alon = [alon dT ];
  arad = [arad dT ];
  parlabel{n+1}='TempI';
  parunit{n+1}='mm/daK';
  n=n+1;
end
parindex(2)=n;

% Set up the harmonics model

for i=1:length(harmonic)
   f=1/harmonic(i);
   alat=[alat sin(2*pi*f*dt) cos(2*pi*f*dt) ];
   alon=[alon sin(2*pi*f*dt) cos(2*pi*f*dt) ];
   arad=[arad sin(2*pi*f*dt) cos(2*pi*f*dt) ];
   parlabel{n+1}=sprintf('s(%3d)',round(harmonic(i)*365));
   parlabel{n+2}=sprintf('c(%3d)',round(harmonic(i)*365));
   parunit{n+1}='mm';
   parunit{n+2}='mm';
   n=n+2;   
end
parindex(3)=n;

% Set up the jump model

for i=1:length(tjump)
   alat=[alat year > tjump(i)];
   alon=[alon year > tjump(i)];
   arad=[arad year > tjump(i)];
   parlabel{n+1}=sprintf('Jump%1d',i);
   parunit{n+1}='mm';
   n=n+1;   
end
parindex(4)=n;

% Solve the equations

alatt= alat'./(ones(n,1)*siglat.^2');qxlat=inv(alatt*alat);xlat=qxlat*alatt*lat;
alont= alon'./(ones(n,1)*siglon.^2');qxlon=inv(alont*alon);xlon=qxlon*alont*lon;
aradt= arad'./(ones(n,1)*sigrad.^2');qxrad=inv(aradt*arad);xrad=qxrad*aradt*rad;

rlat=lat-alat*xlat;
rlon=lon-alon*xlon;
rrad=rad-arad*xrad;

rmslat=sqrt(rlat'*rlat/(m-n));
rmslon=sqrt(rlon'*rlon/(m-n));
rmsrad=sqrt(rrad'*rrad/(m-n));

omtlat=rlat'*(rlat./siglat.^2)/(m-n);
omtlon=rlon'*(rlon./siglon.^2)/(m-n);
omtrad=rrad'*(rrad./sigrad.^2)/(m-n);


% Save results in a structure

tseries.tfirst=tfirst;
tseries.tlast=tlast;
tseries.t0=t0;

tseries.year=year;
tseries.neu=neu;
tseries.sneu=sneu;

tseries.model=model;
tseries.method=method;
tseries.pdegree=pdegree;
tseries.env=env;
tseries.env0=[P0 T0];
tseries.harmonic=harmonic;
tseries.tjump=tjump;

tseries.parindex=parindex;
tseries.parlabel=parlabel;
tseries.parunit=parunit;

tseries.xlat=xlat;
tseries.xlon=xlon;
tseries.xrad=xrad;

tseries.qxlat=qxlat;
tseries.qxlon=qxlon;
tseries.qxrad=qxrad;

tseries.meanstd=mean(sneu);
tseries.empstd=std(diff([ rlat rlon rrad ]))./sqrt(2);
tseries.rms=[rmslat rmslon rmsrad ];
tseries.omt=[omtlat omtlon omtrad ];

if pdegree == 0
  npar=parindex(1);
  sp=spmak(knots,[ xlat(1:npar)'; xlon(1:npar)'; xrad(1:npar)']);
  spvel=fnder(sp);  
  tseries.knots=knots;
  tseries.pspline=pspline;
  tseries.sp=sp;
  tseries.spvel=spvel;
  tseries.vel=fnval(spvel,0)';
  npar2=floor(npar/2);
  tseries.svel=[sqrt(qxlat(npar2,npar2))*sqrt(omtlat) sqrt(qxlon(npar2,npar2))*sqrt(omtlon) sqrt(qxrad(npar2,npar2))*sqrt(omtrad) ];
  %tseries.svel=[ 0.0 0.0  0.0 ];
else
  tseries.knots=[];
  tseries.pspline=[];
  tseries.sp=[];
  tseries.spvel=[];
  tseries.vel=[xlat(2) xlon(2) xrad(2) ];
  tseries.svel=[sqrt(qxlat(2,2))*sqrt(omtlat) sqrt(qxlon(2,2))*sqrt(omtlon) sqrt(qxrad(2,2))*sqrt(omtrad) ];
end
tseries.lattext=sprintf('Rate %5.2f +- %4.2f mm/year  Repeatability %3.1f mm',tseries.vel(1)*1000,tseries.svel(1)*1000,rmslat*1000);
tseries.lontext=sprintf('Rate %5.2f +- %4.2f mm/year  Repeatability %3.1f mm',tseries.vel(2)*1000,tseries.svel(2)*1000,rmslon*1000);
tseries.radtext=sprintf('Rate %5.2f +- %4.2f mm/year  Repeatability %3.1f mm',tseries.vel(3)*1000,tseries.svel(3)*1000,rmsrad*1000);

amplitude=zeros(3,length(harmonic));
npar=parindex(2);
for i=1:length(harmonic)
   amplitude(1,i)=sqrt(xlat(npar+2*(i-1)+1)^2+xlat(npar+2*(i-1)+2)^2);
   amplitude(2,i)=sqrt(xlon(npar+2*(i-1)+1)^2+xlon(npar+2*(i-1)+2)^2);
   amplitude(3,i)=sqrt(xrad(npar+2*(i-1)+1)^2+xrad(npar+2*(i-1)+2)^2);
end
tseries.amplitude=amplitude;

% Time series residuals (fit removed)

tseries.neuresidual=[rlat rlon rrad];

% Time series fit (with jumps)

tseries.neufit=      [ alat*xlat  alon*xlon arad*xrad];

% Time series components

npar=parindex(1);
tseries.neutrend=    [ alat(:,1:npar)*xlat(1:npar) ...
                       alon(:,1:npar)*xlon(1:npar) ...
                       arad(:,1:npar)*xrad(1:npar) ];
                    
nenv=npar;
if ~isempty(P)
   nenv=nenv+1;
   tseries.neuatmld=  [ alat(:,nenv)*xlat(nenv) ...
                        alon(:,nenv)*xlon(nenv) ...
                        arad(:,nenv)*xrad(nenv) ];
else
  tseries.neuatmld=zeros(m,3);
end
if ~isempty(T)
   nenv=nenv+1;
   tseries.neutempi=  [ alat(:,nenv)*xlat(nenv) ...
                        alon(:,nenv)*xlon(nenv) ...
                        arad(:,nenv)*xrad(nenv) ];
else
  tseries.neutempi=zeros(m,3);
end

nenv=parindex(2);
nharmonic=parindex(3);
if nharmonic > nenv
   tseries.neuharmonic= [ alat(:,nenv+1:nharmonic)*xlat(nenv+1:nharmonic) ...
                          alon(:,nenv+1:nharmonic)*xlon(nenv+1:nharmonic) ...
                          arad(:,nenv+1:nharmonic)*xrad(nenv+1:nharmonic) ];
else
  tseries.neuharmonic=zeros(m,3);
end

njump=parindex(4);
if njump > nharmonic
   tseries.neujumps=    [ alat(:,nharmonic+1:njump)*xlat(nharmonic+1:njump) ...
                          alon(:,nharmonic+1:njump)*xlon(nharmonic+1:njump) ...
                          arad(:,nharmonic+1:njump)*xrad(nharmonic+1:njump) ];
else
  tseries.neujumps=zeros(m,3);
end

% Print results

tseriesprint(tseries)

end

function tseriesprint(fid,tseries)
%TSERIESPRINT  Print GPS timeseries analysis results.
%   TSERIESPRINT(TSERIES) prints the results of a timeseries analysis
%   structure TSERIES produced by the function TSERIESANALYSIS.
%
%   TSERIESPRINT(FID,TSERIES) prints the results of a timeseries analysis
%   TSEREIS to file identifier FID. 
%
%   (c) Hans van der Marel, Delft University of Technology, 2004-2015.

if nargin == 1
   tseries=fid;
   fid=1;
elseif nargin < 1 || nargin > 2
   error('Must have one or two arguments')
end
if ~isstruct(tseries)
   error('first or second arguments must be a time series structure')
end

% Print the tseries results

%fprintf(fid,'\nTimeseries for %s (%8.3f - %8.3f)\n\n',upper(tseries.station),tseries.tfirst,tseries.tlast);
fprintf(fid,'\nTimeseries (%8.3f - %8.3f)\n\n',tseries.tfirst,tseries.tlast);

fprintf(fid,'       ');for k=1:length(tseries.parlabel), fprintf('%6s ',tseries.parlabel{k}); end;fprintf('   rms    omt\n') 
fprintf(fid,'       ');for k=1:length(tseries.parunit), fprintf('%6s ',tseries.parunit{k}); end;fprintf('    mm\n') 

fprintf(fid,'Latitude\n')
fprintf(fid,'       ');fprintf('%6.2f ',tseries.xlat*1000);fprintf('%6.2f %6.2f\n',tseries.rms(1)*1000,tseries.omt(1));
fprintf(fid,'       ');fprintf('%6.2f ',sqrt(diag(tseries.qxlat))*1000*sqrt(tseries.omt(1)));fprintf('\n');
fprintf(fid,'Longitude\n')
fprintf(fid,'       ');fprintf('%6.2f ',tseries.xlon*1000);fprintf('%6.2f %6.2f\n',tseries.rms(2)*1000,tseries.omt(2));
fprintf(fid,'       ');fprintf('%6.2f ',sqrt(diag(tseries.qxlon))*1000*sqrt(tseries.omt(2)));fprintf('\n');
fprintf(fid,'Height\n')
fprintf(fid,'       ');fprintf('%6.2f ',tseries.xrad*1000);fprintf('%6.2f %6.2f\n',tseries.rms(3)*1000,tseries.omt(3));
fprintf(fid,'       ');fprintf('%6.2f ',sqrt(diag(tseries.qxrad))*1000*sqrt(tseries.omt(3)));fprintf('\n\n');

for i=1:length(tseries.tjump);
   doy=rem(tseries.tjump(i),1)*365;
   fprintf('Jump%1d %8.3f (%s, doy %5.1f)\n',i,tseries.tjump(i),datestr(datenum(floor(tseries.tjump(i)),0,doy,0,0,0)),doy);
end


% empstd=std(diff(dneu))./sqrt(2)*1000;
% empstd2=std(diff(dneu2))./sqrt(2)*1000;
% meanstd=mean(sneu)*1000;
% 
% f=empstd2./meanstd;
% 
% fprintf('\n                    North [mm]  East [mm]    Up [mm]\n')
% fprintf('Emperical St.Dev.:  %10.3f %10.3f %10.3f\n',empstd);
% fprintf('Idem, after debias: %10.3f %10.3f %10.3f\n',empstd2);
% fprintf('Formal St.Dev.:     %10.3f %10.3f %10.3f\n',meanstd);
% fprintf('Factor (Estimated): %10.3f %10.3f %10.3f\n',f);
% 
% f=[ 1 1 1];
% 
% fprintf('Factor (Applied):   %10.3f %10.3f %10.3f\n\n',f);



end